"""Núcleo numérico — Projeto de Sistemas de Controle.

Sem dependência de Streamlit: tudo aqui pode ser testado/importado isoladamente.
Portado do app HTML (script.js) mantendo as mesmas convenções.
"""
from __future__ import annotations

import ast
import math
import re

import numpy as np

# ----------------------------------------------------------------------------
# Álgebra linear básica
# ----------------------------------------------------------------------------


def eps_for(A) -> float:
    A = np.asarray(A, float)
    m = float(np.max(np.abs(A))) if A.size else 0.0
    return max(1e-12, m * 1e-9)


def mat_inverse(A):
    """Inversa ou None se (numericamente) singular."""
    A = np.asarray(A, float)
    try:
        s = np.linalg.svd(A, compute_uv=False)
        if s[-1] <= max(1e-12, s[0] * 1e-11):
            return None
        return np.linalg.inv(A)
    except np.linalg.LinAlgError:
        return None


def det(A) -> float:
    return float(np.linalg.det(np.asarray(A, float)))


def mat_rank(A) -> int:
    A = np.asarray(A, float)
    return int(np.linalg.matrix_rank(A, tol=eps_for(A)))


def controllability_matrix(A, B):
    A, B = np.asarray(A, float), np.asarray(B, float)
    n = A.shape[0]
    cols, term = [B], B
    for _ in range(1, n):
        term = A @ term
        cols.append(term)
    return np.hstack(cols)


def observability_matrix(A, C):
    A, C = np.asarray(A, float), np.asarray(C, float)
    n = A.shape[0]
    rows, term = [C], C
    for _ in range(1, n):
        term = term @ A
        rows.append(term)
    return np.vstack(rows)


# ----------------------------------------------------------------------------
# Polinômios
# ----------------------------------------------------------------------------


def faddeev_leverrier(A):
    """Coeficientes de det(sI-A) (mônico) + passos intermediários."""
    A = np.asarray(A, float)
    n = A.shape[0]
    I = np.eye(n)
    M = I.copy()
    c = [1.0]
    steps = []
    for k in range(1, n + 1):
        AM = A @ M
        tr = float(np.trace(AM))
        ck = -tr / k
        c.append(ck)
        Mn = AM + ck * I
        steps.append(dict(k=k, Mprev=M.copy(), AM=AM, tr=tr, ck=ck, M=Mn))
        M = Mn
    return np.array(c), steps


def poly_roots(coeffs):
    coeffs = np.asarray(coeffs, float)
    if len(coeffs) <= 1:
        return np.array([], complex)
    return np.roots(coeffs)


def poly_from_roots(roots):
    return np.real(np.poly(np.asarray(roots, complex)))


def eval_poly_matrix(coeffs, A):
    """Horner matricial: coeffs[0]*A^n + ... + coeffs[n]*I."""
    A = np.asarray(A, float)
    n = A.shape[0]
    I = np.eye(n)
    R = coeffs[0] * I
    for c in coeffs[1:]:
        R = R @ A + c * I
    return R


def unpaired_complex_poles(roots, tol=1e-6):
    roots = list(roots)
    used = [False] * len(roots)
    unpaired = []
    for i, r in enumerate(roots):
        if used[i] or abs(r.imag) < 1e-7:
            used[i] = True
            continue
        found = -1
        for j, q in enumerate(roots):
            if used[j] or j == i:
                continue
            if abs(q.real - r.real) < tol and abs(q.imag + r.imag) < tol:
                found = j
                break
        if found == -1:
            unpaired.append(r)
        else:
            used[i] = used[found] = True
    return unpaired


# ----------------------------------------------------------------------------
# Ackermann
# ----------------------------------------------------------------------------


def ackermann_k(A, B, desired):
    """K = -[0..0 1] Wc^-1 q_c(A)   (lei u = Kx, malha A+BK)."""
    A, B = np.asarray(A, float), np.asarray(B, float)
    n = A.shape[0]
    qc = poly_from_roots(desired)
    Wc = controllability_matrix(A, B)
    Wc_inv = mat_inverse(Wc)
    if Wc_inv is None:
        return dict(error="Sistema não controlável (det(W_c) = 0): não é possível alocar "
                          "os polos arbitrariamente por realimentação de estados.")
    qcA = eval_poly_matrix(qc, A)
    en = np.zeros((1, n))
    en[0, -1] = 1
    K = -(en @ Wc_inv @ qcA)
    return dict(K=K, qc=qc, Wc=Wc, Wc_inv=Wc_inv, qcA=qcA)


def ackermann_l(A, C, desired):
    """L = q_o(A) Wo^-1 [0..0 1]^T   (observador A-LC)."""
    A, C = np.asarray(A, float), np.asarray(C, float)
    n = A.shape[0]
    qo = poly_from_roots(desired)
    Wo = observability_matrix(A, C)
    Wo_inv = mat_inverse(Wo)
    if Wo_inv is None:
        return dict(error="Sistema não observável (det(W_o) = 0): não é possível alocar "
                          "os polos do observador arbitrariamente.")
    qoA = eval_poly_matrix(qo, A)
    en = np.zeros((n, 1))
    en[-1, 0] = 1
    L = qoA @ Wo_inv @ en
    return dict(L=L, qo=qo, Wo=Wo, Wo_inv=Wo_inv, qoA=qoA)


# ----------------------------------------------------------------------------
# Discretização
# ----------------------------------------------------------------------------


def matrix_exp(M):
    """exp(M) por escalonamento-e-quadrados + série de Taylor (30 termos)."""
    M = np.asarray(M, float)
    n = M.shape[0]
    norm = float(np.max(np.sum(np.abs(M), axis=1))) if M.size else 0.0
    s = 0
    while norm > 0.5:
        norm /= 2
        s += 1
    Ms = M / (2.0 ** s)
    term = np.eye(n)
    result = np.eye(n)
    for k in range(1, 31):
        term = term @ Ms / k
        result = result + term
    for _ in range(s):
        result = result @ result
    return result


def discretize(A, B, T):
    A, B = np.asarray(A, float), np.asarray(B, float)
    n, m = A.shape[0], B.shape[1]
    M = np.block([[A * T, B * T], [np.zeros((m, n)), np.zeros((m, m))]])
    E = matrix_exp(M)
    return dict(M=M, E=E, G=E[:n, :n], H=E[:n, n:n + m])


# ----------------------------------------------------------------------------
# Formas canônicas
# ----------------------------------------------------------------------------


def controllable_canonical(a, b):
    """a = [α1..αn] (den. mônico), b = [β1..βn]."""
    n = len(a)
    Ac = np.zeros((n, n))
    for i in range(n - 1):
        Ac[i, i + 1] = 1
    for j in range(n):
        Ac[n - 1, j] = -a[n - 1 - j]
    Bc = np.zeros((n, 1))
    Bc[n - 1, 0] = 1
    Cc = np.array([list(b)[::-1]], float)
    return Ac, Bc, Cc


def observable_canonical(a, b):
    n = len(a)
    Ao = np.zeros((n, n))
    for i in range(n - 1):
        Ao[i + 1, i] = 1
    for i in range(n):
        Ao[i, n - 1] = -a[n - 1 - i]
    Bo = np.array([[b[n - 1 - i]] for i in range(n)], float)
    Co = np.zeros((1, n))
    Co[0, -1] = 1
    return Ao, Bo, Co


def canonical_controllable_from_ss(A, B, C):
    """Transformação de similaridade T = (Wc·Ŵ)^-1 (Ogata)."""
    A, B, C = (np.asarray(x, float) for x in (A, B, C))
    n = A.shape[0]
    coeffs, _ = faddeev_leverrier(A)
    a = coeffs[1:]
    Wc = controllability_matrix(A, B)
    seq = list(a[: n - 1][::-1]) + [1.0]
    What = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            idx = i + j
            What[i, j] = seq[idx] if idx < n else 0.0
    T = mat_inverse(Wc @ What)
    if T is None:
        return dict(error="Sistema não controlável: não é possível obter a forma canônica "
                          "controlável por transformação de similaridade.")
    Tinv = mat_inverse(T)
    return dict(Ac=T @ A @ Tinv, Bc=T @ B, Cc=C @ Tinv, T=T, a=a, Wc=Wc, What=What)


# ----------------------------------------------------------------------------
# Seguidor de referência (servo)
#
# Convenções (e = y − r; ξ = integrador do erro):
#   contínuo  : ξ̇ = e = Cx − r            A_a = [[0, C],[0, A]]     B_a = [[0],[B]]
#   discreto  : ξ(k+1) = ξ(k) + e(k)       A_a = [[1, C],[0, G]]     B_a = [[0],[H]]
#   "alt" (Ogata)     : K_a = −[0…1] W_a⁻¹ q_c(A_a);   u = k₁ξ + K₂x        (A_a + B_a K_a)
#   "UFRN" contínuo   : K_a = +[0…1] W_a⁻¹ q_c(A_a);   u = −k₁ξ − K₂x       (A_a − B_a K_a)
#   "UFRN" discreto   : método Ĝ/Ĥ/M do formulário; u(k) = k₁x(k) + k₂ξ(k), ξ(k+1) = ξ(k) + y(k+1) − r(k+1)
# Em todos os casos a malha fechada REAL é reconstruída e seus polos verificados.
# ----------------------------------------------------------------------------


def servo_augmented(A, B, C, discrete=False):
    A, B, C = (np.asarray(x, float) for x in (A, B, C))
    n = A.shape[0]
    Aa = np.block([[np.full((1, 1), 1.0 if discrete else 0.0), C], [np.zeros((n, 1)), A]])
    Ba = np.vstack([np.zeros((1, 1)), B])
    return Aa, Ba


def servo_alt(A, B, C, desired, discrete=False):
    """Ogata: K_a = −[0…1] W_a⁻¹ q_c(A_a);  u = k₁ξ + K₂x  (malha A_a + B_a K_a)."""
    Aa, Ba = servo_augmented(A, B, C, discrete)
    res = ackermann_k(Aa, Ba, desired)
    if "error" in res:
        return res
    K = res["K"]
    return dict(Aa=Aa, Ba=Ba, K=K, k1=float(K[0, 0]), k2=K[:, 1:], Acl=Aa + Ba @ K,
                **{k: res[k] for k in ("qc", "Wc", "Wc_inv", "qcA")})


def servo_continuous_ufrn(A, B, C, desired):
    """UFRN: K_a = +[0…1] W_a⁻¹ q_c(A_a);  u = −k₁ξ − K₂x  (malha A_a − B_a K_a)."""
    base = servo_alt(A, B, C, desired, discrete=False)
    if "error" in base:
        return base
    out = dict(base)
    out["K"] = -base["K"]
    out["k1"] = -base["k1"]
    out["k2"] = -base["k2"]
    out["Acl"] = base["Aa"] - base["Ba"] @ out["K"]
    return out


def servo_discrete_ufrn(G, H, C, desired):
    """Método UFRN/Ogata discreto.

    Ĝ = [[G, H],[0, 0]],  Ĥ = [0…0 1]ᵀ,  K̂ = +[0…1] W_c(Ĝ)⁻¹ q_c(Ĝ)   (malha Ĝ − Ĥ K̂)
    M = [[G − I, H],[CG, CH]],   [K | −k_I] = (K̂ + [0…0 1]) M⁻¹
    u(k) = k₁ x(k) + k₂ ξ(k)  com  k₁ = −K,  k₂ = k_I.
    """
    G, H, C = (np.asarray(x, float) for x in (G, H, C))
    n = G.shape[0]
    Ghat = np.block([[G, H], [np.zeros((1, n + 1))]])
    Hhat = np.zeros((n + 1, 1))
    Hhat[n, 0] = 1
    res = ackermann_k(Ghat, Hhat, desired)
    if "error" in res:
        return dict(error=res["error"])
    Khat = -res["K"]                      # K̂ positivo (Ackermann clássico)
    M = np.block([[G - np.eye(n), H], [C @ G, C @ H]])
    Minv = mat_inverse(M)
    if Minv is None:
        return dict(error="Matriz M = [[G−I, H],[CG, CH]] singular — verifique se (G,H) é controlável "
                          "e se C·H ≠ 0.")
    e_last = np.zeros((1, n + 1))
    e_last[0, n] = 1
    Kp = (Khat + e_last) @ Minv           # [K | −k_I]
    K = Kp[:, :n]
    kI = -float(Kp[0, n])
    k1, k2 = -K, kI
    A11 = G + H @ k1
    Acl = np.block([[A11, H * k2], [C @ A11, np.array([[1.0 + float((C @ H)[0, 0]) * k2]])]])
    return dict(Ghat=Ghat, Hhat=Hhat, Khat=Khat, M=M, Minv=Minv, Kp=Kp, K=K, kI=kI,
                k1=k1, k2=k2, Acl=Acl,
                **{k: res[k] for k in ("qc", "Wc", "Wc_inv", "qcA")})


# ----------------------------------------------------------------------------
# Parsers de entrada
# ----------------------------------------------------------------------------


def normalize_minus(s: str) -> str:
    return re.sub(r"[\u2212\u2013\u2014]", "-", s)


def parse_complex(token: str):
    t = normalize_minus(token).strip().lower().replace(" ", "")
    if not t:
        return None
    t = t.replace("i", "j")
    if "j" not in t:
        try:
            return complex(float(t), 0.0)
        except ValueError:
            return None
    try:
        return complex(t)
    except ValueError:
        # casos como "+j" / "-j" / "2+j" tratados acima; o resto é inválido
        return None


def parse_complex_list(text: str):
    """Lista separada por vírgula/ponto-e-vírgula. Aceita 'a±bi' (par conjugado)."""
    out = []
    for tok in re.split(r"[,;]", normalize_minus(text)):
        tok = tok.strip()
        if not tok:
            continue
        if "±" in tok:
            re_part, im_part = tok.split("±", 1)
            z = parse_complex(f"{re_part}+{im_part}" if im_part.strip().lower().endswith(("i", "j"))
                              else f"{re_part}+{im_part}i")
            if z is None:
                continue
            out += [z, z.conjugate()] if abs(z.imag) > 1e-12 else [z]
            continue
        z = parse_complex(tok)
        if z is not None:
            out.append(z)
    return out


def parse_number_list(text: str):
    vals = []
    for tok in re.split(r"[,;\s]+", normalize_minus(text).strip()):
        if tok == "":
            continue
        try:
            vals.append(safe_eval(tok))
        except ValueError:
            continue
    return vals


_FUNCS = {"sqrt": math.sqrt, "exp": math.exp, "sin": math.sin, "cos": math.cos,
          "tan": math.tan, "log": math.log, "ln": math.log, "abs": abs}
_CONSTS = {"pi": math.pi, "e": math.e}


def _ev(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
        return float(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        v = _ev(node.operand)
        return v if isinstance(node.op, ast.UAdd) else -v
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)):
        a, b = _ev(node.left), _ev(node.right)
        if isinstance(node.op, ast.Add):
            return a + b
        if isinstance(node.op, ast.Sub):
            return a - b
        if isinstance(node.op, ast.Mult):
            return a * b
        if isinstance(node.op, ast.Div):
            return a / b
        if abs(b) > 100:
            raise ValueError("expoente muito grande")
        return a ** b
    if isinstance(node, ast.Name) and node.id in _CONSTS:
        return _CONSTS[node.id]
    if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id in _FUNCS and len(node.args) == 1 and not node.keywords):
        return float(_FUNCS[node.func.id](_ev(node.args[0])))
    raise ValueError("expressão não permitida")


def safe_eval(text) -> float:
    """Avalia número/expressão simples: 0.5, 0,5, -1/3, 2*pi, sqrt(2), 1e-3, 2^3."""
    if isinstance(text, (int, float, np.floating, np.integer)):
        return float(text)
    s = normalize_minus(str(text)).strip()
    if s == "":
        raise ValueError("vazio")
    s = re.sub(r"(?<=\d),(?=\d)", ".", s)   # vírgula decimal
    s = s.replace("^", "**")
    try:
        val = _ev(ast.parse(s, mode="eval").body)
    except (SyntaxError, ZeroDivisionError, OverflowError) as exc:
        raise ValueError(str(exc)) from exc
    if not math.isfinite(val):
        raise ValueError("valor não finito")
    return float(val)


def has_bad_number(*mats) -> bool:
    return any(not np.all(np.isfinite(np.asarray(m, float))) for m in mats)
