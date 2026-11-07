"""Relatórios passo a passo (LaTeX) para cada módulo da bancada."""
from __future__ import annotations

import numpy as np

from latex import Report, Step, cnum, frac, mat, num, poly, roots_list, smat
from nucleo import (
    ackermann_k, ackermann_l, canonical_controllable_from_ss, controllability_matrix,
    controllable_canonical, det, discretize, faddeev_leverrier, has_bad_number, mat_inverse,
    mat_rank, observability_matrix, observable_canonical, poly_from_roots, poly_roots,
    servo_alt, servo_continuous_ufrn, servo_discrete_ufrn,
)

EN_ROW = r"\begin{bmatrix} 0 & \cdots & 0 & 1 \end{bmatrix}"
EN_COL = r"\begin{bmatrix} 0 \\ \vdots \\ 0 \\ 1 \end{bmatrix}"


# ----------------------------------------------------------------------------
# Utilidades
# ----------------------------------------------------------------------------


def stable(poles, discrete: bool) -> bool:
    return all(abs(p) < 1 - 1e-9 for p in poles) if discrete else all(p.real < -1e-9 for p in poles)


def same_poles(a, b, tol=1e-3) -> bool:
    a, b = list(np.asarray(a, complex)), list(np.asarray(b, complex))
    if len(a) != len(b):
        return False
    for x in a:
        j = min(range(len(b)), key=lambda k: abs(b[k] - x), default=None)
        if j is None or abs(b[j] - x) > tol * max(1.0, abs(x)):
            return False
        b.pop(j)
    return True


def factor_tex(r, var="s") -> str:
    r = complex(r)
    if abs(r.imag) < 1e-9:
        if abs(r.real) < 1e-9:
            return var
        return f"({var} + {num(-r.real)})" if r.real < 0 else f"({var} - {num(r.real)})"
    return f"({var} - ({cnum(r)}))"


def char_poly_name(var):  # det(sI - A)
    return rf"\det({var}I - A)"


def _term_label(c, p, a="A") -> str:
    base = "I" if p == 0 else (a if p == 1 else f"{a}^{{{p}}}")
    if p == 0:
        return rf"{num(c)}\,I"
    return base if abs(c - 1) < 1e-9 else rf"{num(c)}\,{base}"


def detail_qA(st: Step, coeffs, A, q: str, a: str = "A"):
    """q(A): fórmula, potências de A, cada termo aₖ·Aᵖ e a soma."""
    coeffs = [float(c) for c in coeffs]
    deg = len(coeffs) - 1
    A = np.asarray(A, float)
    n = A.shape[0]
    P = [np.eye(n)]
    for _ in range(deg):
        P.append(P[-1] @ A)
    st.tex(rf"{q}({a}) = {poly(coeffs, a, ident='I')}")
    if deg >= 2:
        st.text("**Potências de** $" + a + "$:")
        st.mats(*[(f"{a}^{{{k}}}", P[k]) for k in range(2, deg + 1)])
    st.text("**Cada termo** (coeficiente × potência):")
    total = np.zeros((n, n))
    pairs = []
    for i, c in enumerate(coeffs):
        p = deg - i
        T = c * P[p]
        total += T
        pairs.append((_term_label(c, p, a), T))
    st.mats(*pairs)
    st.text("**Soma dos termos:**")
    st.mats((f"{q}({a})", total))
    return total


def detail_inverse(st: Step, W, Winv, w: str, pick: str | None = None):
    """W, det, W⁻¹, conferência W·W⁻¹ = I e (opcional) a linha/coluna usada pelo Ackermann."""
    W = np.asarray(W, float)
    n = W.shape[0]
    I = W @ Winv
    err = float(np.max(np.abs(I - np.eye(n))))
    st.mats((w, W))
    st.tex(rf"\det({w}) = {num(det(W), 5)} \neq 0 \;\Rightarrow\; {w}\ \text{{é inversível}}")
    st.mats((f"{w}^{{-1}}", Winv))
    st.tex(rf"{w}\,{w}^{{-1}} = I \qquad (\text{{erro máximo}} = {err:.1e})")
    if pick == "row":
        en = np.zeros((1, n)); en[0, -1] = 1
        st.text("O Ackermann usa só a **última linha** de $" + w + "^{-1}$:")
        st.mats((rf"{EN_ROW}\,{w}^{{-1}}", en @ Winv))
    elif pick == "col":
        en = np.zeros((n, 1)); en[-1, 0] = 1
        st.text("O Ackermann usa só a **última coluna** de $" + w + "^{-1}$:")
        st.mats((rf"{w}^{{-1}}\,{EN_COL}", Winv @ en))


def desired_poly_step(rep: Report, var, desired, q="q_c", title=None):
    qc = poly_from_roots(desired)
    st = rep.step(title or f"Polinômio característico desejado ${q}({var})$")
    st.tex(rf"{q}({var}) = " + "".join(factor_tex(r, var) for r in desired) + rf" = {poly(qc, var)}")
    return qc


def compare_poly_table(rep: Report, A, var, desired, qc, who="desejado"):
    cur, _ = faddeev_leverrier(A)
    cur_roots = poly_roots(cur)
    st = rep.step(f"Polinômio atual de $A$ × polinômio {who}",
                  why="A realimentação desloca os autovalores de $A$ para os polos escolhidos, "
                      "sem mudar a estrutura física da planta.")
    st.tex(
        r"\begin{array}{l|l|l}"
        r" & \text{Polinômio} & \text{Raízes} \\ \hline "
        rf"\text{{Atual}} & {poly(cur, var)} & {roots_list(cur_roots)} \\ "
        rf"\text{{{who.capitalize()}}} & {poly(qc, var)} & {roots_list(desired)}"
        r"\end{array}")


def ctrl_step(rep: Report, A, B, name="W_c", title="Controlabilidade", sysname="(A,B)", a="A", b="B"):
    n = np.asarray(A).shape[0]
    Wc = controllability_matrix(A, B)
    r = mat_rank(Wc)
    st = rep.step(title, why=f"O par {sysname} é controlável se, e somente se, ${name}$ tem posto completo $n$.")
    st.tex(rf"{name} = \begin{{bmatrix}} {b} & {a}{b} & {a}^2{b} & \cdots & {a}^{{n-1}}{b} \end{{bmatrix}}")
    st.mats((name, Wc))
    st.tex(rf"\det({name}) = {num(det(Wc), 5)} \qquad \operatorname{{posto}}({name}) = {r}\ \text{{de}}\ {n}")
    st.status(r == n, "CONTROLÁVEL", "NÃO CONTROLÁVEL")
    return Wc, r == n


def obs_step(rep: Report, A, C, name="W_o", title="Observabilidade"):
    n = np.asarray(A).shape[0]
    Wo = observability_matrix(A, C)
    r = mat_rank(Wo)
    st = rep.step(title, why=f"O par $(A,C)$ é observável se, e somente se, ${name}$ tem posto completo $n$.")
    st.tex(rf"{name} = \begin{{bmatrix}} C \\ CA \\ CA^2 \\ \vdots \\ CA^{{n-1}} \end{{bmatrix}}")
    st.mats((name, Wo))
    st.tex(rf"\det({name}) = {num(det(Wo), 5)} \qquad \operatorname{{posto}}({name}) = {r}\ \text{{de}}\ {n}")
    st.status(r == n, "OBSERVÁVEL", "NÃO OBSERVÁVEL")
    return Wo, r == n


# ----------------------------------------------------------------------------
# 1. Análise
# ----------------------------------------------------------------------------


def rep_analise(A, B, C, discrete=False) -> Report:
    n = A.shape[0]
    var = "z" if discrete else "s"
    rep = Report(f"Polos, Estabilidade, Controlabilidade e Observabilidade (n = {n})")

    st = rep.step("Sistema informado")
    st.mats(("A", A), ("B", B), ("C", C))
    st.tex(r"x(k+1) = A\,x(k) + B\,u(k)" if discrete else r"\dot{x}(t) = A\,x(t) + B\,u(t)")
    st.tex(r"y(k) = C\,x(k)" if discrete else r"y(t) = C\,x(t)")

    coeffs, steps = faddeev_leverrier(A)
    st = rep.step("Polinômio característico (Faddeev–LeVerrier)",
                  why="O método obtém os coeficientes sem calcular nenhum determinante: em cada "
                      "iteração o traço de $AM_{k-1}$ fornece o próximo coeficiente.")
    st.tex(r"M_0 = I, \qquad c_k = -\dfrac{\operatorname{tr}(A\,M_{k-1})}{k}, \qquad M_k = A\,M_{k-1} + c_k\,I")
    for s in steps:
        k = s["k"]
        st.text(f"**Iteração $k = {k}$**")
        st.mats((f"A\\,M_{{{k-1}}}", s["AM"]))
        st.tex(rf"\operatorname{{tr}}(A\,M_{{{k-1}}}) = {num(s['tr'])} \;\Rightarrow\; "
               rf"c_{{{k}}} = -\dfrac{{{num(s['tr'])}}}{{{k}}} = {num(s['ck'])}")
        st.mats((f"M_{{{k}}} = A\\,M_{{{k-1}}} + c_{{{k}}}I", s["M"]))
    st.tex(rf"\det({var}I - A) = {poly(coeffs, var)}")

    roots = poly_roots(coeffs)
    st = rep.step("Polos e estabilidade")
    crit = (r"|z| < 1 \quad (\text{interior do círculo unitário})" if discrete
            else r"\mathrm{Re}(s) < 0 \quad (\text{semiplano esquerdo})")
    st.tex(rf"\text{{Critério de estabilidade: }} {crit}")
    colname = r"|z|" if discrete else r"\mathrm{Re}(s)"
    rows = []
    for i, r in enumerate(roots, 1):
        v = abs(r) if discrete else r.real
        ok = abs(r) < 1 - 1e-9 if discrete else r.real < -1e-9
        rows.append(rf"p_{{{i}}} & {cnum(r)} & {num(v)} & \text{{{'estável' if ok else 'instável'}}}")
    st.tex(r"\begin{array}{c|l|l|l} & \text{Polo} & " + colname + r" & \text{Situação} \\ \hline "
           + r" \\ ".join(rows) + r"\end{array}")
    st.status(stable(roots, discrete), "SISTEMA ESTÁVEL", "SISTEMA INSTÁVEL")

    ctrl_step(rep, A, B, "W_c", "Matriz de controlabilidade")
    obs_step(rep, A, C, "W_o", "Matriz de observabilidade")
    return rep


# ----------------------------------------------------------------------------
# 2. Realimentação de estados
# ----------------------------------------------------------------------------


def rep_feedback(A, B, desired, discrete=False) -> Report:
    n = A.shape[0]
    var = "z" if discrete else "s"
    rep = Report(f"Realimentação de Estados — Fórmula de Ackermann (n = {n})")

    st = rep.step("Sistema e polos desejados de malha fechada")
    st.mats(("A", A), ("B", B))
    st.tex(rf"\{{p_1,\dots,p_{n}\}} = \{{{roots_list(desired)}\}}")

    Wc, ok = ctrl_step(rep, A, B, "W_c", "Verificação de controlabilidade")
    res = ackermann_k(A, B, desired)
    if "error" in res:
        return rep.fail(res["error"])

    qc = desired_poly_step(rep, var, desired, "q_c")
    compare_poly_table(rep, A, var, desired, qc)

    st = rep.step("Avaliação de $q_c(A)$",
                  why="Substitui-se $s$ por $A$ no polinômio desejado (o termo constante vira $a_n I$).")
    detail_qA(st, qc, A, "q_c", "A")

    st = rep.step("Inversa de $W_c$",
                  why="Só o último vetor-linha de $W_c^{-1}$ entra na fórmula de Ackermann.")
    detail_inverse(st, Wc, res["Wc_inv"], "W_c", "row")

    K = res["K"]
    if has_bad_number(K):
        return rep.fail("Cálculo com valores não finitos — sistema mal-condicionado ($W_c$ quase singular).")
    en = np.zeros((1, n)); en[0, -1] = 1
    st = rep.step("Fórmula de Ackermann")
    st.tex(rf"K = -{EN_ROW}\,W_c^{{-1}}\,q_c(A)")
    st.mats((rf"{EN_ROW}\,W_c^{{-1}}", en @ res["Wc_inv"]), ("q_c(A)", res["qcA"]))
    st.mats(("K", K))
    st.text("Lei de controle: $u = Kx$ (o sinal negativo já está embutido em $K$), "
            "portanto a malha fechada é $A + BK$.")

    Acl = A + B @ K
    poles = np.roots(faddeev_leverrier(Acl)[0])
    st = rep.step("Resultado final", final=True)
    st.mats(("K", K))
    st.text("**Verificação** — matriz de malha fechada $A+BK$:")
    st.mats(("A + BK", Acl))
    st.tex(rf"\text{{polos de }} A+BK:\ {roots_list(poles)}")
    st.status(same_poles(poles, desired), "Os polos de malha fechada coincidem com os desejados.",
              "Os polos de malha fechada NÃO coincidem com os desejados (verifique os dados).")
    st.tex(r"u(k) = K\,x(k)" if discrete else r"u(t) = K\,x(t)")
    return rep


# ----------------------------------------------------------------------------
# 3. Observador de estados
# ----------------------------------------------------------------------------


def rep_observer(A, C, desired, discrete=False) -> Report:
    n = A.shape[0]
    var = "z" if discrete else "s"
    rep = Report(f"Observador de Estados — Ackermann dual (n = {n})")

    st = rep.step("Sistema e polos desejados do observador")
    st.mats(("A", A), ("C", C))
    st.tex(rf"\{{p_1,\dots,p_{n}\}} = \{{{roots_list(desired)}\}}")

    Wo, ok = obs_step(rep, A, C, "W_o", "Verificação de observabilidade")
    res = ackermann_l(A, C, desired)
    if "error" in res:
        return rep.fail(res["error"])

    qo = desired_poly_step(rep, var, desired, "q_o", f"Polinômio característico desejado do observador $q_o({var})$")
    compare_poly_table(rep, A, var, desired, qo, who="do observador")
    rep.steps[-1].why = ("Por convenção, os polos do observador ficam 3–10× mais rápidos que os do "
                         "controlador, para o erro de estimação decair antes de o controlador agir.")

    st = rep.step("Avaliação de $q_o(A)$")
    detail_qA(st, qo, A, "q_o", "A")

    st = rep.step("Inversa de $W_o$", why="Dual do controlador: aqui entra a última **coluna** de $W_o^{-1}$.")
    detail_inverse(st, Wo, res["Wo_inv"], "W_o", "col")

    L = res["L"]
    if has_bad_number(L):
        return rep.fail("Cálculo com valores não finitos — sistema mal-condicionado ($W_o$ quase singular).")
    en = np.zeros((n, 1)); en[-1, 0] = 1
    st = rep.step("Fórmula de Ackermann (dual)",
                  why="É a fórmula da realimentação de estados aplicada ao par $(A^T, C^T)$ e transposta.")
    st.tex(rf"L = q_o(A)\,W_o^{{-1}}\,{EN_COL}")
    st.mats(("q_o(A)", res["qoA"]), (rf"W_o^{{-1}}{EN_COL}", res["Wo_inv"] @ en))
    st.mats(("L", L))

    Aobs = A - L @ C
    poles = np.roots(faddeev_leverrier(Aobs)[0])
    st = rep.step("Resultado final", final=True)
    st.mats(("L", L))
    st.text("**Verificação** — matriz do observador $A-LC$:")
    st.mats(("A - LC", Aobs))
    st.tex(rf"\text{{polos de }} A-LC:\ {roots_list(poles)}")
    st.status(same_poles(poles, desired), "Os polos do observador coincidem com os desejados.",
              "Os polos do observador NÃO coincidem com os desejados (verifique os dados).")
    st.tex(r"\hat{x}(k+1) = A\hat{x}(k) + B u(k) + L\,[\,y(k) - C\hat{x}(k)\,]" if discrete
           else r"\dot{\hat{x}} = A\hat{x} + B u + L\,[\,y - C\hat{x}\,]")
    st.text("**Princípio da separação:** com $K$ e $L$ projetados separadamente e $u = K\\hat{x}$, os polos do "
            "sistema completo são a união dos polos de $A+BK$ com os de $A-LC$.")
    return rep


# ----------------------------------------------------------------------------
# 4. Discretização
# ----------------------------------------------------------------------------


def rep_discretize(A, B, T) -> Report:
    n = A.shape[0]
    rep = Report(f"Discretização com retentor de ordem zero (T = {num(T)} s, n = {n})")
    if T <= 0:
        rep.warnings.append("T ≤ 0 — o resultado será degenerado (G = I, H = 0). Use um período positivo.")

    st = rep.step("Sistema contínuo e período de amostragem")
    st.mats(("A", A), ("B", B))
    st.tex(rf"T = {num(T)}\ \mathrm{{s}}")
    st.tex(r"G = e^{AT}, \qquad H = \left(\int_0^T e^{A\tau}\,d\tau\right) B, \qquad x(k+1) = G\,x(k) + H\,u(k)")

    d = discretize(A, B, T)
    st = rep.step("Matriz aumentada",
                  why="Truque para obter $G$ e $H$ de uma só vez com uma exponencial de matriz; funciona "
                      "mesmo quando $A$ é singular, ao contrário de $H = A^{-1}(G-I)B$.")
    st.tex(r"M = \begin{bmatrix} AT & BT \\ 0 & 0 \end{bmatrix}, \qquad "
           r"e^{M} = \begin{bmatrix} G & H \\ 0 & I \end{bmatrix}")
    st.mats(("M", d["M"]))

    if has_bad_number(d["E"]):
        return rep.fail("A exponencial não convergiu (valores não finitos) — verifique se $AT$ não é grande demais.")
    st = rep.step("Exponencial $e^{M}$",
                  why="Calculada por escalonamento e quadrados com série de Taylor (30 termos).")
    st.mats(("e^{M}", d["E"]))

    st = rep.step("Resultado final", final=True)
    st.mats(("G", d["G"]), ("H", d["H"]))
    st.tex(r"x(k+1) = G\,x(k) + H\,u(k)")
    return rep


# ----------------------------------------------------------------------------
# 5. Formas canônicas
# ----------------------------------------------------------------------------


def _sym_forms(n):
    z, o = "0", "1"
    Ac = [[o if j == i + 1 else z for j in range(n)] for i in range(n - 1)]
    Ac.append([rf"-\alpha_{{{n - j}}}" for j in range(n)])
    Bc = [[z]] * (n - 1) + [[o]]
    Cc = [[rf"\beta_{{{n - j}}}" for j in range(n)]]
    Ao = []
    for i in range(n):
        row = [o if j == i - 1 else z for j in range(n)]
        row[n - 1] = rf"-\alpha_{{{n - i}}}"
        Ao.append(row)
    Bo = [[rf"\beta_{{{n - i}}}"] for i in range(n)]
    Co = [[z] * (n - 1) + [o]]
    return Ac, Bc, Cc, Ao, Bo, Co


def _canonical_common(rep: Report, a, b, D=0.0):
    n = len(a)
    Ac, Bc, Cc, Ao, Bo, Co = _sym_forms(n)
    Ac_n, Bc_n, Cc_n = controllable_canonical(a, b)
    Ao_n, Bo_n, Co_n = observable_canonical(a, b)

    st = rep.step("Coeficientes $\\alpha$ (denominador) e $\\beta$ (numerador)")
    st.tex(rf"G(s) = \dfrac{{\beta_1 s^{{{n-1}}} + \cdots + \beta_{{{n}}}}}{{s^{{{n}}} + \alpha_1 s^{{{n-1}}} + \cdots + \alpha_{{{n}}}}}")
    st.tex(r"\alpha = \begin{bmatrix} " + " & ".join(num(v) for v in a) + r" \end{bmatrix}, \qquad "
           r"\beta = \begin{bmatrix} " + " & ".join(num(v) for v in b) + r" \end{bmatrix}")

    for title, (sA, sB, sC), (nA, nB, nC) in (
        ("Forma Canônica Controlável", (Ac, Bc, Cc), (Ac_n, Bc_n, Cc_n)),
        ("Forma Canônica Observável", (Ao, Bo, Co), (Ao_n, Bo_n, Co_n)),
    ):
        st = rep.step(title, final=True)
        st.text("**Estrutura geral:**")
        st.tex(r" \qquad ".join([smat(sA, "A"), smat(sB, "B"), smat(sC, "C")]))
        st.text("**Valores numéricos:**")
        st.mats(("A", nA), ("B", nB), ("C", nC))
        if abs(D) > 1e-12:
            st.tex(rf"y = Cx + Du, \qquad D = {num(D)}")


def rep_canonical_tf(den, numer) -> Report:
    rep = Report("Formas Canônicas a partir da função de transferência")
    den, numer = list(den), list(numer)
    if len(den) < 2:
        return rep.fail("Informe um denominador de grau ≥ 1 (ao menos 2 coeficientes, ex.: `1, 2` para $s+2$).")
    if abs(den[0]) < 1e-12:
        return rep.fail("O coeficiente líder do denominador não pode ser zero.")
    num_raw, den_raw = list(numer if numer else [0.0]), list(den)
    note = None
    if abs(den[0] - 1) > 1e-9:
        lead = den[0]
        den = [v / lead for v in den]
        numer = [v / lead for v in numer]
        note = f"Todos os coeficientes foram divididos por {num(lead)} para tornar o denominador mônico."
    n = len(den) - 1
    a = den[1:]
    bw = list(numer)
    D = 0.0
    if len(bw) == n + 1:
        D = bw[0]
        bw = [v - D * den[i] for i, v in enumerate(bw)][1:]
    elif len(bw) > n + 1:
        return rep.fail("Função de transferência imprópria: o numerador tem grau maior que o denominador.")
    bw = [0.0] * (n - len(bw)) + bw
    b = bw

    st = rep.step("Função de transferência informada")
    st.tex(rf"G(s) = \dfrac{{{poly(num_raw)}}}{{{poly(den_raw)}}}")
    if note:
        st.text(note)
    st.text(f"Ordem do sistema: $n = {n}$.")
    if abs(D) > 1e-12:
        st.text(f"Grau do numerador = grau do denominador: o termo de transmissão direta $D = {num(D)}$ "
                "foi extraído por divisão polinomial.")
    _canonical_common(rep, a, b, D)
    return rep


def rep_canonical_ss(A, B, C) -> Report:
    n = A.shape[0]
    rep = Report(f"Formas Canônicas a partir das variáveis de estado (n = {n})")
    st = rep.step("Sistema em variáveis de estado informado")
    st.mats(("A", A), ("B", B), ("C", C))

    ctrl = canonical_controllable_from_ss(A, B, C)
    if "error" in ctrl:
        return rep.fail(ctrl["error"])
    if has_bad_number(ctrl["T"], ctrl["Ac"]):
        return rep.fail("A transformação resultou em valores não finitos — sistema quase não controlável.")
    a = ctrl["a"]
    b = list(ctrl["Cc"][0])[::-1]

    st = rep.step("Polinômio característico de $A$",
                  why="Calculado por Faddeev–LeVerrier; fornece $\\alpha_1,\\dots,\\alpha_n$.")
    st.tex(rf"\det(sI - A) = {poly([1.0] + list(a), 's')}")

    st = rep.step("Transformação de similaridade $T$ (Ogata)",
                  why="Nova variável de estado $z = Tx$: preserva os autovalores de $A$ e expõe a estrutura canônica.")
    st.tex(r"T = (W_c\,\widehat{W})^{-1}, \qquad \widehat{W} = \text{matriz de Toeplitz dos coeficientes } \alpha")
    st.mats(("W_c", ctrl["Wc"]), (r"\widehat{W}", ctrl["What"]), ("T", ctrl["T"]))

    st = rep.step(r"Sistema transformado: $\tilde A = TAT^{-1}$, $\tilde B = TB$, $\tilde C = CT^{-1}$",
                  why=r"Os coeficientes $\beta$ do numerador aparecem em $\tilde C = [\beta_n \cdots \beta_1]$.")
    st.mats((r"\tilde A", ctrl["Ac"]), (r"\tilde B", ctrl["Bc"]), (r"\tilde C", ctrl["Cc"]))

    _canonical_common(rep, list(a), b)
    return rep


# ----------------------------------------------------------------------------
# 6. Seguidor de referência
# ----------------------------------------------------------------------------


def _augmented_step(rep, A, B, C, discrete):
    n = A.shape[0]
    if discrete:
        st = rep.step("Sistema aumentado com integrador do erro",
                      why="O integrador $\\xi$ acumula o erro $e = y - r$ e garante erro nulo em regime "
                          "permanente para referência em degrau.")
        st.tex(r"\xi(k+1) = \xi(k) + y(k) - r(k), \qquad "
               r"\begin{bmatrix}\xi(k+1)\\x(k+1)\end{bmatrix} = A_a \begin{bmatrix}\xi(k)\\x(k)\end{bmatrix} + B_a\,u(k) + \cdots")
        st.tex(r"A_a = \begin{bmatrix} 1 & C \\ 0 & G \end{bmatrix}, \qquad B_a = \begin{bmatrix} 0 \\ H \end{bmatrix}")
    else:
        st = rep.step("Sistema aumentado com integrador do erro",
                      why="O integrador $\\xi$ acumula o erro $e = y - r$ e garante erro nulo em regime "
                          "permanente para referência em degrau.")
        st.tex(r"\dot\xi = e = Cx - r, \qquad "
               r"\begin{bmatrix}\dot\xi\\\dot x\end{bmatrix} = A_a \begin{bmatrix}\xi\\x\end{bmatrix} + B_a\,u + \cdots")
        st.tex(r"A_a = \begin{bmatrix} 0 & C \\ 0 & A \end{bmatrix}, \qquad B_a = \begin{bmatrix} 0 \\ B \end{bmatrix}")
    return st


def rep_servo(A, B, C, desired, discrete=False, convention="ufrn") -> Report:
    n = A.shape[0]
    var = "z" if discrete else "s"
    ufrn = convention == "ufrn"
    Am, Bm = ("G", "H") if discrete else ("A", "B")
    conv_label = ("UFRN (formulário oficial)" if ufrn else
                  ("Clássico (Ogata/Franklin)" if discrete else "Alternativo (Ogata)"))
    rep = Report(f"Seguidor de Referência — {'Discreto' if discrete else 'Contínuo'} — {conv_label} (n = {n})")

    st = rep.step("Sistema original e polos desejados do seguidor")
    st.mats((Am, A), (Bm, B), ("C", C))
    st.tex(rf"\{{p_1,\dots,p_{n+1}\}} = \{{{roots_list(desired)}\}} \qquad (n+1 = {n+1}\ \text{{polos}})")
    st.text(f"Convenção: **{conv_label}**.")

    # ---------------- discreto UFRN: método Ĝ/Ĥ/M ----------------
    if discrete and ufrn:
        st = rep.step("Sistema auxiliar (formulário UFRN)",
                      why="O sistema auxiliar embute $u$ na dinâmica e permite usar o Ackermann padrão sobre $(\\hat G,\\hat H)$.")
        st.tex(r"\hat G = \begin{bmatrix} G & H \\ 0 & 0 \end{bmatrix}, \qquad \hat H = \begin{bmatrix} 0 \\ \vdots \\ 0 \\ 1 \end{bmatrix}")
        st.tex(r"\xi(k+1) = \xi(k) + y(k+1) - r(k+1), \qquad u(k) = k_1\,x(k) + k_2\,\xi(k)")
        res = servo_discrete_ufrn(A, B, C, desired)
        if "error" in res:
            return rep.fail(res["error"])
        if has_bad_number(res["Ka"] if "Ka" in res else res["Kp"]):
            return rep.fail("Cálculo com valores não finitos — sistema mal-condicionado.")

        st = rep.step("Matrizes do sistema auxiliar")
        st.mats((r"\hat G", res["Ghat"]), (r"\hat H", res["Hhat"]))

        Wc_hat, _ = ctrl_step(rep, res["Ghat"], res["Hhat"], r"W_c(\hat G)",
                              "Controlabilidade do par $(\\hat G, \\hat H)$", sysname=r"$(\hat G,\hat H)$",
                              a=r"\hat G", b=r"\hat H")

        qc = desired_poly_step(rep, var, desired, "q_c")
        st = rep.step(r"Avaliação de $q_c(\hat G)$")
        detail_qA(st, qc, res["Ghat"], "q_c", r"\hat G")
        st = rep.step(r"Inversa de $W_c(\hat G)$")
        detail_inverse(st, Wc_hat, res["Wc_inv"], r"W_c(\hat G)", "row")

        st = rep.step(r"Ganho auxiliar $\hat K$",
                      why="Neste método $\\hat K$ usa o sinal positivo do Ackermann clássico (malha $\\hat G - \\hat H\\hat K$).")
        st.tex(rf"\hat K = {EN_ROW}\,W_c(\hat G)^{{-1}}\,q_c(\hat G)")
        st.mats((r"\hat K", res["Khat"]))

        st = rep.step("Matriz de conversão $M$ e sua inversa",
                      why="$M$ converte os ganhos do espaço auxiliar $(\\hat G)$ para o espaço físico $(x,\\xi)$.")
        st.tex(r"M = \begin{bmatrix} G - I & H \\ CG & CH \end{bmatrix}")
        st.mats(("M", res["M"]), ("M^{-1}", res["Minv"]))

        st = rep.step("Ganhos finais")
        st.tex(r"\begin{bmatrix} K & -k_I \end{bmatrix} = \left(\hat K + \begin{bmatrix} 0 & \cdots & 0 & 1 \end{bmatrix}\right) M^{-1}"
               r", \qquad k_1 = -K, \quad k_2 = k_I")
        st.mats((r"\begin{bmatrix} K & -k_I \end{bmatrix}", res["Kp"]))
        st.tex(rf"k_1 = {mat(res['k1'])} \qquad k_2 = {num(res['k2'])}")

        poles = np.linalg.eigvals(res["Acl"])
        st = rep.step("Resultado final — servo discreto (método UFRN)", final=True)
        st.tex(rf"k_1\ (\text{{ganho dos estados}}) = {mat(res['k1'])}")
        st.tex(rf"k_2\ (\text{{ganho do integrador }}\xi) = {num(res['k2'])}")
        st.text("**Verificação** — malha fechada real nas variáveis $(x,\\xi)$:")
        st.tex(r"\begin{bmatrix} x(k+1) \\ \xi(k+1) \end{bmatrix} = "
               r"\begin{bmatrix} G + Hk_1 & Hk_2 \\ C(G+Hk_1) & 1 + CHk_2 \end{bmatrix}"
               r"\begin{bmatrix} x(k) \\ \xi(k) \end{bmatrix} + \cdots")
        st.mats(("A_{cl}", res["Acl"]))
        st.tex(rf"\text{{polos de malha fechada:}}\ {roots_list(poles)}")
        st.status(same_poles(poles, desired), "Os polos de malha fechada coincidem com os desejados.",
                  "Os polos de malha fechada NÃO coincidem com os desejados.")
        st.tex(r"u(k) = k_1\,x(k) + k_2\,\xi(k), \qquad \xi(k+1) = \xi(k) + y(k+1) - r(k+1)")
        return rep

    # ---------------- contínuo (UFRN/alt) e discreto clássico ----------------
    _augmented_step(rep, A, B, C, discrete)
    res = (servo_continuous_ufrn(A, B, C, desired) if (ufrn and not discrete)
           else servo_alt(A, B, C, desired, discrete))
    if "error" in res:
        # ainda mostra a controlabilidade para o aluno ver o motivo
        from nucleo import servo_augmented
        Aa, Ba = servo_augmented(A, B, C, discrete)
        ctrl_step(rep, Aa, Ba, "W_a", "Controlabilidade do sistema aumentado", sysname=r"$(A_a,B_a)$",
                  a="A_a", b="B_a")
        return rep.fail(res["error"])
    Aa, Ba = res["Aa"], res["Ba"]
    st = rep.step("Matrizes do sistema aumentado")
    st.mats(("A_a", Aa), ("B_a", Ba))
    Wa, _ = ctrl_step(rep, Aa, Ba, "W_a", "Controlabilidade do sistema aumentado", sysname=r"$(A_a,B_a)$",
                      a="A_a", b="B_a")
    if has_bad_number(res["K"]):
        return rep.fail("Cálculo com valores não finitos — sistema mal-condicionado.")

    qc = desired_poly_step(rep, var, desired, "q_c")
    st = rep.step(r"Avaliação de $q_c(A_a)$")
    detail_qA(st, qc, Aa, "q_c", "A_a")
    st = rep.step(r"Inversa de $W_a$", why="Só a última linha de $W_a^{-1}$ participa do Ackermann.")
    detail_inverse(st, Wa, res["Wc_inv"], "W_a", "row")

    st = rep.step(f"Ganho aumentado $K_a$ — convenção {conv_label}")
    if ufrn and not discrete:
        st.tex(rf"K_a = +{EN_ROW}\,W_a^{{-1}}\,q_c(A_a) = \begin{{bmatrix}} k_1 & K_2 \end{{bmatrix}}")
        st.text("Sem o sinal negativo do Ackermann clássico: o sinal aparece na lei de controle.")
    else:
        st.tex(rf"K_a = -{EN_ROW}\,W_a^{{-1}}\,q_c(A_a) = \begin{{bmatrix}} k_1 & K_2 \end{{bmatrix}}")
    st.mats(("K_a", res["K"]))
    st.tex(rf"k_1 = {num(res['k1'])} \qquad K_2 = {mat(res['k2'])}")

    poles = np.linalg.eigvals(res["Acl"])
    title = ("Resultado final — convenção UFRN" if ufrn else "Resultado final — convenção clássica")
    st = rep.step(title, final=True)
    st.tex(rf"k_1\ (\text{{integrador}}) = {num(res['k1'])}")
    st.tex(rf"K_2\ (\text{{estados}}) = {mat(res['k2'])}")
    if discrete:
        cl = r"A_a + B_a K_a"
        st.text("**Verificação** — malha fechada $A_a + B_aK_a$:")
    elif ufrn:
        cl = r"A_a - B_a K_a"
        st.text("**Verificação** — malha fechada $A_a - B_aK_a$ (a lei UFRN tem sinais negativos):")
    else:
        cl = r"A_a + B_a K_a"
        st.text("**Verificação** — malha fechada $A_a + B_aK_a$:")
    st.mats((cl, res["Acl"]))
    st.tex(rf"\text{{polos de malha fechada:}}\ {roots_list(poles)}")
    st.status(same_poles(poles, desired), "Os polos de malha fechada coincidem com os desejados.",
              "Os polos de malha fechada NÃO coincidem com os desejados.")
    if discrete:
        st.tex(r"u(k) = k_1\,\xi(k) + K_2\,x(k), \qquad \xi(k+1) = \xi(k) + y(k) - r(k)")
    elif ufrn:
        st.tex(r"u(t) = -k_1 \int_0^t e(\tau)\,d\tau - K_2\,x(t), \qquad e = y - r")
    else:
        st.tex(r"u(t) = k_1 \int_0^t e(\tau)\,d\tau + K_2\,x(t), \qquad e = y - r")
        st.text("Na notação original de Ogata, com $\\int(r-y)\\,dt$, o ganho integral troca de sinal: $k_I = -k_1$.")
    return rep
