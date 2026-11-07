"""Bancada de Projeto de Sistemas de Controle — app Streamlit.

Executar:  streamlit run app.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

from exemplos import EXEMPLOS
from latex import Report, to_latex, to_markdown
from modulos import (
    rep_analise, rep_canonical_ss, rep_canonical_tf, rep_discretize, rep_feedback,
    rep_observer, rep_servo,
)
from nucleo import (
    parse_complex_list, parse_number_list, safe_eval, unpaired_complex_poles,
)
from latex import cnum

st.set_page_config(page_title="Bancada de Controle", page_icon="🎛️", layout="wide")

PAGES = {
    "analise": "1 · Análise (polos, controlabilidade, observabilidade)",
    "feedback": "2 · Realimentação de estados (Ackermann)",
    "observer": "3 · Observador de estados",
    "discretize": "4 · Discretização",
    "canonical": "5 · Formas canônicas",
    "servo": "6 · Seguidor de referência",
    "exemplos": "📚 Exemplos das listas",
}
PREFIX = {"analise": "an", "feedback": "fb", "observer": "ob", "discretize": "dz",
          "canonical": "cn", "servo": "sv"}
DOM_C, DOM_D = "Contínuo (s)", "Discreto (z)"
MODE_TF, MODE_SS = "Função de transferência", "Variáveis de estado"
CONV_UFRN, CONV_ALT = "ufrn", "alt"

ss = st.session_state


# ----------------------------------------------------------------------------
# Estado e exemplos
# ----------------------------------------------------------------------------


def apply_example(ex: dict):
    m = ex["module"]
    p = PREFIX[m]
    init = {k: ex[k] for k in ("A", "B", "C") if k in ex}
    ss[f"{p}_init"] = init
    ss[f"{p}_cells"] = {}
    ss[f"{p}_ver"] = ss.get(f"{p}_ver", 0) + 1
    if "n" in ex:
        ss[f"{p}_n"] = ex["n"]
    if m == "analise":
        ss["an_dom"] = DOM_D if ex.get("domain") == "d" else DOM_C
    elif m in ("feedback", "observer"):
        ss[f"{p}_poles"] = ex["poles"]
        ss[f"{p}_dom"] = DOM_D if ex.get("domain") == "d" else DOM_C
    elif m == "discretize":
        ss["dz_T"] = str(ex["T"])
    elif m == "canonical":
        ss["cn_mode"] = MODE_TF if ex["mode"] == "tf" else MODE_SS
        if ex["mode"] == "tf":
            ss["cn_den"], ss["cn_num"] = ex["den"], ex["num"]
    elif m == "servo":
        ss["sv_poles"] = ex["poles"]
        ss["sv_dom"] = DOM_D if ex.get("servoDomain") == "d" else DOM_C
        ss["sv_conv"] = ex.get("servoConvention", CONV_UFRN)


def load_example(i: int):
    ex = EXEMPLOS[i]
    apply_example(ex)
    ss["pagina"] = PAGES[ex["module"]]


def boot():
    """Na primeira execução, cada módulo já abre com o primeiro exemplo correspondente."""
    if ss.get("_boot"):
        return
    seen = set()
    for ex in EXEMPLOS:
        if ex["module"] not in seen:
            seen.add(ex["module"])
            apply_example(ex)
    ss["_boot"] = True
    ss.setdefault("pagina", PAGES["analise"])


# ----------------------------------------------------------------------------
# Entrada de dados
# ----------------------------------------------------------------------------


class Inputs:
    """Coleta erros/avisos de leitura sem interromper a montagem dos widgets."""

    def __init__(self):
        self.errors: list[str] = []
        self.warnings: list[str] = []

    @property
    def ok(self):
        return not self.errors


def _default(name: str, rows: int, cols: int):
    M = np.zeros((rows, cols))
    if name == "A":
        np.fill_diagonal(M, 1)
    elif name == "B":
        M[0, 0] = 1
    elif name == "C":
        M[0, 0] = 1
    return M


def _initial_cells(p: str, name: str, rows: int, cols: int):
    """Células (strings) para montar o editor: edições salvas > exemplo carregado > padrão."""
    saved = ss.get(f"{p}_cells", {}).get(name)
    if saved is not None and np.shape(saved) == (rows, cols):
        return [list(r) for r in saved]
    init = ss.get(f"{p}_init", {}).get(name)
    data = np.array(init, float) if (init is not None and np.shape(init) == (rows, cols)) else _default(name, rows, cols)
    return [[("%g" % v) for v in row] for row in data]


def matrix_input(ctx: Inputs, p: str, name: str, rows: int, cols: int, tex: str):
    ver = ss.get(f"{p}_ver", 0)
    wkey = f"{p}_{name}_{rows}x{cols}_{ver}"
    dfkey = wkey + "__df"
    # O Streamlit descarta o estado de widgets não renderizados (ao trocar de módulo). Por isso o
    # DataFrame inicial é guardado à parte e só é refeito quando o editor é (re)montado.
    if wkey not in ss or dfkey not in ss:
        ss[dfkey] = _initial_cells(p, name, rows, cols)
    df = pd.DataFrame(ss[dfkey], columns=[str(j + 1) for j in range(cols)],
                      index=[str(i + 1) for i in range(rows)])
    st.markdown(f"${tex}$")
    edited = st.data_editor(
        df, key=wkey, num_rows="fixed",
        column_config={c: st.column_config.TextColumn(c) for c in df.columns},
    )
    M = np.zeros((rows, cols))
    empty = False
    cells = []
    for i in range(rows):
        row_cells = []
        for j in range(cols):
            raw = edited.iat[i, j]
            row_cells.append("" if raw is None else str(raw))
            if raw is None or str(raw).strip() == "":
                empty = True
                continue
            try:
                M[i, j] = safe_eval(raw)
            except ValueError:
                ctx.errors.append(f"Célula inválida em ${tex}$ (linha {i + 1}, coluna {j + 1}): `{raw}`.")
        cells.append(row_cells)
    ss.setdefault(f"{p}_cells", {})[name] = cells
    if empty:
        ctx.warnings.append(f"Há células vazias em ${tex}$; foram tratadas como 0.")
    return M


def order_input(p: str, default: int, label="Ordem do sistema $n$"):
    ss.setdefault(f"{p}_n", default)
    return int(st.number_input(label, min_value=1, max_value=6, step=1, key=f"{p}_n"))


def poles_input(ctx: Inputs, key: str, expected: int, label: str):
    ss.setdefault(key, "")
    raw = st.text_input(label, key=key,
                        help="Separe por vírgula. Complexos: `-2+1i`. Atalho: `-2±1i` gera o par conjugado.")
    poles = parse_complex_list(raw)
    if len(poles) != expected:
        ctx.errors.append(f"Foram reconhecidos {len(poles)} polo(s), mas são necessários exatamente "
                          f"**{expected}**. Exemplo: `-0.5+0.5i, -0.5-0.5i, 0.5`.")
    else:
        un = unpaired_complex_poles(poles)
        if un:
            ctx.errors.append("Polo(s) complexo(s) sem o conjugado: "
                              + ", ".join(f"${cnum(u)}$" for u in un)
                              + ". Para que o ganho seja real, todo polo complexo precisa do seu conjugado.")
    return poles


def domain_radio(key: str, label="Domínio"):
    ss.setdefault(key, DOM_C)
    return st.radio(label, [DOM_C, DOM_D], horizontal=True, key=key) == DOM_D


# ----------------------------------------------------------------------------
# Exibição do relatório
# ----------------------------------------------------------------------------


def _render_step(step, idx=None):
    with st.container(border=True):
        if step.final:
            st.markdown(f"#### ✅ {step.title}")
        else:
            st.markdown(f"**{idx}. {step.title}**")
        for it in step.items:
            if it[0] == "text":
                st.markdown(it[1])
            elif it[0] == "tex":
                st.latex(it[1])
            elif it[0] == "status":
                (st.success if it[1] else st.error)(it[2])
        if step.why:
            st.caption("💡 " + step.why)


def show_report(rep: Report, filename: str):
    for w in rep.warnings:
        st.warning(w)
    finals = [s for s in rep.steps if s.final]
    if finals:
        st.markdown("### Resultado")
        for s in finals:
            _render_step(s)
    st.markdown("### Passo a passo")
    idx = 0
    for s in rep.steps:
        if not s.final:
            idx += 1
            _render_step(s, idx)
    if rep.error:
        st.error("**Não foi possível concluir:** " + rep.error)
    c1, c2, _ = st.columns([1, 1, 3])
    c1.download_button("⬇️ Baixar LaTeX (.tex)", to_latex(rep), file_name=f"{filename}.tex",
                       mime="application/x-tex", key=f"dl_tex_{filename}")
    c2.download_button("⬇️ Baixar Markdown (.md)", to_markdown(rep), file_name=f"{filename}.md",
                       mime="text/markdown", key=f"dl_md_{filename}")


def solve(ctx: Inputs, fn, filename: str):
    for w in ctx.warnings:
        st.warning(w)
    if not ctx.ok:
        for e in ctx.errors:
            st.error(e)
        return
    try:
        rep = fn()
    except Exception as exc:  # noqa: BLE001 — mostra o erro ao aluno em vez de quebrar a página
        st.error(f"Erro no cálculo: {exc}")
        return
    st.divider()
    show_report(rep, filename)


# ----------------------------------------------------------------------------
# Páginas
# ----------------------------------------------------------------------------


def page_analise():
    st.header(PAGES["analise"])
    st.latex(r"\dot x = Ax + Bu,\quad y = Cx \qquad\text{ou}\qquad x(k+1) = Ax(k) + Bu(k),\quad y(k) = Cx(k)")
    ctx = Inputs()
    n = order_input("an", 3)
    discrete = domain_radio("an_dom")
    A = matrix_input(ctx, "an", "A", n, n, "A")
    c1, c2 = st.columns(2)
    with c1:
        B = matrix_input(ctx, "an", "B", n, 1, "B")
    with c2:
        C = matrix_input(ctx, "an", "C", 1, n, "C")
    solve(ctx, lambda: rep_analise(A, B, C, discrete), "analise")


def page_feedback():
    st.header(PAGES["feedback"])
    st.latex(r"K = -\begin{bmatrix}0 & \cdots & 0 & 1\end{bmatrix} W_c^{-1}\, q_c(A), \qquad u = Kx")
    ctx = Inputs()
    n = order_input("fb", 3)
    discrete = domain_radio("fb_dom", "Domínio (só altera a notação $s$/$z$ e o critério de estabilidade)")
    A = matrix_input(ctx, "fb", "A", n, n, "A")
    B = matrix_input(ctx, "fb", "B", n, 1, "B")
    desired = poles_input(ctx, "fb_poles", n, f"Polos desejados de malha fechada ({n})")
    solve(ctx, lambda: rep_feedback(A, B, np.array(desired), discrete), "realimentacao")


def page_observer():
    st.header(PAGES["observer"])
    st.latex(r"L = q_o(A)\, W_o^{-1} \begin{bmatrix}0 \\ \vdots \\ 1\end{bmatrix}, \qquad \dot{\hat x} = A\hat x + Bu + L(y - C\hat x)")
    ctx = Inputs()
    n = order_input("ob", 3)
    discrete = domain_radio("ob_dom", "Domínio (só altera a notação $s$/$z$)")
    A = matrix_input(ctx, "ob", "A", n, n, "A")
    C = matrix_input(ctx, "ob", "C", 1, n, "C")
    desired = poles_input(ctx, "ob_poles", n, f"Polos desejados do observador ({n})")
    solve(ctx, lambda: rep_observer(A, C, np.array(desired), discrete), "observador")


def page_discretize():
    st.header(PAGES["discretize"])
    st.latex(r"G = e^{AT}, \quad H = \int_0^T e^{A\tau}d\tau\, B \quad\Longleftarrow\quad "
             r"\exp\!\begin{bmatrix} AT & BT \\ 0 & 0\end{bmatrix} = \begin{bmatrix} G & H \\ 0 & I\end{bmatrix}")
    ctx = Inputs()
    n = order_input("dz", 2)
    ss.setdefault("dz_T", "1")
    T_raw = st.text_input("Período de amostragem $T$ (s)", key="dz_T")
    try:
        T = safe_eval(T_raw)
    except ValueError:
        T = None
        ctx.errors.append("O período de amostragem $T$ não é um número válido.")
    A = matrix_input(ctx, "dz", "A", n, n, "A")
    B = matrix_input(ctx, "dz", "B", n, 1, "B")
    solve(ctx, lambda: rep_discretize(A, B, T), "discretizacao")


def page_canonical():
    st.header(PAGES["canonical"])
    ss.setdefault("cn_mode", MODE_TF)
    mode = st.radio("Entrada", [MODE_TF, MODE_SS], horizontal=True, key="cn_mode")
    ctx = Inputs()
    if mode == MODE_TF:
        st.latex(r"G(s) = \dfrac{\beta_1 s^{n-1} + \cdots + \beta_n}{s^n + \alpha_1 s^{n-1} + \cdots + \alpha_n}")
        ss.setdefault("cn_den", "1, 3, 7, 5")
        ss.setdefault("cn_num", "3, 4, -2")
        c1, c2 = st.columns(2)
        den_raw = c1.text_input("Denominador (potências decrescentes de $s$)", key="cn_den")
        num_raw = c2.text_input("Numerador (potências decrescentes de $s$)", key="cn_num")
        den, numer = parse_number_list(den_raw), parse_number_list(num_raw)
        if not den:
            ctx.errors.append("Informe os coeficientes do denominador, ex.: `1, 3, 7, 5`.")
        solve(ctx, lambda: rep_canonical_tf(den, numer), "formas_canonicas")
    else:
        n = order_input("cn", 3)
        A = matrix_input(ctx, "cn", "A", n, n, "A")
        c1, c2 = st.columns(2)
        with c1:
            B = matrix_input(ctx, "cn", "B", n, 1, "B")
        with c2:
            C = matrix_input(ctx, "cn", "C", 1, n, "C")
        solve(ctx, lambda: rep_canonical_ss(A, B, C), "formas_canonicas")


def page_servo():
    st.header(PAGES["servo"])
    ctx = Inputs()
    ss.setdefault("sv_dom", DOM_C)
    c1, c2 = st.columns(2)
    discrete = c1.radio("Domínio", [DOM_C, DOM_D], horizontal=True, key="sv_dom") == DOM_D
    labels = ({CONV_UFRN: "UFRN (formulário oficial): método $\\hat G/\\hat H/M$",
               CONV_ALT: "Clássico (Ogata/Franklin): sistema aumentado $A_a/B_a$"} if discrete else
              {CONV_UFRN: "UFRN (formulário oficial): $K_a = +[0\\cdots1]W_a^{-1}q_c(A_a)$",
               CONV_ALT: "Alternativo (Ogata): $K_a = -[0\\cdots1]W_a^{-1}q_c(A_a)$"})
    ss.setdefault("sv_conv", CONV_UFRN)
    conv = c2.radio("Convenção", [CONV_UFRN, CONV_ALT], key="sv_conv", format_func=lambda k: labels[k])
    n = order_input("sv", 2)
    am, bm = ("G", "H") if discrete else ("A", "B")
    A = matrix_input(ctx, "sv", "A", n, n, am)
    d1, d2 = st.columns(2)
    with d1:
        B = matrix_input(ctx, "sv", "B", n, 1, bm)
    with d2:
        C = matrix_input(ctx, "sv", "C", 1, n, "C")
    desired = poles_input(ctx, "sv_poles", n + 1,
                          f"Polos desejados de malha fechada ({n + 1} = ordem do sistema aumentado)")
    solve(ctx, lambda: rep_servo(A, B, C, np.array(desired), discrete, conv), "seguidor")


def page_exemplos():
    st.header(PAGES["exemplos"])
    st.caption("Carregar um exemplo abre o módulo correspondente já preenchido.")
    groups: dict[str, list[int]] = {}
    for i, ex in enumerate(EXEMPLOS):
        groups.setdefault(ex["lista"], []).append(i)
    for lista, idxs in groups.items():
        st.subheader(lista)
        for i in idxs:
            ex = EXEMPLOS[i]
            with st.container(border=True):
                c1, c2 = st.columns([6, 1])
                c1.markdown(f"`{ex['tag']}` &nbsp; **{ex['title']}**  \n{ex['desc']}")
                c2.button("Carregar →", key=f"ex_{i}", on_click=load_example, args=(i,))


def sidebar_help():
    with st.sidebar.expander("Como digitar os valores"):
        st.markdown(
            "- **Matrizes:** clique na célula e digite. Aceita expressões: `0.5`, `0,5`, `1/3`, "
            "`-2*pi`, `sqrt(2)`, `1e-3`, `2^3`.\n"
            "- **Polos:** separados por vírgula. Complexos como `-2+1i` (ou `-2+1j`); "
            "o atalho `-2±1i` gera o par conjugado.\n"
            "- **Polos complexos** devem vir sempre com o conjugado.\n"
            "- **Exportar:** cada resolução pode ser baixada em `.tex` (compila em qualquer LaTeX) "
            "ou `.md`."
        )


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

boot()
# Mantém vivo o estado dos widgets simples de módulos que não estão na tela (ver nota em matrix_input).
for _k in (["an_dom", "dz_T", "cn_mode", "cn_den", "cn_num", "sv_dom", "sv_conv"]
           + [f"{p}_n" for p in PREFIX.values()] + [f"{p}_poles" for p in ("fb", "ob", "sv")]
           + [f"{p}_dom" for p in ("fb", "ob")]):
    if _k in ss:
        ss[_k] = ss[_k]
st.sidebar.title("🎛️ Bancada de Controle")
st.sidebar.caption("Projeto de Sistemas de Controle")
page = st.sidebar.radio("Módulo", list(PAGES.values()), key="pagina")
sidebar_help()

{
    PAGES["analise"]: page_analise,
    PAGES["feedback"]: page_feedback,
    PAGES["observer"]: page_observer,
    PAGES["discretize"]: page_discretize,
    PAGES["canonical"]: page_canonical,
    PAGES["servo"]: page_servo,
    PAGES["exemplos"]: page_exemplos,
}[page]()
