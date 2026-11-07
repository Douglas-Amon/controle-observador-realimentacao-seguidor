"""Formatação LaTeX + estrutura de relatório (passos) + exportação.

Tudo que o app mostra é gerado como LaTeX aqui; o mesmo relatório é exibido no
Streamlit (st.latex / st.markdown) e exportado como documento .tex ou .md.
"""
from __future__ import annotations

import re
from datetime import date

import numpy as np

# ----------------------------------------------------------------------------
# Números, matrizes, polinômios
# ----------------------------------------------------------------------------


def num(x, d: int = 4) -> str:
    """Número em LaTeX: inteiros sem casas, até d decimais, notação científica nos extremos."""
    x = float(np.real(x))
    if not np.isfinite(x):
        return r"\text{NaN}" if np.isnan(x) else (r"\infty" if x > 0 else r"-\infty")
    if abs(x) < 1e-9:
        return "0"
    if abs(x) >= 1e7 or abs(x) < 1e-4:
        m, e = f"{x:.3e}".split("e")
        m = m.rstrip("0").rstrip(".")
        return rf"{m}\times 10^{{{int(e)}}}"
    r = round(x, d)
    if abs(r - round(r)) < 1e-9:
        return str(int(round(r)))
    return f"{r:.{d}f}".rstrip("0").rstrip(".")


def cnum(z, d: int = 4) -> str:
    z = complex(z)
    if abs(z.imag) < 1e-6:
        return num(z.real, d)
    sign = "+" if z.imag >= 0 else "-"
    im = num(abs(z.imag), d)
    im = "" if im == "1" else im
    return f"{num(z.real, d)} {sign} {im}\\mathrm{{i}}"


def mat(M, label: str | None = None, bracket: str = "b") -> str:
    """Matriz numérica -> \\begin{bmatrix}...\\end{bmatrix}, com rótulo opcional ('A', 'W_c^{-1}', ...)."""
    M = np.atleast_2d(np.asarray(M, float))
    rows = r" \\ ".join(" & ".join(num(v) for v in row) for row in M)
    body = rf"\begin{{{bracket}matrix}} {rows} \end{{{bracket}matrix}}"
    return f"{label} = {body}" if label else body


def smat(rows, label: str | None = None) -> str:
    """Matriz simbólica (lista de listas de strings LaTeX)."""
    body = r" \\ ".join(" & ".join(r) for r in rows)
    body = rf"\begin{{bmatrix}} {body} \end{{bmatrix}}"
    return f"{label} = {body}" if label else body


def poly(coeffs, var: str = "s", ident: str | None = None) -> str:
    """Polinômio em LaTeX. Se ident='I', o termo constante vira c·I (polinômio matricial)."""
    coeffs = [float(c) for c in coeffs]
    n = len(coeffs) - 1
    parts = []
    for i, c in enumerate(coeffs):
        p = n - i
        if abs(c) < 1e-9 and not (p == 0 and not parts):
            continue
        a = abs(c)
        if p == 0:
            if ident:
                t = ident if abs(a - 1) < 1e-9 else rf"{num(a)}\,{ident}"
            else:
                t = num(a)
        else:
            v = var if p == 1 else f"{var}^{{{p}}}"
            t = v if abs(a - 1) < 1e-9 else rf"{num(a)}\,{v}"
        parts.append((c < 0, t))
    if not parts:
        return "0"
    out = ("-" if parts[0][0] else "") + parts[0][1]
    for neg, t in parts[1:]:
        out += (" - " if neg else " + ") + t
    return out


def roots_list(roots) -> str:
    return r",\; ".join(cnum(r) for r in roots)


def frac(n, d) -> str:
    return rf"\dfrac{{{n}}}{{{d}}}"


# ----------------------------------------------------------------------------
# Relatório
# ----------------------------------------------------------------------------


class Step:
    def __init__(self, title: str, why: str | None = None, final: bool = False):
        self.title = title
        self.why = why
        self.final = final
        self.items: list[tuple] = []

    def text(self, s: str):
        self.items.append(("text", s))
        return self

    def tex(self, *ss: str):
        for s in ss:
            self.items.append(("tex", s))
        return self

    def mats(self, *pairs, max_cols: int = 9):
        """pairs = (label, matriz). Agrupa matrizes pequenas na mesma linha."""
        groups, cur, cur_cols = [], [], 0
        for label, M in pairs:
            M = np.atleast_2d(np.asarray(M, float))
            c = M.shape[1]
            if cur and cur_cols + c > max_cols:
                groups.append(cur)
                cur, cur_cols = [], 0
            cur.append(mat(M, label))
            cur_cols += c
        if cur:
            groups.append(cur)
        for g in groups:
            self.items.append(("tex", r" \qquad ".join(g)))
        return self

    def status(self, ok: bool, msg_ok: str, msg_bad: str):
        self.items.append(("status", bool(ok), msg_ok if ok else msg_bad))
        return self


class Report:
    def __init__(self, title: str):
        self.title = title
        self.steps: list[Step] = []
        self.warnings: list[str] = []
        self.error: str | None = None

    def step(self, title: str, why: str | None = None, final: bool = False) -> Step:
        s = Step(title, why, final)
        self.steps.append(s)
        return s

    def fail(self, msg: str):
        self.error = msg
        return self


# ----------------------------------------------------------------------------
# Exportação
# ----------------------------------------------------------------------------

_TEX_ESC = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "#": r"\#", "_": r"\_",
            "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\^{}"}


def _esc(s: str) -> str:
    return "".join(_TEX_ESC.get(ch, ch) for ch in s)


def md_to_tex(s: str) -> str:
    """Markdown mínimo ($math$, **negrito**, `código`) -> LaTeX."""
    parts = s.split("$")
    out = []
    for i, p in enumerate(parts):
        if i % 2 == 1:
            out.append(f"${p}$")
        else:
            p = _esc(p)
            p = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", p)
            p = re.sub(r"`(.+?)`", r"\\texttt{\1}", p)
            out.append(p)
    return "".join(out)


def _numbered(report: Report):
    idx = 0
    for st in report.steps:
        if st.final:
            yield None, st
        else:
            idx += 1
            yield idx, st


def to_latex(report: Report) -> str:
    L = [
        r"\documentclass[11pt,a4paper]{article}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage{amsmath,amssymb}",
        r"\usepackage[margin=2cm]{geometry}",
        r"\usepackage{xcolor}",
        r"\setlength{\parindent}{0pt}\setlength{\parskip}{6pt}",
        rf"\title{{{md_to_tex(report.title)}}}",
        r"\author{Projeto de Sistemas de Controle}",
        rf"\date{{{date.today().strftime('%d/%m/%Y')}}}",
        r"\begin{document}",
        r"\maketitle",
    ]
    for w in report.warnings:
        L.append(rf"\textcolor{{orange}}{{\textbf{{Aviso:}} {md_to_tex(w)}}}" + "\n")
    for idx, st in _numbered(report):
        head = md_to_tex(st.title)
        L.append(rf"\subsection*{{{idx}. {head}}}" if idx else rf"\subsection*{{\textcolor{{teal}}{{{head}}}}}")
        for it in st.items:
            if it[0] == "text":
                L.append(md_to_tex(it[1]) + "\n")
            elif it[0] == "tex":
                L.append(r"\[" + "\n" + it[1] + "\n" + r"\]")
            elif it[0] == "status":
                L.append(rf"\textbf{{{md_to_tex(it[2])}}}" + "\n")
        if st.why:
            L.append(rf"\emph{{{md_to_tex(st.why)}}}" + "\n")
    if report.error:
        L.append(rf"\textcolor{{red}}{{\textbf{{Erro:}} {md_to_tex(report.error)}}}")
    L.append(r"\end{document}")
    return "\n".join(L) + "\n"


def to_markdown(report: Report) -> str:
    L = [f"# {report.title}", ""]
    for w in report.warnings:
        L += [f"> ⚠️ {w}", ""]
    for idx, st in _numbered(report):
        L += [f"## {idx}. {st.title}" if idx else f"## ✅ {st.title}", ""]
        for it in st.items:
            if it[0] == "text":
                L += [it[1], ""]
            elif it[0] == "tex":
                L += ["$$", it[1], "$$", ""]
            elif it[0] == "status":
                L += [f"**{it[2]}**", ""]
        if st.why:
            L += [f"*{st.why}*", ""]
    if report.error:
        L += [f"> ❌ **Erro:** {report.error}", ""]
    return "\n".join(L)
