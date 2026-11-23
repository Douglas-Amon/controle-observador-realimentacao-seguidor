"""Formatação LaTeX + estrutura de relatório (passos) + exportação.

Tudo que o app mostra é gerado como LaTeX aqui; o mesmo relatório é exibido no
Streamlit (st.latex / st.markdown) e exportado como documento .tex ou .md.
"""
from __future__ import annotations

import html as _html
import re
import unicodedata
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


# ----------------------------------------------------------------------------
# Exportação HTML (com o nome do aluno)
# ----------------------------------------------------------------------------

KATEX_VER = "0.16.11"


def slug(nome: str) -> str:
    """Nome -> parte de nome de arquivo (sem acentos/espacos): 'João da Silva' -> 'Joao_da_Silva'."""
    t = unicodedata.normalize("NFKD", nome).encode("ascii", "ignore").decode()
    t = re.sub(r"[^A-Za-z0-9]+", "_", t).strip("_")
    return t or "aluno"


def md_to_html(s: str) -> str:
    """Markdown mínimo ($math$, **negrito**, `código`) -> HTML com delimitadores \\( \\) para o KaTeX."""
    out = []
    for i, p in enumerate(s.split("$")):
        if i % 2 == 1:
            out.append(r"\(" + _html.escape(p, quote=False) + r"\)")
        else:
            p = _html.escape(p, quote=False)
            p = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", p)
            p = re.sub(r"`(.+?)`", r"<code>\1</code>", p)
            out.append(p)
    return "".join(out)


_HTML_CSS = """
:root{--bg:#f4f6f8;--card:#fff;--ink:#1c2430;--mut:#5b6675;--line:#d9dee5;--acc:#0b6b6f;
--ok:#14703a;--okbg:#e6f4ea;--bad:#a12626;--badbg:#fdeaea;--warn:#8a5a00;--warnbg:#fff4de}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.55 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:900px;margin:0 auto;padding:24px 16px 56px}
header.top{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px 22px;margin-bottom:18px}
header.top .inst{font-size:.85rem;color:var(--mut);margin:0}
header.top h1{font-size:1.35rem;margin:.35rem 0 .8rem}
.meta{display:flex;flex-wrap:wrap;gap:6px 28px;font-size:.95rem}
.meta b{color:var(--acc)}
.step{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:14px 20px;margin:14px 0;break-inside:avoid}
.step h2{font-size:1.05rem;margin:.1rem 0 .6rem}
.step.final{border-color:var(--acc);border-width:2px}
.step.final h2{color:var(--acc)}
.math{overflow-x:auto;padding:2px 0}
.why{color:var(--mut);font-size:.9rem;font-style:italic;margin:.5rem 0 0}
.status{border-radius:8px;padding:8px 12px;font-weight:600;margin:8px 0}
.status.ok{background:var(--okbg);color:var(--ok)}
.status.bad{background:var(--badbg);color:var(--bad)}
.warn{background:var(--warnbg);color:var(--warn);border-radius:8px;padding:8px 12px;margin:10px 0}
.err{background:var(--badbg);color:var(--bad);border-radius:8px;padding:10px 14px;margin:14px 0;font-weight:600}
code{background:#eef1f4;padding:1px 5px;border-radius:4px}
#offline{background:var(--warnbg);color:var(--warn);border-radius:8px;padding:8px 12px;margin:10px 0}
button.print{margin:0 0 14px;padding:8px 14px;border:1px solid var(--acc);background:var(--acc);color:#fff;border-radius:8px;cursor:pointer;font:inherit}
footer{color:var(--mut);font-size:.8rem;margin-top:26px;text-align:center}
@media print{body{background:#fff}button.print,#offline{display:none!important}
.step,header.top{border-color:#999}main{padding:0}}
"""


def to_html(report: Report, aluno: str) -> str:
    """Documento HTML autônomo com o nome do aluno e a solução (matemática via KaTeX/CDN)."""
    nome = _html.escape(aluno.strip(), quote=False) or "&nbsp;" * 30
    titulo = _html.escape(report.title, quote=False)
    hoje = date.today().strftime("%d/%m/%Y")
    L = [
        "<!DOCTYPE html>",
        '<html lang="pt-BR"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{titulo} — {_html.escape(aluno.strip(), quote=False) or 'Resolução'}</title>",
        f'<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@{KATEX_VER}/dist/katex.min.css">',
        f"<style>{_HTML_CSS}</style></head><body><main>",
        '<button class="print" onclick="window.print()">🖨️ Imprimir / salvar em PDF</button>',
        '<div id="offline" hidden>Não foi possível carregar o renderizador de fórmulas (KaTeX). '
        "Conecte-se à internet e recarregue a página.</div>",
        '<header class="top">',
        '<p class="inst">Universidade Federal do Rio Grande do Norte · Centro de Tecnologia · '
        "Departamento de Engenharia de Computação e Automação</p>",
        f"<h1>{titulo}</h1>",
        '<div class="meta">',
        f"<span><b>Aluno(a):</b> {nome}</span>",
        f"<span><b>Data:</b> {hoje}</span>",
        "<span><b>Disciplina:</b> DCA-3701.0 — Projeto de Sistemas de Controle</span>",
        "</div></header>",
    ]
    for w in report.warnings:
        L.append(f'<div class="warn"><strong>Aviso:</strong> {md_to_html(w)}</div>')
    for idx, st in _numbered(report):
        head = md_to_html(st.title)
        L.append(f'<section class="step{" final" if st.final else ""}">')
        L.append(f"<h2>{'✅ ' + head if st.final else f'{idx}. {head}'}</h2>")
        for it in st.items:
            if it[0] == "text":
                L.append(f"<p>{md_to_html(it[1])}</p>")
            elif it[0] == "tex":
                L.append('<div class="math">\\[' + _html.escape(it[1], quote=False) + "\\]</div>")
            elif it[0] == "status":
                L.append(f'<div class="status {"ok" if it[1] else "bad"}">{md_to_html(it[2])}</div>')
        if st.why:
            L.append(f'<p class="why">{md_to_html(st.why)}</p>')
        L.append("</section>")
    if report.error:
        L.append(f'<div class="err">Erro: {md_to_html(report.error)}</div>')
    L += [
        f"<footer>Resolução gerada para {nome} em {hoje} — Bancada de Controle</footer>",
        f'<script defer src="https://cdn.jsdelivr.net/npm/katex@{KATEX_VER}/dist/katex.min.js"></script>',
        f'<script defer src="https://cdn.jsdelivr.net/npm/katex@{KATEX_VER}/dist/contrib/auto-render.min.js" '
        'onload="renderMathInElement(document.body,{delimiters:['
        "{left:'\\\\[',right:'\\\\]',display:true},{left:'\\\\(',right:'\\\\)',display:false}],"
        'throwOnError:false})"></script>',
        "<script>window.addEventListener('load',function(){if(!window.katex)"
        "document.getElementById('offline').hidden=false})</script>",
        "</main></body></html>",
    ]
    return "\n".join(L) + "\n"
