# 🎛️ Bancada de Controle — Streamlit

App para a disciplina **Projeto de Sistemas de Controle**: seis módulos que resolvem cada exercício
**passo a passo**, com todas as matrizes e fórmulas em LaTeX.

## Como rodar

```bash
pip install -r requirements.txt
streamlit run app.py
```

Para publicar de graça: suba esta pasta para um repositório no GitHub e crie o app no
[Streamlit Community Cloud](https://streamlit.io/cloud) apontando para `app.py`.

## Módulos

| # | Módulo | O que faz |
|---|--------|-----------|
| 1 | Análise | Faddeev–LeVerrier (cada iteração), polos, estabilidade (s ou z), $W_c$ e $W_o$ com det e posto |
| 2 | Realimentação de estados | Ackermann: $q_c(s)$, $q_c(A)$ detalhado, $W_c^{-1}$ (det, inversa, conferência), $K$, verificação de $A+BK$ |
| 3 | Observador | Ackermann dual: $q_o(A)$, $W_o^{-1}$, $L$, verificação de $A-LC$ |
| 4 | Discretização | $G$ e $H$ pela exponencial da matriz aumentada |
| 5 | Formas canônicas | A partir de $G(s)$ ou de $(A,B,C)$ (matriz $T$ de Ogata); forma geral simbólica + valores |
| 6 | Seguidor de referência | Contínuo/discreto, convenção UFRN ou Ogata/Franklin |

Cada resolução pode ser baixada em **`.tex`** (compila com `pdflatex`, sem pacotes exóticos) ou **`.md`**.
A aba **Exemplos** carrega as questões das Listas 1–3 já preenchidas.

## Entrada de dados

- Matrizes: células editáveis; aceitam `0.5`, `0,5`, `1/3`, `-2*pi`, `sqrt(2)`, `1e-3`, `2^3`.
- Polos: `-2+1i, -2-1i, -10`. O atalho `-2±1i` gera o par conjugado. Complexos precisam do conjugado.

## Convenções do seguidor (e = y − r, ξ = integrador do erro)

| Caso | Sistema aumentado | Ganho | Lei de controle | Malha fechada |
|------|-------------------|-------|-----------------|---------------|
| Contínuo, UFRN | $A_a=\begin{bmatrix}0&C\\0&A\end{bmatrix}$ | $K_a=+[0\cdots1]W_a^{-1}q_c(A_a)$ | $u=-k_1\!\int e-K_2x$ | $A_a-B_aK_a$ |
| Contínuo, Ogata | idem | $K_a=-[0\cdots1]W_a^{-1}q_c(A_a)$ | $u=k_1\!\int e+K_2x$ | $A_a+B_aK_a$ |
| Discreto, UFRN | $\hat G=\begin{bmatrix}G&H\\0&0\end{bmatrix}$, $\hat H=[0\cdots1]^T$ | $[K\,|\,-k_I]=(\hat K+e_n^T)M^{-1}$, $\hat K=+[0\cdots1]W_c^{-1}q_c(\hat G)$ | $u=k_1x+k_2\xi$, $k_1=-K$, $k_2=k_I$ | verificada na malha $(x,\xi)$ |
| Discreto, clássico | $A_a=\begin{bmatrix}1&C\\0&G\end{bmatrix}$ | $K_a=-[0\cdots1]W_a^{-1}q_c(A_a)$ | $u=k_1\xi+K_2x$ | $A_a+B_aK_a$ |

**Em todos os casos o app reconstrói a malha fechada real e confere que os polos coincidem com os desejados.**

## Arquivos

- `app.py` — interface Streamlit
- `nucleo.py` — álgebra linear, Ackermann, discretização, canônicas, servo, parsers (sem Streamlit)
- `modulos.py` — monta o passo a passo de cada módulo
- `latex.py` — formatação LaTeX, estrutura do relatório e exportação
- `exemplos.py` — exercícios das listas
- `tests/test_nucleo.py` — testes (`python -m pytest -q`)
