
# 🎛️ Bancada de Controle — Streamlit

Aplicativo desenvolvido para a disciplina **Projeto de Sistemas de Controle**, com seis módulos de resolução de exercícios de sistemas de controle.

O aplicativo apresenta as resoluções **passo a passo**, incluindo cálculos intermediários, matrizes, equações matemáticas em LaTeX e verificações dos resultados.

## 🚀 Acessar o aplicativo

**[Abrir a Bancada de Controle no Streamlit](https://SEU-APP.streamlit.app)**

> Substitua `https://SEU-APP.streamlit.app` pelo endereço público do aplicativo.

## 📚 Módulos disponíveis

| Nº | Módulo | Funcionalidades |
|---|---|---|
| 1 | Análise de sistemas | Faddeev–LeVerrier, polos, estabilidade, controlabilidade e observabilidade |
| 2 | Realimentação de estados | Fórmula de Ackermann, cálculo de ganhos e verificação da malha fechada |
| 3 | Observador de estados | Ackermann dual, ganho do observador e verificação dos polos |
| 4 | Discretização | Matrizes discretas pela exponencial da matriz aumentada |
| 5 | Formas canônicas | Representações canônicas e transformações de similaridade |
| 6 | Seguidor de referência | Projeto contínuo e discreto pelas convenções UFRN e Ogata/Franklin |

## 🔬 1. Análise de sistemas

O módulo realiza:

- Cálculo do polinômio característico pelo método de Faddeev–LeVerrier, mostrando cada iteração.
- Determinação dos autovalores (polos).
- Análise de estabilidade no plano contínuo (s) ou discreto (z).
- Construção da matriz de controlabilidade.
- Construção da matriz de observabilidade.
- Cálculo do determinante e do posto das matrizes.

Matriz de controlabilidade:

$$
W_c =
\begin{bmatrix}
B & AB & A^2B & \cdots & A^{n-1}B
\end{bmatrix}
$$

Matriz de observabilidade:

$$
W_o =
\begin{bmatrix}
C \\
CA \\
CA^2 \\
\vdots \\
CA^{n-1}
\end{bmatrix}
$$

Para sistemas SISO, a controlabilidade e a observabilidade completas são verificadas pelas condições:

$$
\operatorname{posto}(W_c)=n
$$

$$
\operatorname{posto}(W_o)=n
$$

## 🎯 2. Realimentação de estados

O módulo calcula o ganho de realimentação utilizando a fórmula de Ackermann.

São apresentados:

- Polinômio característico desejado.
- Cálculo detalhado da matriz polinomial.
- Matriz de controlabilidade.
- Determinante e inversa da matriz de controlabilidade.
- Verificação da inversa.
- Cálculo do ganho de realimentação.
- Verificação dos polos da malha fechada.

Para a convenção:

$$
u=Kx
$$

A matriz de malha fechada é:

$$
A_{mf}=A+BK
$$

O ganho correspondente é:

$$
K=-e_n^T W_c^{-1}q_c(A)
$$

em que:

$$
e_n^T=
\begin{bmatrix}
0 & 0 & \cdots & 1
\end{bmatrix}
$$

e $q_c(A)$ é o polinômio característico desejado avaliado na matriz $A$.

## 👁️ 3. Observador de estados

O módulo realiza o projeto de observadores pelo princípio da dualidade.

São apresentados:

- Polinômio característico desejado do observador.
- Matriz de observabilidade.
- Inversa da matriz de observabilidade.
- Cálculo do ganho do observador.
- Verificação dos autovalores da dinâmica do erro.

A dinâmica do erro de estimação é:

$$
\dot{\tilde{x}}=(A-LC)\tilde{x}
$$

A matriz do observador é:

$$
A_o=A-LC
$$

Para um sistema SISO observável, o ganho pode ser calculado por:

$$
L=q_o(A)W_o^{-1}e_n
$$

em que:

$$
e_n=
\begin{bmatrix}
0\\
0\\
\vdots\\
1
\end{bmatrix}
$$

## ⏱️ 4. Discretização

O módulo converte modelos contínuos para discretos usando a exponencial da matriz aumentada.

Modelo contínuo:

$$
\dot{x}(t)=Ax(t)+Bu(t)
$$

Modelo discreto:

$$
x[k+1]=Gx[k]+Hu[k]
$$

As matrizes discretas são:

$$
G=e^{AT}
$$

$$
H=\int_0^T e^{A\tau}B\,d\tau
$$

O cálculo é realizado por:

$$
e^{
\begin{bmatrix}
A & B\\
0 & 0
\end{bmatrix}T
}
=
\begin{bmatrix}
G & H\\
0 & I
\end{bmatrix}
$$

em que $T$ é o período de amostragem.

## 🔄 5. Formas canônicas

O módulo permite obter representações em espaço de estados a partir de:

- Função de transferência.
- Matrizes de espaço de estados.
- Transformações de similaridade.
- Procedimentos de Ogata.

Para a transformação:

$$
x=Tz
$$

as matrizes transformadas são:

$$
\bar{A}=T^{-1}AT
$$

$$
\bar{B}=T^{-1}B
$$

$$
\bar{C}=CT
$$

O aplicativo apresenta as expressões simbólicas gerais e os valores numéricos calculados.

## 🎛️ 6. Seguidor de referência

O módulo realiza projetos de seguidores de referência para sistemas contínuos e discretos.

São consideradas diferentes convenções de sinal utilizadas na disciplina e em referências clássicas.

### 6.1. Seguidor contínuo — Convenção UFRN

Define-se o erro:

$$
e=y-r
$$

A matriz aumentada é:

$$
A_a=
\begin{bmatrix}
0 & C\\
0 & A
\end{bmatrix}
$$

Para a ordem de estados aumentados $(\xi,x)$:

$$
B_a=
\begin{bmatrix}
0\\
B
\end{bmatrix}
$$

O ganho é calculado por:

$$
K_a=
\begin{bmatrix}
0 & \cdots & 1
\end{bmatrix}
W_a^{-1}q_c(A_a)
$$

A lei de controle é:

$$
u=-k_1\int e\,dt-K_2x
$$

A matriz de malha fechada é:

$$
A_{mf}=A_a-B_aK_a
$$

### 6.2. Seguidor contínuo — Convenção Ogata

A matriz aumentada é:

$$
A_a=
\begin{bmatrix}
0 & C\\
0 & A
\end{bmatrix}
$$

O ganho é:

$$
K_a=-
\begin{bmatrix}
0 & \cdots & 1
\end{bmatrix}
W_a^{-1}q_c(A_a)
$$

A lei de controle é:

$$
u=k_1\int e\,dt+K_2x
$$

A matriz de malha fechada é:

$$
A_{mf}=A_a+B_aK_a
$$

### 6.3. Seguidor discreto — Convenção UFRN

A matriz aumentada é:

$$
\hat{G}=
\begin{bmatrix}
G & H\\
0 & 0
\end{bmatrix}
$$

O vetor de entrada aumentado é:

$$
\hat{H}=
\begin{bmatrix}
0\\
\vdots\\
0\\
1
\end{bmatrix}
$$

O ganho auxiliar é:

$$
\hat{K}=
\begin{bmatrix}
0 & \cdots & 1
\end{bmatrix}
W_c^{-1}q_c(\hat{G})
$$

A transformação dos ganhos é:

$$
\begin{bmatrix}
K & -k_I
\end{bmatrix}
=
(\hat{K}+e_n^T)M^{-1}
$$

A lei de controle é:

$$
u=k_1x+k_2\xi
$$

com:

$$
k_1=-K
$$

$$
k_2=k_I
$$

A verificação dos polos é realizada utilizando a matriz real da malha fechada nas variáveis $(x,\xi)$.

### 6.4. Seguidor discreto — Convenção clássica

A matriz aumentada é:

$$
A_a=
\begin{bmatrix}
1 & C\\
0 & G
\end{bmatrix}
$$

O ganho é:

$$
K_a=-
\begin{bmatrix}
0 & \cdots & 1
\end{bmatrix}
W_a^{-1}q_c(A_a)
$$

A lei de controle é:

$$
u=k_1\xi+K_2x
$$

A matriz de malha fechada é:

$$
A_{mf}=A_a+B_aK_a
$$

### Verificação dos resultados

Em todos os casos, o aplicativo reconstrói a matriz real de malha fechada e verifica se os polos obtidos coincidem com os polos desejados.

## 📝 Entrada de dados

### Matrizes

As matrizes são preenchidas em células editáveis.

São aceitas expressões como:

```text
0.5
0,5
1/3
-2*pi
sqrt(2)
1e-3
2^3
```

### Polos desejados

Exemplo de polos reais e complexos:

```text
-2+1i, -2-1i, -10
```

Também é possível utilizar:

```text
-2±1i
```

O aplicativo interpreta automaticamente o par conjugado.

Os polos complexos devem aparecer em pares conjugados para projetos com matrizes e ganhos reais.

## 📂 Exemplos resolvidos

A aba **Exemplos** disponibiliza questões das Listas 1, 2 e 3 da disciplina.

Os dados podem ser carregados diretamente nos módulos para acompanhar o desenvolvimento dos cálculos.

## 📄 Exportação dos resultados

Cada resolução pode ser exportada nos formatos:

- `.tex` — documento LaTeX.
- `.md` — documento Markdown.

Os arquivos LaTeX podem ser compilados utilizando:

```bash
pdflatex arquivo.tex
```

## 💻 Como executar localmente

### 1. Instalar as dependências

```bash
pip install -r requirements.txt
```

### 2. Executar o aplicativo

```bash
streamlit run app.py
```

### 3. Acessar no navegador

Normalmente, o Streamlit disponibiliza o aplicativo em:

```text
http://localhost:8501
```

## ☁️ Publicação no Streamlit Community Cloud

1. Envie os arquivos do projeto para um repositório no GitHub.
2. Acesse [Streamlit Community Cloud](https://streamlit.io/cloud).
3. Crie um novo aplicativo.
4. Selecione o repositório do projeto.
5. Informe `app.py` como arquivo principal.
6. Realize a publicação.

Os arquivos Python devem permanecer na mesma pasta para que as importações funcionem corretamente.

## 🗂️ Estrutura do projeto

```text
controle_app/
├── app.py
├── exemplos.py
├── latex.py
├── modulos.py
├── nucleo.py
├── requirements.txt
├── README.md
└── tests/
    └── test_nucleo.py
```

### Descrição dos arquivos

| Arquivo | Função |
|---|---|
| `app.py` | Interface gráfica Streamlit |
| `nucleo.py` | Álgebra linear, Ackermann, discretização, formas canônicas, seguidor e parsers |
| `modulos.py` | Montagem das resoluções passo a passo |
| `latex.py` | Formatação matemática e exportação |
| `exemplos.py` | Exercícios e exemplos das listas |
| `requirements.txt` | Dependências Python |
| `tests/test_nucleo.py` | Testes automatizados |

## 🧪 Testes

Para executar os testes do núcleo:

```bash
python -m pytest -q
```

## 🎓 Finalidade acadêmica

Este aplicativo foi desenvolvido como ferramenta de apoio ao estudo de **Projeto de Sistemas de Controle**, permitindo visualizar os procedimentos matemáticos envolvidos na análise e no projeto de sistemas de controle.

O objetivo é facilitar a compreensão dos métodos por meio de resoluções detalhadas, cálculos intermediários e verificações numéricas.
