
# 🎛️ Bancada de Controle — Streamlit

Aplicativo desenvolvido para a disciplina **Projeto de Sistemas de Controle**, com seis módulos de resolução de exercícios de sistemas de controle.

O aplicativo apresenta as resoluções **passo a passo**, incluindo cálculos intermediários, matrizes, equações matemáticas e verificações dos resultados.

## 🚀 Acessar o aplicativo

**[Abrir a Bancada de Controle no Streamlit](https://SEU-APP.streamlit.app)**

> Substitua o endereço acima pelo link público do seu aplicativo.

## 📚 Módulos disponíveis

| Nº | Módulo | Funcionalidades |
|---|---|---|
| 1 | Análise de sistemas | Faddeev–LeVerrier, polos, estabilidade, controlabilidade e observabilidade |
| 2 | Realimentação de estados | Fórmula de Ackermann, cálculo de ganhos e verificação da malha fechada |
| 3 | Observador de estados | Ackermann dual, ganho do observador e verificação dos polos |
| 4 | Discretização | Matrizes discretas pela exponencial da matriz aumentada |
| 5 | Formas canônicas | Representações canônicas e transformações de similaridade |
| 6 | Seguidor de referência | Projeto contínuo e discreto pelas convenções UFRN e Ogata/Franklin |

---

## 🔬 1. Análise de sistemas

O módulo realiza:

- Cálculo do polinômio característico pelo método de Faddeev–LeVerrier, mostrando cada iteração.
- Determinação dos autovalores (polos).
- Análise de estabilidade no plano contínuo (s) ou discreto (z).
- Construção da matriz de controlabilidade.
- Construção da matriz de observabilidade.
- Cálculo do determinante e do posto das matrizes.

### Matriz de controlabilidade

```text
Wc = [B  AB  A²B  ...  A⁽ⁿ⁻¹⁾B]
```

### Matriz de observabilidade

```text
     [ C       ]
Wo = [ CA      ]
     [ CA²     ]
     [ ...     ]
     [ CA⁽ⁿ⁻¹⁾ ]
```

### Condições de controlabilidade e observabilidade

Um sistema de ordem `n` é completamente controlável quando:

```text
posto(Wc) = n
```

É completamente observável quando:

```text
posto(Wo) = n
```

Para matrizes quadradas, também é possível verificar o determinante:

```text
det(Wc) != 0
det(Wo) != 0
```

---

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

### Convenção de controle

```text
u = Kx
```

### Matriz de malha fechada

```text
Amf = A + BK
```

### Fórmula de Ackermann

```text
K = -enᵀ * inv(Wc) * qc(A)
```

Em que:

```text
enᵀ = [0  0  ...  1]
```

- `Wc`: matriz de controlabilidade.
- `inv(Wc)`: inversa da matriz de controlabilidade.
- `qc(A)`: polinômio característico desejado avaliado na matriz A.
- `enᵀ`: vetor linha que seleciona a última linha da matriz resultante.

A verificação é realizada comparando os autovalores de `Amf` com os polos desejados.

---

## 👁️ 3. Observador de estados

O módulo realiza o projeto de observadores pelo princípio da dualidade.

São apresentados:

- Polinômio característico desejado do observador.
- Matriz de observabilidade.
- Inversa da matriz de observabilidade.
- Cálculo do ganho do observador.
- Verificação dos autovalores da dinâmica do erro.

### Dinâmica do erro de estimação

```text
d(x_til)/dt = (A - LC) * x_til
```

Em que `x_til` representa o erro de estimação.

### Matriz do observador

```text
Ao = A - LC
```

### Ganho do observador

Para um sistema SISO observável:

```text
L = qo(A) * inv(Wo) * en
```

Em que:

```text
     [0]
en = [0]
     [.]
     [1]
```

- `qo(A)`: polinômio característico desejado do observador avaliado em A.
- `Wo`: matriz de observabilidade.
- `en`: último vetor da base canônica.

A verificação é realizada pelos autovalores de:

```text
A - LC
```

---

## ⏱️ 4. Discretização

O módulo converte modelos contínuos para discretos usando a exponencial da matriz aumentada.

### Modelo contínuo

```text
dx(t)/dt = Ax(t) + Bu(t)
```

### Modelo discreto

```text
x[k+1] = Gx[k] + Hu[k]
```

### Matrizes discretas

```text
G = exp(A*T)
```

```text
H = integral de 0 até T de exp(A*tau)*B dtau
```

Em que:

- `T`: período de amostragem.
- `G`: matriz de estados discreta.
- `H`: matriz de entrada discreta.

### Exponencial da matriz aumentada

O cálculo é realizado por:

```text
         [ A  B ]
M = exp( [      ] * T )
         [ 0  0 ]

Resultado:

    [ G  H ]
M = [      ]
    [ 0  I ]
```

Esse procedimento permite calcular simultaneamente as matrizes G e H, inclusive quando A é singular.

---

## 🔄 5. Formas canônicas

O módulo permite obter representações em espaço de estados a partir de:

- Função de transferência.
- Matrizes de espaço de estados.
- Transformações de similaridade.
- Procedimentos de Ogata.

### Transformação de estados

```text
x = Tz
```

### Matrizes transformadas

```text
A_barra = inv(T) * A * T
```

```text
B_barra = inv(T) * B
```

```text
C_barra = C * T
```

Em que:

- `T`: matriz de transformação não singular.
- `inv(T)`: inversa da matriz de transformação.
- `A_barra`: matriz de estados transformada.
- `B_barra`: matriz de entrada transformada.
- `C_barra`: matriz de saída transformada.

O aplicativo apresenta as expressões simbólicas gerais e os valores numéricos calculados.

---

## 🎛️ 6. Seguidor de referência

O módulo realiza projetos de seguidores de referência para sistemas contínuos e discretos.

São consideradas diferentes convenções de sinal utilizadas na disciplina e em referências clássicas.

### 6.1. Seguidor contínuo — Convenção UFRN

Define-se o erro:

```text
e = y - r
```

### Matriz aumentada

Para a ordem de estados aumentados `(xi, x)`:

```text
     [ 0  C ]
Aa = [      ]
     [ 0  A ]
```

### Vetor de entrada aumentado

```text
     [ 0 ]
Ba = [   ]
     [ B ]
```

### Ganho aumentado

```text
Ka = [0  ...  1] * inv(Wa) * qc(Aa)
```

Em que:

- `Wa`: matriz de controlabilidade do sistema aumentado.
- `qc(Aa)`: polinômio desejado avaliado na matriz aumentada.

### Lei de controle

```text
u = -k1 * integral(e dt) - K2*x
```

### Matriz de malha fechada

```text
Amf = Aa - Ba*Ka
```

---

### 6.2. Seguidor contínuo — Convenção Ogata

### Matriz aumentada

```text
     [ 0  C ]
Aa = [      ]
     [ 0  A ]
```

### Ganho aumentado

```text
Ka = -[0  ...  1] * inv(Wa) * qc(Aa)
```

### Lei de controle

```text
u = k1 * integral(e dt) + K2*x
```

### Matriz de malha fechada

```text
Amf = Aa + Ba*Ka
```

A principal diferença em relação à convenção UFRN está no sinal adotado para o ganho e na lei de controle.

---

### 6.3. Seguidor discreto — Convenção UFRN

### Matriz aumentada

```text
       [ G  H ]
G_hat = [      ]
       [ 0  0 ]
```

### Vetor de entrada aumentado

```text
        [ 0 ]
        [ . ]
H_hat = [ . ]
        [ 0 ]
        [ 1 ]
```

### Ganho auxiliar

```text
K_hat = [0  ...  1] * inv(Wc) * qc(G_hat)
```

### Transformação dos ganhos

```text
[K  -kI] = (K_hat + enᵀ) * inv(M)
```

Em que:

- `M`: matriz de transformação utilizada no procedimento.
- `K_hat`: ganho auxiliar.
- `kI`: ganho integral.
- `enᵀ`: vetor linha da base canônica.

### Lei de controle

```text
u = k1*x + k2*xi
```

Com:

```text
k1 = -K
k2 = kI
```

A verificação dos polos é realizada utilizando a matriz real da malha fechada nas variáveis `(x, xi)`.

---

### 6.4. Seguidor discreto — Convenção clássica

### Matriz aumentada

```text
     [ 1  C ]
Aa = [      ]
     [ 0  G ]
```

### Ganho aumentado

```text
Ka = -[0  ...  1] * inv(Wa) * qc(Aa)
```

### Lei de controle

```text
u = k1*xi + K2*x
```

### Matriz de malha fechada

```text
Amf = Aa + Ba*Ka
```

### Verificação dos resultados

Em todos os casos, o aplicativo reconstrói a matriz real de malha fechada e verifica se os polos obtidos coincidem com os polos desejados.

---

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

---

## 📂 Exemplos resolvidos

A aba **Exemplos** disponibiliza questões das Listas 1, 2 e 3 da disciplina.

Os dados podem ser carregados diretamente nos módulos para acompanhar o desenvolvimento dos cálculos.

---

## 📄 Exportação dos resultados

Cada resolução pode ser exportada nos formatos:

- `.tex` — documento LaTeX.
- `.md` — documento Markdown.

Os arquivos LaTeX podem ser compilados utilizando:

```bash
pdflatex arquivo.tex
```

---

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

---

## ☁️ Publicação no Streamlit Community Cloud

1. Envie os arquivos do projeto para um repositório no GitHub.
2. Acesse [Streamlit Community Cloud](https://streamlit.io/cloud).
3. Crie um novo aplicativo.
4. Selecione o repositório do projeto.
5. Informe `app.py` como arquivo principal.
6. Realize a publicação.

Os arquivos Python devem permanecer na mesma pasta para que as importações funcionem corretamente.

---

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

---

## 🧪 Testes

Para executar os testes do núcleo:

```bash
python -m pytest -q
```

---

## 🎓 Finalidade acadêmica

Este aplicativo foi desenvolvido como ferramenta de apoio ao estudo de **Projeto de Sistemas de Controle**, permitindo visualizar os procedimentos matemáticos envolvidos na análise e no projeto de sistemas de controle.

O objetivo é facilitar a compreensão dos métodos por meio de resoluções detalhadas, cálculos intermediários e verificações numéricas.
