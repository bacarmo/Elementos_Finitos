### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 2f9d63db-0ffe-43b1-a29f-b4c92305a472
using PlutoUI,PlutoTeachingTools, Plots, LaTeXStrings, GaussQuadrature, SparseArrays, DataFrames

# ╔═╡ 0cc5dd29-8119-4fe5-9047-b2ee62eb6893
# ╠═╡ skip_as_script = true
#=╠═╡
PlutoUI.TableOfContents(title="Índice",depth=5,include_definitions=false)
  ╠═╡ =#

# ╔═╡ 354a7b37-db60-405d-a93b-a6d02a561cc0
# ╠═╡ skip_as_script = true
#=╠═╡
ChooseDisplayMode()
  ╠═╡ =#

# ╔═╡ 34ae963b-4aed-4504-9cd2-b6bfd5237833
md"# Aspectos teóricos"

# ╔═╡ c2e3404b-1584-407b-8f37-11f0b72a204c
md"## Formulação forte"

# ╔═╡ 57338d26-7c9a-471b-9005-3a2695753bf8
md"""
Dada uma função ``f:\overline{\Omega}\rightarrow\mathbb{R}`` e constantes reais ``\alpha>0`` e ``\beta\geq 0``, determine ``u:\overline{\Omega}\rightarrow\mathbb{R}`` tal que

```math
\begin{align}
\left\{
\begin{aligned}
&-\alpha \Delta u(x) + \beta u(x) = f(x),\quad x\in\Omega,
\\[4pt]
&u(x) = 0, \quad x\in\Gamma,
\end{aligned}\right.
\end{align}
```

em que, ``\Omega`` é um subconjunto do ``\mathbb{R}^2``, ``\Gamma`` é a fronteira de ``\Omega`` e ``\overline{\Omega}=\Omega\cup\Gamma``. 
"""

# ╔═╡ ae6e9c26-3f0a-4d49-be01-63ebe7750630
Foldable("Notação",
md"""
``\displaystyle \Delta u(x) = \frac{\partial^2 u}{\partial x_1^2}(x) + \frac{\partial^2 u}{\partial x_2^2}(x)``, sendo ``x=(x_1^{},x_2^{})`` um ponto do espaço ``\mathbb{R}^2``.
"""
)

# ╔═╡ cbd7811e-f1c5-4815-9fa5-d9e7f406734b
function exemplo1()
	α = 1.0
	β = 1.0
	f(x₁,x₂) = (2*α*π^2+β) * sin(π*x₁) * sin(π*x₂)
	u(x₁,x₂) = sin(π*x₁) * sin(π*x₂)

	return α, β, f, u
end

# ╔═╡ 41ed9d0d-819b-439a-8c3e-b58e09005924
md"""## Transição entre a formulação forte e fraca"""

# ╔═╡ db2d618a-6a27-4f19-93f9-522cf7fa9ad7
md"""
Considere ``V`` como o espaço formado por funções ``v:\overline{\Omega}\rightarrow\mathbb{R}``, com ``v`` suficientemente regular e satisfazendo a condição ``v(x)=0,\; \forall x\in\Gamma``.

Dado ``v \in V``, ao multiplicar a equação diferencial da formulação forte por ``v`` e integrar sobre ``\Omega``, obtemos

```math
-\alpha\int_\Omega u_{x_1^{}x_1^{}}^{}(x)v(x)d\Omega
-\alpha\int_\Omega u_{x_2^{}x_2^{}}^{}(x)v(x)d\Omega
+\beta\int_\Omega u(x)v(x)d\Omega
= 
  \int_\Omega f(x)v(x)d\Omega . 
```

Com isso, dado que
```math
- \int_\Omega u_{x_i^{}x_i^{}}^{}(x)v(x)d\Omega =
  \int_\Omega u_{x_i^{}}^{}(x)v_{x_i^{}}^{}(x)d\Omega,
```

temos que

```math
 \alpha\int_\Omega u_{x_1^{}}^{}(x)v_{x_1^{}}^{}(x)d\Omega
+\alpha\int_\Omega u_{x_2^{}}^{}(x)v_{x_2^{}}^{}(x)d\Omega
+\beta\int_\Omega u(x)v(x)d\Omega
= 
  \int_\Omega f(x)v(x)d\Omega. 
```
"""

# ╔═╡ 50781f95-ee8f-4d17-ace3-8eb44d8ae727
Foldable(
"Notação",
md"""
Veja que a expressão acima pode ser reescrita como

```math
 \alpha\int_\Omega \nabla u(x)\cdot\nabla v(x)d\Omega
+\beta\int_\Omega u(x)v(x)d\Omega
= 
  \int_\Omega f(x)v(x)d\Omega. 
```

Para simplificar ainda mais a escrita, vamos reescrever essa expressão em função de um operador bilinear ``\kappa`` e do produto interno de ``L^2(\Omega)``.
Considere

```math
\begin{align}
\kappa:V\times V\rightarrow\mathbb{R},\;(w,v)\mapsto \kappa(w,v)=\alpha\int_\Omega \nabla w(x)\cdot\nabla v(x)d\Omega + \beta  \int_\Omega w(x)v(x)d\Omega
\end{align}
```

e

```math
\begin{align*}
\big(w,v\big) = \int_\Omega w(x) v(x)d\Omega.
\end{align*}
```

Com isso, a equação integral pode ser reescrita como

```math
\kappa(u,v) = (f,v).
```
""")

# ╔═╡ c8177a78-1849-4a77-b81f-91d4849df06f
Foldable(
"Extra",
md"""
Demonstração da identidade
``\displaystyle - \int_\Omega u_{x_i^{}x_i^{}}^{}(x)v(x)d\Omega =  \int_\Omega u_{x_i^{}}^{}(x)v_{x_i^{}}^{}(x)d\Omega.``

---

Considere ``v\in V``. A derivada em relação a componente ``x_i^{}`` do produto ``u_{x_i^{}}(x)v(x)`` é dada por

```math
\frac{\partial}{\partial x_i^{}} \Big(u_{x_i^{}}^{}(x)v(x)\Big)
= 
u_{x_i^{}x_i^{}}^{}(x)v(x) + u_{x_i^{}}^{}(x)v_{x_i^{}}^{}(x). 
```

A partir dessa igualdade, tomando a integral em ``\Omega``, obtemos

```math
\int_\Omega\frac{\partial}{\partial x_i^{}} \Big(u_{x_i^{}}^{}(x)v(x)\Big)d\Omega
= 
  \int_\Omega u_{x_i^{} x_i^{}}^{}(x)v(x)d\Omega + \int_\Omega u_{x_i^{}}^{}(x)v_{x_i^{}}^{}(x)d\Omega . 
```

Para o lado esquerdo da equação acima, temos via o Teorema da Divergência que

```math
\int_\Omega
\frac{\partial}{\partial x_i^{}} \Big(u_{x_i^{}}^{}(x)v(x)\Big) d\Omega
= 
\int_\Gamma u_{x_i^{}}^{}(x)v(x) \eta_i^{} d\Gamma,  
```

onde ``\eta_i^{}`` é a ``i``-ésima componente do vetor normal unitário exterior ``\eta`` em ``x\in\Gamma``.
Com isso, fazendo uso dessa identidade,

```math
\int_\Gamma u_{x_i^{}}^{}(x)v(x) \eta_i^{} d\Gamma
= 
  \int_\Omega u_{x_i^{} x_i^{}}^{}(x)v(x)d\Omega + \int_\Omega u_{x_i^{}}^{}(x)v_{x_i^{}}^{}(x)d\Omega . 
```

Por fim, sendo ``v(x)=0,\,\forall x\in\Gamma``, a integral sobre a fronteira ``\Gamma`` é nula, resultando em 

```math
  \int_\Omega u_{x_i^{} x_i^{}}^{}(x)v(x)d\Omega + \int_\Omega u_{x_i^{}}^{}(x)v_{x_i^{}}^{}(x)d\Omega = 0. 
```
"""
)

# ╔═╡ bf596662-c956-4ed4-8711-89e7870a04db
md"## Formulação fraca"

# ╔═╡ 2f25b167-978a-4e5b-a22f-a63a0fc6e593
md"""
Dada uma função ``f:\overline{\Omega}\rightarrow\mathbb{R}`` e constantes reais ``\alpha>0`` e ``\beta\geq 0``, determine ``u\in V`` tal que

```math
\kappa(u,v) = (f,v),\quad \forall v\in V.
```
"""

# ╔═╡ 156276ab-1928-49c6-8cff-51ff7d987212
md"## Problema aproximado - Galerkin"

# ╔═╡ aedd68e1-c955-43dd-834b-dc8719713569
md"""
Dada uma função ``f:\overline{\Omega}\rightarrow\mathbb{R}`` e constantes reais ``\alpha>0`` e ``\beta\geq 0``, determine ``u_h^{}\in V_m`` tal que

```math
\kappa(u_h^{},v) = (f,v),\quad \forall v\in V_m,
```

sendo ``V_m`` um subespaço vetorial de ``V`` com dimensão finita ``m``.
"""

# ╔═╡ 203da506-72ef-4ca1-9c10-d615ebad7715
md"### Transição entre o problema aproximado e sua forma matriz-vetor"

# ╔═╡ 8e205232-794f-475a-912d-3cb822bc986f
md"""
Sendo ``V_m`` um subespaço vetorial de ``V`` com dimensão finita ``m``, ele possui uma base, que denotaremos por ``\varphi_1^{},\varphi_2^{},\dots,\varphi_m^{}``.
Como buscamos ``u_h^{}\in V_m``, o mesmo pode ser descrito como uma combinação linear dos elementos da base de ``V_m``, ou seja,

```math
u_h^{}(x) = \sum_{j=1}^{m} c_j^{}\varphi_j^{}(x).
```

Assim, determinar ``u_h^{}`` se resume a encontrar os coeficientes ``c_j^{}``, para ``j=1,2,\dots,m``.
Como há ``m`` incógnitas, são necessárias ``m`` equações, que podem ser obtidas ao tomar no problema aproximado ``v=\varphi_i^{}`` para ``i=1,2,\dots,m``.
Procedendo dessa forma, obtemos

```math
\begin{align*}
\left\{
\begin{aligned}
& \kappa\big(\sum_{j=1}^mc_j\varphi_j,\varphi_1 \big) = (f,\varphi_1), 
\\
& \kappa\big(\sum_{j=1}^mc_j\varphi_j,\varphi_2 \big) = (f,\varphi_2),
\\
& \vdots
\\
& \kappa\big(\sum_{j=1}^mc_j\varphi_j,\varphi_m \big) = (f,\varphi_m).
\end{aligned}\right.
\end{align*}
```

"""

# ╔═╡ 2c045f17-9dd7-4759-b19a-55933ee31728
md"""
Prosseguindo, dado que o operador $\kappa$ é linear em cada componente, segue que
```math
\begin{align*}
\left\{
\begin{aligned}
& \kappa\big(\varphi_1,\varphi_1\big)c_1 + \kappa\big(\varphi_2,\varphi_1\big)c_2 + \dots + \kappa\big(\varphi_m,\varphi_1\big)c_m = (f,\varphi_1), 
\\
& \kappa\big(\varphi_1,\varphi_2\big)c_1 + \kappa\big(\varphi_2,\varphi_2\big)c_2 + \dots + \kappa\big(\varphi_m,\varphi_2\big)c_m = (f,\varphi_2), 
\\
& \vdots
\\
& \kappa\big(\varphi_1,\varphi_m\big)c_1 + \kappa\big(\varphi_2,\varphi_m\big)c_2 + \dots + \kappa\big(\varphi_m,\varphi_m\big)c_m = (f,\varphi_m).
\end{aligned}\right.
\end{align*}
``` 
"""

# ╔═╡ a497d172-0374-4471-bb14-0ae8f081bf5a
md"""
Note que esse sistema pode ser reescrito na seguinte forma matricial:
```math
\begin{align*}
\begin{bmatrix}
\kappa\big(\varphi_1,\varphi_1\big)&
\kappa\big(\varphi_2,\varphi_1\big)&
\dots&
\kappa\big(\varphi_m,\varphi_1\big)
\\
\kappa\big(\varphi_1,\varphi_2\big)&
\kappa\big(\varphi_2,\varphi_2\big)&
\dots&
\kappa\big(\varphi_m,\varphi_2\big)
\\
\vdots&
\vdots&
\ddots&
\vdots  
\\
\kappa\big(\varphi_1,\varphi_m\big)&
\kappa\big(\varphi_2,\varphi_m\big)&
\dots&
\kappa\big(\varphi_m,\varphi_m\big)
\\
\end{bmatrix}
\begin{bmatrix}
c_1\\
c_2\\
\vdots\\
c_m
\end{bmatrix}
=
\begin{bmatrix}
(f,\varphi_1)\\
(f,\varphi_2)\\
\vdots\\
(f,\varphi_m)\\
\end{bmatrix}.
\end{align*}
```
"""

# ╔═╡ b9f73af1-9be5-4142-8292-bccd18af7978
md"### Formulação matricial"

# ╔═╡ 45fcb91c-0721-45f6-a2bc-514858791a34
md"""
Dada a matriz ``K`` e o vetor ``F``, determine o vetor ``c\in\mathbb{R}^m`` tal que 
```math
Kc=F,
```
onde 
```math
K_{i,j} = \kappa\big(\varphi_j,\varphi_i\big)
\quad\text{e}\quad
F_i = (f,\varphi_i), 
\quad\text{com}\quad
i,j\in\{1,2,\dots,m\}.
```
"""

# ╔═╡ e8a0bbc1-b3ca-49e0-9fd3-49cd883aa07f
md"# Aspectos numéricos"

# ╔═╡ 738c986b-e524-4cc1-a6f8-ff5650861a55
md"""
Nesta seção, consideraremos que cada elemento finito da malha, resultante da partição do domínio ``\Omega``, é um quadrilátero. 
Além disso, abordaremos apenas polinômios lineares por partes como base do subespaço aproximado ``V_m``.
"""

# ╔═╡ 93547f8d-e0f4-415f-abbe-800e52349746
md"## Definições de ``\phi_a``, ``\displaystyle\frac{\partial\phi_a}{\partial\xi_1}`` e ``\displaystyle\frac{\partial\phi_a}{\partial\xi_2}``"

# ╔═╡ 49c12023-3b0e-4246-9d28-33ebbbbbe0df
function ϕ(ξ₁,ξ₂,a)
	if a == 1
		return (1-ξ₁)*(1-ξ₂)/4
	elseif a == 2
		return (1+ξ₁)*(1-ξ₂)/4
	elseif a == 3
		return (1+ξ₁)*(1+ξ₂)/4
	elseif a == 4
		return (1-ξ₁)*(1+ξ₂)/4
	else
		error("a deve ser 1, 2, 3 ou 4.")
	end
end

# ╔═╡ 26f290b8-ce09-44d8-8d60-aad4fa04367a
function ∂ϕ_∂ξ₁(ξ₁,ξ₂,a)
	if a == 1
		return -(1-ξ₂)/4
	elseif a == 2
		return  (1-ξ₂)/4
	elseif a == 3
		return  (1+ξ₂)/4
	elseif a == 4
		return -(1+ξ₂)/4
	else
		error("a deve ser 1, 2, 3 ou 4.")
	end
end

# ╔═╡ 6f269a14-472a-4cd3-b5a3-99a61d9a3f23
function ∂ϕ_∂ξ₂(ξ₁,ξ₂,a)
	if a == 1
		return -(1-ξ₁)/4
	elseif a == 2
		return -(1+ξ₁)/4
	elseif a == 3
		return  (1+ξ₁)/4
	elseif a == 4
		return  (1-ξ₁)/4
	else
		error("a deve ser 1, 2, 3 ou 4.")
	end
end

# ╔═╡ 417b8cd2-bcda-4354-9e8d-d9c2492f8d7a
# ╠═╡ skip_as_script = true
#=╠═╡
function plot_base()
	# Partição uniforme de [-1, 1] no eixo ξ₁
	N₁ = 2^8
	ξ₁ = -1:2/N₁:1
	
	# Partição uniforme de [-1, 1] no eixo ξ₂
	N₂ = 2^8
	ξ₂ = -1:2/N₂:1
	
	# Gráficos das superfícies correspondentes às funções base ϕ₁, ϕ₂, ϕ₃ e ϕ₄
	plots = [plot(ξ₁, ξ₂, (ξ₁, ξ₂) -> ϕ(ξ₁, ξ₂, a), seriestype=:surface,
	              title=L"\phi_%$a", xlabel="ξ₁", ylabel="ξ₂", 
	              xticks=[-1, 0, 1], yticks=[-1, 0, 1], zticks=[0, 1],  # Definir ticks dos eixos
	              colorbar=false, color=:viridis)  # Personalizar cor da superfície (color=:diverging_bwr_55_98_c37_n256; color=:Purples_5)
	         for a in 1:4]  # Iterar para cada função base ϕ
	
	# Exibir os 4 gráficos lado a lado, organizados em uma única linha
	plot(plots[1],plots[2],plots[3],plots[4], layout=(1, 4), size=(1200, 300))
end
  ╠═╡ =#

# ╔═╡ e1df5449-2dc0-44da-883c-05cefaa286d2
# ╠═╡ skip_as_script = true
#=╠═╡
plot_base()
  ╠═╡ =#

# ╔═╡ 61fccc23-8111-42ec-9fe0-2e3d4ef993b3
md"## Caso 1: Domínios ``\Omega`` e ``\Omega^e`` com lados paralelos aos eixos cartesianos"

# ╔═╡ f6cfbb5a-933b-48da-a488-0e21972ca15d
md"""
Por simplicidade, consideremos ``\displaystyle \Omega = ]0,1[\times]0,1[``.
A partição de ``\Omega`` será realizada em função de dois parâmetros, 
```math
N_{x_1^{}}
\quad\text{e}\quad
N_{x_2^{}},
```
que representam, respectivamente, o número de intervalos nos eixos ``x_1^{}`` e ``x_2^{}``.

Para simplificar a notação, faremos uso de dois parâmetros adicionais, ``n_{x_1^{}}`` e ``n_{x_2^{}}``, definidos em função de ``N_{x_1^{}}`` e ``N_{x_2^{}}`` e do grau dos polinômios utilizados como base do espaço aproximado ``V_m``.
Dado que estamos considerando polinômios lineares por partes como base do espaço aproximado,
```math
n_{x_1^{}} = N_{x_1^{}} + 1
\quad\text{e}\quad
n_{x_2^{}} = N_{x_2^{}} + 1.
```

Na figura exibida a seguir, apresentamos a partição do domínio, com a numeração de cada elemento finito e a numeração de cada função base global.
"""

# ╔═╡ 8c7c8177-b8e8-4257-82b1-ba2125dc4d4a
html"""
<center>
<img src="https://raw.githubusercontent.com/bacarmo/Elementos_Finitos/main/imagens/malha2D.jpeg">
</center>
"""

# ╔═╡ 2450b17f-982b-4af8-bb68-7cca4f1944d6
md"### Mudança de variável"

# ╔═╡ 2f8e0a58-3476-4c04-aebe-4dec4f0f1db8
md"""
A seguir, definimos uma aplicação que mapeia cada ponto ``\xi\in\mathcal{R}=]-1,1[\times]-1,1[`` a um ponto ``x\in\Omega^e``.
Considere
```math
x(\xi) = \Big(
x_1(\xi),
x_2(\xi)
\Big),
\quad\hbox{onde}\quad
\begin{align}
\left\{
\begin{aligned}
x_1(\xi) = \frac{h_1}{2}(\xi_1+1) + p_1^e,
\\[5pt]
x_2(\xi) = \frac{h_2}{2}(\xi_2+1) + p_2^e.
\end{aligned}\right.
\end{align}
```

Observações:

``\bullet``  ``p^e = (p_1^e,p_2^e)`` são as coordenadas do ponto inferior esquerdo do elemento finito retangular ``\Omega^e``.

``\bullet``  ``h_1`` e ``h_2`` representam a base e altura do retângulo ``\Omega^e``.

``\bullet``  O diâmetro do retângulo ``\Omega^e`` é o comprimento da hipotenusa, dado por ``h=\sqrt{(h_1)^2+(h_1)^2}``.

``\bullet``  O determinante Jacobiano dessa aplicação é dado por 
```math
J
=
\frac{\partial(x_1,x_2)}{\partial(\xi_1,\xi_2)}
=
\begin{vmatrix} 
\displaystyle \frac{\partial x_1}{\partial \xi_1}(\xi) 
& 
\displaystyle \frac{\partial x_1}{\partial \xi_2}(\xi)
\\
\displaystyle \frac{\partial x_2}{\partial \xi_1}(\xi) 
&
\displaystyle \frac{\partial x_2}{\partial \xi_2}(\xi)
\end{vmatrix}
=
\begin{vmatrix} 
\displaystyle \frac{h_1}{2}
& 
\displaystyle 0
\\
\displaystyle 0
&
\displaystyle \frac{h_2}{2}
\end{vmatrix}
=
\frac{h_1 h_2}{4}.
```
"""

# ╔═╡ 0f85f3a1-beaa-4321-abb2-c671e75d0b15
md"""
``\bullet`` Temos que ``\varphi_a^e\big(x(\xi)\big)=\phi_a(\xi)``. Desta relação, precisamos descrever os termos 
```math
\frac{\partial\varphi_a^e}{\partial x_1}\big(x(\xi)\big)
\quad\hbox{e}\quad
\frac{\partial\varphi_a^e}{\partial x_2}\big(x(\xi)\big)
``` 
em função dos termos
```math
\frac{\partial\phi_a}{\partial \xi_1}(\xi)
\quad\hbox{e}\quad
\frac{\partial\phi_a}{\partial \xi_2}(\xi).
``` 
"""

# ╔═╡ e0da177e-18da-4c57-bca4-bd432ee19a68
md"""
Nesse sentido, derivando a identidade  ``\varphi_a^e\big(x(\xi)\big)=\phi_a(\xi)`` em relação ``\xi_1`` e em relação a ``\xi_2``, temos
```math
\begin{align}
(\bigstar)
\left\{
\begin{aligned}
\frac{\partial\varphi_a^e}{\partial x_1}\big(x(\xi)\big) 
\frac{\partial x_1}{\partial\xi_1}(\xi) + 
\frac{\partial\varphi_a^e}{\partial x_2}\big(x(\xi)\big) 
\frac{\partial x_2}{\partial\xi_1}(\xi)
=
\frac{\partial\phi_a}{\partial \xi_1}(\xi),
\\[10pt]
\frac{\partial\varphi_a^e}{\partial x_1}\big(x(\xi)\big) 
\frac{\partial x_1}{\partial\xi_2}(\xi) + 
\frac{\partial\varphi_a^e}{\partial x_2}\big(x(\xi)\big) 
\frac{\partial x_2}{\partial\xi_2}(\xi)
=
\frac{\partial\phi_a}{\partial \xi_2}(\xi).
\end{aligned}
\right.
\end{align}
```
Daí, dado que
```math
\frac{\partial x_1}{\partial\xi_1}(\xi) = \frac{h_1}{2},
\quad
\frac{\partial x_2}{\partial\xi_1}(\xi) = 0,
\quad
\frac{\partial x_1}{\partial\xi_2}(\xi) = 0
\quad\hbox{e}\quad
\frac{\partial x_2}{\partial\xi_2}(\xi) = \frac{h_2}{2},
```
obtemos
```math
\begin{align}
\left\{
\begin{aligned}
\frac{\partial\varphi_a^e}{\partial x_1}\big(x(\xi)\big) 
=\frac{2}{h_1}\frac{\partial\phi_a}{\partial \xi_1}(\xi),
\\[10pt]
\frac{\partial\varphi_a^e}{\partial x_2}\big(x(\xi)\big)
=\frac{2}{h_2}\frac{\partial\phi_a}{\partial \xi_2}(\xi).
\end{aligned}
\right.
\end{align}
```
"""

# ╔═╡ dd0daf35-6cdf-4884-8f0c-2d5e0387db78
Foldable("Spoiler",
md"""
O sistema em ``(\bigstar)`` pode ser reescrito como 
```math
\begin{bmatrix}
\displaystyle\frac{\partial x_1}{\partial\xi_1}(\xi)
&
\displaystyle\frac{\partial x_2}{\partial\xi_1}(\xi)
\\
\displaystyle\frac{\partial x_1}{\partial\xi_2}(\xi)
&
\displaystyle\frac{\partial x_2}{\partial\xi_2}(\xi)
\end{bmatrix}
\begin{bmatrix}
\displaystyle
\frac{\partial\varphi_a^e}{\partial x_1}\big(x(\xi)\big)
\\
\displaystyle
\frac{\partial\varphi_a^e}{\partial x_2}\big(x(\xi)\big)
\end{bmatrix}
=
\begin{bmatrix}
\displaystyle
\frac{\partial\phi_a}{\partial \xi_1}(\xi)
\\
\displaystyle
\frac{\partial\phi_a}{\partial \xi_2}(\xi)
\end{bmatrix}.
```
Note que o determinante da matriz ``2\times 2`` é igual a ``J``, o determinante Jacobiano da aplicação que realiza a mudança de variável de ``\xi`` para ``x``. 
Portanto, sempre que ``J`` for diferente de zero, podemos tomar a inversa, obtendo:
```math
\begin{bmatrix}
\displaystyle
\frac{\partial\varphi_a^e}{\partial x_1}\big(x(\xi)\big)
\\
\displaystyle
\frac{\partial\varphi_a^e}{\partial x_2}\big(x(\xi)\big)
\end{bmatrix}
=
\frac{1}{J}
\begin{bmatrix}
  \displaystyle\frac{\partial x_2}{\partial\xi_2}(\xi)
&
- \displaystyle\frac{\partial x_2}{\partial\xi_1}(\xi)
\\
- \displaystyle\frac{\partial x_1}{\partial\xi_2}(\xi)
&
  \displaystyle\frac{\partial x_1}{\partial\xi_1}(\xi)
\end{bmatrix}
\begin{bmatrix}
\displaystyle
\frac{\partial\phi_a}{\partial \xi_1}(\xi)
\\
\displaystyle
\frac{\partial\phi_a}{\partial \xi_2}(\xi)
\end{bmatrix}.
```
"""
)

# ╔═╡ a23031f6-073c-4d90-bbf3-280cd7d74bf5
@doc raw"""
	x₁_de_ξ(ξ₁::Float64, h₁::Float64, p₁::Vector{Float64}) -> Float64

Calcula a primeira componente da mudança de variável de um ponto ``\xi`` pertencente ao domínio de referência ``\mathcal{R}=]-1,1[ \times]-1,1[`` para o ponto correspondente `x` no elemento finito retangular ``\Omega^e``. 
```math
x_1(\xi) = \frac{h_1}{2}(\xi_1+1) + p_1.
```

# Parâmetros
- `ξ₁::Float64`: Primeira componente do ponto ``\xi\in\mathcal{R}=]-1,1[ \times]-1,1[``.
- `h₁::Float64`: Comprimento da base do elemento finito retangular ``\Omega^e``..
- `p₁::Float64`: Primeira componente do ponto inferior esquerdo do retângulo ``\Omega^e``..

# Retorno
- `x₁::Float64`: Primeira coordenada de ``x \in \Omega^e`` resultante da mudança de variável em `ξ`.
"""
function x₁_de_ξ(ξ₁::Float64, h₁::Float64, p₁::Float64)::Float64
    x₁ = p₁ + (h₁ / 2) * (ξ₁ + 1)
    return x₁
end

# ╔═╡ 3e07b842-9d3e-4424-b187-a1a4a2ae5ddf
@doc raw"""
	x₂_de_ξ(ξ₂::Float64, h₂::Float64, p₂::Vector{Float64}) -> Float64

Calcula a segunda componente da mudança de variável de um ponto ``\xi`` pertencente ao domínio de referência ``\mathcal{R}=]-1,1[ \times]-1,1[`` para o ponto correspondente `x` no elemento finito retangular ``\Omega^e``.
```math
x_2(\xi) = \frac{h_2}{2}(\xi_2+1) + p_2.
```

# Parâmetros
- `ξ₂::Float64`: Segunda componente do ponto ``\xi\in\mathcal{R}=]-1,1[ \times]-1,1[``.
- `h₂::Float64`: Comprimento da altura do elemento finito retangular ``\Omega^e``.
- `p₂::Float64`: Segunda componente do ponto inferior esquerdo do retângulo ``\Omega^e``.

# Retorno
- `x₂::Float64`: Segunda coordenada de ``x \in \Omega^e`` resultante da mudança de variável em `ξ`.
"""
function x₂_de_ξ(ξ₂::Float64, h₂::Float64, p₂::Float64)::Float64
    x₂ = p₂ + (h₂ / 2) * (ξ₂ + 1)
    return x₂
end

# ╔═╡ a952e33c-e0ea-4d37-9a5b-3e3aec41ab63
md"### Monta locais"

# ╔═╡ 9158f0d7-78cd-4a5e-a24c-d648273eec52
md"#### Cálculo do vetor local ``F^e``"

# ╔═╡ cdd4766a-b661-4bfb-9e2a-57708ce90438
md"""
```math
\begin{align}
F_a^e 
=& \int_{\Omega^e}f(x)\varphi_a^e(x)d\Omega
\\[5pt]
=& \int_{-1}^1\int_{-1}^1 f\big(x(\xi)\big)\varphi_a^e\big(x(\xi)\big)
\left|\frac{\partial(x_1,x_2)}{\partial(\xi_1,\xi_2)}\right| d\xi_1\,d\xi_2
\\[5pt]
=& \frac{h_1h_2}{4}\int_{-1}^1\int_{-1}^1 f\big(x(\xi)\big) \phi_a(\xi) d\xi_1\,d\xi_2.
\end{align}
```

"""

# ╔═╡ 3661a97c-6d06-41d0-b7df-71a8c1f70adc
@doc raw"""
    monta_Fᵉ!(Fᵉ::Vector{Float64}, f::Function, h₁::Float64, h₂::Float64, p₁::Float64, p₂::Float64, P::Vector{Float64}, W::Vector{Float64})

Na entrada `a` do vetor local `Fᵉ` é armazenado o valor aproximado da expressão 
```math
\frac{h₁h₂}{4} \int_{-1}^1\int_{-1}^1 f\Big(x_1(ξ_1,ξ_2),x_2(ξ_1,ξ_2)\Big)\phi_a(ξ_1,ξ_2)dξ_1dξ_2
``` 
utilizando quadratura gaussiana.

# Parâmetros
- `Fᵉ::Vector{Float64}`: Vetor força local a ser modificado, de tamanho 4.
- `f::Function`: Função ``f(x_1,x_2)`` fornecida como dado de entrada da edp.
- `h₁::Float64`: Comprimento da base do elemento finito retangular `Ωᵉ`.
- `h₂::Float64`: Comprimento da altura do elemento finito retangular `Ωᵉ`.
- `p₁::Float64`: Primeira componente do ponto inferior esquerdo do retângulo `Ωᵉ`.
- `p₂::Float64`: Segunda componente do ponto inferior esquerdo do retângulo `Ωᵉ`.
- `P::Vector{Float64}`: Pontos de quadratura no intervalo padrão `[-1, 1]`.
- `W::Vector{Float64}`: Pesos de quadratura associados a `P`.
"""
function monta_Fᵉ!(Fᵉ::Vector{Float64}, f::Function, h₁::Float64, h₂::Float64, p₁::Float64, p₂::Float64, P::Vector{Float64}, W::Vector{Float64})
    # Zera as entradas do vetor local Fᵉ
    fill!(Fᵉ, 0.0)

	# Pré-calcula o fator constante da integral
	cst = h₁ * h₂ / 4

    # Loop sobre as entradas do vetor local Fᵉ
    for a in 1:4
        # Inicializa a contribuição da quadratura para a entrada a
        soma = 0.0
        
        # Integração via quadratura gaussiana dupla
        for i in 1:length(P)          # Loop sobre os pontos de quadratura no eixo ξ₁
			ξ₁ = P[i]
            x₁ = x₁_de_ξ(ξ₁, h₁, p₁)
            for j in 1:length(P)      # Loop sobre os pontos de quadratura no eixo ξ₂
				ξ₂ = P[j]
                x₂ = x₂_de_ξ(ξ₂, h₂, p₂)

                # Acumula a contribuição da quadratura
                soma += W[i] * W[j] * 
                    f(x₁, x₂) * ϕ(ξ₁, ξ₂, a) * cst
            end
        end
        
        # Armazena a contribuição no vetor Fᵉ
        Fᵉ[a] = soma
    end
end

# ╔═╡ 4e4a5e81-ee38-4aad-a699-c44cb9f93178
md"""
__Verificação da montagem do vetor local ``F^e``__

__Teste 1.__ Observando a expressão a ser calculada, note que, se considerarmos
```math
f(x_1,x_2) = \frac{4}{h_1h_2},\quad
h_1 = h_2 = 1/4\quad\hbox{e}\quad
p_1 = p_2 = 0.0,
```
devemos obter
```math
F_a^e = \int_{-1}^1\int_{-1}^1 \phi_a(ξ_1,ξ_2)dξ_1dξ_2 = 1,
\quad\hbox{para}\; a\in\{1,2,3,4\}.
```
__Teste 2.__ Se considerarmos
```math
f(x_1,x_2) = \frac{16}{(h_1h_2)^2} 9x_1 x_2,\quad
h_1 = h_2 = 1/4\quad\hbox{e}\quad
p_1 = p_2 = 0.0,
```
devemos obter
```math
F_a^e = 9\int_{-1}^1\int_{-1}^1 (\xi_1+1)(\xi_2+1)\phi_a(ξ_1,ξ_2)dξ_1dξ_2,
\quad\hbox{para}\; a\in\{1,2,3,4\}.
```
Fazendo as contas analiticamente, segue que
```math
F^e = [4,\, 8,\, 16, \, 8].
```
"""

# ╔═╡ 2e5bae39-727b-4401-b425-23e39f1d982f
# ╠═╡ skip_as_script = true
#=╠═╡
function teste_monta_Fᵉ()
	Fᵉ = zeros(4)
	
	P, W = legendre(5)
	
	h₁ = 1/4
	h₂ = 1/4
	
	p₁ = 0.0
	p₂ = 0.0
	
	monta_Fᵉ!(Fᵉ, (x₁,x₂) -> 4/(h₁*h₂), h₁, h₂, p₁, p₂, P, W)
    display("Fᵉ - Teste 1")
	display(Fᵉ)
	
	monta_Fᵉ!(Fᵉ, (x₁,x₂) -> (16*9*x₁*x₂)/((h₁*h₂)^2), h₁, h₂, p₁, p₂, P, W)
    display("Fᵉ - Teste 2")
	display(Fᵉ)
end
  ╠═╡ =#

# ╔═╡ 3e78c4c9-d906-40f1-aa73-221529730b02
# ╠═╡ skip_as_script = true
#=╠═╡
teste_monta_Fᵉ()
  ╠═╡ =#

# ╔═╡ 0c0a81e8-cdc3-41ed-a8ad-a06042fdf3bb
md"#### Cálculo da matriz local ``K^e``"

# ╔═╡ cf54c8f7-f28d-4870-a306-227e25365b63
md"""
```math
\begin{align}
K_{a,b}^e 
=&
 \alpha\int_{\Omega^e} 
\frac{\partial\varphi_b^e}{\partial x_1}(x)
\frac{\partial\varphi_a^e}{\partial x_1}(x)
d\Omega
+\alpha\int_{\Omega^e} 
\frac{\partial\varphi_b^e}{\partial x_2}(x)
\frac{\partial\varphi_a^e}{\partial x_2}(x)
d\Omega
+\beta \int_{\Omega^e} 
\varphi_b^e(x)
\varphi_a^e(x)
d\Omega
\\[10pt]
=&
\alpha
\frac{h_2}{h_1}
\int_{-1}^1\int_{-1}^1
\frac{\partial\phi_b}{\partial \xi_1}(\xi)
\frac{\partial\phi_a}{\partial \xi_1}(\xi)
d\xi_1\,d\xi_2
\\[10pt]
& +
\alpha
\frac{h_1}{h_2}
\int_{-1}^1\int_{-1}^1
\frac{\partial\phi_b}{\partial \xi_2}(\xi)
\frac{\partial\phi_a}{\partial \xi_2}(\xi)
d\xi_1\,d\xi_2
\\[10pt]
& +
\beta
\frac{h_1h_2}{4} 
\int_{-1}^1\int_{-1}^1
\phi_b(\xi)
\phi_a(\xi)
d\xi_1\,d\xi_2.
\end{align}
```
"""

# ╔═╡ 399a4af6-593d-4a91-b0f9-b661ad2fba74
Foldable("Passo a passo",
md"""
``` math
\begin{align}
K_{a,b}^e 
=&
 \alpha\int_{\Omega^e} 
\frac{\partial\varphi_b^e}{\partial x_1}(x)
\frac{\partial\varphi_a^e}{\partial x_1}(x)
d\Omega
+\alpha\int_{\Omega^e} 
\frac{\partial\varphi_b^e}{\partial x_2}(x)
\frac{\partial\varphi_a^e}{\partial x_2}(x)
d\Omega
+\beta \int_{\Omega^e} 
\varphi_b^e(x)
\varphi_a^e(x)
d\Omega

\\[20pt]
= &
\alpha\int_{-1}^1\int_{-1}^1
\frac{\partial\varphi_b^e}{\partial x_1}\big(x(\xi)\big)
\frac{\partial\varphi_a^e}{\partial x_1}\big(x(\xi)\big)
\left|\frac{\partial(x_1,x_2)}{\partial(\xi_1,\xi_2)}\right| d\xi_1\,d\xi_2
\\[10pt]
& +
\alpha\int_{-1}^1\int_{-1}^1
\frac{\partial\varphi_b^e}{\partial x_2}\big(x(\xi)\big)
\frac{\partial\varphi_a^e}{\partial x_2}\big(x(\xi)\big)
\left|\frac{\partial(x_1,x_2)}{\partial(\xi_1,\xi_2)}\right| d\xi_1\,d\xi_2
\\[10pt]
& +
\beta \int_{-1}^1\int_{-1}^1
\varphi_b^e\big(x(\xi)\big)
\varphi_a^e\big(x(\xi)\big)
\left|\frac{\partial(x_1,x_2)}{\partial(\xi_1,\xi_2)}\right| d\xi_1\,d\xi_2

\\[20pt]
=&
\alpha\int_{-1}^1\int_{-1}^1
\frac{\partial\phi_b}{\partial \xi_1}(\xi)\frac{2}{h_1}
\frac{\partial\phi_a}{\partial \xi_1}(\xi)\frac{2}{h_1}
\frac{h_1h_2}{4} d\xi_1\,d\xi_2
\\[10pt]
& +
\alpha\int_{-1}^1\int_{-1}^1
\frac{\partial\phi_b}{\partial \xi_2}(\xi)\frac{2}{h_2}
\frac{\partial\phi_a}{\partial \xi_2}(\xi)\frac{2}{h_2}
\frac{h_1h_2}{4} d\xi_1\,d\xi_2
\\[10pt]
& +
\beta \int_{-1}^1\int_{-1}^1
\phi_b(\xi)
\phi_a(\xi)
\frac{h_1h_2}{4} d\xi_1\,d\xi_2

\\[20pt]
=&
\alpha
\frac{h_2}{h_1}
\int_{-1}^1\int_{-1}^1
\frac{\partial\phi_b}{\partial \xi_1}(\xi)
\frac{\partial\phi_a}{\partial \xi_1}(\xi)
d\xi_1\,d\xi_2
\\[10pt]
& +
\alpha
\frac{h_1}{h_2}
\int_{-1}^1\int_{-1}^1
\frac{\partial\phi_b}{\partial \xi_2}(\xi)
\frac{\partial\phi_a}{\partial \xi_2}(\xi)
d\xi_1\,d\xi_2
\\[10pt]
& +
\beta
\frac{h_1h_2}{4} 
\int_{-1}^1\int_{-1}^1
\phi_b(\xi)
\phi_a(\xi)
d\xi_1\,d\xi_2.
\end{align}
```
"""
)

# ╔═╡ 90787447-7ddf-41dc-9b76-3261f4168e5e
@doc raw"""
    monta_Kᵉ(α::Float64, β::Float64, h₁::Float64, h₂::Float64, P::Vector{Float64}, W::Vector{Float64}) -> Matrix{Float64}

Na entrada `[a,b]` da matriz local `Kᵉ` é armazenado o valor da expressão 
```math
\begin{align}
\alpha\frac{h_2}{h_1}
\int_{\mathcal{R}}
\frac{\partial\phi_b}{\partial \xi_1}(\xi)
\frac{\partial\phi_a}{\partial \xi_1}(\xi)
d\xi
+
\alpha\frac{h_1}{h_2}
\int_{\mathcal{R}}
\frac{\partial\phi_b}{\partial \xi_2}(\xi)
\frac{\partial\phi_a}{\partial \xi_2}(\xi)
d\xi
+
\beta\frac{h_1h_2}{4}
\int_{\mathcal{R}}
\phi_b(\xi)
\phi_a(\xi)
d\xi
\end{align}
```
utilizando quadratura gaussiana.

# Parâmetros
- `α::Float64`: Constante fornecida como dado de entrada da edp.
- `β::Float64`: Constante fornecida como dado de entrada da edp.
- `h₁::Float64`: Comprimento da base do elemento finito retangular `Ωᵉ`.
- `h₂::Float64`: Comprimento da altura do elemento finito retangular `Ωᵉ`.
- `P::Vector{Float64}`: Pontos de quadratura no intervalo padrão `[-1, 1]`.
- `W::Vector{Float64}`: Pesos de quadratura associados a `P`.

# Retorno
- `Kᵉ::Matrix{Float64}`: Matriz local `Kᵉ` de tamanho ``4\times 4``.
"""
function monta_Kᵉ(α::Float64, β::Float64, h₁::Float64, h₂::Float64, P::Vector{Float64}, W::Vector{Float64}) ::Matrix{Float64}
    # Inicializa a matriz local Kᵉ
    Kᵉ = zeros(4, 4)

    # Pré-calcula os fatores constantes das integrais
    cst1 = α * (h₂ / h₁)
    cst2 = α * (h₁ / h₂)
    cst3 = β * (h₁ * h₂ / 4)

    # Loop sobre as colunas (b) e linhas (a) da matriz local Kᵉ
    for b in 1:4
        for a in 1:4
            # Inicializa a contribuição da quadratura para a entrada (a, b)
            soma = 0.0

            # Integração via quadratura gaussiana dupla
            for i in 1:length(P)     # Loop sobre os pontos de quadratura no eixo ξ₁
				ξ₁ = P[i]
                for j in 1:length(P) # Loop sobre os pontos de quadratura no eixo ξ₂
                    ξ₂ = P[j]

                    # Acumula a contribuição da quadratura
                    soma += W[i] * W[j] * (
                        cst1 * ∂ϕ_∂ξ₁(ξ₁, ξ₂, b) * ∂ϕ_∂ξ₁(ξ₁, ξ₂, a) + 
                        cst2 * ∂ϕ_∂ξ₂(ξ₁, ξ₂, b) * ∂ϕ_∂ξ₂(ξ₁, ξ₂, a) + 
                        cst3 *      ϕ(ξ₁, ξ₂, b) *      ϕ(ξ₁, ξ₂, a)    
                    )
                end
            end

            # Armazena a contribuição acumulada na entrada [a, b] da matriz Kᵉ
            Kᵉ[a, b] = soma
        end
    end

    return Kᵉ
end

# ╔═╡ 074adb78-415f-405e-a8cf-3f4988019819
md"""
__Verificação da montagem da matriz local ``K^e``__

```math
K^e = 
\frac{\alpha h_2}{6h_1}
\begin{bmatrix}
\phantom{-}2&          -2&          -1&\phantom{-}1\\
          -2&\phantom{-}2&\phantom{-}1&          -1\\
          -1&\phantom{-}1&\phantom{-}2&          -2\\
\phantom{-}1&          -1&          -2&\phantom{-}2
\end{bmatrix}
+
\frac{\alpha h_1}{6h_2}
\begin{bmatrix}
\phantom{-}2&\phantom{-}1&          -1&          -2\\
\phantom{-}1&\phantom{-}2&          -2&          -1\\
          -1&          -2&\phantom{-}2&\phantom{-}1\\
          -2&          -1&\phantom{-}1&\phantom{-}2
\end{bmatrix}
+
\frac{\beta h_1h_2}{9\cdot 4}
\begin{bmatrix}
4& 2& 1& 2\\
2& 4& 2& 1\\
1& 2& 4& 2\\
2& 1& 2& 4
\end{bmatrix}
```

__Teste 1.__ Se considerarmos
```math
\alpha = 6.0,\quad
\beta  = 0.0 \quad\hbox{e}\quad
h_1 = h_2 = 1/4,
```
devemos obter
```math
K^e
=
\begin{bmatrix}
\phantom{-}4&          -1&          -2&          -1\\
          -1&\phantom{-}4&          -1&          -2\\
          -2&          -1&\phantom{-}4&          -1\\
          -1&          -2&          -1&\phantom{-}4
\end{bmatrix}
```

__Teste 2.__ Se considerarmos
```math
\alpha = 0.0,\quad
\beta  = \frac{9\cdot 4}{h_1h_2} \quad\hbox{e}\quad
h_1 = h_2 = 1/4,
```
devemos obter
```math
K^e
=
\begin{bmatrix}
4& 2& 1& 2\\
2& 4& 2& 1\\
1& 2& 4& 2\\
2& 1& 2& 4
\end{bmatrix}
```
"""

# ╔═╡ 1e6a2be3-ea18-476b-ad15-b4a4584f4ff5
# ╠═╡ skip_as_script = true
#=╠═╡
function teste_monta_Kᵉ()
	h₁ = 1/4
	h₂ = 1/4

	P, W = legendre(2)

	# Teste 1
	α = 6.0
	β = 0.0
	Kᵉ = monta_Kᵉ(α, β, h₁, h₂, P, W)
	display("Kᵉ - Teste 1")
	display(Kᵉ)

	# Teste 2
	α = 0.0
	β = (9*4)/(h₁*h₂)
	Kᵉ = monta_Kᵉ(α, β, h₁, h₂, P, W)
	display("Kᵉ - Teste 2")
	display(Kᵉ)	
end
  ╠═╡ =#

# ╔═╡ 958e69fc-12ec-432e-b4e9-babdc9a86c33
# ╠═╡ skip_as_script = true
#=╠═╡
teste_monta_Kᵉ()
  ╠═╡ =#

# ╔═╡ 52992179-f61b-4a18-bd57-fa96e571acf6
md"### Monta LG e EQ"

# ╔═╡ 2e274ee4-76c2-48f6-b980-34bf8a086849
md"#### Construção da matriz LG"

# ╔═╡ 83e1a70c-9307-41fb-9e1e-0d49276ab03f
html"""
<center>
<img src="https://raw.githubusercontent.com/bacarmo/Elementos_Finitos/main/imagens/malha2D.jpeg">
</center>
"""

# ╔═╡ 98cb3a5e-9628-4e31-9869-bdb319393d60
md"""
```math
LG = 
\begin{bmatrix}
%%%%%%%%%%%%%%% Bloco 1
\begin{bmatrix}
\\
\hbox{Bloco }1\\
\\
\end{bmatrix}_{4\times N_{x_1}},
%%%%%%%%%%%%%%% Bloco 2
\begin{bmatrix}
\\
\hbox{Bloco }2\\
\\
\end{bmatrix}_{4\times N_{x_1}},
\dots,
%%%%%%%%%%%%%%% Bloco N_{x_2}
\begin{bmatrix}
\\
\hbox{Bloco }N_{x_2}\\
\\
\end{bmatrix}_{4\times N_{x_1}}
\end{bmatrix}
```
"""

# ╔═╡ 956a6fcb-ee3a-4bb2-aaf5-109236982854
md"""
```math
\begin{align}
%%%%%%%%%%%%%%%%%%%%%%%% Primeiro bloco, i.e., e=1,..., N_{x_1}
\hbox{Bloco }1 = &
\begin{bmatrix}
1         & 2         & \dots &  n_{x_1}-1 \\
2         & 3         & \dots &  n_{x_1}   \\
n_{x_1}+2 & n_{x_1}+3 & \dots & 2n_{x_1}   \\
n_{x_1}+1 & n_{x_1}+2 & \dots & 2n_{x_1}-1
\end{bmatrix}_{4\times N_{x_1}},
\\[10pt]
%%%%%%%%%%%%%%%%%%%%%%%% Segundo bloco, i.e., e=N_{x_1}+1,..., 2N_{x_1}
\hbox{Bloco }2 = &
\begin{bmatrix}
 n_{x_1}+1 &  n_{x_1}+2 & \dots & 2n_{x_1}-1 \\
 n_{x_1}+2 &  n_{x_1}+3 & \dots & 2n_{x_1}   \\
2n_{x_1}+2 & 2n_{x_1}+3 & \dots & 3n_{x_1}   \\
2n_{x_1}+1 & 2n_{x_1}+2 & \dots & 3n_{x_1}-1
\end{bmatrix}_{4\times N_{x_1}},
\\[10pt]
%%%%%%%%%%%%%%%%%%%%%%%% 
\vdots
\\[10pt]
%%%%%%%%%%%%%%%%%%%%%%%% Último bloco, i.e., e=N_{x_2}-1)N_{x_1}+1,...,N_{x_2}N_{x_1}
\hbox{Bloco }N_{x_2} = &
\begin{bmatrix}
(n_{x_2}-2)n_{x_1}+1 & (n_{x_2}-2)n_{x_1}+2 & \dots & (n_{x_2}-1)n_{x_1}-1 \\
(n_{x_2}-2)n_{x_1}+2 & (n_{x_2}-2)n_{x_1}+3 & \dots & (n_{x_2}-1)n_{x_1}   \\
(n_{x_2}-1)n_{x_1}+2 & (n_{x_2}-1)n_{x_1}+3 & \dots &  n_{x_2}   n_{x_1}   \\
(n_{x_2}-1)n_{x_1}+1 & (n_{x_2}-1)n_{x_1}+2 & \dots &  n_{x_2}   n_{x_1}-1 \\
\end{bmatrix}_{4\times N_{x_1}}.
\end{align}
```


"""

# ╔═╡ dcbca563-105a-48cf-95a5-bca36795eed7
@doc raw"""
    monta_LG(Nx1::Int64, Nx2::Int64) -> Matrix{Int64}

Matriz de conectividade local/global (LG). Relaciona a numeração local e global das funções ``\varphi``.

# Parâmetros
- `Nx1::Int64`: Número de subdivisões ao longo do eixo x₁.
- `Nx2::Int64`: Número de subdivisões ao longo do eixo x₂.

# Retorno
- `LG::Matrix{Int64}`: Matriz de tamanho 4 x (Nx1 * Nx2), onde LG(a,e) fornece a numeração global da função local ``\varphi_a^e``.

# Exemplos
```jldoctest
julia> monta_LG(4,3)
4×12 Matrix{Int64}
  1  2  3   4   6   7   8   9  11  12  13  14
  2  3  4   5   7   8   9  10  12  13  14  15
  7  8  9  10  12  13  14  15  17  18  19  20
  6  7  8   9  11  12  13  14  16  17  18  19
```
```jldoctest
julia> monta_LG(4,4)
4×16 Matrix{Int64}:
  1  2  3   4   6   7   8   9  11  12  13  14  16  17  18  19
  2  3  4   5   7   8   9  10  12  13  14  15  17  18  19  20
  7  8  9  10  12  13  14  15  17  18  19  20  22  23  24  25
  6  7  8   9  11  12  13  14  16  17  18  19  21  22  23  24
```
"""
function monta_LG(Nx1::Int64, Nx2::Int64)::Matrix{Int64}
    # Define o número de funções φ no eixo x₁ e x₂
    nx1 = Nx1 + 1
    nx2 = Nx2 + 1

    # M[:,j] contém a numeração da primeira linha do "Bloco j"
    M = (1:nx1-1) .+ (0:nx1:(nx2-2)*nx1)'

    # LG[1,:] contém a numeração global da primeira função local de cada elemento 
    linha1 = reshape(M, 1, :)

    # Constrói a matriz LG
    LG = vcat(linha1, linha1 .+ 1, linha1 .+ (nx1+1), linha1 .+ nx1)

    return LG
end

# ╔═╡ 973ed431-0e2d-43f5-a407-fe73f868ab40
md"#### Construção do vetor EQ"

# ╔═╡ dc66e8b9-2e82-4e3b-bca3-9e316510d5a3
"""
    monta_EQ(Nx1::Int64, Nx2::Int64) -> Tuple{Int64, Vector{Int64}}

Calcula o valor de `m`, que representa a dimensão do subespaço aproximado `Vₘ`, e gera o vetor `EQ`, que fornece a reenumeração das funções globais `φ`.

# Parâmetros
- `Nx1::Int64`: Número de subdivisões ao longo do eixo x₁.
- `Nx2::Int64`: Número de subdivisões ao longo do eixo x₂.

# Retorno
Uma tupla contendo:
- `m::Int64`: Número de funções globais `φ` que compõem a base do espaço `Vₘ`. Para a malha retangular considerada, tomando-se uma base linear e com condições de contorno prescritas, temos `m = (Nx1-1)*(Nx2-1)`.
- `EQ::Vector{Int64}`: Vetor de inteiros de tamanho `(Nx1 + 1) * (Nx2 + 1)` que mapeia cada função global `φ` para sua nova numeração. As funções globais que fazem parte da base recebem a numeração de `1` a `m`, enquanto as que não fazem parte recebem o valor `m + 1`.

# Exemplos
```jldoctest
julia> m, EQ = monta_EQ(4, 3)
(6, [7, 7, 7, 7, 7,
     7, 1, 2, 3, 7,
     7, 4, 5, 6, 7,
     7, 7, 7, 7, 7])
```
```jldoctest
julia> m, EQ = monta_EQ(4, 4)
(9, [10, 10, 10, 10, 10, 
     10,  1,  2,  3, 10, 
     10,  4,  5,  6, 10, 
     10,  7,  8,  9, 10, 
     10, 10, 10, 10, 10])
```
"""
function monta_EQ(Nx1::Int64, Nx2::Int64) :: Tuple{Int64,Vector{Int64}}
    # Define o número de funções φ no eixo x₁ e x₂
    nx1 = Nx1 + 1
    nx2 = Nx2 + 1
    
    # Calcula o número de funções globais φ que compõem a base do espaço Vₘ
    m = (nx1-2) * (nx2-2)

    # Inicializa o vetor EQ preenchido com m+1
    EQ = fill(m+1, nx1 * nx2)

    # Vetor contendo os índices das funções globais φ que compõem a base do espaço Vₘ
    L = reshape( (0:nx1-3) .+ (nx1+2:nx1:(nx2-2)*nx1+2)' , :,1)

    # Atribui os valores de 1 até m as funções globais φ que compõem a base do espaço Vₘ
    EQ[L] = 1:m

    return m, EQ
end

# ╔═╡ 73081276-d764-401a-8f92-cc6b09620613
md"### Monta globais"

# ╔═╡ 301719ff-80d1-44ec-bb02-f92795d14e6a
md"#### Monta o vetor global ``F``"

# ╔═╡ 38e32ff9-a9c2-4018-b69b-8442563aa133
"""
    monta_F(f::Function, Nx1::Int64, Nx2::Int64, m::Int64, EQoLG::Matrix{Int64}) -> Vector{Float64}

Constrói o vetor global `F` de tamanho `m`, agregando os vetores locais `Fᵉ` associados a cada elemento finito na malha.

# Parâmetros
- `f::Function`: Função ``f(x_1,x_2)`` fornecida como dado de entrada da EDP.
- `Nx1::Int64`: Número de subdivisões ao longo do eixo x₁.
- `Nx2::Int64`: Número de subdivisões ao longo do eixo x₂.
- `m::Int64`: Dimensão do espaço aproximado `Vₘ`.
- `EQoLG::Matrix{Int64}`: EQ[LG].

# Retorno
- `F::Vector{Float64}`: Vetor global de tamanho `m`.
"""
function monta_F(f::Function, Nx1::Int64, Nx2::Int64, m::Int64, EQoLG::Matrix{Int64}) :: Vector{Float64}
    # Comprimento da base (h₁) e altura (h₂) de cada elemento retangular Ωᵉ
    h₁ = 1 / Nx1
    h₂ = 1 / Nx2

    # Número total de elementos finitos na malha
    ne = Nx1 * Nx2

    # P: Pontos de quadratura de Gauss-Legendre (ordem 5)
    # W: Pesos de quadratura de Gauss-Legendre
    P, W = legendre(5)

    # Inicializa o vetor local Fᵉ
    Fᵉ = zeros(4)

    # Inicializa o vetor global F com tamanho (m+1)
    F = zeros(m+1)
    
    # Loop sobre os elementos Ωᵉ (percorrendo cada subdivisão ao longo de x₂ e x₁)
    for j = 1:Nx2
        # Segunda componente do ponto inferior esquerdo do retângulo `Ωᵉ`.
        p₂ = (j-1)*h₂ 
        for i = 1:Nx1
            # Primeira componente do ponto inferior esquerdo do retângulo `Ωᵉ`.
            p₁ = (i-1)*h₁
            # Numeração do elemento finito atual (e) na malha
            e = (j-1)*Nx1 + i 

            # Calcula o vetor local Fᵉ para o elemento e
            monta_Fᵉ!(Fᵉ, f, h₁, h₂, p₁, p₂, P, W) 

            # Acumula Fᵉ em F usando a matriz de conectividade EQoLG
            for a = 1:4
                F[EQoLG[a,e]] += Fᵉ[a]
            end
        end
    end

    # Retorna o vetor global F com tamanho `m`, excluindo a última entrada adicional
    return F[1:m]
end

# ╔═╡ 8fd613d4-ee9f-4bde-aeb2-ea7950f822f6
function teste_monta_F()
	# Teste 1
	Nx1 = 4; Nx2 = 3; h₁ = 1/Nx1; h₂ = 1/Nx2;
	
	m, EQ = monta_EQ(Nx1,Nx2)
	LG = monta_LG(Nx1,Nx2)
	EQoLG = EQ[LG]
	
	F = monta_F((x₁,x₂) -> 4.0/(h₁*h₂), Nx1, Nx2, m, EQoLG)
    display("F - Teste 1")
	display("Nx1 = 4; Nx2 = 3; h₁ = 1/Nx1; h₂ = 1/Nx2; f(x₁,x₂) = 4/(h₁*h₂);")
	display(F)
	
	# Teste 2
	Nx1 = 4; Nx2 = 4; h₁ = 1/Nx1; h₂ = 1/Nx2;
	
	m, EQ = monta_EQ(Nx1,Nx2)
	LG = monta_LG(Nx1,Nx2)
	EQoLG = EQ[LG]
	
	F = monta_F((x₁,x₂) -> (16*9*x₁*x₂)/((h₁*h₂)^2), Nx1, Nx2, m, EQoLG)
    display("F - Teste 2")
	display("Nx1 = Nx2 = 4; h₁ = 1/Nx1; h₂ = 1/Nx2; f(x₁,x₂) = (16*9*x₁*x₂)/((h₁*h₂)^2);")
	display(F)
end

# ╔═╡ 12d6be03-ef98-4174-b1cc-2743b1198cf8
teste_monta_F()

# ╔═╡ bfd24f64-a567-478c-a70b-addc2632e4bf
md"#### Monta a matriz global ``K``"

# ╔═╡ 14427308-9ef8-46e8-b53d-7c1deb5d4226
"""
    monta_K(α::Float64, β::Float64, Nx1::Int64, Nx2::Int64, m::Int64, EQoLG::Matrix{Int64}) -> SparseMatrixCSC{Float64, Int64}

Gera a matriz esparsa `K` de tamanho `m x m`, montada a partir das matrizes locais `Kᵉ`.

# Parâmetros
- `α::Float64`: Parâmetro fornecido como dado de entrada da EDP.
- `β::Float64`: Parâmetro fornecido como dado de entrada da EDP.
- `Nx1::Int64`: Número de subdivisões ao longo do eixo x₁.
- `Nx2::Int64`: Número de subdivisões ao longo do eixo x₂.
- `m::Int64`: Dimensão do espaço aproximado `Vₘ`.
- `EQoLG::Matrix{Int64}`: ≡ EQ[LG].

# Retorno
- `K::SparseMatrixCSC{Float64, Int64}`: Matriz esparsa `K` de tamanho `m x m`.
"""
function monta_K(α::Float64, β::Float64, Nx1::Int64, Nx2::Int64, m::Int64, EQoLG::Matrix{Int64}) :: SparseMatrixCSC{Float64, Int64}
    h₁ = 1 / Nx1 # Comprimento da base do elemento finito retangular Ωᵉ
    h₂ = 1 / Nx2 # Comprimento da altura do elemento finito retangular Ωᵉ

    # Número de elementos na malha
    ne = Nx1 * Nx2

    # Pontos e pesos de quadratura de Gauss-Legendre
    P, W = legendre(2)

    # Calcula a matriz local Kᵉ
    Kᵉ = monta_Kᵉ(α, β, h₁, h₂, P, W)

    # Inicializa a matriz esparsa K com tamanho (m+1) x (m+1)
    K = spzeros(m+1, m+1)

    # Loop sobre os elementos
    for e = 1:ne
        # Loop sobre as colunas (b) e linhas (a) da matriz local Kᵉ
        for b = 1:4
            j = EQoLG[b,e]
            for a = 1:4
                i = EQoLG[a,e]
                K[i,j] += Kᵉ[a,b]
            end
        end
    end
   
    # Remove a última linha e coluna da matriz K
    return K[1:m, 1:m]
end

# ╔═╡ 8c2529d9-02e4-4d2c-a03a-cad4298b6f16
function teste_monta_K()
	α = 1.0; β = 1.0; Nx1 = 4; Nx2 = 3
	
	m, EQ = monta_EQ(Nx1,Nx2)
	LG = monta_LG(Nx1,Nx2)
	EQoLG = EQ[LG]
	
	K = monta_K(α, β, Nx1, Nx2, m, EQoLG)

	display("Parâmetros de entrada: α = 1.0; β = 1.0; Nx1 = 4; Nx2 = 3")
	display(K)
end

# ╔═╡ 3c066c9c-f0ac-46ff-82ca-7d4b7d316f4c
teste_monta_K()

# ╔═╡ 6c1557e9-3192-4cf6-87d3-3a807e0b1169
md"### Cálculo do erro"

# ╔═╡ 8543aafc-33f5-4e9b-afb2-8df58d12de81
md"""
Sejam ``u(x)`` a solução exata e ``\displaystyle u_h(x)=\sum_{j=1}^m c_j\varphi_j(x)`` a solução aproximada da EDP em estudo.
Denotando por ``\bar{c}=[c;0]`` o vetor de coeficientes da solução aproximada estendido com um zero na entrada ``m+1``, o erro na norma do espaço ``L^2(\Omega)`` é dado por:
```math
\begin{align}
\|u-u_h\|^2
= &
\int_\Omega |u(x)-u_h(x)|^2dx
\\[5pt]
=&
\sum_{e=1}^{n_e} \int_{\Omega^e} |u(x)-u_h(x)|^2dx
\\[5pt]
=&
\sum_{e=1}^{n_e} \int_{-1}^{1} \int_{-1}^{1}
\big|u\big(x(\xi)\big)
-\sum_{a=1}^{4}\bar{c}_{EQoLG(a,e)}\phi_a(\xi)\big|^2|J|d\xi
\\[5pt]
=&
\frac{h_1h_2}{4}\sum_{e=1}^{n_e} \int_{-1}^{1} \int_{-1}^{1}
\big|u\big(x(\xi)\big)
-\sum_{a=1}^{4}\bar{c}_{EQoLG(a,e)}\phi_a(\xi)\big|^2d\xi.
\end{align}
```
"""

# ╔═╡ 29097174-9c4e-4a47-afad-66a07c088e70
@doc raw"""
    erro_norma_L2(u::Function, c̄::Vector{Float64}, Nx1::Int64, Nx2::Int64, EQoLG::Matrix{Int64}) -> Float64

Calcula o erro na norma ``L^2(\Omega)`` entre a solução exata `u` e a solução aproximada representada pelos coeficientes `c̄` em uma malha de elementos retangulares.

# Parâmetros
- `u::Function`: Função u(x₁,x₂) que representa a solução exata da EDP.
- `c̄::Vector{Float64}`: Vetor com os coeficientes da solução aproximada acrescido de um zero, i.e., `c̄ = [c; 0]`.
- `Nx1::Int64`: Número de subdivisões ao longo do eixo x₁.
- `Nx2::Int64`: Número de subdivisões ao longo do eixo x₂.
- `EQoLG::Matrix{Int}`: EQ[LG].

# Retorna
- `Float64`: O valor do erro na norma ``L^2(\Omega)`` entre a solução exata e aproximada, calculada pela integração em todos os elementos da malha utilizando quadratura de Gauss-Legendre.
"""
function erro_norma_L2(u::Function, c̄::Vector{Float64}, Nx1::Int64, Nx2::Int64, EQoLG::Matrix{Int64}) :: Float64
    # Comprimentos da base (h₁) e altura (h₂) de cada elemento retangular Ωᵉ
    h₁ = 1 / Nx1
    h₂ = 1 / Nx2

    # Inicializa o erro
    erro = 0.0

    # Define o número de pontos de quadratura (Npg) e os pontos (P) e pesos (W)
    Npg = 5
    P, W = legendre(Npg)

    # Loop sobre os elementos Ωᵉ (percorrendo cada subdivisão ao longo de x₂ e x₁)
    for j = 1:Nx2
        # Calcula a segunda componente do ponto inferior esquerdo do retângulo `Ωᵉ`.
        p₂ = (j - 1) * h₂ 
        for i = 1:Nx1
            # Calcula a primeira componente do ponto inferior esquerdo de `Ωᵉ`.
            p₁ = (i - 1) * h₁
            # Identifica o número do elemento finito atual (e) na malha
            e = (j - 1) * Nx1 + i

            # Obtém os coeficientes `c` da solução aproximada no elemento `e` 
            c1e = c̄[EQoLG[1, e]]
            c2e = c̄[EQoLG[2, e]]
            c3e = c̄[EQoLG[3, e]]
            c4e = c̄[EQoLG[4, e]]

            # Realiza a integração dupla usando quadratura de Gauss-Legendre
            for b = 1:Npg
                # Coordenada ξ₂ no espaço de referência e x₂ no espaço físico
                ξ₂ = P[b]
                x₂ = x₂_de_ξ(ξ₂, h₂, p₂)
                for a = 1:Npg
                    # Coordenada ξ₁ no espaço de referência e x₁ no espaço físico
                    ξ₁ = P[a]
                    x₁ = x₁_de_ξ(ξ₁, h₁, p₁)

                    # Acumula a contribuição da quadratura
                    erro += W[a] * W[b] * (u(x₁, x₂) - 
                        c1e * ϕ(ξ₁, ξ₂, 1) -
                        c2e * ϕ(ξ₁, ξ₂, 2) -
                        c3e * ϕ(ξ₁, ξ₂, 3) -
                        c4e * ϕ(ξ₁, ξ₂, 4))^2
                end
            end
        end
    end

    return sqrt(erro * h₁ * h₂ / 4)
end

# ╔═╡ fa49b560-e000-47cc-ba85-16c111d89acd
md"# Simulações numéricas"

# ╔═╡ 8bef9e7a-6def-4751-bae8-d53343e94a3f
md"## Solução aproximada vs solução exata"

# ╔═╡ de6cac6c-1f5c-4c44-9e1d-78c6712d3d45
function plot_solução_aproximada(c̄::Vector{Float64}, Nx1::Int64, Nx2::Int64, EQoLG::Matrix{Int64})
    # Comprimentos da base (h₁) e altura (h₂) de cada elemento retangular Ωᵉ
    h₁ = 1 / Nx1
    h₂ = 1 / Nx2

    # Define uma discretização do intervalo de referência [-1, 1] nos eixos ξ₁ e ξ₂
    P = collect(-1:0.1:1)

    # Inicializa um objeto de gráfico que acumulará as superfícies
    plt = plot(seriestype = :surface, title="Solução Aproximada")
	
    # Loop nos elementos da malha (percorrendo cada subdivisão ao longo de x₂ e x₁)
    for j = 1:Nx2
        # Define a segunda coordenada do ponto inferior esquerdo do retângulo `Ωᵉ`.
        p₂ = (j - 1) * h₂
        # Calcula os valores correspondentes no eixo x₂ para os pontos P no eixo ξ₂
        x₂ = x₂_de_ξ.(P, h₂, p₂)

        for i = 1:Nx1
            # Primeira coordenada do ponto inferior esquerdo do retângulo `Ωᵉ`.
            p₁ = (i - 1) * h₁
            # Valores correspondentes no eixo x₁ para os pontos P no eixo ξ₁
            x₁ = x₁_de_ξ.(P, h₁, p₁)

            # Determina a numeração do elemento finito atual (`e`) na malha
            e = (j - 1) * Nx1 + i

            # Obtém os coeficientes `c` da solução aproximada no elemento `e`
            c1e = c̄[EQoLG[1, e]]
            c2e = c̄[EQoLG[2, e]]
            c3e = c̄[EQoLG[3, e]]
            c4e = c̄[EQoLG[4, e]]

            # Gera o gráfico da solução aproximada sobre o elemento `e`
            plot!(plt, x₁, x₂, (x₁, x₂) -> 
				  c1e*ϕ(x₁,x₂,1) 
				+ c2e*ϕ(x₁,x₂,2) 
				+ c3e*ϕ(x₁, x₂, 3) 
				+ c4e * ϕ(x₁, x₂, 4),
                seriestype = :surface, colorbar=false, color=:viridis,
				xticks=[0, 0.5, 1],yticks=[0, 0.5, 1],zticks=[0, 1],
				zlims=(0,1)
            )
        end
    end

	return plt
end

# ╔═╡ ea18dc94-c885-4b67-a4bd-5ab43258b42c
function solução_aproximada_vs_exata()
    # Carrega os parâmetros de entrada da EDP
    α, β, f, u = exemplo1()
	
    # Define o número de subdivisões ao longo dos eixos x₁ e x₂
    Nx1 = 4
    Nx2 = 4

    # Exibe os parâmetros de entrada
    println("Parâmetros de entrada:")
    println("Exemplo 1 & Nx1 = $Nx1, Nx2 = $Nx2")
    
    # Define o comprimento da base (h₁) e altura (h₂) de cada elemento retangular Ωᵉ
    h₁ = 1/Nx1
    h₂ = 1/Nx2

    # Gera a matriz de conectividade local/global (LG)
    LG = monta_LG(Nx1, Nx2)

    # Gera o valor de `m` e o vetor `EQ`
    m, EQ = monta_EQ(Nx1, Nx2)

    # Define a matriz de conectividade EQoLG
    EQoLG = EQ[LG]

    # Monta a matriz esparsa `K`
    K = monta_K(α, β, Nx1, Nx2, m, EQoLG)

    # Monta o vetor global `F`
    F = monta_F(f, Nx1, Nx2, m, EQoLG)

    # Resolve o sistema linear Kc = F
    c = K \ F

    # Calcula a solução exata nos nós internos da malha
    c_exato = [u(x₁,x₂) for x₂ in h₂:h₂:1-h₂, x₁ in h₁:h₁:1-h₁]
    
    # Exibe a solução aproximada (vetor c) e a solução exata (vetor c_exato)
    println("Solução aproximada:")
    println(c)
    println("Solução exata:")
    println(c_exato)

	# # Discretização no eixo x₁
	# X = 0:0.01:1
	
	# plt = plot_solução_aproximada([c;0], Nx1, Nx2, EQoLG)
	# plot(
	# plot(X,X,(x₁, x₂) -> u(x₁, x₂),
	# 	title="Solução Exata",
	# 	seriestype = :surface, colorbar=false, color=:viridis, 
	# 	xticks=[0, 0.5, 1],yticks=[0, 0.5, 1],zticks=[0, 1]),
	# plt,
	# layout=(1, 2))
end

# ╔═╡ 13fa4b92-ce6a-4c9f-bb88-6b91024cff4b
solução_aproximada_vs_exata()

# ╔═╡ 8ff37696-bb35-445b-8ac0-575f2c74b9db
md"## Estudo da convergência do erro"

# ╔═╡ b766d178-ce59-4055-bd43-035b66e43920
function estudo_do_erro()
	# Carrega os dados de entrada da EDP
    α, β, f, u = exemplo1()

    # Define os valores para Nx1 (subdivisões em x₁) e Nx2 (subdivisões em x₂)
    vec_Nx1 = [2^i for i in 2:8]
    vec_Nx2 = vec_Nx1

    # Calcula os valores de `h` para cada combinação de vec_Nx1[i] e vec_Nx2[i]
    vec_h = sqrt.( (1 ./ vec_Nx1).^2 + (1 ./ vec_Nx2).^2 )

    # Inicializa o vetor para armazenar os erros da norma L2
    vec_erro = zeros(length(vec_Nx1))

    # Loop sobre os valores de Nx1 e Nx2
    for i = 1:length(vec_Nx1)
        Nx1 = vec_Nx1[i]
        Nx2 = vec_Nx2[i]

        # Gera a matriz de conectividade local/global (LG)
        LG = monta_LG(Nx1, Nx2)

        # Gera o valor de `m` e o vetor `EQ`
        m, EQ = monta_EQ(Nx1, Nx2)

        # Define a matriz de conectividade EQoLG
        EQoLG = EQ[LG]

        # Monta a matriz esparsa `K`
        K = monta_K(α, β, Nx1, Nx2, m, EQoLG)

        # Monta o vetor global `F`
        F = monta_F(f, Nx1, Nx2, m, EQoLG)

        # Resolve o sistema linear Kc = F
        c = K \ F

        # Calcula o erro na norma L2 entre a solução exata `u` e a solução aproximada
        vec_erro[i] = erro_norma_L2(u, [c;0], Nx1, Nx2, EQoLG)
    end

    return vec_h, vec_erro
end

# ╔═╡ a173a245-a285-47a6-8ed6-8c3b540ed440
function plot_estudo_do_erro()
	# Realiza o estudo do erro e mede o tempo de execução
	@time begin
		vec_h, vec_erro = estudo_do_erro()
	end

    # Cria o gráfico do erro na norma L2 em função de h, diâmetro de cada Ωᵉ
    plt = plot(
        vec_h, vec_erro, lw=3, linestyle=:solid, markershape=:circle,
        label="Erro", title="Estudo do erro",
        xscale=:log10, yscale=:log10, legend=:topleft
    )

    # Adiciona a curva teórica de h² ao gráfico
    plot!(plt, vec_h, vec_h.^2, lw=3, linestyle=:solid, label="h²")

    # Adiciona rótulos aos eixos
    xlabel!("h")
    ylabel!("Erro")

    # Exibe uma tabela com os valores de h e erro
    display("Tabela com os valores de h e erro:")
    display(DataFrame(h=vec_h, erro=vec_erro))

	# Retorna o gráfico 
	return plt
end

# ╔═╡ 08ccecca-f127-4adb-bd0e-0fdb2979b62a
plt = plot_estudo_do_erro()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
GaussQuadrature = "d54b0c1a-921d-58e0-8e36-89d8069c0969"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
DataFrames = "~1.6.1"
GaussQuadrature = "~0.5.8"
LaTeXStrings = "~1.3.1"
Plots = "~1.40.8"
PlutoTeachingTools = "~0.2.15"
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.0"
manifest_format = "2.0"
project_hash = "17a0daaec47b389b8bd75fd155ceb6a0dfcd3261"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "7eee164f122511d3e4e1ebadb7956939ea7e1c77"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.6"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "ee28ddcd5517d54e417182fec3886e7412d3926f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.8"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f31929b9e67066bee48eec8b03c0df47d31a74b3"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.8+0"

[[deps.GaussQuadrature]]
deps = ["SpecialFunctions"]
git-tree-sha1 = "eb6f1f48aa994f3018cbd029a17863c6535a266d"
uuid = "d54b0c1a-921d-58e0-8e36-89d8069c0969"
version = "0.5.8"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "674ff0db93fffcd11a3573986e550d66cd4fd71f"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "39d64b09147620f5ffbf6b2d3255be3c901bec63"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.8"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "2984284a8abcfcc4784d95a9e2ea4e352dd8ede7"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.36"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "36bdbc52f13a7d1dcb0f3cd694e01677a515655b"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "b404131d06f7886402758c9ce2214b636eb4d54a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "96d2a4a668f5c098fb8a26ce7da53cde3e462a80"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "45470145863035bb124ca51b320ed35d071cc6c2"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.8"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "5d9ab1a4faf25a62bb9d07ef0003396ac258ef1c"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.15"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "2d4e5de3ac1c348fd39ddf8adbef82aa56b65576"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.6.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ff11acffdb082493657550959d4feb4b6149e73a"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.5"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a6b1675a536c5ad1a60e5a5153e1fee12eb146e3"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "936081b536ae4aa65415d869287d43ef3cb576b2"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.53.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╠═2f9d63db-0ffe-43b1-a29f-b4c92305a472
# ╠═0cc5dd29-8119-4fe5-9047-b2ee62eb6893
# ╟─354a7b37-db60-405d-a93b-a6d02a561cc0
# ╟─34ae963b-4aed-4504-9cd2-b6bfd5237833
# ╟─c2e3404b-1584-407b-8f37-11f0b72a204c
# ╟─57338d26-7c9a-471b-9005-3a2695753bf8
# ╟─ae6e9c26-3f0a-4d49-be01-63ebe7750630
# ╠═cbd7811e-f1c5-4815-9fa5-d9e7f406734b
# ╟─41ed9d0d-819b-439a-8c3e-b58e09005924
# ╟─db2d618a-6a27-4f19-93f9-522cf7fa9ad7
# ╟─50781f95-ee8f-4d17-ace3-8eb44d8ae727
# ╟─c8177a78-1849-4a77-b81f-91d4849df06f
# ╟─bf596662-c956-4ed4-8711-89e7870a04db
# ╟─2f25b167-978a-4e5b-a22f-a63a0fc6e593
# ╟─156276ab-1928-49c6-8cff-51ff7d987212
# ╟─aedd68e1-c955-43dd-834b-dc8719713569
# ╟─203da506-72ef-4ca1-9c10-d615ebad7715
# ╟─8e205232-794f-475a-912d-3cb822bc986f
# ╟─2c045f17-9dd7-4759-b19a-55933ee31728
# ╟─a497d172-0374-4471-bb14-0ae8f081bf5a
# ╟─b9f73af1-9be5-4142-8292-bccd18af7978
# ╟─45fcb91c-0721-45f6-a2bc-514858791a34
# ╟─e8a0bbc1-b3ca-49e0-9fd3-49cd883aa07f
# ╟─738c986b-e524-4cc1-a6f8-ff5650861a55
# ╟─93547f8d-e0f4-415f-abbe-800e52349746
# ╟─49c12023-3b0e-4246-9d28-33ebbbbbe0df
# ╟─26f290b8-ce09-44d8-8d60-aad4fa04367a
# ╟─6f269a14-472a-4cd3-b5a3-99a61d9a3f23
# ╟─417b8cd2-bcda-4354-9e8d-d9c2492f8d7a
# ╟─e1df5449-2dc0-44da-883c-05cefaa286d2
# ╟─61fccc23-8111-42ec-9fe0-2e3d4ef993b3
# ╟─f6cfbb5a-933b-48da-a488-0e21972ca15d
# ╟─8c7c8177-b8e8-4257-82b1-ba2125dc4d4a
# ╟─2450b17f-982b-4af8-bb68-7cca4f1944d6
# ╟─2f8e0a58-3476-4c04-aebe-4dec4f0f1db8
# ╟─0f85f3a1-beaa-4321-abb2-c671e75d0b15
# ╟─e0da177e-18da-4c57-bca4-bd432ee19a68
# ╟─dd0daf35-6cdf-4884-8f0c-2d5e0387db78
# ╟─a23031f6-073c-4d90-bbf3-280cd7d74bf5
# ╟─3e07b842-9d3e-4424-b187-a1a4a2ae5ddf
# ╟─a952e33c-e0ea-4d37-9a5b-3e3aec41ab63
# ╟─9158f0d7-78cd-4a5e-a24c-d648273eec52
# ╟─cdd4766a-b661-4bfb-9e2a-57708ce90438
# ╟─3661a97c-6d06-41d0-b7df-71a8c1f70adc
# ╟─4e4a5e81-ee38-4aad-a699-c44cb9f93178
# ╟─2e5bae39-727b-4401-b425-23e39f1d982f
# ╟─3e78c4c9-d906-40f1-aa73-221529730b02
# ╟─0c0a81e8-cdc3-41ed-a8ad-a06042fdf3bb
# ╟─cf54c8f7-f28d-4870-a306-227e25365b63
# ╟─399a4af6-593d-4a91-b0f9-b661ad2fba74
# ╟─90787447-7ddf-41dc-9b76-3261f4168e5e
# ╟─074adb78-415f-405e-a8cf-3f4988019819
# ╟─1e6a2be3-ea18-476b-ad15-b4a4584f4ff5
# ╟─958e69fc-12ec-432e-b4e9-babdc9a86c33
# ╟─52992179-f61b-4a18-bd57-fa96e571acf6
# ╟─2e274ee4-76c2-48f6-b980-34bf8a086849
# ╟─83e1a70c-9307-41fb-9e1e-0d49276ab03f
# ╟─98cb3a5e-9628-4e31-9869-bdb319393d60
# ╟─956a6fcb-ee3a-4bb2-aaf5-109236982854
# ╟─dcbca563-105a-48cf-95a5-bca36795eed7
# ╟─973ed431-0e2d-43f5-a407-fe73f868ab40
# ╟─dc66e8b9-2e82-4e3b-bca3-9e316510d5a3
# ╟─73081276-d764-401a-8f92-cc6b09620613
# ╟─301719ff-80d1-44ec-bb02-f92795d14e6a
# ╟─38e32ff9-a9c2-4018-b69b-8442563aa133
# ╟─8fd613d4-ee9f-4bde-aeb2-ea7950f822f6
# ╟─12d6be03-ef98-4174-b1cc-2743b1198cf8
# ╟─bfd24f64-a567-478c-a70b-addc2632e4bf
# ╟─14427308-9ef8-46e8-b53d-7c1deb5d4226
# ╟─8c2529d9-02e4-4d2c-a03a-cad4298b6f16
# ╟─3c066c9c-f0ac-46ff-82ca-7d4b7d316f4c
# ╟─6c1557e9-3192-4cf6-87d3-3a807e0b1169
# ╟─8543aafc-33f5-4e9b-afb2-8df58d12de81
# ╟─29097174-9c4e-4a47-afad-66a07c088e70
# ╟─fa49b560-e000-47cc-ba85-16c111d89acd
# ╟─8bef9e7a-6def-4751-bae8-d53343e94a3f
# ╟─de6cac6c-1f5c-4c44-9e1d-78c6712d3d45
# ╟─ea18dc94-c885-4b67-a4bd-5ab43258b42c
# ╟─13fa4b92-ce6a-4c9f-bb88-6b91024cff4b
# ╟─8ff37696-bb35-445b-8ac0-575f2c74b9db
# ╟─b766d178-ce59-4055-bd43-035b66e43920
# ╟─a173a245-a285-47a6-8ed6-8c3b540ed440
# ╟─08ccecca-f127-4adb-bd0e-0fdb2979b62a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
