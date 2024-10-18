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
