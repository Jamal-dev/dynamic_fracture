# Installing deal.II

[Deal.II](https://www.dealii.org/) can be download from the offical websit.

For easy installation please follow instructions for [candi](https://github.com/dealii/candi).
For linux distribution

``` bash
git clone https://github.com/dealii/candi.git
cd candi
./candi.sh
```

```latex {cmd = True}
\documentclass{standalone}
\usepackage{bm}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{mathtools}
\begin{document}
  $$
  \bm{a}
  $$
\end{document}```

# Dynamic Phase field equations

Find a solution $\{\bm{u}, \varphi \} \in \{ \bm{u}_D + V\} \times W$ such that:

$$
\begin{align*}
\rho \left( \partial_{tt} \bm{u}, \bm{w} \right)
+
\left( \left( (1-\kappa) \varphi^2 +\kappa \right) {\bm{\sigma}}^+ , e(\bm{w}) \right)
+
\left( {\bm{\sigma}}^- , e(\bm{w}) \right)
& =0 & & \forall \bm{w} \in V
\\
(1-\kappa) \left( \varphi \left(\bm{\sigma}^+ : e(\bm{u}) \right) , \psi - \varphi \right)
+
G_c
\left( \frac{1}{\epsilon} (1- \varphi) , \psi - \varphi \right)
+
G_c
\left(    \nabla \varphi   , \nabla \psi - \nabla \varphi \right)
& \geq 0 & & \forall \psi \in W_{in}
\end{align*}
$$

* $\bm{u}$ is the displacement $(\bm{u} \coloneqq \Omega \longrightarrow \R^2)$
* $\bm{\sigma}$ is the stress tensor
* $e(\bm{u})$ is the symmetric strain tensor $(e(\bm{u}) \coloneqq \frac{1}{2} \left( \nabla{\bm{u}} + \nabla{\bm{u}}^T \right))$
* $G_c$ is the energy release rate, and it is strictly positive
* $\kappa$ and $\epsilon$ are the regularization parameters
* ${\bm{\sigma}}^-$ is the compressive stress tensor
* ${\bm{\sigma}}^+$ is the tensile stress tensor

By introducing another variable velocity $\bm{v} \coloneq \R^2 \longrightarrow \R^2$, we can write the formulation as:

#### Formulation

Find a solution $\{\bm{u}, \bm{v}, \varphi \} \in \{ \bm{u}_D + V\} \times \{ \bm{v}_D + L\} \times W$ such that:

$$
\begin{align*}
\rho \left( \partial_{t} \bm{v}, \bm{w} \right)
+
\left( \left( (1-\kappa) \varphi^2 +\kappa \right) {\bm{\sigma}}^+ , e(\bm{w}) \right)
+
\left( {\bm{\sigma}}^- , e(\bm{w}) \right)
& =0 & & \forall \bm{w} \in V
\\
(1-\kappa) \left( \varphi \left(\bm{\sigma}^+ : e(\bm{u}) \right) , \psi - \varphi \right)
+
G_c
\left( \frac{1}{\epsilon} (1- \varphi) , \psi - \varphi \right)
+
G_c
\left(    \nabla \varphi   , \nabla \psi - \nabla \varphi \right)
& \geq 0 & & \forall \psi \in W_{in}
\\
\rho \left( \partial_t \bm{u} - \bm{v}, \bm{w} \right)
& =0 & & \forall \bm{w} \in L
\end{align*}
$$

## Temporal Discritization

### One-step-$\theta$ scheme

### Formulation

Let $\bm{u}^{n-1}$, and $\bm{v}^{n-1}$ be the displacement and velocity at the previous time step. Find $\bm{u}, \bm{v} \in \{ \bm{u}_D + V\} \times \{ \bm{v}_D + L\}$ such that:

$$
\begin{align*}
\rho \left( \frac{ \bm{v} - {\bm{v}}^{n-1} }{\delta t} , \bm{w} \right)
& +
\theta 
\left( \left( (1-\kappa) \varphi^2 +\kappa \right) {{\bm{\sigma}}^+}({\bm{u}}) , e(\bm{w}) \right)
& & &
\\
&+
(1-\theta) 
\left( \left( (1-\kappa) \varphi^2 +\kappa \right) 
{{\bm{\sigma}}^+} ({\bm{u}}^{n-1})
, e(\bm{w}) 
\right)
& & &
\\
&+
\theta 
    \left( {\bm{\sigma}}^-({\bm{u}}) , e(\bm{w}) \right)
& & &
\\
&+
(1-\theta) 
    \left( {\bm{\sigma}}^-({\bm{u}}^{n-1}) , e(\bm{w}) \right)
& =0 & & \forall \bm{w} \in V
\\
% second equation
(1-\kappa) \left( \varphi \left(\bm{\sigma}^+ : e(\bm{u}) \right) , \psi - \varphi \right)
& +
G_c
\left( \frac{1}{\epsilon} (1- \varphi) , \psi - \varphi \right)
+
G_c
\left(    \nabla \varphi   , \nabla \psi - \nabla \varphi \right)
& \geq 0 & & \forall \psi \in W_{in}
\\
% third equation
\rho \left( \frac{\bm{u} - {\bm{u}}^{n-1} }{\delta t}, \bm{w} \right) &
- (\theta \rho) \left(  \bm{v}, \bm{w} \right)
- \left( (1-\theta) \rho \right) \left(  {\bm{v}}^{n-1}, \bm{w} \right) 
& =0 & & \forall \bm{w} \in L
\end{align*}
$$

# Stress $\bm{\sigma}$

The stress tensor is defined as:

$$ 
\bm{\sigma}(\bm{u}) \coloneqq 2 \mu \ e(\bm{u}) + \lambda \ tr(e(\bm{u})) \ \bm{I} 
$$

where $\mu$ and $\lambda$ are the Lame parameters, and $\bm{I}$ is the second order identity tensor.

# Stress decomposition

The stress tensor can be decomposed into two parts as suggested by Miehe et all (2010):

$$
{\bm{\sigma}} \coloneqq 
                          \left( 
        \left(1-\kappa \right) \varphi^2 + \kappa
                          \right)
                          {\bm{\sigma}}^+
                          + {\bm{\sigma}}^-
$$

$$
\begin{align*}
    {\bm{\sigma}}^{+} &= 2 \mu \ e(\bm{u}^{+}) + \lambda \ tr(e(\bm{u}^{+})) \ \bm{I} \ ,\\
    {\bm{\sigma}}^{-} &= 2 \mu \ \left( e(\bm{u})-e(\bm{u}^{+}) \right)
    + \lambda \ \left( tr(e(\bm{u}))- tr(e(\bm{u}^{+})) \right) \ \bm{I} \ .
\end{align*}
$$

This e^+ and e^- can be obtained by eigen value decomposition of the strain tensor.

Let $P$ be matrix of eigen vectors, and $\lambda_i$ be a eigen value coresponding to its column. By modal decomposition of $e$

$$            
    \begin{align*}
                e &\coloneqq P \Lambda P^T
                \\
                e^+ &\coloneqq P \Lambda^+ P^T
                \\
                e^- &\coloneqq P \Lambda^- P^T
                \\
                \Lambda^+ &= diag\left(
                            [<\lambda_i>]_{i=1}^2
                            \right)
                \\
                \Lambda^- &= \Lambda - \Lambda^+
                \\
                <x> &\coloneqq 
                \begin{dcases}
                            x, & \text{ if } x>0
                            \\
                            0, & \text{ if } x \leq 0
                \end{dcases}
    \end{align*}
$$

# Acknowledgement
This code is initially developed by ![Prof. Thomas Wick](https://thomaswick.org/). Then modified by Monika Senftl in her master thesis. Her thesis pdf is included in the docs folder. The code is further modified by me to make it more general and easy to use.
