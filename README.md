# Installing deal.II

[Deal.II](https://www.dealii.org/) can be download from the offical websit.

For easy installation please follow instructions for [candi](https://github.com/dealii/candi).
For linux distribution

``` bash
git clone https://github.com/dealii/candi.git
cd candi
./candi.sh
```

# Dynamic Phase field equations

Find a solution $\{\boldsymbol{u}, \varphi \} \in \{ \boldsymbol{u}_D + V\} \times W$ such that:

$$
\begin{align*}
\rho \left( \partial_{tt} \boldsymbol{u}, \boldsymbol{w} \right)
+
\left( \left( (1-\kappa) \varphi^2 +\kappa \right) {\boldsymbol{\sigma}}^+ , e(\boldsymbol{w}) \right)
+
\left( {\boldsymbol{\sigma}}^- , e(\boldsymbol{w}) \right)
& =0 & & \forall \boldsymbol{w} \in V
\\
(1-\kappa) \left( \varphi \left(\boldsymbol{\sigma}^+ : e(\boldsymbol{u}) \right) , \psi - \varphi \right)
+
G_c
\left( \frac{1}{\epsilon} (1- \varphi) , \psi - \varphi \right)
+
G_c
\left(    \nabla \varphi   , \nabla \psi - \nabla \varphi \right)
& \geq 0 & & \forall \psi \in W_{in}
\end{align*}
$$

* $\boldsymbol{u}$ is the displacement $(\boldsymbol{u} \coloneqq \Omega \longrightarrow \R^2)$
* $\boldsymbol{\sigma}$ is the stress tensor
* $e(\boldsymbol{u})$ is the symmetric strain tensor $(e(\boldsymbol{u}) \coloneqq \frac{1}{2} \left( \nabla{\boldsymbol{u}} + \nabla{\boldsymbol{u}}^T \right))$
* $G_c$ is the energy release rate, and it is strictly positive
* $\kappa$ and $\epsilon$ are the regularization parameters
* ${\boldsymbol{\sigma}}^-$ is the compressive stress tensor
* ${\boldsymbol{\sigma}}^+$ is the tensile stress tensor

By introducing another variable velocity $\boldsymbol{v} \coloneq \R^2 \longrightarrow \R^2$, we can write the formulation as:

#### Formulation

Find a solution $\{\boldsymbol{u}, \boldsymbol{v}, \varphi \} \in \{ \boldsymbol{u}_D + V\} \times \{ \boldsymbol{v}_D + L\} \times W$ such that:

$$
\begin{align*}
\rho \left( \partial_{t} \boldsymbol{v}, \boldsymbol{w} \right)
+
\left( \left( (1-\kappa) \varphi^2 +\kappa \right) {\boldsymbol{\sigma}}^+ , e(\boldsymbol{w}) \right)
+
\left( {\boldsymbol{\sigma}}^- , e(\boldsymbol{w}) \right)
& =0 & & \forall \boldsymbol{w} \in V
\\
(1-\kappa) \left( \varphi \left(\boldsymbol{\sigma}^+ : e(\boldsymbol{u}) \right) , \psi - \varphi \right)
+
G_c
\left( \frac{1}{\epsilon} (1- \varphi) , \psi - \varphi \right)
+
G_c
\left(    \nabla \varphi   , \nabla \psi - \nabla \varphi \right)
& \geq 0 & & \forall \psi \in W_{in}
\\
\rho \left( \partial_t \boldsymbol{u} - \boldsymbol{v}, \boldsymbol{w} \right)
& =0 & & \forall \boldsymbol{w} \in L
\end{align*}
$$

## Temporal Discritization

### One-step-$\theta$ scheme

### Formulation

Let $\boldsymbol{u}^{n-1}$, and $\boldsymbol{v}^{n-1}$ be the displacement and velocity at the previous time step. Find $\boldsymbol{u}, \boldsymbol{v} \in \{ \boldsymbol{u}_D + V\} \times \{ \boldsymbol{v}_D + L\}$ such that:

$$
\begin{align*}
\rho \left( \frac{ \boldsymbol{v} - {\boldsymbol{v}}^{n-1} }{\delta t} , \boldsymbol{w} \right)
& +
\theta 
\left( \left( (1-\kappa) \varphi^2 +\kappa \right) {{\boldsymbol{\sigma}}^+}({\boldsymbol{u}}) , e(\boldsymbol{w}) \right)
& & &
\\
&+
(1-\theta) 
\left( \left( (1-\kappa) \varphi^2 +\kappa \right) 
{{\boldsymbol{\sigma}}^+} ({\boldsymbol{u}}^{n-1})
, e(\boldsymbol{w}) 
\right)
& & &
\\
&+
\theta 
    \left( {\boldsymbol{\sigma}}^-({\boldsymbol{u}}) , e(\boldsymbol{w}) \right)
& & &
\\
&+
(1-\theta) 
    \left( {\boldsymbol{\sigma}}^-({\boldsymbol{u}}^{n-1}) , e(\boldsymbol{w}) \right)
& =0 & & \forall \boldsymbol{w} \in V
\\
% second equation
(1-\kappa) \left( \varphi \left(\boldsymbol{\sigma}^+ : e(\boldsymbol{u}) \right) , \psi - \varphi \right)
& +
G_c
\left( \frac{1}{\epsilon} (1- \varphi) , \psi - \varphi \right)
+
G_c
\left(    \nabla \varphi   , \nabla \psi - \nabla \varphi \right)
& \geq 0 & & \forall \psi \in W_{in}
\\
% third equation
\rho \left( \frac{\boldsymbol{u} - {\boldsymbol{u}}^{n-1} }{\delta t}, \boldsymbol{w} \right) &
- (\theta \rho) \left(  \boldsymbol{v}, \boldsymbol{w} \right)
- \left( (1-\theta) \rho \right) \left(  {\boldsymbol{v}}^{n-1}, \boldsymbol{w} \right) 
& =0 & & \forall \boldsymbol{w} \in L
\end{align*}
$$

# Stress $\boldsymbol{\sigma}$

The stress tensor is defined as:

$$ 
\boldsymbol{\sigma}(\boldsymbol{u}) \coloneqq 2 \mu \ e(\boldsymbol{u}) + \lambda \ tr(e(\boldsymbol{u})) \ \boldsymbol{I} 
$$

where $\mu$ and $\lambda$ are the Lame parameters, and $\boldsymbol{I}$ is the second order identity tensor.

# Stress decomposition

The stress tensor can be decomposed into two parts as suggested by Miehe et all (2010):

$$
{\boldsymbol{\sigma}} \coloneqq 
                          \left( 
        \left(1-\kappa \right) \varphi^2 + \kappa
                          \right)
                          {\boldsymbol{\sigma}}^+
                          + {\boldsymbol{\sigma}}^-
$$

$$
\begin{align*}
    {\boldsymbol{\sigma}}^{+} &= 2 \mu \ e(\boldsymbol{u}^{+}) + \lambda \ tr(e(\boldsymbol{u}^{+})) \ \boldsymbol{I} \ ,\\
    {\boldsymbol{\sigma}}^{-} &= 2 \mu \ \left( e(\boldsymbol{u})-e(\boldsymbol{u}^{+}) \right)
    + \lambda \ \left( tr(e(\boldsymbol{u}))- tr(e(\boldsymbol{u}^{+})) \right) \ \boldsymbol{I} \ .
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
