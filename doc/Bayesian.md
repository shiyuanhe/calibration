
## Nonparametric Functional Calibration of Computer Models

The proposed model
$$
\begin{align}
y | \theta^{(x)}, \lambda_y &\sim \mathcal{N}_N(\eta(\theta^{(x)}, 
\lambda^{-1} I) \\
\lambda_y &\sim \text{Ga}(a_y, b_y),\ a_y > 0,\ b_y > 0 \\
\theta_1 | \lambda_\theta, \rho_\theta &\sim
\mathcal{GP}(\mu_\theta, \lambda_\theta^{-1} 
R_{\rho\theta}(\cdot, \cdot)) \\
\theta_2 &\sim\text{ Unif}(0,1) \\
\lambda_\theta &\sim \text{Ga}(a_\theta, b_\theta),
a_\theta, b_\theta > 0 \\
\rho_\theta &\sim \text{Beta}(1, b_\theta), b_\theta >0.
\end{align}
$$
The joint posterior is
$$
\begin{align}
& \pi(\theta_1^{(x)}, \theta_2, \rho_\theta, \lambda_\theta,
\lambda_y | y) \\
\propto\ &\lambda_y^{N/2}\exp\left\{
-\frac{\lambda_y}{2} (y - \eta(\theta_1^{(x)}) )^T
(y - \eta(\theta_1^{(x)}) )\right\}
\times \lambda_y^{a_y-1} \exp(-b_y \lambda_y) \\
&\times \lambda_\theta^{N/2} \vert R_\rho \vert^{-1/2}
\exp\left\{ 
-\frac{\lambda_\theta}{2} (\theta_1^{(x)} - \mu_\theta 1)^T
R^{-1}_\rho (\theta_1^{(x)} - \mu_\theta 1) 
\right\} \\
&\times \rho_\theta^{a_\theta -1}
\exp(-b_\theta \lambda_\theta) (1-\rho_\theta)^{b_\theta -1}
\end{align}
$$
with parameterization $\nu = \log(-\log(\rho_\theta))$.

> **Remark:** When it is swifted to the multivariate setting, $\theta$ becomes a matrix. Each column represents one multivariate component. Both $\nu$ and $\lambda_\theta$ become vectors and each vector element corresponds to one multivariate component. 



The correlation $\rho = \exp(-\exp(\nu)) \sim \mathrm{Beta}(1, b_\theta)$. 
For observation $x_i, x_j$, their correlation is 
$$
\rho^{\Vert x_i - x_j\Vert^2} =
\exp\left(-\exp(\nu)\Vert x_i - x_j\Vert^2\right) = 
e^{\beta_1}  \exp\left(-\Vert x_i - x_j\Vert^2\right / e^{\beta_2})
$$
which corresponds to squared exponential kernel with $\beta_1 = 0$ and $\beta_2 = -\nu$.


Define notation:
$$
\begin{align}
\mathrm{expS1} &= -\frac{1}{2}[y - \eta(\theta(x), x )]^T [y - \eta(\theta(x), x)] \\
\mathrm{expS2} &= -\frac{1}{2}(\theta - \mu_\theta 1)^TR_\nu^{-1} (\theta - \mu_\theta 1)
\end{align}
$$
Note we have  decomposition $R_\nu = U\Lambda U^T$.


* In the following, the gamma distribution has shape parameter $a$ and rate parameter $1/s$. 
* In the coding, the gamma distribution has shape parameter $a$ and scaling parameter $s$.

### Update $\nu$:
$$
\pi(\nu|\theta(x), \lambda_\theta) \propto
\exp\left(
-0.5\log\det R_\nu + \lambda_\theta \times\mathrm{expS2} + \nu - e^\nu
+ (b_\theta - 1) (1 - \exp\{-e^\nu\})
\right)
$$
Propose new value by $\nu^+ = \nu + C_1 z$.

> This is specific to each multivariate component. $R_v$ and $v$

### Update $\theta$:
$$
\pi(\theta | \lambda_\theta, \lambda_y, \nu) \propto
\exp\left(\lambda_y\times\mathrm{expS1} + \lambda_\theta\times\mathrm{expS2}\right)
$$
Propose new value by $\theta^+ = \theta + C_2 U\Lambda^{1/2} z$, with  decomposition $R_\nu = U\Lambda U^T$.

> This is specific to each multivariate component. $R_\nu = U\Lambda U^T$
> and $\mathrm{expS2}$


### Update $\lambda_\theta$ and $\lambda_y$:
$$
\pi(\lambda_y | \lambda_\theta, \theta, \nu) \propto
\mathrm{Ga}(a_y + N/2, b_y - \mathrm{expS1})
$$
> This is general for all components.


$$
\pi(\lambda_\theta | \lambda_y, \theta, \nu) \propto
\mathrm{Ga}(a_\theta + N/2, b_\theta - \mathrm{expS2})
$$
> This is specific to each multivariate component. $\mathrm{expS2}$



## Modification

- For MCMC, the outputed storage has changed to row order. Change this in R code!
- Unified framework. One component as a special case? 
- methodBase, one more init para, numComp
- How about the function `computeLikS_nu`?
- Code testing? 
- Constructor now has four parameters


## Parameters

