# Calibration Model Confidence Bands

## Basic Parameters

### Univariate



### Multivariate 

```r
problemIndex = 2 # 0: GP process; 1: simu3; 2: simu4
linkIndex = 2 # 1: identity link 2: logit link
```


## RKHS Model

### $\theta$ class
```theta_fun.hpp```

* Constant theta $\theta(x) = \theta_0$.
* Parametric theta $\theta(x) = \gamma_0 exp(-\gamma_1 x)$.
* Parametric theta $\theta(x) = \gamma_0 + \gamma_1 x + \gamma_2 x^2$.
* Linear, cubic and squared exp reproducing kernel.

Observed $X$ is scaled to $[0,1]$. 

```thetaRKSH_cubic_v2.hpp``` class
The first-order soblov space for $x\in[0,1]$
$$ \theta(x) = \alpha_0 + \alpha_1 x + \sum_{i=1}^n \beta_i k(x, x_i) $$ with $k(x, y) = \min(x,y)$. All the parameter $\gamma = (\alpha_0, \alpha_1, \beta_1, \cdots, \beta_n)$.


First order derivative `theta_deriv(gamma)`. Output: each column for one observation, each row for one parameter. 
$$
\left(\frac{\partial \theta(x_1)}{\partial \gamma},
\frac{\partial \theta(x_2)}{\partial \gamma}, \cdots, 
\frac{\partial \theta(x_n)}{\partial \gamma}\right)
$$
at training data set.
Second order derivative `theta_deriv_2O(gamma, coef)`
$$
\sum_{i=1}^n c_i\frac{\partial^2 \theta(x_i)}{\partial \gamma\partial \gamma^T}
$$

### RKHS Class

Data X should be scaled to $[0,1]$ when computing kerenels.

### Link class

```link_class.hpp```

- **Identity link** $l(\theta) = \theta$
- **Logit link** $l(\theta) = \frac{U - L}{1 + \exp(-\theta)} + L\in[L,U]$.

The function $l(\cdot)$, $l'(\cdot)$ and 
$l''(\cdot)$ are implemented in `theta_trans`, 
`theta_trans_deriv` and `theta_trans_deriv2O`, 
respectively. 

### Main method

```methodRKHS.hpp```

Objective function $$ \frac{1}{2}\sum_{i=1}^n \left[y_i - y_m(l\circ\theta(x_i), x_i) \right]^2 + n\lambda \Vert \theta\Vert_\mathcal{H},$$

Gradient
$$
- \sum_{i=1}^n \left[y_i - y_m(l\circ\theta(x_i), x_i) \right] \frac{\partial y_m}{\partial \theta} 
\frac{\partial l}{\partial \theta}
\frac{\partial \theta(x_i)}{\partial \gamma}
+ 2 \lambda n K\gamma.
$$ with matrix $K = (k(x_i, x_j))$

Hessian
$$
\begin{align}
& \sum_{i=1}^n \left(\frac{\partial y_m}{\partial \theta} 
\frac{\partial l}{\partial \theta} \right)^2
\frac{\partial \theta(x_i)}{\partial \gamma} 
\frac{\partial \theta(x_i)}{\partial \gamma^T}\\
& - \sum_{i=1}^n \left[y_i - y_m(l\circ\theta(x_i), x_i) \right] 
\frac{\partial^2 y_m}{\partial \theta^2} 
\left(\frac{\partial l}{\partial \theta}\right)^2
\frac{\partial \theta(x_i)}{\partial \gamma}
\frac{\partial \theta(x_i)}{\partial \gamma^T}\\
& - \sum_{i=1}^n \left[y_i - y_m(l\circ\theta(x_i), x_i) \right] 
\frac{\partial y_m}{\partial \theta} 
\frac{\partial^2 l}{\partial \theta^2}
\frac{\partial \theta(x_i)}{\partial \gamma}
\frac{\partial \theta(x_i)}{\partial \gamma^T}\\
& - \sum_{i=1}^n \left[y_i - y_m(l\circ\theta(x_i), x_i) \right] \frac{\partial y_m}{\partial \theta} 
\frac{\partial l}{\partial \theta}
\frac{\partial^2 \theta(x_i)}{\partial \gamma \partial \gamma^T}
+ 2 \lambda n K.\\
= & - \sum_{i=1}^n td_i
\frac{\partial \theta(x_i)}{\partial \gamma}
\frac{\partial \theta(x_i)}{\partial \gamma^T}\\
& - \sum_{i=1}^n \left[y_i - y_m(l\circ\theta(x_i), x_i) \right] \frac{\partial y_m}{\partial \theta} 
\frac{\partial l}{\partial \theta}
\frac{\partial^2 \theta(x_i)}{\partial \gamma \partial \gamma^T}
+ 2 \lambda n K.
\end{align}
$$
with ```derivTmp```
$$
td_i = \left(\frac{\partial y_m}{\partial \theta}  \frac{\partial l}{\partial \theta} \right)^2 - \left[y_i - y_m(l\circ\theta(x_i), x_i) \right] 
\frac{\partial^2 y_m}{\partial \theta^2} 
\left(\frac{\partial l}{\partial \theta}\right)^2- 
\left[y_i - y_m(l\circ\theta(x_i), x_i) \right] 
\frac{\partial y_m}{\partial \theta} 
\frac{\partial^2 l}{\partial \theta^2}
$$
```thetaCore->theta_deriv_2O(gamma, coef)```compute the value of  $\sum_{i=1}^n c_i
\frac{\partial^2 \theta(x_i)}{\partial \gamma \partial \gamma^T} $. 

## GP for expensive code

- `squaredExp.*pp`
- `gaussianProcess.*pp`

### Kernel

> Input observation is column-wise.

The kernel `scalingParam` default value is 1.

For kernel
$$
\theta_i = \exp(\beta_i),\ i = 1,2
$$

$$
\theta_1 \exp\left(- \Vert X_1 - X_2 \Vert_2^2/\theta_2\right)
$$

First order derivative `kernelCov_partial_cj(X1, X2, cj)` with respect to the `j`th component of `X1`
$$
\theta_1 \exp\left(- \Vert X_1 - X_2 \Vert_2^2/\theta_2\right)
\times \frac{-2(X_{1j} - X_{2j})}{\theta_2}
$$
Second order derivative `kernelCov_partial2O_cj(X1, X2, cj)` with respect to the `j`th component of `X1`
$$
\theta_1 \exp\left(- \Vert X_1 - X_2 \Vert_2^2/\theta_2\right)
\times \frac{4(X_{1j} - X_{2j})^2}{\theta_2^2}
- \frac{2\theta_1}{\theta_2}\exp\left(- \Vert X_1 - X_2 \Vert_2^2/\theta_2\right)
$$
Second order derivative with respect to `j`,`k`
$$
\theta_1 \exp\left(- \Vert X_1 - X_2 \Vert_2^2/\theta_2\right)
\times \frac{2(X_{1j} - X_{2j})}{\theta_2}
\times \frac{2(X_{1k} - X_{2k})}{\theta_2}
$$
implemented in `kernelCov_partial2OCross`.

### GP

> Input observation is row-wise.
> Internal computation is column-wise. 

Given $\beta$ and data, the $\alpha_i$'s are fixed training parameters. 

For emulator, predict at new value $$
y_e( \theta, x) = \sum_{i=1}^n \alpha_i e^{\beta_1}
\exp\left(-\frac{1}{ e^{\beta_2}}[(x - x_i)^2 + (\theta - \theta_i)^2]\right )
$$
Partial derivative with respect to $\theta$
$$
\frac{\partial}{\partial \theta}\, y_e( \theta, x) = \sum_{i=1}^n \alpha_i e^{\beta_1}
\exp\left(-\frac{1}{ e^{\beta_2}}[(x - x_i)^2 + (\theta - \theta_i)^2]\right )
\times \frac{-2(\theta - \theta_i)}{ e^{\beta_2}}
$$
which is implemented by the function ```kernelCov_partialX1_cj``` in the file ```squaredExp.hpp```. 


## Multivariate Version


### Parametric Model



Objective function 

$$ \frac{1}{2}\sum_{i=1}^n 
\left[y_i - y_m(l_1 \circ\theta_1(x_i), 
\cdots, l_q \circ\theta_q(x_i), x_i) \right]^2,$$

Gradient with respect to $\gamma_{jk}$ the $k$th parameter in group $j$. $\gamma_{j} = (\gamma_{j1},\gamma_{j2},\cdots)^T$
$$
- \sum_{i=1}^n \left[y_i - y_m(l_1 \circ\theta_1(x_i), 
\cdots, l_q \circ\theta_q(x_i), x_i) \right] \frac{\partial y_m}{\partial \theta_j} 
\frac{\partial l_j}{\partial \theta_j}
\frac{\partial \theta_j(x_i)}{\partial \gamma_{j}}.
$$ 

Hessian, diagonal blocks
$$
\begin{align}
& \sum_{i=1}^n \left(\frac{\partial y_m}{\partial \theta_j} 
\frac{\partial l_j}{\partial \theta_j} \right)^2
\frac{\partial \theta_j(x_i)}{\partial \gamma_j} 
\frac{\partial \theta_j(x_i)}{\partial \gamma^T_j}\\
& - \sum_{i=1}^n \left[y_i - y_m( x_i) \right] 
\frac{\partial^2 y_m}{\partial \theta_j^2} 
\left(\frac{\partial l_j}{\partial \theta_j}\right)^2
\frac{\partial \theta_j(x_i)}{\partial \gamma_j}
\frac{\partial \theta_j(x_i)}{\partial \gamma_j^T}\\
& - \sum_{i=1}^n \left[y_i - y_m( x_i) \right] 
\frac{\partial y_m}{\partial \theta_j} 
\frac{\partial^2 l_j}{\partial \theta_j^2}
\frac{\partial \theta_j(x_i)}{\partial \gamma_j}
\frac{\partial \theta_j(x_i)}{\partial \gamma_j^T}\\
& - \sum_{i=1}^n \left[y_i - y_m(x_i) \right] \frac{\partial y_m}{\partial \theta_j} 
\frac{\partial l_j}{\partial \theta_j}
\frac{\partial^2 \theta_j(x_i)}{\partial \gamma_j \partial \gamma_j^T}\\
= & - \sum_{i=1}^n td_i
\frac{\partial \theta_j(x_i)}{\partial \gamma_j}
\frac{\partial \theta_j(x_i)}{\partial \gamma_j^T}\\
& - \sum_{i=1}^n \left[y_i - y_m( x_i) \right] \frac{\partial y_m}{\partial \theta_j} 
\frac{\partial l_j}{\partial \theta_j}
\frac{\partial^2 \theta_j(x_i)}{\partial \gamma_j \partial \gamma_j^T}.
\end{align}
$$
with ```derivTmp```
$$
td_i = \left(\frac{\partial y_m}{\partial \theta_j}  \frac{\partial l_j}{\partial \theta_j} \right)^2 
- \left[y_i - y_m(x_i) \right] 
\frac{\partial^2 y_m}{\partial \theta_j^2} 
\left(\frac{\partial l_j}{\partial \theta_j}\right)^2- 
\left[y_i - y_m( x_i) \right] 
\frac{\partial y_m}{\partial \theta_j} 
\frac{\partial^2 l_j}{\partial \theta_j^2}
$$

Off-diagonal blocks for $j\neq k$ is
$$
\begin{align}
& \sum_{i=1}^n \left(\frac{\partial y_m}{\partial \theta_j} 
\frac{\partial l_j}{\partial \theta_j} \right)
\left(\frac{\partial y_m}{\partial \theta_k} 
\frac{\partial l_k}{\partial \theta_k} \right) \times
\frac{\partial \theta_j(x_i)}{\partial \gamma_j} 
\frac{\partial \theta_k(x_i)}{\partial \gamma^T_k}\\
& - \sum_{i=1}^n \left[y_i - y_m( x_i) \right] 
\frac{\partial^2 y_m}{\partial \theta_j\partial \theta_k} 
\left(\frac{\partial l_j}{\partial \theta_j} \frac{\partial l_k}{\partial \theta_k}\right)
\frac{\partial \theta_j(x_i)}{\partial \gamma_j}
\frac{\partial \theta_k(x_i)}{\partial \gamma_k^T}
\end{align}
$$



## Version One


The objective function for the observed data $\mathcal{D} = \{(x_i, y_i)\}_{i=1}^n$, 
$$ f(\gamma) = \frac{1}{2}\sum_{i=1}^n \left[y_i - y_m(l\circ\theta(x_i), x_i) \right]^2 + n\lambda \Vert \theta\Vert_\mathcal{H},$$

 * Use function $l(\theta)$ to constraint the range of $\theta(x)$. *Identity link* $l(\theta) = \theta$. *Logit link* $l(\theta) = \frac{U - L}{1 + \exp(-\theta)} + L\in[L,U]$.
 * The cubic splines for $x\in[0,1]$ $ \theta(x) = \alpha_0 + \alpha_1 x + \sum_{i=1}^n \beta_i k(x, x_i) $ with $k(x, y)$ being the corresponding kernel. All the parameter $\gamma = (\alpha_0, \alpha_1, \beta_1, \cdots, \beta_n)$.
 * The objective function is optimized for $\gamma$ as parameters.

The posterior is proportion to 
$$
p(\gamma | \mathcal{D}) \propto \exp\left\{-\frac{f(\gamma)}{\hat{\sigma}^2}\right\}
$$
The $\hat{\sigma}^2$ is estimated from the residuals $r_i = y_i - \hat{y}_i$. At the optimal $\hat{\gamma}$, the first order derivative $\frac{\partial}{\partial \gamma} f(\hat{\gamma}) = 0$, and the posterior is approximated by $$
p(\gamma | \mathcal{D}) \propto \exp\left\{-\frac{1}{\hat{\sigma}^2} (\gamma-\hat{\gamma}\,)^T\frac{\partial^2 f(\hat{\gamma}\,)}{\partial \gamma\partial\gamma^T}(\gamma-\hat{\gamma}\,)\right\}
$$ The posterior of $\gamma | \mathcal{D}$ is approximately normal with mean $\hat{\gamma}$, and the covariance matrix $$\Sigma = \mathrm{Var}(\gamma - \hat{\gamma}) = \hat{\sigma}^2\left(\frac{\partial^2 f(\hat{\gamma}\,)}{\partial \gamma\partial\gamma^T}\right)^{-1}\, .$$

At an abitratry location $x$, let $\xi = (1, x, k(x,x_1), \cdots, k(x, x_n))^T$. The predicted value of $\theta(x)$ is $\hat{\theta}(x) = \xi^T\hat{\gamma}$. Its uncertainty is denoted by $e(x)$, and computed as
$$
e^2(x) = \mathrm{Var}\left(\theta(x) - \hat{\theta}(x)\right)
= \mathrm{Var}\left(\xi^T(\gamma- \hat{\gamma})\right)
= \xi^T \Sigma \xi
$$

The final *calibration parameter of interest* is $l(\hat{\theta})$, and its uncertainty is computed as
$$
\begin{align}
d^2(x) = \mathrm{Var}\left(l(\theta) - l(\hat{\theta})\right) 
&\approx \mathrm{Var}\left\{l'(\hat{\theta})(\theta - \hat{\theta})
+0.5 \times l''(\hat{\theta})(\theta - \hat{\theta})^2
\right\} \\
& = [l'(\hat{\theta})]^2\mathrm{Var}(\theta - \hat{\theta})
+\frac{1}{4} [l''(\hat{\theta})]^2\mathrm{Var}(\theta - \hat{\theta})^2  \\
& = [l'(\hat{\theta})]^2 e^2(x)
+\frac{1}{4} [l''(\hat{\theta})]^2 [3e^4(x) - e^4(x)] \\ 
& = [l'(\hat{\theta})]^2 e^2(x)
+\frac{1}{2} [l''(\hat{\theta})]^2 e^4(x) \\ 
\end{align} 
$$
Finally, the $(1-\alpha)$ pointwise coverage for $l(\hat{\theta})$ is
$$
\left[ l(\hat{\theta}) - z_{\alpha/2}d(x), l(\hat{\theta}) + z_{\alpha/2}d(x)\right]
$$

**Note**

The gradient of the objective function
$$ - \sum_{i=1}^n \left[y_i - y_m(l\circ\theta(x_i), x_i) \right] \frac{\partial y_m}{\partial \theta} 
\frac{\partial l}{\partial \theta}
\frac{\partial \theta(x_i)}{\partial \gamma}
+ 2 \lambda n K\gamma. $$ 
with matrix $K = \left(k(x_i, x_j)\right)_{n\times n}$. 
The hessian of the objective function
$$
\begin{align}
& \sum_{i=1}^n \left(\frac{\partial y_m}{\partial \theta} 
\frac{\partial l}{\partial \theta} \right)^2
\frac{\partial \theta(x_i)}{\partial \gamma} 
\frac{\partial \theta(x_i)}{\partial \gamma^T}\\
& - \sum_{i=1}^n \left[y_i - y_m(l\circ\theta(x_i), x_i) \right] 
\frac{\partial^2 y_m}{\partial \theta^2} 
\left(\frac{\partial l}{\partial \theta}\right)^2
\frac{\partial \theta(x_i)}{\partial \gamma}
\frac{\partial \theta(x_i)}{\partial \gamma^T}\\
& - \sum_{i=1}^n \left[y_i - y_m(l\circ\theta(x_i), x_i) \right] 
\frac{\partial y_m}{\partial \theta} 
\frac{\partial^2 l}{\partial \theta^2}
\frac{\partial \theta(x_i)}{\partial \gamma}
\frac{\partial \theta(x_i)}{\partial \gamma^T}\\
& - \sum_{i=1}^n \left[y_i - y_m(l\circ\theta(x_i), x_i) \right] \frac{\partial y_m}{\partial \theta} 
\frac{\partial l}{\partial \theta}
\frac{\partial^2 \theta(x_i)}{\partial \gamma \partial \gamma^T}
+ 2 \lambda n K.
\end{align}
$$
