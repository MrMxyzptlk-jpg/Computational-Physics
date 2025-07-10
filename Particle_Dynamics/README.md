# Particle Dynamics
Currently only one short-range potential (Lennard-Jones) and one long-range potential (Coulomb) are implemented. The former supports three formalisms (MD, MC and BM) while the latter only two (MD and BM).

## Lennar-Jones potential


## Coulomb potential
We consider the summation of the particles in the reference super-cell to get the real space contribution (short-range contribution) of the potential and forces. The long-range term is taken into account by summing over a ball of k-vectors in reciprocal space. Thus we get:
$$ V = V_r + V_k - V_s + V_0 $$

$$ V_r = \frac{1}{2} \sum_i^N \sum_j^N q_i q_j \left( \sum_{\mathbf{R}} \frac{\text{erfc}(\frac{|\mathbf{r}_{ij}+ \mathbf{R}|}{\sigma})}{|\mathbf{r}_{ij}+ \mathbf{R}|}\right)  $$

$$ V_k = \frac{2\pi}{V} \sum_{\mathbf{k} \neq 0} \hat{G}(k)\hat{\rho}^q\bf(k) \hat{\rho}^q\bf(-k) $$

$$ V_s = \sum_{i} \frac{q_i^2}{\sqrt{\pi}\sigma} $$

$$ V_0 = \frac{2\pi}{3V} |\sum_i q_i\mathbf{r}_i |^2 $$

$$ \hat{G}(k) = 4\pi e^{-(\frac{k\sigma}{2})^2} $$

$$ \hat{\rho}^q(\mathbf{k}) = \sum_i q_i e^{-i\mathbf{k \cdot r}}  $$

where $V_r$ is calculated in real space, $V_k$ in reciprocal space, $V_s$ is the self interaction term and $V_0$ is a dipole term for charges in a vacuum. $\hat{G}(k)$ and $\hat{\rho}^q\bf(k)$ are used to simplify the equations. Notice the former and the self interaction term need only be calculated once at the beginning of the program. The $V_0$ is needed to go from a conducting medium to a vacuum by the following expression:

$$ V^{qq}(\epsilon_s=\infty) = V^{qq}(\epsilon_s=0) - V^{qq}_0 $$

The supra-index qq is used to denote that we are dealing with charge-charge interactions.

In out implementation, the real space summation is reduced to the reference super-cell specified in the input, thus the summation over all lattice vectors $\mathbf{R}$ is obviated.

Using the standard force formula $f^{qq}_i = -\nabla _{\mathbf{r}_i}V^{qq}$ one finds:

$$ f = f_r + f_k + f_0 $$

$$ f_r = q_i \sum_{j \neq i} q_j \left( \frac{\text{erfc}(\frac{r_{ij}}{\sigma})}{r_{ij}}  +  \frac{2}{\sqrt{\pi}\sigma}e^{-(\frac{r_{ij}}{\sigma})^2}\right) \frac{\mathbf{r}_{ij}}{r^2_{ij}}  $$

$$ f_k = - \frac{q_i}{V} \sum_{\mathbf{k} \neq 0} \hat{G}(k) Im(\hat{\rho}^q\bf(-k) e^{-i\mathbf{k \cdot r_i}}) \mathbf{k}  $$

$$ f_0 = - \frac{4\pi q_i}{3V} \sum_j q_j\mathbf{r}_j  $$

where the $f_0$ term only applies in the case of a vacuum.

