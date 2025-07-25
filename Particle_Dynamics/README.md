# Particle Dynamics
Currently only one short-range potential (Lennard-Jones) and one long-range potential (Coulomb) are implemented. The former supports three formalisms (MD, MC and BM) while the latter only two (MD and BM).

## Lennar-Jones potential
The Lennard-Jones potential is given by:

$$U^{LJ}(\mathbf{r}) = 4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right] $$

### Potential implementation
We consider a cut-off radius ($r_{cut}$) to reduce needless calculation of insignificant interactions (due to the short-range nature of the interaction). Thus, the potential is truncated and displaces, resulting in:

optimize_kSpace

$$
U(r) =
\begin{cases}
    U^{LJ}(r) - U^{LJ}(r_\text{cut}) & r \le r_\text{cut} \\
    0 & r > r_\text{cut}
\end{cases}
$$

### Forces
Using the standard force formula $\mathbf{f}_i = -\nabla _{\mathbf{r}_i}U$ one finds the force between two particles is given by:

$$ \mathbf{f}_{ij} = 48\epsilon \left[\left(\frac{\sigma}{r_{ij}}\right)^{12} - \frac{1}{2}\left(\frac{\sigma}{r_{ij}}\right)^6 \right] \frac{\mathbf{r}_{ij}}{r^2_{ij}} $$

### Pressure
It can be shown that the potential contribution to the pressure is:

$$P = \delta\cdot K_bT + \frac{1}{3V}\frac{1}{2} \left( \sum_i^N \sum_j^N \mathbf{f}_{ij} \cdot \mathbf{r}_{ij} \right)$$

Notice there is no correction due to the displacement of the potential. The truncation, however, is implicitly considered in the summation, as no forces between particles further than $r_{cut}$ is even calculated.

## Coulomb interactions
The Coulomb potential between two charged particles ($q_i$, $q_j$)is given by:

$$U(\mathbf{r}_{ij}) = \frac{\epsilon q_i q_j}{r_{ij}}$$

### Potential implementation

We consider the summation of the particles in the reference super-cell to get the real space contribution (short-range contribution) of the potential and forces. The long-range term is taken into account by summing over a ball of k-vectors in reciprocal space. Thus, setting $\epsilon=4\pi\epsilon_0$, we get:
$$
U = U_r + U_k - U_s - U_q - U_0
$$
$$ U_r = \epsilon \frac{1}{2} \sum_i^N \sum_j^N q_i q_j \left( \sum_{\mathbf{R}} \frac{\text{erfc}(\frac{|\mathbf{r}_{ij}+ \mathbf{R}|}{\sigma})}{|\mathbf{r}_{ij}+ \mathbf{R}|}\right)  $$

$$ U_k = \epsilon  \frac{1}{2V} \sum_{\mathbf{k} \neq 0} G(k) \rho^q\bf(k) \rho^q\bf(-k) $$

$$ U_s = \epsilon \sum_{i} \frac{q_i^2}{\sqrt{\pi}\sigma} $$

$$ U_{q} = \epsilon \frac{\pi\sigma^2}{2V} |\sum_i q_i|^2 $$

$$ U_0 = \epsilon \frac{2\pi}{3V} |\sum_i q_i\mathbf{r}_i |^2 $$

$$ G(k) = 4\pi \frac{e^{-(\frac{k\sigma}{2})^2}}{k^2} $$

$$ \rho^q(\mathbf{k}) = \sum_i q_i e^{-i\mathbf{k \cdot r_i}}  $$

where $U_r$ is calculated in real space, $U_k$ in reciprocal space, $U_s$ is the self interaction term, $U_q$ is a term to neutralize the system charge (necessary for Ewald sums) and $U_0$ is a dipole term for charges in a vacuum. $G(k)$ and $\rho^q\bf(k)$ are used to simplify the equations. Notice the former and the self interaction term need only be calculated once at the beginning of the program. The $U_0$ is needed to go from a conducting medium to a vacuum by the following expression:

$$ U(\epsilon_s=\infty) = U(\epsilon_s=0) - U_0 $$

In out implementation, the real space summation is reduced to the reference super-cell specified in the input, thus the summation over all lattice vectors $\mathbf{R}$ is obviated.

The potential energy difference due to the translation $\mathbf{r} \rightarrow\widetilde{\mathbf{r}}$ is simply $\Delta U = \Delta U_r + \Delta U_k$, each given by:

$$\Delta U_r = \epsilon \frac{q_i}{2} \sum_j^N q_j \left(\frac{\text{erfc}(\frac{\mathbf{\widetilde{r}}_{ij}}{\sigma})}{\widetilde{r}_{ij}} - \frac{\text{erfc}(\frac{\mathbf{r}_{ij}}{\sigma})}{r_{ij}} \right)  $$


$$\Delta  U_k = \epsilon \frac{1}{2V} \sum_{\mathbf{k} \neq 0} G(k) \left[2\mathcal{Re}(\rho^q (\mathbf{k}) \Delta\rho^q(\mathbf{-k})) + |\Delta\rho^q(\mathbf{-k})|^2 \right]$$

$$ \Delta\rho^q(\mathbf{k}) = q_i(e^{-i\mathbf{k \cdot \widetilde{r}_i}} - e^{-i\mathbf{k \cdot r_i}}) $$

Notice that for the MC run, only the variation of the reciprocal charge ($\Delta\rho^q(\mathbf{k})$) needs to be calculated, thus we store the full reciprocal charge ($\rho^q(\mathbf{k})$) and update it only if the trial is accepted. That is, once a trial is accepted, the full reciprocal charge is recalculated. If we simply attempt to add the variation of the reciprocal charge to the full reciprocal charge, a huge error is accumulated throughout the run. Even a recalculation every MC step is not satisfactory enough (though it improves the results). Thus, the full reciprocal charge is recalculated from scratch after $\textit{every accepted MC trial}$.


### Forces

Using the standard force formula $ \mathbf{f}_i = -\nabla _{\mathbf{r}_i}U $ one finds:

$$ \mathbf{f} = \mathbf{f}_r + \mathbf{f}_k + \mathbf{f}_0 $$

$$ \mathbf{f}_r = q_i \sum_{j \neq i} q_j \left( \frac{\text{erfc}(\frac{r_{ij}}{\sigma})}{r_{ij}}  +  \frac{2}{\sqrt{\pi}\sigma}e^{-(\frac{r_{ij}}{\sigma})^2}\right) \frac{\mathbf{r}_{ij}}{r^2_{ij}}  $$

$$ \mathbf{f}_k = - \frac{q_i}{V} \sum_{\mathbf{k} \neq 0} G(k) \mathcal{Im} (\rho^q\bf(-k) e^{-i\mathbf{k \cdot r_i}}) \mathbf{k}  $$

$$ \mathbf{f}_0 = - \frac{4\pi q_i}{3V} \sum_j q_j\mathbf{r}_j  $$

where the $\mathbf{f}_0$ term only applies in the case of a vacuum.

### Pressure
It can be shown that the potential contribution to the pressure is simply $\frac{U}{3}$. Considering the contribution due to the reference temperature, we get:
$$P = \delta\cdot K_bT_0 + \frac{1}{3V}\frac{U}{3}$$
