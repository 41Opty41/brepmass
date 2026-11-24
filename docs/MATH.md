# Mathematical Formulations in `brepmass`

## 1. Closed Solid Volume & Moments (Divergence Theorem)

For a watertight oriented triangle mesh with outward face winding and triangles $(\mathbf{a}_k, \mathbf{b}_k, \mathbf{c}_k)$, decompose each triangle into a signed tetrahedron with the origin $(0, \mathbf{a}, \mathbf{b}, \mathbf{c})$:

$$v_k = \frac{1}{6}\, \mathbf{a} \cdot (\mathbf{b} \times \mathbf{c})$$

Total signed volume $V = \sum_k v_k$. Let $\mathbf{S} = \mathbf{a} + \mathbf{b} + \mathbf{c}$. Then:

- **First moment** $\int \mathbf{x}\, dV$:
  $$\mathbf{M}_1 = \sum_k \frac{v_k}{4} \mathbf{S}, \qquad \mathbf{c} = \frac{\mathbf{M}_1}{V}$$
- **Second moment (covariance) tensor** $C_{ij} = \int \rho\, x_i x_j\, dV$:
  $$C_{ij} = \rho \sum_k \frac{v_k}{20}\left(a_i a_j + b_i b_j + c_i c_j + S_i S_j\right)$$

## 2. Inertia Tensors & Principal Axes

Given mass $M = \rho V$, centroid $\mathbf{c}$, and mass covariance tensor $C$:

- **Inertia about the origin**:
  $$I_0 = \operatorname{tr}(C)\, I_3 - C$$
- **Inertia about the centroid** (Parallel Axis Theorem):
  $$I_{\text{cm}} = I_0 - M \left(\|\mathbf{c}\|^2 I_3 - \mathbf{c} \otimes \mathbf{c}\right)$$
- **Principal moments & axes**: Computed via cyclic Jacobi rotations on the real symmetric $3 \times 3$ matrix $I_{\text{cm}}$, sorted in ascending order $I_1 \le I_2 \le I_3$.

## 3. Open Sheet Integrals (`dim = 2`)

When half-edge pairing indicates an open boundary, `brepmass` evaluates surface-area integrals per triangle with area $A_k = \frac{1}{2}\|(\mathbf{b}-\mathbf{a})\times(\mathbf{c}-\mathbf{a})\|$:

$$\int \mathbf{x}\, dA = \sum_k \frac{A_k}{3}\mathbf{S}, \qquad \int x_i x_j\, dA = \sum_k \frac{A_k}{12}\left(a_i a_j + b_i b_j + c_i c_j + S_i S_j\right)$$

## 4. Conical Frustum Reference

For $z \in [0, h]$ with bottom radius $r_1$ and top radius $r_2$:

- $V = \frac{\pi h}{3}(r_1^2 + r_1 r_2 + r_2^2)$
- $z_{\text{cm}} = \frac{h}{4}\frac{r_1^2 + 2r_1 r_2 + 3r_2^2}{r_1^2 + r_1 r_2 + r_2^2}$
- $I_{zz} = \rho\frac{\pi h}{10}(r_1^4 + r_1^3 r_2 + r_1^2 r_2^2 + r_1 r_2^3 + r_2^4)$ (singularity-free when $r_1 = r_2$)
- $I_{xx,0} = I_{yy,0} = \frac{1}{2}I_{zz} + \rho\pi h^3 \left(\frac{r_1^2}{30} + \frac{r_1 r_2}{10} + \frac{r_2^2}{5}\right)$
