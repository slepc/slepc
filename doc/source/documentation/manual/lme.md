(ch:lme)=
# LME: Linear Matrix Equation

The Linear Matrix Equation (`LME`) solver object encapsulates functionality intended for the solution of linear matrix equations such as Lyapunov, Sylvester, and their various generalizations. A matrix equation is an equation where the unknown and the coefficients are all matrices, and it is linear if there are no nonlinear terms involving the unknown. This excludes some types of equations such as the Algebraic Riccati Equation or other quadratic matrix equations, which are not considered here.

The general form of a linear matrix equation is
```{math}
:label: eq:lme

AXE+DXB=C,
```
where $X$ is the solution.

Since SLEPc is mainly concerned with iterative methods, so is the `LME` module. This implies that `LME` can address only the case where the solution $X$ has low rank. In many situations, $X$ does not have low rank, which means that the `LME` solver, if it converges, will truncate the solution to an artificially small rank, with a large approximation error.

:::{warning}
The `LME` module should be considered experimental. As explained below, the currently implemented functionality is quite limited, and may be extended in the future.
:::

## Current Functionality

The user interface of the `LME` module is prepared for different types of equations, see a list in table [](#tab:eqtype). However, currently there is only basic support for continuous-time Lyapunov equations, and the rest are just a wish list.

:::{table} Problem types considered in `LME`
:name: tab:eqtype

  Problem Type             | Equation         |`LMEProblemType`   | $A$ | $B$ | $D$ | $E$
  -------------------------|------------------|-------------------|-----|-----|-----|-----
  Continuous-Time Lyapunov | $AX+XA^*=-C$     |`LME_LYAPUNOV`     | yes |$A^*$|  -  |  -
  Sylvester                | $AX+XB=C$        |`LME_SYLVESTER`    | yes | yes |  -  |  -
  Generalized Lyapunov     | $AXD^*+DXA^*=-C$ |`LME_GEN_LYAPUNOV` | yes |$A^*$| yes |$D^*$
  Generalized Sylvester    | $AXE+DXB=C$      |`LME_GEN_SYLVESTER`| yes | yes | yes | yes
  Discrete-Time Lyapunov   | $AXA^*-X=-C$     |`LME_DT_LYAPUNOV`  | yes |  -  |  -  |$A^*$
  Stein                    | $AXE-X=-C$       |`LME_STEIN`        | yes |  -  |  -  | yes

:::

### Continuous-Time Lyapunov Equation

Given two matrices $A,C\in\mathbb{C}^{n\times n}$, with $C$ Hermitian positive-definite, the continuous-time Lyapunov equation, assuming that the right-hand side matrix $C$ has low rank, can be expressed as
```{math}
:label: eq:clyap

AX+XA^*=-C_1C_1^*,
```
where we have written the right-hand side in factored form, $C=C_1C_1^*$, $C_1\in\mathbb{C}^{n\times r}$ with $r = \operatorname{rank}(C)$. In general, $X$ will not have low rank, but our intention is to compute a low-rank approximation $\tilde X=X_1X_1^*$, $X_1\in\mathbb{C}^{n\times p}$ for a given value $p$. Note that  the equation has a unique solution if and only if $A$ is stable, i.e., all of its eigenvalues have strictly negative real part.

The reason why SLEPc provides a solver for equation {math:numref}`eq:clyap` is that it is required in the eigensolver `EPSLYAPII`, which implements the Lyapunov inverse iteration {cite:p}`Mee10,Elm13`.

:::{note}
The `LME` solvers currently implemented in SLEPc are very basic. Users interested in applying `LME` to matrix equations appearing in their applications are encouraged to contact SLEPc developers describing their use case.
:::

### Krylov Solver

Currently, the solver for the continuous-time Lyapunov equation available in SLEPc is based on the following procedure for each column $c$ of $C_1$:

  1. Build an _Arnoldi factorization_ for the Krylov subspace $\mathcal{K}_m(A,c)$,
     $AV_m=V_mH_m+h_{m+1,m}v_{m+1}e_m^*$.
  2. Solve the _compressed Lyapunov equation_,
     $H_mY+YH_m^*=-\tilde c\,\tilde c^* \quad\mathrm{with}\quad\tilde c=V_m^*c$.
  3. Set the approximate solution to $X_m=V_mYV_m^*$.

Furthermore, our Krylov solver incorporates an Eiermann-Ernst-type restart as proposed by {cite:t}`Kre08`. The restart will discard the Arnoldi basis $V_m$ but keep $H_m$, then continue the Arnoldi recurrence from $v_{m+1}$. The upper Hessenberg matrices generated at each restart are glued together. This implies that at each restart the size of the compressed Lyapunov equation grows. After convergence, once the full $Y$ is available, we perform a second pass to reconstruct the successive $V_m$ bases.

The restarted algorithm is as follows:

>  **for** $k=1,2,\ldots$
>
>  - Run Arnoldi $AV_m^{(k)}=V_m^{(k)}H_m^{(k)}+h_{m+1,m}^{(k)}v_{m+1}^{(k)}e_m^*$.
>  - Set $H_{km}=\left[\begin{array}{cc}H_{(k-1)m}&0\\h_{m+1,m}^{(k-1)}e_1e_{(k-1)m}^*&H_m^{(k)}\end{array}\right]$.
>  - Solve compressed equation $H_{km}Y+YH_{km}^*=-\|c\|_2^2\,e_1e_1^*$.
>  - Compute residual norm $\|R_k\|_2=h_{m+1,m}^{(k)}\|e_{km}^*Y\|_2$.
>
>  **end**

  The residual norm is used for the stopping criterion.

### Projected Problem

As mentioned in the previous section, the projected problem is a smaller Lyapunov equation, which must be solved in each outer iteration of the method (restart). The solution of this small dense matrix equation is also implemented within the `LME` module (there is no auxiliary class for this).

If $A$ is stable, then $X$ is symmetric positive semidefinite. Hence, looking at step 3 of the algorithm in the previous section, it is better to compute the Cholesky factor $Y=LL^*$ to obtain the low rank factor $X_1=V_mL$. The method of Hammarling will compute $L$ directly, without computing $Y$ first.

If SLEPc has been configured with the option `--download-slicot` (see [](#sec:wrap)), it will be possible to use SLICOT subroutines to apply Hammarling's method, otherwise a more basic algorithm will be used.

:::{note}
Most of the SLICOT subroutines are implemented for real arithmetic only, so it is not possible to use them in a complex build of PETSc/SLEPc.
:::

In the restarted algorithm described above, in the second pass to recalculate the successive $V_m$'s, the update is based on the truncated SVD of the Cholesky factor
  $$L=\begin{bmatrix}Q_1 & Q_2\end{bmatrix}\begin{bmatrix}\Sigma_1 & 0 \\ 0 & \Sigma_2\end{bmatrix}\begin{bmatrix}Z_1 & Z_2\end{bmatrix}^*\quad\mathrm{with}\quad\|\Sigma_2\|_2<\varepsilon,$$
so that the final $X_1$ has the appropriate rank.

## Basic Usage

Once the solver object has been created with `LMECreate()`, the user has to define the matrix equation to be solved. First, the type of equation must be set with `LMESetProblemType()`, choosing one from table [](#tab:eqtype).

The coefficient matrices $A$, $B$, $D$, $E$ must be provided via `LMESetCoefficients()`, but some of them are optional depending on the matrix equation. For Lyapunov equations, only $A$ must be set, which is normally a large, sparse matrix stored in {external:doc}`MATAIJ` format. Note that in table [](#tab:eqtype), the notation $A^*$ means that this matrix need not be passed, but the user may choose to pass an explicit transpose of matrix $A$ (for improved efficiency). Also note that some of the equation types impose restrictions on the properties of the coefficient matrices (e.g., $A$ stable in Lyapunov equations) and possibly on the right-hand side $C$.

To conclude the definition of the equation, we must pass the right-hand side $C$ with `LMESetRHS()`. Note that $C$ is expressed as the outer product of a tall-skinny matrix, $C=C_1C_1^*$. This can be represented in PETSc using a special type of matrix, {external:doc}`MATLRC`, that can be created with {external:doc}`MatCreateLRC` passing a dense matrix $C_1$.

The call to `LMESolve()` will run the solver to compute the solution $X$. The rank of $X$ may be prescribed by the user or selected dynamically by the solver. Next, we discuss these two options:

  * To prescribe a fixed rank for $X$, we must preallocate it, creating a {external:doc}`MATLRC` matrix, and then pass it via `LMESetSolution()` _before_ calling `LMESolve()`. In this way, the solver will be restricted to a rank equal to the number of columns provided in $X_1$.

  * If we do not call `LMESetSolution()` beforehand, the solver will select the rank dynamically. Then after `LMESolve()` we must call `LMEGetSolution()` to retrieve a {external:doc}`MATLRC` matrix representing $X$.

```{only} html
<p class="rubric">References</p>
```
```{bibliography}
:filter: docname in docnames
```
