# -----------------------------------------------------------------------------

cdef class Util:
    """Util."""

    @classmethod
    def createMatBSE(cls, Mat R: petsc4py.PETSc.Mat, Mat C: petsc4py.PETSc.Mat) -> petsc4py.PETSc.Mat:
        """
        Create a matrix that can be used to define a BSE type problem.

        Collective.

        Create a matrix that can be used to define a structured eigenvalue
        problem of type BSE (Bethe-Salpeter Equation).

        Parameters
        ----------
        R
            The matrix for the diagonal block (resonant).
        C
            The matrix for the off-diagonal block (coupling).

        Returns
        -------
        petsc4py.PETSc.Mat
            The matrix with the block form :math:`H = [ R\; C; {-C}^H\; {-R}^T ]`.
        """
        cdef Mat H = Mat()
        CHKERR( MatCreateBSE(R.mat, C.mat, &H.mat) )
        return H

# -----------------------------------------------------------------------------
