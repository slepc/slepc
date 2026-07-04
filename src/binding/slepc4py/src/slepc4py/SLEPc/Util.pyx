# -----------------------------------------------------------------------------

cdef class Util:
    """
    Other utilities such as the creation of structured matrices.
    """

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
            The matrix with the block form :math:`H = [ R\; C; {-C}^*\; {-R}^T ]`.

        See Also
        --------
        slepc.MatCreateBSE
        """
        cdef Mat H = Mat()
        CHKERR( MatCreateBSE(R.mat, C.mat, &H.mat) )
        return H

    @classmethod
    def createMatHamiltonian(cls, Mat A: petsc4py.PETSc.Mat, Mat B: petsc4py.PETSc.Mat, Mat C: petsc4py.PETSc.Mat) -> petsc4py.PETSc.Mat:
        """
        Create matrix to be used for a structured Hamiltonian eigenproblem.

        Collective.

        Parameters
        ----------
        A
            The matrix for (0,0) block.
        B
            The matrix for (0,1) block, must be real symmetric or Hermitian.
        C
            The matrix for (1,0) block, must be real symmetric or Hermitian.

        Returns
        -------
        petsc4py.PETSc.Mat
            The matrix with the block form :math:`H = [ A\; B; C\; -A^* ]`.

        See Also
        --------
        slepc.MatCreateHamiltonian
        """
        cdef Mat H = Mat()
        CHKERR( MatCreateHamiltonian(A.mat, B.mat, C.mat, &H.mat) )
        return H

    @classmethod
    def createMatLREP(cls, Mat AK: petsc4py.PETSc.Mat, Mat BM: petsc4py.PETSc.Mat, red: bool = False) -> petsc4py.PETSc.Mat:
        """
        Create a matrix that can be used to define a LREP type problem.

        Collective.

        Create a matrix that can be used to define a structured Linear
        Response eigenvalue problem.

        Parameters
        ----------
        AK
            The matrix for the diagonal block or the top block.
        BM
            The matrix for the off-diagonal block or the bottom block.
        red
            Whether the reduced form should be built.

        Returns
        -------
        petsc4py.PETSc.Mat
            The matrix with the block form :math:`H = [ A\; B; -B\; -A ]`
            (non-reduced) or :math:`H = [ 0\; K; M\; 0 ]` (reduced).

        See Also
        --------
        slepc.MatCreateLREP
        """
        cdef Mat H = Mat()
        cdef PetscBool tval = asBool(red)
        CHKERR( MatCreateLREP(AK.mat, BM.mat, tval, &H.mat) )
        return H

# -----------------------------------------------------------------------------
