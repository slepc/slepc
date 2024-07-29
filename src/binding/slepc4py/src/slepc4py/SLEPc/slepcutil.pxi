cdef extern from * nogil:
    PetscErrorCode MatCreateBSE(PetscMat,PetscMat,PetscMat*)
    PetscErrorCode MatCreateHamiltonian(PetscMat,PetscMat,PetscMat,PetscMat*)

