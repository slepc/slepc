# --------------------------------------------------------------------

cdef extern from "Python.h":
    enum: PY_SSIZE_T_MAX
    void *PyMem_Malloc(size_t)
    void *PyMem_Realloc(void*, size_t)
    void PyMem_Free(void*)

# @cython.final
# @cython.internal
cdef class _p_mem:
    cdef void *buf

    def __cinit__(self):
        self.buf = NULL

    def __dealloc__(self):
        PyMem_Free(self.buf)

cdef inline object allocate(size_t n, void **buf):
    cdef _p_mem ob = <_p_mem>_p_mem.__new__(_p_mem)
    ob.buf = PyMem_Malloc(<size_t>n)
    if ob.buf == NULL: raise MemoryError
    if buf != NULL: buf[0] = ob.buf
    return ob

# --------------------------------------------------------------------
