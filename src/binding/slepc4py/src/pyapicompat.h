/* Author:  Lisandro Dalcin   */
/* Contact: dalcinl@gmail.com */

/* -------------------------------------------------------------------------- */

#if defined(Py_LIMITED_API) && Py_LIMITED_API+0 < 0x030B0000

#define Py_bf_getbuffer 1
#define Py_bf_releasebuffer 2

typedef struct {
  void *buf;
  PyObject *obj;
  Py_ssize_t len;
  Py_ssize_t itemsize;
  int readonly;
  int ndim;
  char *format;
  Py_ssize_t *shape;
  Py_ssize_t *strides;
  Py_ssize_t *suboffsets;
  void *internal;
} Py_buffer;

#define PyBUF_SIMPLE 0
#define PyBUF_WRITABLE 0x0001

#define PyBUF_FORMAT 0x0004
#define PyBUF_ND 0x0008
#define PyBUF_STRIDES (0x0010 | PyBUF_ND)
#define PyBUF_C_CONTIGUOUS (0x0020 | PyBUF_STRIDES)
#define PyBUF_F_CONTIGUOUS (0x0040 | PyBUF_STRIDES)
#define PyBUF_ANY_CONTIGUOUS (0x0080 | PyBUF_STRIDES)
#define PyBUF_INDIRECT (0x0100 | PyBUF_STRIDES)

#define PyBUF_CONTIG (PyBUF_ND | PyBUF_WRITABLE)
#define PyBUF_CONTIG_RO (PyBUF_ND)

#define PyBUF_STRIDED (PyBUF_STRIDES | PyBUF_WRITABLE)
#define PyBUF_STRIDED_RO (PyBUF_STRIDES)

#define PyBUF_RECORDS (PyBUF_STRIDES | PyBUF_WRITABLE | PyBUF_FORMAT)
#define PyBUF_RECORDS_RO (PyBUF_STRIDES | PyBUF_FORMAT)

#define PyBUF_FULL (PyBUF_INDIRECT | PyBUF_WRITABLE | PyBUF_FORMAT)
#define PyBUF_FULL_RO (PyBUF_INDIRECT | PyBUF_FORMAT)

#define PyBUF_READ  0x100
#define PyBUF_WRITE 0x200

PyAPI_FUNC(int)  PyObject_CheckBuffer(PyObject *);
PyAPI_FUNC(int)  PyObject_GetBuffer(PyObject *, Py_buffer *, int);
PyAPI_FUNC(void) PyBuffer_Release(Py_buffer *);
PyAPI_FUNC(int)  PyBuffer_FillInfo(Py_buffer *, PyObject *,
                                   void *, Py_ssize_t, int, int);

#endif

/* -------------------------------------------------------------------------- */

#if defined(Py_LIMITED_API) && Py_LIMITED_API+0 < 0x030D0000

/* https://github.com/cython/cython/pull/6914 */

#define PyDict_GetItemStringRef PyDict_GetItemStringRef_313
static inline int PyDict_GetItemStringRef(PyObject *d,
                                          const char *key,
                                          PyObject **result)
{
  PyObject *key_obj = PyUnicode_FromString(key);
  if (key_obj == NULL) return -1;
  *result = PyDict_GetItemWithError(d, key_obj);
  Py_DECREF(key_obj);
  if (*result == NULL) return PyErr_Occurred() ? -1 : 0;
  Py_INCREF(*result);
  return 1;
}

#endif

/* -------------------------------------------------------------------------- */
