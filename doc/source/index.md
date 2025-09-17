# SLEPc

## Scalable Library for Eigenvalue Problem Computations

SLEPc is a software library for the solution of large scale eigenvalue problems on parallel computers. It is an extension of [PETSc](https://petsc.org) and can be used for several related linear algebra problems such as:

- _Linear eigenvalue problems_ in either standard or generalized form, with real or complex arithmetic.
- _Partial SVD_ of a large, sparse, rectangular matrix, as well as other types of singular value problems.
- _Nonlinear eigenvalue problems_, either polynomial or general.
- The computation of the action of a _matrix function_ on a vector.

SLEPc relies on PETSc data structures and linear solvers. The model of parallelism employs the [MPI standard](https://www.mpi-forum.org) for message-passing communication, with additional support for GPU as is done in PETSc. SLEPc has bindings for C, Fortran, and Python (via slepc4py). It is being developed mainly by researchers from Universitat Polit&egrave;cnica de Val&egrave;ncia (Spain).

The current version of SLEPc is {{env.config.version}}, released in {{release_date}}.

## Latest News and Forthcoming Events

<script src="_static/js/news.js"></script>

<div class="pst-scrollable-table-container">
<table class="table" id="newstable">
<thead></thead>
<tbody>
<script>
var newsitems = 5;
for (var i=0;i<newsitems;i++) {
    if (i%2) {
        document.write("<tr class='row-even'>");
    } else {
        document.write("<tr class='row-odd'>");
    }
    document.write("<td style='width:15%'><p><strong>"+news[i][0]+"</strong></p></td>");
    document.write("<td style='width:85%'><p>"+news[i][1]+"</p></td>");
    document.write("</tr>");
}
</script>
</tbody>
</table>
</div>

<div style="float: right"><a class="reference internal nav-link" onclick="batch=4;expandTable(newsitems,batch);newsitems+=batch;">[Show older news]</a></div>

***
<br/>

```{image} https://www.upv.es/imagenes/svg/logo-upv.svg
:alt: logo-upv
:align: left
:width: 300px
:target: https://www.upv.es/
```

```{image} _static/images/logo-dsic.svg
:alt: logo-dsic
:align: right
:width: 180px
:target: https://www.dsic.upv.es/
```

***

```{toctree}
:maxdepth: 1
:hidden:

about/index
installation/index
documentation/index
manualpages/index
slepc4py/index
material/index
contact/index
```


