# SLEPc

## Scalable Library for Eigenvalue Problem Computation

This is the home page of **SLEPc, the Scalable Library for Eigenvalue Problem Computations**. SLEPc is a software library for the solution of large scale sparse eigenvalue problems on parallel computers. It is an extension of [PETSc](https://petsc.org) and can be used for _linear eigenvalue problems_ in either standard or generalized form, with real or complex arithmetic. It can also be used for computing a _partial SVD_ of a large, sparse, rectangular matrix, and to solve _nonlinear eigenvalue problems_ (polynomial or general).  Additionally, SLEPc provides solvers for the computation of the action of a _matrix function_ on a vector.

The current version of SLEPc is {{env.config.version}}, released in {{release_date}}.

SLEPc is based on the PETSc data structures and it employs the [MPI standard](https://www.mpi-forum.org) for message-passing communication. It is being developed by researchers from Universitat Politecnica de Valencia (Spain).

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

<div class="hideLink"><a style="float: right" onclick="batch=4;expandTable(newsitems,batch);newsitems+=batch;">[Show older news]</a></div>

***
<br/>

```{image} https://www.upv.es/imagenes/svg/logo-upv.svg
:alt: logo-upv
:align: left
:width: 340px
:target: https://www.upv.es/
```

```{image} _static/images/logo-dsic.svg
:alt: logo-dsic
:align: right
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


