/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test RG interface functions.\n\n";

#include <slepcrg.h>

#define NPOINTS 10
#define NVERTEX 7

int main(int argc,char **argv)
{
  RG             rg;
  PetscInt       i,inside,nv;
  PetscBool      triv;
  PetscReal      re,im,a,b,c,d;
  PetscScalar    ar,ai,cr[NPOINTS],ci[NPOINTS],vr[NVERTEX],vi[NVERTEX],*pr,*pi;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(RGCreate(PETSC_COMM_WORLD,&rg));

  /* ellipse */
  PetscCall(RGSetType(rg,RGELLIPSE));
  PetscCall(RGIsTrivial(rg,&triv));
  PetscCheck(triv,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Region should be trivial before setting parameters");
  PetscCall(RGEllipseSetParameters(rg,1.1,2,0.1));
  PetscCall(RGSetFromOptions(rg));
  PetscCall(RGIsTrivial(rg,&triv));
  PetscCheck(!triv,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Region should be non-trivial after setting parameters");
  PetscCall(RGView(rg,NULL));
  PetscCall(RGViewFromOptions(rg,NULL,"-rg_ellipse_view"));
  re = 0.1; im = 0.3;
#if PetscDefined(USE_COMPLEX)
  ar = PetscCMPLX(re,im);
#else
  ar = re; ai = im;
#endif
  PetscCall(RGCheckInside(rg,1,&ar,&ai,&inside));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Point (%g,%g) is %s the region\n",(double)re,(double)im,(inside>=0)?"inside":"outside"));

  PetscCall(RGComputeBoundingBox(rg,&a,&b,&c,&d));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"The bounding box is [%g,%g]x[%g,%g]\n",(double)a,(double)b,(double)c,(double)d));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Contour points: "));
  PetscCall(RGComputeContour(rg,NPOINTS,cr,ci));
  for (i=0;i<NPOINTS;i++) {
#if PetscDefined(USE_COMPLEX)
    re = PetscRealPart(cr[i]);
    im = PetscImaginaryPart(cr[i]);
#else
    re = cr[i];
    im = ci[i];
#endif
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"(%.3g,%.3g) ",(double)re,(double)im));
  }
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));

  /* interval */
  PetscCall(RGSetType(rg,RGINTERVAL));
  PetscCall(RGIsTrivial(rg,&triv));
  PetscCheck(triv,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Region should be trivial before setting parameters");
  PetscCall(RGIntervalSetEndpoints(rg,-1,1,-0.1,0.1));
  PetscCall(RGSetFromOptions(rg));
  PetscCall(RGIsTrivial(rg,&triv));
  PetscCheck(!triv,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Region should be non-trivial after setting parameters");
  PetscCall(RGView(rg,NULL));
  PetscCall(RGViewFromOptions(rg,NULL,"-rg_interval_view"));
  re = 0.2; im = 0;
#if PetscDefined(USE_COMPLEX)
  ar = PetscCMPLX(re,im);
#else
  ar = re; ai = im;
#endif
  PetscCall(RGCheckInside(rg,1,&ar,&ai,&inside));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Point (%g,%g) is %s the region\n",(double)re,(double)im,(inside>=0)?"inside":"outside"));

  PetscCall(RGComputeBoundingBox(rg,&a,&b,&c,&d));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"The bounding box is [%g,%g]x[%g,%g]\n",(double)a,(double)b,(double)c,(double)d));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Contour points: "));
  PetscCall(RGComputeContour(rg,NPOINTS,cr,ci));
  for (i=0;i<NPOINTS;i++) {
#if PetscDefined(USE_COMPLEX)
    re = PetscRealPart(cr[i]);
    im = PetscImaginaryPart(cr[i]);
#else
    re = cr[i];
    im = ci[i];
#endif
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"(%.3g,%.3g) ",(double)re,(double)im));
  }
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));

  /* polygon */
#if PetscDefined(USE_COMPLEX)
  vr[0] = PetscCMPLX(0.0,2.0);
  vr[1] = PetscCMPLX(1.0,4.0);
  vr[2] = PetscCMPLX(2.0,5.0);
  vr[3] = PetscCMPLX(4.0,3.0);
  vr[4] = PetscCMPLX(5.0,4.0);
  vr[5] = PetscCMPLX(6.0,1.0);
  vr[6] = PetscCMPLX(2.0,0.0);
#else
  vr[0] = 0.0; vi[0] = 1.0;
  vr[1] = 0.0; vi[1] = -1.0;
  vr[2] = 0.6; vi[2] = -0.8;
  vr[3] = 1.0; vi[3] = -1.0;
  vr[4] = 2.0; vi[4] = 0.0;
  vr[5] = 1.0; vi[5] = 1.0;
  vr[6] = 0.6; vi[6] = 0.8;
#endif
  PetscCall(RGSetType(rg,RGPOLYGON));
  PetscCall(RGIsTrivial(rg,&triv));
  PetscCheck(triv,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Region should be trivial before setting parameters");
  PetscCall(RGPolygonSetVertices(rg,NVERTEX,vr,vi));
  PetscCall(RGSetFromOptions(rg));
  PetscCall(RGIsTrivial(rg,&triv));
  PetscCheck(!triv,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Region should be non-trivial after setting parameters");
  PetscCall(RGView(rg,NULL));
  PetscCall(RGViewFromOptions(rg,NULL,"-rg_polygon_view"));
  re = 5; im = 0.9;
#if PetscDefined(USE_COMPLEX)
  ar = PetscCMPLX(re,im);
#else
  ar = re; ai = im;
#endif
  PetscCall(RGCheckInside(rg,1,&ar,&ai,&inside));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Point (%g,%g) is %s the region\n",(double)re,(double)im,(inside>=0)?"inside":"outside"));

  PetscCall(RGComputeBoundingBox(rg,&a,&b,&c,&d));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"The bounding box is [%g,%g]x[%g,%g]\n",(double)a,(double)b,(double)c,(double)d));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Contour points: "));
  PetscCall(RGComputeContour(rg,NPOINTS,cr,ci));
  for (i=0;i<NPOINTS;i++) {
#if PetscDefined(USE_COMPLEX)
    re = PetscRealPart(cr[i]);
    im = PetscImaginaryPart(cr[i]);
#else
    re = cr[i];
    im = ci[i];
#endif
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"(%.3g,%.3g) ",(double)re,(double)im));
  }
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));

  /* check vertices */
  PetscCall(RGPolygonGetVertices(rg,&nv,&pr,&pi));
  PetscCheck(nv==NVERTEX,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Wrong number of vertices: %" PetscInt_FMT,nv);
  for (i=0;i<nv;i++) {
#if !PetscDefined(USE_COMPLEX)
    if (pr[i]!=vr[i] || pi[i]!=vi[i])
#else
    if (pr[i]!=vr[i])
#endif
       SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Vertex number %" PetscInt_FMT " does not match",i);
  }

  PetscCall(PetscFree(pr));
#if !PetscDefined(USE_COMPLEX)
  PetscCall(PetscFree(pi));
#endif
  PetscCall(RGDestroy(&rg));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      requires: !complex

   test:
      suffix: 1_complex
      requires: complex

   test:
      suffix: 2
      requires: !complex
      args: -rg_ellipse_view draw:tikz:ellipse.tikz -rg_interval_view draw:tikz:interval.tikz -rg_polygon_view draw:tikz:polygon.tikz
      filter: cat - ellipse.tikz interval.tikz polygon.tikz
      requires: !single

   test:
      suffix: 2_complex
      requires: complex !single
      args: -rg_ellipse_view draw:tikz:ellipse.tikz -rg_interval_view draw:tikz:interval.tikz -rg_polygon_view draw:tikz:polygon.tikz
      filter: cat - ellipse.tikz interval.tikz polygon.tikz

TEST*/
