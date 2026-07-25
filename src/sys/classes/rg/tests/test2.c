/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test the ring region.\n\n";

#include <slepcrg.h>

#define NPOINTS 11

PetscErrorCode CheckPoint(RG rg,PetscReal re,PetscReal im)
{
  PetscInt       inside;
  PetscScalar    ar,ai;

  PetscFunctionBeginUser;
#if PetscDefined(USE_COMPLEX)
  ar = PetscCMPLX(re,im);
#else
  ar = re; ai = im;
#endif
  PetscCall(RGCheckInside(rg,1,&ar,&ai,&inside));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Point (%g,%g) is %s the region\n",(double)re,(double)im,(inside>=0)?"inside":"outside"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  RG             rg;
  RGType         rtype;
  PetscInt       i;
  PetscBool      triv;
  PetscReal      re,im,radius,vscale,start_ang,end_ang,width,a,b,c,d;
  PetscScalar    center,cr[NPOINTS],ci[NPOINTS];

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(RGCreate(PETSC_COMM_WORLD,&rg));

  PetscCall(RGSetType(rg,RGRING));
  PetscCall(RGIsTrivial(rg,&triv));
  PetscCheck(triv,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Region should be trivial before setting parameters");
  PetscCall(RGRingSetParameters(rg,2,PETSC_DETERMINE,0.5,0.25,0.75,0.1));
  PetscCall(RGSetFromOptions(rg));
  PetscCall(RGIsTrivial(rg,&triv));
  PetscCheck(!triv,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Region should be non-trivial after setting parameters");
  PetscCall(RGView(rg,NULL));
  PetscCall(RGViewFromOptions(rg,NULL,"-rg_view"));

  PetscCall(RGGetType(rg,&rtype));
  PetscCall(RGRingGetParameters(rg,&center,&radius,&vscale,&start_ang,&end_ang,&width));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"%s region: \n  center=%g, radius=%g, vscale=%g\n  start angle=%g, end angle=%g, width=%g\n\n",rtype,(double)PetscRealPart(center),(double)radius,(double)vscale,(double)start_ang,(double)end_ang,(double)width));

  PetscCall(CheckPoint(rg,3.0,0.3));
  PetscCall(CheckPoint(rg,1.1747,0.28253));

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

  PetscCall(RGDestroy(&rg));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      args: -rg_ring_width 0.015

   test:
      suffix: 2
      args: -rg_ring_width 0.015 -rg_scale 1.5

   test:
      suffix: 3
      args: -rg_view draw:tikz:test2_3_ring.tikz
      filter: cat - test2_3_ring.tikz
      requires: !single

   test:
      suffix: 4
      args: -rg_ring_width 0.015 -rg_ring_center 3 -rg_ring_radius 0.3 -rg_ring_vscale 1
      requires: !single

   test:
      suffix: 5
      args: -rg_ring_width 0.1 -rg_ring_center 0.35 -rg_ring_radius 0.86 -rg_ring_vscale 1 -rg_ring_startangle 0.75 -rg_ring_endangle 0.25

TEST*/
