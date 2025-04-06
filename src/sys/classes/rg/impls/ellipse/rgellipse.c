/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Region enclosed in an ellipse (aligned with the coordinate axes)
*/

#include <slepc/private/rgimpl.h>      /*I "slepcrg.h" I*/
#include <petscdraw.h>

typedef struct {
  PetscScalar center;     /* center of the ellipse */
  PetscReal   radius;     /* radius of the ellipse */
  PetscReal   vscale;     /* vertical scale of the ellipse */
} RG_ELLIPSE;

static PetscErrorCode RGEllipseSetParameters_Ellipse(RG rg,PetscScalar center,PetscReal radius,PetscReal vscale)
{
  RG_ELLIPSE *ctx = (RG_ELLIPSE*)rg->data;

  PetscFunctionBegin;
  ctx->center = center;
  if (radius == (PetscReal)PETSC_DETERMINE) {
    ctx->radius = 1.0;
  } else if (radius != (PetscReal)PETSC_CURRENT) {
    PetscCheck(radius>0.0,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The radius argument must be > 0.0");
    ctx->radius = radius;
  }
  if (vscale == (PetscReal)PETSC_DETERMINE) {
    ctx->vscale = 1.0;
  } else if (vscale != (PetscReal)PETSC_CURRENT) {
    PetscCheck(vscale>0.0,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The vscale argument must be > 0.0");
    ctx->vscale = vscale;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGEllipseSetParameters - Sets the parameters defining the ellipse region.

   Logically Collective

   Input Parameters:
+  rg     - the region context
.  center - center of the ellipse
.  radius - radius of the ellipse
-  vscale - vertical scale of the ellipse

   Options Database Keys:
+  -rg_ellipse_center - Sets the center
.  -rg_ellipse_radius - Sets the radius
-  -rg_ellipse_vscale - Sets the vertical scale

   Notes:
   In the case of complex scalars, a complex center can be provided in the
   command line with [+/-][realnumber][+/-]realnumberi with no spaces, e.g.
   -rg_ellipse_center 1.0+2.0i

   When PETSc is built with real scalars, the center is restricted to a real value.

   For radius and vscale, you can use PETSC_CURRENT to keep the current value, and
   PETSC_DETERMINE to set them to a default of 1.

   Level: advanced

.seealso: RGEllipseGetParameters()
@*/
PetscErrorCode RGEllipseSetParameters(RG rg,PetscScalar center,PetscReal radius,PetscReal vscale)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidLogicalCollectiveScalar(rg,center,2);
  PetscValidLogicalCollectiveReal(rg,radius,3);
  PetscValidLogicalCollectiveReal(rg,vscale,4);
  PetscTryMethod(rg,"RGEllipseSetParameters_C",(RG,PetscScalar,PetscReal,PetscReal),(rg,center,radius,vscale));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGEllipseGetParameters_Ellipse(RG rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale)
{
  RG_ELLIPSE *ctx = (RG_ELLIPSE*)rg->data;

  PetscFunctionBegin;
  if (center) *center = ctx->center;
  if (radius) *radius = ctx->radius;
  if (vscale) *vscale = ctx->vscale;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGEllipseGetParameters - Gets the parameters that define the ellipse region.

   Not Collective

   Input Parameter:
.  rg     - the region context

   Output Parameters:
+  center - center of the region
.  radius - radius of the region
-  vscale - vertical scale of the region

   Level: advanced

.seealso: RGEllipseSetParameters()
@*/
PetscErrorCode RGEllipseGetParameters(RG rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscUseMethod(rg,"RGEllipseGetParameters_C",(RG,PetscScalar*,PetscReal*,PetscReal*),(rg,center,radius,vscale));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGView_Ellipse(RG rg,PetscViewer viewer)
{
  RG_ELLIPSE     *ctx = (RG_ELLIPSE*)rg->data;
  PetscBool      isdraw,isascii;
  int            winw,winh;
  PetscDraw      draw;
  PetscDrawAxis  axis;
  PetscReal      cx,cy,r,ab,cd,lx,ly,w,scale=1.2;
  char           str[50];

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw));
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->center,PETSC_FALSE));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  center: %s, radius: %g, vscale: %g\n",str,RGShowReal(ctx->radius),RGShowReal(ctx->vscale)));
  } else if (isdraw) {
    PetscCall(PetscViewerDrawGetDraw(viewer,0,&draw));
    PetscCall(PetscDrawCheckResizedWindow(draw));
    PetscCall(PetscDrawGetWindowSize(draw,&winw,&winh));
    winw = PetscMax(winw,1); winh = PetscMax(winh,1);
    PetscCall(PetscDrawClear(draw));
    PetscCall(PetscDrawSetTitle(draw,"Ellipse region"));
    PetscCall(PetscDrawAxisCreate(draw,&axis));
    cx = PetscRealPart(ctx->center)*rg->sfactor;
    cy = PetscImaginaryPart(ctx->center)*rg->sfactor;
    r  = ctx->radius*rg->sfactor;
    lx = 2*r;
    ly = 2*r*ctx->vscale;
    ab = cx;
    cd = cy;
    w  = scale*PetscMax(lx/winw,ly/winh)/2;
    PetscCall(PetscDrawAxisSetLimits(axis,ab-w*winw,ab+w*winw,cd-w*winh,cd+w*winh));
    PetscCall(PetscDrawAxisDraw(axis));
    PetscCall(PetscDrawAxisDestroy(&axis));
    PetscCall(PetscDrawEllipse(draw,cx,cy,2*r,2*r*ctx->vscale,PETSC_DRAW_RED));
    PetscCall(PetscDrawFlush(draw));
    PetscCall(PetscDrawSave(draw));
    PetscCall(PetscDrawPause(draw));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGIsTrivial_Ellipse(RG rg,PetscBool *trivial)
{
  RG_ELLIPSE *ctx = (RG_ELLIPSE*)rg->data;

  PetscFunctionBegin;
  if (rg->complement) *trivial = PetscNot(ctx->radius);
  else *trivial = PetscNot(ctx->radius<PETSC_MAX_REAL);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGComputeContour_Ellipse(RG rg,PetscInt n,PetscScalar *cr,PetscScalar *ci)
{
  RG_ELLIPSE *ctx = (RG_ELLIPSE*)rg->data;
  PetscReal  theta;
  PetscInt   i;

  PetscFunctionBegin;
  for (i=0;i<n;i++) {
    theta = 2.0*PETSC_PI*(i+0.5)/n;
#if defined(PETSC_USE_COMPLEX)
    cr[i] = ctx->center + ctx->radius*PetscCMPLX(PetscCosReal(theta),ctx->vscale*PetscSinReal(theta));
#else
    if (cr) cr[i] = ctx->center + ctx->radius*PetscCosReal(theta);
    if (ci) ci[i] = ctx->radius*ctx->vscale*PetscSinReal(theta);
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGComputeBoundingBox_Ellipse(RG rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d)
{
  RG_ELLIPSE *ctx = (RG_ELLIPSE*)rg->data;

  PetscFunctionBegin;
  if (a) *a = PetscRealPart(ctx->center) - ctx->radius;
  if (b) *b = PetscRealPart(ctx->center) + ctx->radius;
  if (c) *c = PetscImaginaryPart(ctx->center) - ctx->radius*ctx->vscale;
  if (d) *d = PetscImaginaryPart(ctx->center) + ctx->radius*ctx->vscale;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGComputeQuadrature_Ellipse(RG rg,RGQuadRule quad,PetscInt n,PetscScalar *z,PetscScalar *zn,PetscScalar *w)
{
  RG_ELLIPSE *ctx = (RG_ELLIPSE*)rg->data;
  PetscReal  theta;
  PetscInt   i;

  PetscFunctionBegin;
  for (i=0;i<n;i++) {
#if defined(PETSC_USE_COMPLEX)
    theta = 2.0*PETSC_PI*(i+0.5)/n;
    zn[i] = PetscCMPLX(PetscCosReal(theta),ctx->vscale*PetscSinReal(theta));
    w[i]  = (PetscCMPLX(ctx->vscale*PetscCosReal(theta),PetscSinReal(theta)))*ctx->radius*rg->sfactor/n;
#else
    theta = PETSC_PI*(i+0.5)/n;
    zn[i] = PetscCosReal(theta);
    w[i]  = PetscCosReal((n-1)*theta)/n;
    z[i]  = rg->sfactor*(ctx->center + ctx->radius*zn[i]);
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGCheckInside_Ellipse(RG rg,PetscReal px,PetscReal py,PetscInt *inside)
{
  RG_ELLIPSE *ctx = (RG_ELLIPSE*)rg->data;
  PetscReal  dx,dy,r;

  PetscFunctionBegin;
#if defined(PETSC_USE_COMPLEX)
  dx = (px-PetscRealPart(ctx->center))/ctx->radius;
  dy = (py-PetscImaginaryPart(ctx->center))/ctx->radius;
#else
  dx = (px-ctx->center)/ctx->radius;
  dy = py/ctx->radius;
#endif
  r = 1.0-dx*dx-(dy*dy)/(ctx->vscale*ctx->vscale);
  *inside = PetscSign(r);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGIsAxisymmetric_Ellipse(RG rg,PetscBool vertical,PetscBool *symm)
{
  RG_ELLIPSE *ctx = (RG_ELLIPSE*)rg->data;

  PetscFunctionBegin;
  if (vertical) *symm = (PetscRealPart(ctx->center) == 0.0)? PETSC_TRUE: PETSC_FALSE;
  else *symm = (PetscImaginaryPart(ctx->center) == 0.0)? PETSC_TRUE: PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGSetFromOptions_Ellipse(RG rg,PetscOptionItems PetscOptionsObject)
{
  PetscScalar    s;
  PetscReal      r1,r2;
  PetscBool      flg1,flg2,flg3;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"RG Ellipse Options");

    PetscCall(RGEllipseGetParameters(rg,&s,&r1,&r2));
    PetscCall(PetscOptionsScalar("-rg_ellipse_center","Center of ellipse","RGEllipseSetParameters",s,&s,&flg1));
    PetscCall(PetscOptionsReal("-rg_ellipse_radius","Radius of ellipse","RGEllipseSetParameters",r1,&r1,&flg2));
    PetscCall(PetscOptionsReal("-rg_ellipse_vscale","Vertical scale of ellipse","RGEllipseSetParameters",r2,&r2,&flg3));
    if (flg1 || flg2 || flg3) PetscCall(RGEllipseSetParameters(rg,s,r1,r2));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGDestroy_Ellipse(RG rg)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(rg->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)rg,"RGEllipseSetParameters_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)rg,"RGEllipseGetParameters_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode RGCreate_Ellipse(RG rg)
{
  RG_ELLIPSE     *ellipse;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ellipse));
  ellipse->center = 0.0;
  ellipse->radius = PETSC_MAX_REAL;
  ellipse->vscale = 1.0;
  rg->data = (void*)ellipse;

  rg->ops->istrivial         = RGIsTrivial_Ellipse;
  rg->ops->computecontour    = RGComputeContour_Ellipse;
  rg->ops->computebbox       = RGComputeBoundingBox_Ellipse;
  rg->ops->computequadrature = RGComputeQuadrature_Ellipse;
  rg->ops->checkinside       = RGCheckInside_Ellipse;
  rg->ops->isaxisymmetric    = RGIsAxisymmetric_Ellipse;
  rg->ops->setfromoptions    = RGSetFromOptions_Ellipse;
  rg->ops->view              = RGView_Ellipse;
  rg->ops->destroy           = RGDestroy_Ellipse;
  PetscCall(PetscObjectComposeFunction((PetscObject)rg,"RGEllipseSetParameters_C",RGEllipseSetParameters_Ellipse));
  PetscCall(PetscObjectComposeFunction((PetscObject)rg,"RGEllipseGetParameters_C",RGEllipseGetParameters_Ellipse));
  PetscFunctionReturn(PETSC_SUCCESS);
}
