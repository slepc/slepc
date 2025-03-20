/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Ring region, similar to the ellipse but with a start and end angle,
   together with the width
*/

#include <slepc/private/rgimpl.h>      /*I "slepcrg.h" I*/
#include <petscdraw.h>

typedef struct {
  PetscScalar center;     /* center of the ellipse */
  PetscReal   radius;     /* radius of the ellipse */
  PetscReal   vscale;     /* vertical scale of the ellipse */
  PetscReal   start_ang;  /* start angle */
  PetscReal   end_ang;    /* end angle */
  PetscReal   width;      /* ring width */
} RG_RING;

static PetscErrorCode RGRingSetParameters_Ring(RG rg,PetscScalar center,PetscReal radius,PetscReal vscale,PetscReal start_ang,PetscReal end_ang,PetscReal width)
{
  RG_RING *ctx = (RG_RING*)rg->data;

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
  if (start_ang == (PetscReal)PETSC_DETERMINE) {
    ctx->start_ang = 0.0;
  } else if (start_ang != (PetscReal)PETSC_CURRENT) {
    PetscCheck(start_ang>=0.0,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The right-hand side angle argument must be >= 0.0");
    PetscCheck(start_ang<=1.0,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The right-hand side angle argument must be <= 1.0");
    ctx->start_ang = start_ang;
  }
  if (end_ang == (PetscReal)PETSC_DETERMINE) {
    ctx->end_ang = 1.0;
  } else if (end_ang != (PetscReal)PETSC_CURRENT) {
    PetscCheck(end_ang>=0.0,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The left-hand side angle argument must be >= 0.0");
    PetscCheck(end_ang<=1.0,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The left-hand side angle argument must be <= 1.0");
    ctx->end_ang = end_ang;
  }
#if !defined(PETSC_USE_COMPLEX)
  PetscCheck(ctx->start_ang+ctx->end_ang==1.0,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_WRONG,"In real scalars the region must be symmetric wrt real axis");
#endif
  if (width == (PetscReal)PETSC_DETERMINE) {
    ctx->width = 0.1;
  } else if (width != (PetscReal)PETSC_CURRENT) {
    PetscCheck(width>0.0,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The width argument must be > 0.0");
    ctx->width = width;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGRingSetParameters - Sets the parameters defining the ring region.

   Logically Collective

   Input Parameters:
+  rg        - the region context
.  center    - center of the ellipse
.  radius    - radius of the ellipse
.  vscale    - vertical scale of the ellipse
.  start_ang - the right-hand side angle
.  end_ang   - the left-hand side angle
-  width     - width of the ring

   Options Database Keys:
+  -rg_ring_center     - Sets the center
.  -rg_ring_radius     - Sets the radius
.  -rg_ring_vscale     - Sets the vertical scale
.  -rg_ring_startangle - Sets the right-hand side angle
.  -rg_ring_endangle   - Sets the left-hand side angle
-  -rg_ring_width      - Sets the width of the ring

   Notes:
   The values of center, radius and vscale have the same meaning as in the
   ellipse region. The startangle and endangle define the span of the ring
   (by default it is the whole ring), while the width is the separation
   between the two concentric ellipses (above and below the radius by
   width/2).

   The start and end angles are expressed as a fraction of the circumference.
   The allowed range is [0..1], with 0 corresponding to 0 radians, 0.25 to
   pi/2 radians, and so on. It is allowed to have startangle>endangle, in
   which case the ring region crosses over the zero angle.

   In the case of complex scalars, a complex center can be provided in the
   command line with [+/-][realnumber][+/-]realnumberi with no spaces, e.g.
   -rg_ring_center 1.0+2.0i

   When PETSc is built with real scalars, the center is restricted to a real value,
   and the start and end angles must be such that the region is symmetric with
   respect to the real axis.

   For all arguments except center, you can use PETSC_CURRENT to keep the current
   value, and PETSC_DETERMINE to set them to a default value (1 for radius, vscale,
   end_ang, 0 for start_and, and 0.1 for width).

   Level: advanced

.seealso: RGRingGetParameters()
@*/
PetscErrorCode RGRingSetParameters(RG rg,PetscScalar center,PetscReal radius,PetscReal vscale,PetscReal start_ang,PetscReal end_ang,PetscReal width)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidLogicalCollectiveScalar(rg,center,2);
  PetscValidLogicalCollectiveReal(rg,radius,3);
  PetscValidLogicalCollectiveReal(rg,vscale,4);
  PetscValidLogicalCollectiveReal(rg,start_ang,5);
  PetscValidLogicalCollectiveReal(rg,end_ang,6);
  PetscValidLogicalCollectiveReal(rg,width,7);
  PetscTryMethod(rg,"RGRingSetParameters_C",(RG,PetscScalar,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal),(rg,center,radius,vscale,start_ang,end_ang,width));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGRingGetParameters_Ring(RG rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscReal *start_ang,PetscReal *end_ang,PetscReal *width)
{
  RG_RING *ctx = (RG_RING*)rg->data;

  PetscFunctionBegin;
  if (center)    *center    = ctx->center;
  if (radius)    *radius    = ctx->radius;
  if (vscale)    *vscale    = ctx->vscale;
  if (start_ang) *start_ang = ctx->start_ang;
  if (end_ang)   *end_ang   = ctx->end_ang;
  if (width)     *width     = ctx->width;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGRingGetParameters - Gets the parameters that define the ring region.

   Not Collective

   Input Parameter:
.  rg     - the region context

   Output Parameters:
+  center    - center of the region
.  radius    - radius of the region
.  vscale    - vertical scale of the region
.  start_ang - the right-hand side angle
.  end_ang   - the left-hand side angle
-  width     - width of the ring

   Level: advanced

.seealso: RGRingSetParameters()
@*/
PetscErrorCode RGRingGetParameters(RG rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscReal *start_ang,PetscReal *end_ang,PetscReal *width)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscUseMethod(rg,"RGRingGetParameters_C",(RG,PetscScalar*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*),(rg,center,radius,vscale,start_ang,end_ang,width));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGView_Ring(RG rg,PetscViewer viewer)
{
  RG_RING        *ctx = (RG_RING*)rg->data;
  int            winw,winh;
  PetscBool      isdraw,isascii;
  PetscDraw      draw;
  PetscDrawAxis  axis;
  PetscReal      cx,cy,radius,width,ab,cd,lx,ly,w,end_ang,x1,y1,x2,y2,r,theta,scale=1.2;
  char           str[50];

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw));
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->center,PETSC_FALSE));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  center: %s, radius: %g, vscale: %g, start angle: %g, end angle: %g, ring width: %g\n",str,RGShowReal(ctx->radius),RGShowReal(ctx->vscale),(double)ctx->start_ang,(double)ctx->end_ang,(double)ctx->width));
  } else if (isdraw) {
    PetscCall(PetscViewerDrawGetDraw(viewer,0,&draw));
    PetscCall(PetscDrawCheckResizedWindow(draw));
    PetscCall(PetscDrawGetWindowSize(draw,&winw,&winh));
    winw = PetscMax(winw,1); winh = PetscMax(winh,1);
    PetscCall(PetscDrawClear(draw));
    PetscCall(PetscDrawSetTitle(draw,"Ring region"));
    PetscCall(PetscDrawAxisCreate(draw,&axis));
    cx = PetscRealPart(ctx->center)*rg->sfactor;
    cy = PetscImaginaryPart(ctx->center)*rg->sfactor;
    radius = ctx->radius*rg->sfactor;
    width  = ctx->width*rg->sfactor;
    lx = 2*(radius+width);
    ly = 2*(radius+width)*ctx->vscale;
    ab = cx;
    cd = cy;
    w  = scale*PetscMax(lx/winw,ly/winh)/2;
    PetscCall(PetscDrawAxisSetLimits(axis,ab-w*winw,ab+w*winw,cd-w*winh,cd+w*winh));
    PetscCall(PetscDrawAxisDraw(axis));
    PetscCall(PetscDrawAxisDestroy(&axis));
    /* draw outer ellipse */
    PetscCall(PetscDrawEllipse(draw,cx,cy,2*(radius+width),2*(radius+width)*ctx->vscale,PETSC_DRAW_ORANGE));
    /* remove inner part */
    PetscCall(PetscDrawEllipse(draw,cx,cy,2*(radius-width),2*(radius-width)*ctx->vscale,PETSC_DRAW_WHITE));
    if (ctx->start_ang!=ctx->end_ang) {
      /* remove section from end_ang to start_ang */
      end_ang = (ctx->start_ang<ctx->end_ang)? ctx->end_ang-1: ctx->end_ang;
      theta = end_ang;
      r = scale*(radius+width);
      if (ctx->vscale>1) r *= ctx->vscale;
      x1 = PetscMin(PetscMax(ab+r*PetscCosReal(2.0*PETSC_PI*theta),ab-w*winw),ab+w*winw);
      y1 = PetscMin(PetscMax(cd+r*PetscSinReal(2.0*PETSC_PI*theta),cd-w*winh),cd+w*winh);
      do {
        theta = PetscMin(PetscFloorReal(8*theta+1)/8,ctx->start_ang);
        x2 = PetscMin(PetscMax(ab+r*PetscCosReal(2.0*PETSC_PI*theta),ab-w*winw),ab+w*winw);
        y2 = PetscMin(PetscMax(cd+r*PetscSinReal(2.0*PETSC_PI*theta),cd-w*winh),cd+w*winh);
        PetscCall(PetscDrawTriangle(draw,cx,cy,x1,y1,x2,y2,PETSC_DRAW_WHITE,PETSC_DRAW_WHITE,PETSC_DRAW_WHITE));
        x1 = x2; y1 = y2;
      } while (theta<ctx->start_ang);
    }
    PetscCall(PetscDrawFlush(draw));
    PetscCall(PetscDrawSave(draw));
    PetscCall(PetscDrawPause(draw));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGIsTrivial_Ring(RG rg,PetscBool *trivial)
{
  RG_RING *ctx = (RG_RING*)rg->data;

  PetscFunctionBegin;
  if (rg->complement) *trivial = PetscNot(ctx->radius);
  else *trivial = PetscNot(ctx->radius<PETSC_MAX_REAL);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGComputeContour_Ring(RG rg,PetscInt n,PetscScalar *cr,PetscScalar *ci)
{
  RG_RING   *ctx = (RG_RING*)rg->data;
  PetscReal theta,start_ang;
  PetscInt  i,n2=n/2;

  PetscFunctionBegin;
  start_ang = (ctx->start_ang>ctx->end_ang)? ctx->start_ang-1: ctx->start_ang;
  for (i=0;i<n;i++) {
    if (i < n2) {
      theta = ((ctx->end_ang-start_ang)*i/n2 + start_ang)*2.0*PETSC_PI;
#if defined(PETSC_USE_COMPLEX)
      cr[i] = ctx->center + (ctx->radius+ctx->width/2.0)*PetscCMPLX(PetscCosReal(theta),ctx->vscale*PetscSinReal(theta));
#else
      if (cr) cr[i] = ctx->center + (ctx->radius+ctx->width/2.0)*PetscCosReal(theta);
      if (ci) ci[i] = (ctx->radius+ctx->width/2.0)*ctx->vscale*PetscSinReal(theta);
#endif
    } else {
      theta = ((ctx->end_ang-start_ang)*(n-i)/n2 + start_ang)*2.0*PETSC_PI;
#if defined(PETSC_USE_COMPLEX)
      cr[i] = ctx->center + (ctx->radius-ctx->width/2.0)*PetscCMPLX(PetscCosReal(theta),ctx->vscale*PetscSinReal(theta));
#else
      if (cr) cr[i] = ctx->center + (ctx->radius-ctx->width/2.0)*PetscCosReal(theta);
      if (ci) ci[i] = (ctx->radius-ctx->width/2.0)*ctx->vscale*PetscSinReal(theta);
#endif
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGComputeBoundingBox_Ring(RG rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d)
{
  RG_RING *ctx = (RG_RING*)rg->data;

  PetscFunctionBegin;
  /* current implementation does not return a tight bounding box */
  if (a) *a = PetscRealPart(ctx->center) - (ctx->radius+ctx->width/2.0);
  if (b) *b = PetscRealPart(ctx->center) + (ctx->radius+ctx->width/2.0);
  if (c) *c = PetscImaginaryPart(ctx->center) - (ctx->radius+ctx->width/2.0)*ctx->vscale;
  if (d) *d = PetscImaginaryPart(ctx->center) + (ctx->radius+ctx->width/2.0)*ctx->vscale;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGComputeQuadrature_Ring(RG rg,RGQuadRule quad,PetscInt n,PetscScalar *z,PetscScalar *zn,PetscScalar *w)
{
  RG_RING     *ctx = (RG_RING*)rg->data;
  PetscReal   max_w=0.0;
  PetscScalar tmp,tmp2;
  PetscInt    i,j;

  PetscFunctionBegin;
  if (quad == RG_QUADRULE_CHEBYSHEV) {
#if defined(PETSC_USE_COMPLEX)
    PetscReal theta;
    for (i=0;i<n;i++) {
      theta = PETSC_PI*(i+0.5)/n;
      zn[i] = PetscCosReal(theta);
      w[i]  = PetscCosReal((n-1)*theta)/n;
      theta = (ctx->start_ang*2.0+(ctx->end_ang-ctx->start_ang)*(PetscRealPart(zn[i])+1.0))*PETSC_PI;
      z[i] = rg->sfactor*(ctx->center + ctx->radius*PetscCMPLX(PetscCosReal(theta),ctx->vscale*PetscSinReal(theta)));
    }
#else
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Integration points on a vertical line require complex arithmetic");
#endif
  } else {  /* RG_QUADRULE_TRAPEZOIDAL */
    for (i=0;i<n;i++) {
      zn[i] = (z[i]-rg->sfactor*ctx->center)/(rg->sfactor*ctx->radius);
      tmp = 1.0; tmp2 = 1.0;
      for (j=0;j<n;j++) {
        tmp *= z[j];
        if (i != j) tmp2 *= z[j]-z[i];
      }
      w[i] = tmp/tmp2;
      max_w = PetscMax(PetscAbsScalar(w[i]),max_w);
    }
    for (i=0;i<n;i++) w[i] /= (PetscScalar)max_w;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGCheckInside_Ring(RG rg,PetscReal px,PetscReal py,PetscInt *inside)
{
  RG_RING   *ctx = (RG_RING*)rg->data;
  PetscReal dx,dy,r;

  PetscFunctionBegin;
  /* outer ellipse */
#if defined(PETSC_USE_COMPLEX)
  dx = (px-PetscRealPart(ctx->center))/(ctx->radius+ctx->width/2.0);
  dy = (py-PetscImaginaryPart(ctx->center))/(ctx->radius+ctx->width/2.0);
#else
  dx = (px-ctx->center)/(ctx->radius+ctx->width/2.0);
  dy = py/(ctx->radius+ctx->width/2.0);
#endif
  r = 1.0-dx*dx-(dy*dy)/(ctx->vscale*ctx->vscale);
  *inside = PetscSign(r);
  /* inner ellipse */
#if defined(PETSC_USE_COMPLEX)
  dx = (px-PetscRealPart(ctx->center))/(ctx->radius-ctx->width/2.0);
  dy = (py-PetscImaginaryPart(ctx->center))/(ctx->radius-ctx->width/2.0);
#else
  dx = (px-ctx->center)/(ctx->radius-ctx->width/2.0);
  dy = py/(ctx->radius-ctx->width/2.0);
#endif
  r = -1.0+dx*dx+(dy*dy)/(ctx->vscale*ctx->vscale);
  *inside *= PetscSign(r);
  if (*inside == 1) {  /* check angles */
#if defined(PETSC_USE_COMPLEX)
    dx = (px-PetscRealPart(ctx->center));
    dy = (py-PetscImaginaryPart(ctx->center));
#else
    dx = px-ctx->center;
    dy = py;
#endif
    if (dx == 0) {
      if (dy == 0) r = -1;
      else if (dy > 0) r = 0.25;
      else r = 0.75;
    } else if (dx > 0) {
      r = PetscAtanReal((dy/ctx->vscale)/dx);
      if (dy >= 0) r /= 2*PETSC_PI;
      else r = r/(2*PETSC_PI)+1;
    } else r = PetscAtanReal((dy/ctx->vscale)/dx)/(2*PETSC_PI)+0.5;
    if (ctx->start_ang>ctx->end_ang) {
      if (r>ctx->end_ang && r<ctx->start_ang) *inside = -1;
    } else {
      if (r<ctx->start_ang || r>ctx->end_ang) *inside = -1;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGIsAxisymmetric_Ring(RG rg,PetscBool vertical,PetscBool *symm)
{
  RG_RING *ctx = (RG_RING*)rg->data;

  PetscFunctionBegin;
  if (vertical) *symm = (PetscRealPart(ctx->center) == 0.0 && PetscAbs(ctx->start_ang+ctx->end_ang-PetscRealConstant(1.0)) == 0.5)? PETSC_TRUE: PETSC_FALSE;
  else *symm = (PetscImaginaryPart(ctx->center) == 0.0 && ctx->start_ang+ctx->end_ang == 1.0)? PETSC_TRUE: PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGSetFromOptions_Ring(RG rg,PetscOptionItems PetscOptionsObject)
{
  PetscScalar    s;
  PetscReal      r1,r2,r3,r4,r5;
  PetscBool      flg1,flg2,flg3,flg4,flg5,flg6;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"RG Ring Options");

    PetscCall(RGRingGetParameters(rg,&s,&r1,&r2,&r3,&r4,&r5));
    PetscCall(PetscOptionsScalar("-rg_ring_center","Center of ellipse","RGRingSetParameters",s,&s,&flg1));
    PetscCall(PetscOptionsReal("-rg_ring_radius","Radius of ellipse","RGRingSetParameters",r1,&r1,&flg2));
    PetscCall(PetscOptionsReal("-rg_ring_vscale","Vertical scale of ellipse","RGRingSetParameters",r2,&r2,&flg3));
    PetscCall(PetscOptionsReal("-rg_ring_startangle","Right-hand side angle","RGRingSetParameters",r3,&r3,&flg4));
    PetscCall(PetscOptionsReal("-rg_ring_endangle","Left-hand side angle","RGRingSetParameters",r4,&r4,&flg5));
    PetscCall(PetscOptionsReal("-rg_ring_width","Width of ring","RGRingSetParameters",r5,&r5,&flg6));
    if (flg1 || flg2 || flg3 || flg4 || flg5 || flg6) PetscCall(RGRingSetParameters(rg,s,r1,r2,r3,r4,r5));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGDestroy_Ring(RG rg)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(rg->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)rg,"RGRingSetParameters_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)rg,"RGRingGetParameters_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode RGCreate_Ring(RG rg)
{
  RG_RING        *ring;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ring));
  ring->center    = 0.0;
  ring->radius    = PETSC_MAX_REAL;
  ring->vscale    = 1.0;
  ring->start_ang = 0.0;
  ring->end_ang   = 1.0;
  ring->width     = 0.1;
  rg->data = (void*)ring;

  rg->ops->istrivial         = RGIsTrivial_Ring;
  rg->ops->computecontour    = RGComputeContour_Ring;
  rg->ops->computebbox       = RGComputeBoundingBox_Ring;
  rg->ops->computequadrature = RGComputeQuadrature_Ring;
  rg->ops->checkinside       = RGCheckInside_Ring;
  rg->ops->isaxisymmetric    = RGIsAxisymmetric_Ring;
  rg->ops->setfromoptions    = RGSetFromOptions_Ring;
  rg->ops->view              = RGView_Ring;
  rg->ops->destroy           = RGDestroy_Ring;
  PetscCall(PetscObjectComposeFunction((PetscObject)rg,"RGRingSetParameters_C",RGRingSetParameters_Ring));
  PetscCall(PetscObjectComposeFunction((PetscObject)rg,"RGRingGetParameters_C",RGRingGetParameters_Ring));
  PetscFunctionReturn(PETSC_SUCCESS);
}
