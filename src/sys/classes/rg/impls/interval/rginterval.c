/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Interval of the real axis or more generally a (possibly open) rectangle
   of the complex plane
*/

#include <slepc/private/rgimpl.h>      /*I "slepcrg.h" I*/
#include <petscdraw.h>

typedef struct {
  PetscReal   a,b;     /* interval in the real axis */
  PetscReal   c,d;     /* interval in the imaginary axis */
} RG_INTERVAL;

static PetscErrorCode RGIntervalSetEndpoints_Interval(RG rg,PetscReal a,PetscReal b,PetscReal c,PetscReal d)
{
  RG_INTERVAL *ctx = (RG_INTERVAL*)rg->data;

  PetscFunctionBegin;
  PetscCheck(a || b || c || d,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_WRONG,"At least one argument must be nonzero");
  PetscCheck(a!=b || !a,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_WRONG,"Badly defined interval, endpoints must be distinct (or both zero)");
  PetscCheck(a<=b,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_WRONG,"Badly defined interval, must be a<b");
  PetscCheck(c!=d || !c,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_WRONG,"Badly defined interval, endpoints must be distinct (or both zero)");
  PetscCheck(c<=d,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_WRONG,"Badly defined interval, must be c<d");
#if !defined(PETSC_USE_COMPLEX)
  PetscCheck(c==-d,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_WRONG,"In real scalars the region must be symmetric wrt real axis");
#endif
  ctx->a = a;
  ctx->b = b;
  ctx->c = c;
  ctx->d = d;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGIntervalSetEndpoints - Sets the parameters defining the interval region.

   Logically Collective

   Input Parameters:
+  rg  - the region context
.  a - left endpoint of the interval in the real axis
.  b - right endpoint of the interval in the real axis
.  c - bottom endpoint of the interval in the imaginary axis
-  d - top endpoint of the interval in the imaginary axis

   Options Database Keys:
.  -rg_interval_endpoints - the four endpoints

   Note:
   The region is defined as [a,b]x[c,d]. Particular cases are an interval on
   the real axis (c=d=0), similar for the imaginary axis (a=b=0), the whole
   complex plane (a=-inf,b=inf,c=-inf,d=inf), and so on.

   When PETSc is built with real scalars, the region must be symmetric with
   respect to the real axis.

   Level: advanced

.seealso: RGIntervalGetEndpoints()
@*/
PetscErrorCode RGIntervalSetEndpoints(RG rg,PetscReal a,PetscReal b,PetscReal c,PetscReal d)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidLogicalCollectiveReal(rg,a,2);
  PetscValidLogicalCollectiveReal(rg,b,3);
  PetscValidLogicalCollectiveReal(rg,c,4);
  PetscValidLogicalCollectiveReal(rg,d,5);
  PetscTryMethod(rg,"RGIntervalSetEndpoints_C",(RG,PetscReal,PetscReal,PetscReal,PetscReal),(rg,a,b,c,d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGIntervalGetEndpoints_Interval(RG rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d)
{
  RG_INTERVAL *ctx = (RG_INTERVAL*)rg->data;

  PetscFunctionBegin;
  if (a) *a = ctx->a;
  if (b) *b = ctx->b;
  if (c) *c = ctx->c;
  if (d) *d = ctx->d;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGIntervalGetEndpoints - Gets the parameters that define the interval region.

   Not Collective

   Input Parameter:
.  rg - the region context

   Output Parameters:
+  a - left endpoint of the interval in the real axis
.  b - right endpoint of the interval in the real axis
.  c - bottom endpoint of the interval in the imaginary axis
-  d - top endpoint of the interval in the imaginary axis

   Level: advanced

.seealso: RGIntervalSetEndpoints()
@*/
PetscErrorCode RGIntervalGetEndpoints(RG rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscUseMethod(rg,"RGIntervalGetEndpoints_C",(RG,PetscReal*,PetscReal*,PetscReal*,PetscReal*),(rg,a,b,c,d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGView_Interval(RG rg,PetscViewer viewer)
{
  RG_INTERVAL    *ctx = (RG_INTERVAL*)rg->data;
  PetscBool      isdraw,isascii;
  int            winw,winh;
  PetscDraw      draw;
  PetscDrawAxis  axis;
  PetscReal      a,b,c,d,ab,cd,lx,ly,w,scale=1.2;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw));
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) PetscCall(PetscViewerASCIIPrintf(viewer,"  region: [%g,%g]x[%g,%g]\n",RGShowReal(ctx->a),RGShowReal(ctx->b),RGShowReal(ctx->c),RGShowReal(ctx->d)));
  else if (isdraw) {
    PetscCall(PetscViewerDrawGetDraw(viewer,0,&draw));
    PetscCall(PetscDrawCheckResizedWindow(draw));
    PetscCall(PetscDrawGetWindowSize(draw,&winw,&winh));
    winw = PetscMax(winw,1); winh = PetscMax(winh,1);
    PetscCall(PetscDrawClear(draw));
    PetscCall(PetscDrawSetTitle(draw,"Interval region"));
    PetscCall(PetscDrawAxisCreate(draw,&axis));
    a  = ctx->a*rg->sfactor;
    b  = ctx->b*rg->sfactor;
    c  = ctx->c*rg->sfactor;
    d  = ctx->d*rg->sfactor;
    lx = b-a;
    ly = d-c;
    ab = (a+b)/2;
    cd = (c+d)/2;
    w  = scale*PetscMax(lx/winw,ly/winh)/2;
    PetscCall(PetscDrawAxisSetLimits(axis,ab-w*winw,ab+w*winw,cd-w*winh,cd+w*winh));
    PetscCall(PetscDrawAxisDraw(axis));
    PetscCall(PetscDrawAxisDestroy(&axis));
    PetscCall(PetscDrawRectangle(draw,a,c,b,d,PETSC_DRAW_BLUE,PETSC_DRAW_BLUE,PETSC_DRAW_BLUE,PETSC_DRAW_BLUE));
    PetscCall(PetscDrawFlush(draw));
    PetscCall(PetscDrawSave(draw));
    PetscCall(PetscDrawPause(draw));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGIsTrivial_Interval(RG rg,PetscBool *trivial)
{
  RG_INTERVAL *ctx = (RG_INTERVAL*)rg->data;

  PetscFunctionBegin;
  if (rg->complement) *trivial = (ctx->a==ctx->b && ctx->c==ctx->d)? PETSC_TRUE: PETSC_FALSE;
  else *trivial = (ctx->a<=-PETSC_MAX_REAL && ctx->b>=PETSC_MAX_REAL && ctx->c<=-PETSC_MAX_REAL && ctx->d>=PETSC_MAX_REAL)? PETSC_TRUE: PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGComputeContour_Interval(RG rg,PetscInt n,PetscScalar *cr,PetscScalar *ci)
{
  RG_INTERVAL *ctx = (RG_INTERVAL*)rg->data;
  PetscInt    i,N,Nv,Nh,k1,k0;
  PetscReal   hv,hh,t;

  PetscFunctionBegin;
  PetscCheck(ctx->a>-PETSC_MAX_REAL && ctx->b<PETSC_MAX_REAL && ctx->c>-PETSC_MAX_REAL && ctx->d<PETSC_MAX_REAL,PetscObjectComm((PetscObject)rg),PETSC_ERR_SUP,"Contour not defined in unbounded regions");
  if (ctx->a==ctx->b || ctx->c==ctx->d) {
    if (ctx->a==ctx->b) {hv = (ctx->d-ctx->c)/(n-1); hh = 0.0;}
    else {hh = (ctx->b-ctx->a)/(n-1); hv = 0.0;}
    for (i=0;i<n;i++) {
#if defined(PETSC_USE_COMPLEX)
      cr[i] = PetscCMPLX(ctx->a+hh*i,ctx->c+hv*i);
#else
      if (cr) cr[i] = ctx->a+hh*i;
      if (ci) ci[i] = ctx->c+hv*i;
#endif
    }
  } else {
    PetscCheck(n>3,PetscObjectComm((PetscObject)rg),PETSC_ERR_SUP,"Minimum number of contour points: 4");
    N = n/2;
    t = ((ctx->d-ctx->c)/(ctx->d-ctx->c+ctx->b-ctx->a))*N;
    Nv = t-PetscFloorReal(t)>0.5?PetscCeilReal(t):PetscFloorReal(t);
    if (Nv==0) Nv++;
    else if (Nv==N) Nv--;
    Nh = N-Nv;
    hh = (ctx->b-ctx->a)/Nh;
    hv = (ctx->d-ctx->c)/Nv;
    /* positive imaginary part first */
    k1 = Nv/2+1;
    k0 = Nv-k1;

    for (i=k1;i<Nv;i++) {
#if defined(PETSC_USE_COMPLEX)
      cr[i-k1]   = PetscCMPLX(ctx->b,ctx->c+i*hv);
      cr[i-k1+N] = PetscCMPLX(ctx->a,ctx->d-i*hv);
#else
      if (cr) {cr[i-k1] = ctx->b;      cr[i-k1+N] = ctx->a;}
      if (ci) {ci[i-k1] = ctx->c+i*hv; ci[i-k1+N] = ctx->d-i*hv;}
#endif
    }
    for (i=0;i<Nh;i++) {
#if defined(PETSC_USE_COMPLEX)
      cr[i+k0]   = PetscCMPLX(ctx->b-i*hh,ctx->d);
      cr[i+k0+N] = PetscCMPLX(ctx->a+i*hh,ctx->c);
#else
      if (cr) {cr[i+k0] = ctx->b-i*hh; cr[i+k0+N] = ctx->a+i*hh;}
      if (ci) {ci[i+k0] = ctx->d;      ci[i+k0+N] = ctx->c;}
#endif
    }
    for (i=0;i<k1;i++) {
#if defined(PETSC_USE_COMPLEX)
      cr[i+k0+Nh]   = PetscCMPLX(ctx->a,ctx->d-i*hv);
      cr[i+k0+Nh+N] = PetscCMPLX(ctx->b,ctx->c+i*hv);
#else
      if (cr) {cr[i+k0+Nh] = ctx->a; cr[i+k0+Nh+N] = ctx->b;}
      if (ci) {ci[i+k0+Nh] = ctx->d+i*hv; ci[i+k0+Nh+N] = ctx->c-i*hv;}
#endif
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGComputeBoundingBox_Interval(RG rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d)
{
  RG_INTERVAL *ctx = (RG_INTERVAL*)rg->data;

  PetscFunctionBegin;
  if (a) *a = ctx->a;
  if (b) *b = ctx->b;
  if (c) *c = ctx->c;
  if (d) *d = ctx->d;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGComputeQuadrature_Interval(RG rg,RGQuadRule quad,PetscInt n,PetscScalar *z,PetscScalar *zn,PetscScalar *w)
{
  RG_INTERVAL *ctx = (RG_INTERVAL*)rg->data;
  PetscReal   theta,max_w=0.0,radius=1.0;
  PetscScalar tmp,tmp2,center=0.0;
  PetscInt    i,j;

  PetscFunctionBegin;
  if (quad == RG_QUADRULE_CHEBYSHEV) {
    PetscCheck((ctx->c==ctx->d && ctx->c==0.0) || (ctx->a==ctx->b && ctx->a==0.0),PetscObjectComm((PetscObject)rg),PETSC_ERR_SUP,"Endpoints of the imaginary axis or the real axis must be both zero");
    for (i=0;i<n;i++) {
      theta = PETSC_PI*(i+0.5)/n;
      zn[i] = PetscCosReal(theta);
      w[i]  = PetscCosReal((n-1)*theta)/n;
      if (ctx->c==ctx->d) z[i] = ((zn[i]+1.0)*(ctx->b-ctx->a)/2.0+ctx->a)*rg->sfactor;
      else if (ctx->a==ctx->b) {
#if defined(PETSC_USE_COMPLEX)
        z[i] = ((zn[i]+1.0)*(ctx->d-ctx->c)/2.0+ctx->c)*rg->sfactor*PETSC_i;
#else
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Integration points on a vertical line require complex arithmetic");
#endif
      }
    }
  } else {  /* RG_QUADRULE_TRAPEZOIDAL */
#if defined(PETSC_USE_COMPLEX)
    center = PetscCMPLX(ctx->b+ctx->a,ctx->d+ctx->c)*rg->sfactor/2.0;
#else
    center = (ctx->b+ctx->a)*rg->sfactor/2.0;
#endif
    radius = PetscSqrtReal(PetscPowRealInt(rg->sfactor*(ctx->b-ctx->a)/2.0,2)+PetscPowRealInt(rg->sfactor*(ctx->d-ctx->c)/2.0,2));
    for (i=0;i<n;i++) {
      zn[i] = (z[i]-center)/radius;
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

static PetscErrorCode RGCheckInside_Interval(RG rg,PetscReal dx,PetscReal dy,PetscInt *inside)
{
  RG_INTERVAL *ctx = (RG_INTERVAL*)rg->data;

  PetscFunctionBegin;
  if (dx>ctx->a && dx<ctx->b) *inside = 1;
  else if (dx==ctx->a || dx==ctx->b) *inside = 0;
  else *inside = -1;
  if (*inside>=0) {
    if (dy>ctx->c && dy<ctx->d) ;
    else if (dy==ctx->c || dy==ctx->d) *inside = 0;
    else *inside = -1;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGIsAxisymmetric_Interval(RG rg,PetscBool vertical,PetscBool *symm)
{
  RG_INTERVAL *ctx = (RG_INTERVAL*)rg->data;

  PetscFunctionBegin;
  if (vertical) *symm = (ctx->a == -ctx->b)? PETSC_TRUE: PETSC_FALSE;
  else *symm = (ctx->c == -ctx->d)? PETSC_TRUE: PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGSetFromOptions_Interval(RG rg,PetscOptionItems PetscOptionsObject)
{
  PetscBool      flg;
  PetscInt       k;
  PetscReal      array[4]={0,0,0,0};

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"RG Interval Options");

    k = 4;
    PetscCall(PetscOptionsRealArray("-rg_interval_endpoints","Interval endpoints (two or four real values separated with a comma without spaces)","RGIntervalSetEndpoints",array,&k,&flg));
    if (flg) {
      PetscCheck(k>1,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_SIZ,"Must pass at least two values in -rg_interval_endpoints (comma-separated without spaces)");
      PetscCall(RGIntervalSetEndpoints(rg,array[0],array[1],array[2],array[3]));
    }

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RGDestroy_Interval(RG rg)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(rg->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)rg,"RGIntervalSetEndpoints_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)rg,"RGIntervalGetEndpoints_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode RGCreate_Interval(RG rg)
{
  RG_INTERVAL    *interval;

  PetscFunctionBegin;
  PetscCall(PetscNew(&interval));
  interval->a = -PETSC_MAX_REAL;
  interval->b = PETSC_MAX_REAL;
  interval->c = -PETSC_MAX_REAL;
  interval->d = PETSC_MAX_REAL;
  rg->data = (void*)interval;

  rg->ops->istrivial         = RGIsTrivial_Interval;
  rg->ops->computecontour    = RGComputeContour_Interval;
  rg->ops->computebbox       = RGComputeBoundingBox_Interval;
  rg->ops->computequadrature = RGComputeQuadrature_Interval;
  rg->ops->checkinside       = RGCheckInside_Interval;
  rg->ops->isaxisymmetric    = RGIsAxisymmetric_Interval;
  rg->ops->setfromoptions    = RGSetFromOptions_Interval;
  rg->ops->view              = RGView_Interval;
  rg->ops->destroy           = RGDestroy_Interval;
  PetscCall(PetscObjectComposeFunction((PetscObject)rg,"RGIntervalSetEndpoints_C",RGIntervalSetEndpoints_Interval));
  PetscCall(PetscObjectComposeFunction((PetscObject)rg,"RGIntervalGetEndpoints_C",RGIntervalGetEndpoints_Interval));
  PetscFunctionReturn(PETSC_SUCCESS);
}
