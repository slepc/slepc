/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc matrix function solver: "cram"

   Method: Chebyshev Rational Approximation Method (CRAM)

   Algorithm:

       Evaluates exp(A*t)*n_0 using a rational partial fraction expansion
       of order k (k=32 in quadruple precision, k=16 in double precision,
       k=8 in single precision), solving k/2 shifted linear systems
       of the form (A*t - theta_i*I)*x_i = n_0.

   References:

       [1] M. Pusa, "Rational approximations to the matrix exponential in
           burnup calculations", Nucl. Sci. Eng. 169(2):155-167, 2011.

       [2] C. Fichtlscherer, "High-precision higher-order coefficients for
           the Chebyshev rational approximation method (CRAM)", Nucl. Sci.
           Eng., 2026.
*/

#include <slepc/private/mfnimpl.h>         /*I "slepcmfn.h" I*/

#if PetscDefined(USE_COMPLEX)

typedef struct {
  KSP ksp;
} MFN_CRAM;

#if PetscDefined(USE_REAL_SINGLE)

static const PetscInt k_half = 4; /* 8 poles (4 conjugate pairs) -> 4 linear systems */

/* Poles (theta) for 8th-order CRAM: {Real, Imag} (Fichtlscherer, 2026) */
static const PetscReal poles[4][2] = {
  {-3.408539510e+00F, +8.773034554e+00F},
  {+2.694909809e-01F, +6.082032587e+00F},
  {+3.220945240e+00F, +1.193619605e+00F},
  {+2.292249142e+00F, +3.600771493e+00F}
};

/* Residues (alpha) for 8th-order CRAM: {Real, Imag} (Fichtlscherer, 2026) */
static const PetscReal residues[4][2] = {
  {+2.812975702e-02F, +1.157738421e-02F},
  {-6.325880513e-01F, -4.439230969e-01F},
  {-1.831771706e+00F, -9.525608074e+00F},
  {+2.436240726e+00F, +3.716755614e+00F}
};

/* Limit at infinity (alpha_0), strictly real */
static const PetscReal alpha_0 = +1.172265212e-08F;

#elif PetscDefined(USE_REAL_DOUBLE)

static const PetscInt k_half = 8; /* 16 poles (8 conjugate pairs) -> 8 linear systems */

/* Poles (theta) for 16th-order CRAM: {Real, Imag} (Fichtlscherer, 2026) */
static const PetscReal poles[8][2] = {
  {-1.0843917078696988026e+01, +1.9277446167181652284e+01},
  {+6.4161776990994341857e+00, +1.1941223933701386699e+00},
  {-5.2649713434426468908e+00, +1.6220221473167927305e+01},
  {-1.4139284624888862117e+00, +1.3497725698892745388e+01},
  {+1.4193758971856659905e+00, +1.0925363484496722585e+01},
  {+5.9481522689511774823e+00, +3.5874573620183223162e+00},
  {+3.5091036084149180718e+00, +8.4361989858843750942e+00},
  {+4.9931747377179964192e+00, +5.9968817136039421951e+00}
};

/* Residues (alpha) for 16th-order CRAM: {Real, Imag} (Fichtlscherer, 2026) */
static const PetscReal residues[8][2] = {
  {-5.0901521865224928712e-07, -2.4220017652852287986e-05},
  {-6.4500878025539644564e+01, -2.2459440762652096092e+02},
  {+2.1151742182466031443e-04, +4.3892969647380673895e-03},
  {+4.1023136835410020949e-02, -1.5743466173455468195e-01},
  {-1.4793007113558000013e+00, +1.7686588323782937902e+00},
  {+1.1339775178483930464e+02, +1.0194721704215856386e+02},
  {+1.5059585270023467196e+01, -5.7514052776421820767e+00},
  {-6.2518392463207919933e+01, -1.1190391094283228881e+01}
};

/* Limit at infinity (alpha_0), strictly real */
static const PetscReal alpha_0 = +2.1248537104952237488e-16;

#else /* PETSC_USE_REAL___FLOAT128 */

static const PetscInt k_half = 16; /* 32 poles (16 conjugate pairs) -> 16 linear systems */

/* Poles (theta) for 32nd-order CRAM: {Real, Imag} (Fichtlscherer, 2026) */
static const PetscReal poles[16][2] = {
  {+1.272430588935006438619831346751499107999e+01Q, +1.194256446871873253814182775373845975982e+00Q},
  {-2.734798356101249028927136505131250496230e+01Q, +4.069231229671755061003931347709070433726e+01Q},
  {-1.958248015332671268829049473658241810680e+01Q, +3.710329224494451308018012644195943867764e+01Q},
  {-1.377054550962031303785013257426553756472e+01Q, +3.399214933428828134752235969733162088172e+01Q},
  {+1.248804035352348427522632727331537189658e+01Q, +3.584068696364334819273345217102862931955e+00Q},
  {-9.054506283982660325362485409093544062950e+00Q, +3.111173004087708982894818872949048202132e+01Q},
  {-5.098669890503974413539297099462474752433e+00Q, +2.837060507799691809709678191350235829766e+01Q},
  {+1.201313618215790812041619674962043121768e+01Q, +5.977830703724749868246555267506325605461e+00Q},
  {-1.727168130828051491700852331556141151279e+00Q, +2.572302624565581743656536239564595165169e+01Q},
  {+1.166498099877330444729979805678576482636e+00Q, +2.314228200408636858376649587524819849744e+01Q},
  {+1.129472112530613370369458932045984100204e+01Q, +8.378353071561875050996029275690726473591e+00Q},
  {+3.652646888760183870244619989890222182065e+00Q, +2.061118862125450611473158543516741826938e+01Q},
  {+5.780372991371440163796598281821897283591e+00Q, +1.811790953587012021504390558422512356722e+01Q},
  {+1.032515517992622155082674666967603079800e+01Q, +1.078873149419210022345695806278668984329e+01Q},
  {+7.585235004412852240648130191654563529173e+00Q, +1.565385209991484235798025241201654188964e+01Q},
  {+9.093582328296484346497675301567147762861e+00Q, +1.321250259877625667484986061086199364770e+01Q}
};

/* Residues (alpha) for 32nd-order CRAM: {Real, Imag} (Fichtlscherer, 2026) */
static const PetscReal residues[16][2] = {
  {-4.108255209226283820322653413591444491895e+04Q, -1.210418418929308645605002527539202295314e+05Q},
  {+1.138416890221272716569582224693904653440e-12Q, +1.846814581089479468813662486303093287698e-12Q},
  {-3.186074124953601116509269949438175305764e-09Q, -1.786077444562727737576281059401608901941e-09Q},
  {+8.805808613442816490946594094530109998343e-07Q, +4.637400254039854072442808265525955508572e-07Q},
  {+8.478038897863382471020644246238551309997e+04Q, +5.674642335847666608528498310077385144070e+04Q},
  {-7.364462586735769244840153987931373686307e-05Q, -6.066726191425203694817073500684972816821e-05Q},
  {+2.260249700352534764829664864534364772318e-03Q, +3.773111384528383154253854297689743820342e-03Q},
  {-6.469184596136451936285918531668443634023e+04Q, +4.000759720932547650266072413920725832786e+03Q},
  {-1.225633449829473842630611313426775343246e-02Q, -1.144257713565144519325304771362493677635e-01Q},
  {-7.191074931096019249201049126849476067628e-01Q, +1.749727380850491962756362480303976906946e+00Q},
  {+2.473784861357228354159575507900717889708e+04Q, -2.124017891241354126908046560688841290930e+04Q},
  {+1.678896739588043068167387079497727246642e+01Q, -1.245541197261106504125544765597804325553e+01Q},
  {-1.625640278379788357174529998350408514262e+02Q, +1.009563734411555346376356065466148026895e+01Q},
  {-2.808616041847053963726245323360213641340e+03Q, +1.257087618419213913653894147332835105731e+04Q},
  {+7.951725784483317070133077817095961090145e+02Q, +4.756232053890630195061920905332604572897e+02Q},
  {-1.583891838392792538117705657587952178212e+03Q, -3.619327242725449650665673945499726008997e+03Q}
};

/* Limit at infinity (alpha_0), strictly real */
static const PetscReal alpha_0 = +6.932444346272943657112238792737063129652e-32Q;

#endif

static PetscErrorCode MFNCRAMGetKSP_CRAM(MFN mfn,KSP *ksp)
{
  MFN_CRAM *ctx = (MFN_CRAM*)mfn->data;
  PC       pc;

  PetscFunctionBegin;
  if (!ctx->ksp) {
    PetscCall(KSPCreate(PetscObjectComm((PetscObject)mfn),&ctx->ksp));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)ctx->ksp,(PetscObject)mfn,1));
    PetscCall(KSPSetOptionsPrefix(ctx->ksp,((PetscObject)mfn)->prefix));
    PetscCall(KSPAppendOptionsPrefix(ctx->ksp,"mfn_cram_"));
    PetscCall(PetscObjectSetOptions((PetscObject)ctx->ksp,((PetscObject)mfn)->options));
    PetscCall(KSPSetType(ctx->ksp,KSPPREONLY));
    PetscCall(KSPGetPC(ctx->ksp,&pc));
    PetscCall(PCSetType(pc,PCLU));
  }
  *ksp = ctx->ksp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MFNSetUp_CRAM(MFN mfn)
{
  PetscBool isexp;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)mfn->fn,FNEXP,&isexp));
  PetscCheck(isexp,PetscObjectComm((PetscObject)mfn),PETSC_ERR_SUP,"This solver only supports the exponential function");
  if (mfn->max_it==PETSC_DETERMINE) mfn->max_it = 1;
  if (mfn->ncv==PETSC_DETERMINE) mfn->ncv = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MFNSolve_CRAM(MFN mfn,Vec n_0,Vec y)
{
  MFN_CRAM           *ctx = (MFN_CRAM*)mfn->data;
  Mat                At;
  Vec                aux,n_0_work;
  PetscInt           i;
  PetscScalar        t,sfactor,alpha,theta,theta0=0.0;
  KSPConvergedReason kspreason;

  PetscFunctionBegin;
  PetscCall(FNGetScale(mfn->fn,&t,&sfactor));
  PetscCheck(PetscImaginaryPart(t)==0.0,PetscObjectComm((PetscObject)mfn),PETSC_ERR_SUP,"This solver requires a real scaling factor t");

  /* Define temporary matrix At and vector aux for workspace */
  PetscCall(MatDuplicate(mfn->A,MAT_COPY_VALUES,&At));
  PetscCall(MatScale(At,t));
  PetscCall(VecDuplicate(n_0,&aux));
  if (n_0==y) {
    PetscCall(VecDuplicate(n_0,&n_0_work));
    PetscCall(VecCopy(n_0,n_0_work));
  } else n_0_work = n_0;

  if (!ctx->ksp) PetscCall(MFNCRAMGetKSP(mfn,&ctx->ksp));
  PetscCall(VecSet(y,0.0));
  /* Solve the k_half linear systems: (A*t - theta_i*I) * x_i = n_0 */
  for (i=0;i<k_half;i++) {
    theta = PetscCMPLX(poles[i][0],poles[i][1]);
    alpha = PetscCMPLX(residues[i][0],residues[i][1]);

    PetscCall(MatShift(At,theta0-theta));
    theta0 = theta; /* Save the previous shift to unshift in the next iteration */

    PetscCall(KSPSetOperators(ctx->ksp,At,At));
    if (mfn->transpose_solve) PetscCall(KSPSolveTranspose(ctx->ksp,n_0_work,aux));
    else PetscCall(KSPSolve(ctx->ksp,n_0_work,aux));
    PetscCall(KSPGetConvergedReason(ctx->ksp,&kspreason));

    if (kspreason < 0) {
      PetscCall(PetscInfo(mfn,"Linear solve failed at step %" PetscInt_FMT ", stopping CRAM\n",i));
      mfn->reason = MFN_DIVERGED_BREAKDOWN;
      break;
    }
    PetscCall(VecAXPY(y,alpha,aux));
  }

  /* Reconstruct final real solution: y = 2 * Re(y) + alpha_0 * n_0 */
  if (mfn->reason != MFN_DIVERGED_BREAKDOWN) {
    PetscCall(VecRealPart(y));
    PetscCall(VecAXPBY(y,alpha_0,2.0,n_0_work));
    PetscCall(VecScale(y,sfactor));
    mfn->reason = MFN_CONVERGED_ITS;
  }

  mfn->its = 1; /* CRAM is a one-step method */

  PetscCall(MatDestroy(&At));
  PetscCall(VecDestroy(&aux));
  if (n_0==y) PetscCall(VecDestroy(&n_0_work));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MFNSetFromOptions_CRAM(MFN mfn,PetscOptionItems PetscOptionsObject)
{
  MFN_CRAM *ctx = (MFN_CRAM*)mfn->data;

  PetscFunctionBegin;
  if (!ctx->ksp) PetscCall(MFNCRAMGetKSP(mfn,&ctx->ksp));
  PetscCall(KSPSetFromOptions(ctx->ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MFNView_CRAM(MFN mfn,PetscViewer viewer)
{
  MFN_CRAM  *ctx = (MFN_CRAM*)mfn->data;
  PetscBool isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  CRAM approximation order: %" PetscInt_FMT "\n",2*k_half));
    if (!ctx->ksp) PetscCall(MFNCRAMGetKSP(mfn,&ctx->ksp));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(KSPView(ctx->ksp,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MFNReset_CRAM(MFN mfn)
{
  MFN_CRAM *ctx = (MFN_CRAM*)mfn->data;

  PetscFunctionBegin;
  if (ctx->ksp) PetscCall(KSPReset(ctx->ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MFNDestroy_CRAM(MFN mfn)
{
  MFN_CRAM *ctx = (MFN_CRAM*)mfn->data;

  PetscFunctionBegin;
  PetscCall(KSPDestroy(&ctx->ksp));
  PetscCall(PetscFree(mfn->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)mfn,"MFNCRAMGetKSP_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

/*@
   MFNCRAMGetKSP - Obtain the linear solver object associated to the CRAM solver.

   Collective

   Input Parameter:
.  mfn - the matrix function context

   Output Parameter:
.  ksp - the linear solver object

   Level: advanced

.seealso: [](ch:mfn), `MFNCRAM`
@*/
PetscErrorCode MFNCRAMGetKSP(MFN mfn,KSP *ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscAssertPointer(ksp,2);
  PetscUseMethod(mfn,"MFNCRAMGetKSP_C",(MFN,KSP*),(mfn,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   MFNCRAM - MFNCRAM = "cram" - The Chebyshev Rational Approximation Method
   for the matrix exponential.

   Notes:
   This solver approximates the action of the matrix exponential $\exp(At)n_0$
   by evaluating a rational function in partial fraction form of order $k$,
   as described by {cite:t}`Pus11`, using high-precision coefficients computed
   by {cite:t}`Fic26`. The approximation order is selected automatically
   according to the floating-point precision: $k=32$ (16 complex conjugate pole
   pairs) in quadruple precision (`__float128`), $k=16$ (8 complex conjugate
   pole pairs) in double precision, and $k=8$ (4 complex conjugate pole pairs)
   in single precision.

   The computation is performed in a single step by solving $k/2$ shifted linear
   systems $(At - \theta_i I)x_i = n_0$ and reconstructing the real solution as
   $y = 2 \mathrm{Re}\left( \sum_{i=1}^{k/2} \alpha_i x_i \right) + \alpha_0 n_0$.
   Options for the internal linear solver can be set using the `-mfn_cram_ksp_`
   and `-mfn_cram_pc_` prefixes.

   This method is tailored for real matrices whose spectrum lies on or near the
   negative real axis $(-\infty, 0]$ with real $t$ and $n_0$, and requires complex
   scalars.

   Level: beginner

.seealso: [](ch:mfn), `MFN`, `MFNType`, `MFNSetType()`, `MFNCRAMGetKSP()`
M*/
SLEPC_EXTERN PetscErrorCode MFNCreate_CRAM(MFN mfn)
{
#if PetscDefined(USE_COMPLEX)
  MFN_CRAM *ctx;
#endif

  PetscFunctionBegin;
  PetscCheck(PetscDefined(USE_COMPLEX),PetscObjectComm((PetscObject)mfn),PETSC_ERR_SUP,"The CRAM solver currently requires complex scalars");
#if PetscDefined(USE_COMPLEX)
  PetscCall(PetscNew(&ctx));
  mfn->data = (void*)ctx;

  mfn->ops->solve          = MFNSolve_CRAM;
  mfn->ops->setup          = MFNSetUp_CRAM;
  mfn->ops->setfromoptions = MFNSetFromOptions_CRAM;
  mfn->ops->view           = MFNView_CRAM;
  mfn->ops->reset          = MFNReset_CRAM;
  mfn->ops->destroy        = MFNDestroy_CRAM;

  PetscCall(PetscObjectComposeFunction((PetscObject)mfn,"MFNCRAMGetKSP_C",MFNCRAMGetKSP_CRAM));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}
