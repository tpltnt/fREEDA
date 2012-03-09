/*
 * dsp.h
 */
extern void Dsp_DoublevCopy(int n, doublev_t r, doublev_t x);
extern void Dsp_DoublevLastCopy(int n, int l, doublev_t r, doublev_t x);
extern void Dsp_Doublev2Dcxv(int n, dcxv_t r, doublev_t x);
extern void Dsp_TimeRepeatCopy(int n, int l, doublev_t r, doublev_t x);
extern void Dsp_Double2DoublevCopy(int n, doublev_t r, double x );
extern void Dsp_DoublevAdd(int n, doublev_t r, doublev_t x, doublev_t y);
extern void Dsp_DoublevSub(int n, doublev_t r, doublev_t x, doublev_t y);
extern void Dsp_DoublevMult(int n, doublev_t r, doublev_t x, doublev_t y);
extern void Dsp_DoublevDiv(int n, doublev_t r, doublev_t x, doublev_t y);
extern void Dsp_DoublevRepeatCopy(int n, int l, doublev_t r, doublev_t x);
extern void Dsp_DoublevSclAdd(int n, doublev_t r, doublev_t x, double s);
extern void Dsp_DoublevSclSub(int n, doublev_t r, doublev_t x, double s);
extern void Dsp_DoublevSclMult(int n, doublev_t r, doublev_t x, double s);
extern void Dsp_DoublevSclDiv(int n, doublev_t r, doublev_t x, double s);

extern void Dsp_DoublevRecip(int n, doublev_t r, doublev_t x);
extern void Dsp_DoublevNeg(int n, doublev_t r, doublev_t x);
extern void Dsp_DoublevDB(int n, doublev_t r, doublev_t x);
extern void Dsp_DoublevDB10(int n, doublev_t r, doublev_t x);
extern void Dsp_DoublevRad2Deg(int n, doublev_t r, doublev_t x);
extern void Dsp_DoublevDeg2Rad(int n, doublev_t r, doublev_t x);
extern void Dsp_DoublevAbs(int n, doublev_t r, doublev_t x);
extern void Dsp_DoublevMinLmt(int n, doublev_t r, doublev_t x, double minlmt);
extern void Dsp_DoublevMaxLmt(int n, doublev_t r, doublev_t x, double maxlmt);
extern void Dsp_DoublevLog(int n, doublev_t r, doublev_t x);
extern void Dsp_DoublevLog10(int n, doublev_t r, doublev_t x);
extern void Dsp_DoublevExp(int n, doublev_t r, doublev_t x);
extern void Dsp_DoublevExp10(int n, doublev_t r, doublev_t x);

extern void Dsp_DoublevDiff(int n, doublev_t diff, doublev_t x);
extern void Dsp_DoublevDeriv(int n, doublev_t deriv, doublev_t x);
extern void Dsp_DoublevSum(int n, doublev_t sum, doublev_t x);
extern void Dsp_DoublevInteg(int n, doublev_t integ, doublev_t x);

extern void Dsp_RealInverseFFT(int realPts, doublev_t x, dcxv_t packedSpect);
extern void Dsp_RealForwardFFT(int realPts, dcxv_t packedSpect, doublev_t x);
extern void Dsp_DcxvUnpackSpect(int packedPts, dcxv_t spect);
extern void Dsp_DcxvPackSpect(int packedPts, dcxv_t spect);

extern void Dsp_RealCircConv(int n, doublev_t result, doublev_t s, 
			     doublev_t r);
extern void Dsp_RealSpectConv(int n, doublev_t result, doublev_t s, 
			      dcxv_t packedSpect);

extern void Dsp_DcxvCopy(int n, dcxv_t r, dcxv_t x);
extern void Dsp_DcxvLastCopy(int n, int l, dcxv_t r, dcxv_t x);
extern void Dsp_DcxvRepeatCopy(int n, int l, dcxv_t r, dcxv_t x);
extern void Dsp_Dcx2DcxvCopy(int n, dcxv_t r, dcx_t x);

extern void Dsp_DcxvAdd(int n, dcxv_t r, dcxv_t x, dcxv_t y);
extern void Dsp_DcxvSub(int n, dcxv_t r, dcxv_t x, dcxv_t y);
extern void Dsp_DcxvMult(int n, dcxv_t r, dcxv_t x, dcxv_t y);
extern void Dsp_DcxvDiv(int n, dcxv_t r, dcxv_t x, dcxv_t y);
extern void Dsp_DcxvNeg(int n, dcxv_t r, dcxv_t x);
extern void Dsp_DcxvRecip(int n, dcxv_t r, dcxv_t x);
extern void Dsp_DcxvConj(int n, dcxv_t r, dcxv_t x);

extern void Dsp_DcxvScale(int n, dcxv_t r, double scalar, dcxv_t x);

extern void Dsp_DcxvReal(int n, doublev_t r, dcxv_t x);
extern void Dsp_DcxvImag(int n, doublev_t r, dcxv_t x);
extern void Dsp_DcxvMag(int n, doublev_t r, dcxv_t x);
extern void Dsp_DcxvDB(int n, doublev_t r, dcxv_t x);
extern void Dsp_DcxvDB10(int n, doublev_t r, dcxv_t x);
extern void Dsp_DcxvPrinPhase(int n, doublev_t phase, dcxv_t x);
extern void Dsp_DcxvContPhase(int n, doublev_t phase, dcxv_t x);

extern void Dsp_DoublevFFT(int order, dcxv_t r, doublev_t x);
extern void Dsp_DoublevInvFFT(int order, doublev_t r, dcxv_t x);

extern void Dsp_LPBwthRsp(int n, dcxv_t response, doublev_t frq,
			  double corner, int order);
