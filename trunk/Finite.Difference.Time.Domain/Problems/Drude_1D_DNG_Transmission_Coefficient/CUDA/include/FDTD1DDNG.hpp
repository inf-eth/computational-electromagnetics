#ifndef FDTD1DDNG_H_
#define FDTD1DDNG_H_

// Constants.
#define c0	299792458.
#define PI	3.14159265358979323846
#define PRECISION double

#include <Timer.h>

class CFDTD1DDNG
{
private:
	// Simulation parameters.
	const unsigned int Size;
	const unsigned int MaxTime;
	const unsigned int PulseWidth;
	const unsigned int td;
	const unsigned int SourceLocation;
	const unsigned int SlabLeft;
	const unsigned int SlabRight;
	const unsigned int SnapshotInterval;

	// Choice of source.
	// 1. Gaussian pulse 2. Sine wave 3. Ricker wavelet
	const unsigned int SourceChoice;

	const PRECISION e0;
	const PRECISION u0;
	const PRECISION dt;
	const PRECISION dz;
	const PRECISION Sc;

	// Frequency, wavelength, wave number.
	const PRECISION l;
	const PRECISION f;
	const PRECISION fmax;
	const PRECISION w;
	const PRECISION k0;
	const PRECISION fp; // Ricker wavelet peak frequency.
	const PRECISION dr; // Ricker wavelet delay.

	// Data arrays.
	PRECISION *Ex_;
	PRECISION *Dx_;
	PRECISION *Hy_;
	PRECISION *By_;
	unsigned int frame;

	// Incident and transmitted fields.
	PRECISION *Exi;
	PRECISION *Ext;
	PRECISION *Extt;
	const unsigned int x1;

	// Refractive index.
	const unsigned int Z1;
	const PRECISION z1;
	const unsigned int Z2;
	const PRECISION z2;
	PRECISION *Exz1;
	PRECISION *Exz2;

	// Drude material parameters.
	PRECISION *einf;
	PRECISION *uinf;
	PRECISION *wpesq;
	PRECISION *wpmsq;
	PRECISION *ge;
	PRECISION *gm;

	// Auxiliary field scalars.
	PRECISION *ae0, *ae, *be, *ce, *de, *ee;
	PRECISION *am0, *am, *bm, *cm, *dm, *em;

	// Time indices.
	unsigned int nf, n0, np;

	// Timer variables.
	__int64 tStart;
	__int64 tEnd;

	// ====================== Device arrays ======================
	// Data arrays.
	PRECISION *d_Ex_;
	PRECISION *d_Dx_;
	PRECISION *d_Hy_;
	PRECISION *d_By_;
	// Incident and transmitted fields.
	PRECISION *d_Exi;
	PRECISION *d_Ext;
	PRECISION *d_Extt;
	// Refractive Index.
	PRECISION *d_Exz1;
	PRECISION *d_Exz2;
	// Drude material parameters.
	PRECISION *d_einf;
	PRECISION *d_uinf;
	PRECISION *d_wpesq;
	PRECISION *d_wpmsq;
	PRECISION *d_ge;
	PRECISION *d_gm;
	// Auxiliary field scalars.
	PRECISION *d_ae0, *d_ae, *d_be, *d_ce, *d_de, *d_ee;
	PRECISION *d_am0, *d_am, *d_bm, *d_cm, *d_dm, *d_em;
	// ===========================================================

public:
	CFDTD1DDNG(unsigned int=4U*1024U, unsigned int=10U, unsigned int=16U, unsigned int=1U);

	// Space calculations.
	unsigned long SimSize();
	unsigned long HDDSpace();

	// Memory allocation and initialisation.
	int AllocateMemoryCPU();
	int InitialiseCPU();
	int InitialiseExHyCPU();
	int AllocateMemoryGPU();//int initializeFDTD1DDNGKernel();
	int CopyDataCPUtoGPU();
	int CopyExHyCPUtoGPU();

	// Simulations.
	int DryRunCPU();
	int RunSimulationCPU(bool=true);
	int DryRunGPU();
	int RunSimulationGPU(bool=true);//int runFDTD1DDNGKernels(bool=true);

	// Complete Runs on CPU and GPU.
	int CompleteRunCPU(bool=true);
	int CompleteRunGPU(bool=true);

	// Timing.
	void StartTimer();
	void StopTimer();
	PRECISION GetElapsedTime();

	int SafeCall(int, const char[]=NULL);

	int CleanupCPU();
	int CleanupGPU();
	~CFDTD1DDNG();
};
#endif // #ifndef FDTD1DDNG_H_
