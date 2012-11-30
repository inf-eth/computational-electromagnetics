#ifndef FDTD1DDNG_H_
#define FDTD1DDNG_H_

// Constants.
#define c0	299792458.
#define PI	3.14159265358979323846
#define PRECISION double

#include <Timer.h>
#include <string>
#include <CL/cl.h>

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
	__int64 tDelta;
	bool tPaused;

		// ====================== Device arrays ======================
	// Data arrays.
	cl_mem d_Ex_;
	cl_mem d_Dx_;
	cl_mem d_Hy_;
	cl_mem d_By_;
	// Incident and transmitted fields.
	cl_mem d_Exi;
	cl_mem d_Ext;
	cl_mem d_Extt;
	// Refractive Index.
	cl_mem d_Exz1;
	cl_mem d_Exz2;
	// Drude material parameters.
	cl_mem d_einf;
	cl_mem d_uinf;
	cl_mem d_wpesq;
	cl_mem d_wpmsq;
	cl_mem d_ge;
	cl_mem d_gm;
	// Auxiliary field scalars.
	cl_mem d_ae0, d_ae, d_be, d_ce, d_de, d_ee;
	cl_mem d_am0, d_am, d_bm, d_cm, d_dm, d_em;
	// ===========================================================

	// ============ OPENCL related parameters ===========
	// OPENCL context/device/program
	cl_context context;
	cl_device_id *devices;
	cl_command_queue commandQueue;
	cl_program program;

	// Dry run and simulation kernel handles.
	cl_kernel DryRun_kernel_M;
	cl_kernel DryRun_kernel_E;
	cl_kernel Simulation_kernel_M;
	cl_kernel Simulation_kernel_E;
	// ==================================================

public:
	CFDTD1DDNG(unsigned int=4U*1024U, unsigned int=4U*4096U, unsigned int=10U, unsigned int=16U, unsigned int=1U);

	// Space calculations.
	unsigned long SimSize();
	unsigned long HDDSpace();

	// Memory allocation and initialisation.
	int AllocateMemoryCPU();
	int InitialiseCPU();
	int InitialiseCL();			// Search and allocate a device.
	int InitialiseExHyCPU();
	int AllocateMemoryGPU();
	int InitialiseCLKernelsGPU(); // Build/attach kernels to respective kernel functions and set arguments.
	int ReinitialiseExHyGPU();

	// Simulations.
	int DryRunCPU();
	int RunSimulationCPU(bool=true);
	int DryRunGPU();
	int RunSimulationGPU(bool=true);

	// Complete Runs on CPU and GPU.
	int CompleteRunCPU(bool=true);
	int CompleteRunGPU(bool=true);

	// Timing.
	void StartTimer();
	void StopTimer();
	void ResetTimer();
	PRECISION GetElapsedTime();

	std::string convertToString(const char * filename);

	int SafeCall(cl_int, const char []=NULL);

	int CleanupCPU();
	int CleanupCL();
	int CleanupGPU();
	~CFDTD1DDNG ();
};
#endif // #ifndef FDTD1DDNG_H_
