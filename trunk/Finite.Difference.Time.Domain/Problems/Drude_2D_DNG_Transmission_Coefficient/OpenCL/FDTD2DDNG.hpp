#ifndef FDTD2DDNG_H_
#define FDTD2DDNG_H_

#define PRECISION double

#include <Timer.h>
#include <string>
#include <CL/cl.h>

class CFDTD2DDNG
{
private:
	// Simulation parameters.
	const unsigned int I;
	const unsigned int J;
	const unsigned int PMLw;
	const unsigned int SlabLeft;
	const unsigned int SlabRight;
	const unsigned int MaxTime;
	const unsigned int PulseWidth;
	const unsigned int td;
	const unsigned int SnapshotResolution;
	const unsigned int SnapshotInterval;
	// Choice of source.
	// 1. Gaussian pulse 2. Sine wave 3. Ricker wavelet
	const unsigned int SourceChoice;
	const unsigned int SourcePlane;	// 0. Omni, 1. Plane wave.
	const unsigned int SourceLocationX;
	const unsigned int SourceLocationY;

	const PRECISION c;
	const PRECISION pi;
	const PRECISION e0;
	const PRECISION u0;
	const PRECISION dt;
	const PRECISION delta;
	const PRECISION Sc;

	// Frequency, wavelength, wave number.
	const PRECISION l;
	const PRECISION f;
	const PRECISION fmax;
	const PRECISION w;
	const PRECISION k0;
	const PRECISION fp; // Ricker wavelet peak frequency.
	const PRECISION dr; // Ricker wavelet delay.

	// Data array sizes.
	const unsigned int IEz, JEz;
	const unsigned int IHx, JHx;
	const unsigned int IHy, JHy;

	// Data arrays.
	PRECISION *Ez_;
	PRECISION *Dz_;
	PRECISION *Hx_;
	PRECISION *Bx_;
	PRECISION *Hy_;
	PRECISION *By_;

	// Incident and transmitted fields.
	PRECISION *Ezi;
	PRECISION *Ezt;
	PRECISION *Eztt;
	const unsigned int x1;

	// Refractive index.
	const unsigned int Y1;
	const PRECISION y1;
	const unsigned int Y2;
	const PRECISION y2;
	PRECISION *Ezy1;
	PRECISION *Ezy2;

	// Drude material parameters.
	PRECISION *einf_;
	PRECISION *uinf_;
	PRECISION *wpesq_;
	PRECISION *wpmsq_;
	PRECISION *ge_;
	PRECISION *gm_;

	// Auxiliary field scalars.
	PRECISION *ae0_, *ae_, *be_, *ce_, *de_, *ee_;
	PRECISION *am0_, *am_, *bm_, *cm_, *dm_, *em_;

	// PML arrays.
	PRECISION *PsiEzX_, *PsiEzY_;
	PRECISION *PsiHyX_, *PsiHxY_;

	// PML parameters.
	const PRECISION kappe, kappm;
	const PRECISION kappex, kappey, kappmx, kappmy;
	const PRECISION aex, aey, amx, amy;
	const PRECISION sigex, sigey, sigmx, sigmy;
	const PRECISION bex, bey, bmx, bmy;
	const PRECISION Cex, Cey, Cmx, Cmy;

	// Snapshot frame number.
	unsigned int frame;

	// Time indices.
	unsigned int nf, n0, np;

	// Timer variables.
	__int64 tStart;
	__int64 tEnd;

		// ====================== Device arrays ======================
	// Data arrays.
	cl_mem d_Ez_;
	cl_mem d_Dz_;
	cl_mem d_Hx_;
	cl_mem d_Bx_;
	cl_mem d_Hy_;
	cl_mem d_By_;
	// Incident and transmitted fields.
	cl_mem d_Ezi;
	cl_mem d_Ezt;
	cl_mem d_Eztt;
	// Refractive Index.
	cl_mem d_Ezy1;
	cl_mem d_Ezy2;
	// Drude material parameters.
	cl_mem d_einf_;
	cl_mem d_uinf_;
	cl_mem d_wpesq_;
	cl_mem d_wpmsq_;
	cl_mem d_ge_;
	cl_mem d_gm_;
	// Auxiliary field scalars.
	cl_mem d_ae0_, d_ae_, d_be_, d_ce_, d_de_, d_ee_;
	cl_mem d_am0_, d_am_, d_bm_, d_cm_, d_dm_, d_em_;
	// PML arrays.
	cl_mem d_PsiEzX_, d_PsiEzY_;
	cl_mem d_PsiHyX_, d_PsiHxY_;
	// ===========================================================

	// ========== OPENCL related parameters =========
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
	// ==============================================

public:
	CFDTD2DDNG(
				unsigned int=512,	// I
				unsigned int=512,	// J
				unsigned int=64,	// PMLw
				unsigned int=4*512,	// MaxTime
				unsigned int=1,		// Snapshot resolution
				unsigned int=10,	// Snapshot interval
				unsigned int=1,		// Source choice
				unsigned int=1,		// Source is plane wave?
				unsigned int=256,	// Source location X
				unsigned int=64+5);	// Source location Y

	// Space calculations.
	unsigned long SimSize();
	unsigned long HDDSpace();

	// Memory allocation and initialisation.
	int AllocateMemoryCPU();
	int InitialiseCPU();
	int InitialiseCL(bool=false);		// Search and allocate a device.
	int InitialiseForSimulationCPU();
	int AllocateMemoryGPU();
	int InitialiseCLKernelsGPU();		// Build/attach kernels to respective kernel functions and set arguments.
	int InitialiseForSimulationGPU();

	// Simulations.
	int DryRunCPU();
	int RunSimulationCPU(bool=true);
	int DryRunGPU();
	int RunSimulationGPU(bool=true);

	// Complete Runs on CPU and GPU.
	int CompleteRunCPU(bool=true);
	int CompleteRunGPU(bool=true, bool=false);

	// Timing.
	void StartTimer();
	void StopTimer();
	PRECISION GetElapsedTime();

	std::string convertToString(const char * filename);

	int SafeCall(int, const char []=NULL);

	int CleanupCPU();
	int CleanupCL();
	int CleanupGPU();
	~CFDTD2DDNG ();
};
#endif  /* #ifndef FDTD1DDNG_H_ */
