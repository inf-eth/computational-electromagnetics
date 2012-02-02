#ifndef FDTD2D_H_
#define FDTD2D_H_

#ifndef WIN32
#include <fcntl.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <string>

#include <FDTD2D_Kernels.cu>
#include <cutil.h>

class CFDTD2D
{
private:

	// Generic simulation parameters.
	const unsigned int I;				// Width.
	const unsigned int J;				// Height.
	const float c;				// Speed of light.
	const float delta;			// dx and dy.
	const float dx;			// dx if being used.
	const float dy;			// dy if being used.
	const float dtscalar;		// dt scale factor. dt will be divided by this factor.
	const float dt;			// dt.
	const unsigned int PMLw;			// Width of PML layer.

	const unsigned int NMax;			// Maximum n

	// Constants.
	const float f;				// frequency
	const float pi;			// pi
	const float e0;			// epsilon nought
	const float u0;			// mu nought

	// miscellenaeous
	const float Two_pi_f_deltat;
	const unsigned int NHW;			// Half wave cycle.
	const unsigned int Js;				// J position of plane wave front.
	const unsigned int Is;				// I position of plane wave front.
	unsigned int n, n0, n1, n2;			// past/present/future time indices.
	const unsigned int tResolution;		// Snapshots will be saved after this much time.
	const unsigned int xResolution;		// Resolution of plotted field will be divided by this factor.
	const unsigned int yResolution;

	// TMz parameters.
	const unsigned int IHx;
	const unsigned int JHx;
	const unsigned int IHy;
	const unsigned int JHy;
	const unsigned int IEz;
	const unsigned int JEz;

	// Geometry parameters.
	const unsigned int XCenter;
	const unsigned int YCenter;
	const float ra;
	const float rb;
	const float imp0;			// Impedence of free space

	// =========== Data Arrays ==========
	float *Hx;					// Hx, magnetic field.
	float *Bx;
	float *Hy;					// Hy, magnetic field.
	float *By;

	float *Ez;					// Ez, electric field.
	float *Dz;

	float *Dzx;
	float *Dzy;
	float *EzSnapshots;		// For storing Ez snapshots.

	// ========= Field specific arrays =======

	// Permeability and permittivity.
	float *urHx;
	float *urHy;
	float *erEz;

	// Magnetic and electric conductances.
	float *smHx;
	float *smHy;
	float *sEz;

	// s*dt/(2*er) and sm*dt/(2*ur)
	float *Sc;
	float *ScmHx;
	float *ScmHy;

	// PML conductance arrays.
	float *sex;	// sigma ex
	float *sey;	// sigma ey
	float *smx;	// sigma mx
	float *smy;	// sigma my

	float *Scsx;
	float *Scsy;
	float *ScmsmxHy;
	float *ScmsmyHx;

	// Timing.
	clock_t tStart;
	clock_t tEnd;
	clock_t tElapsed;

	// === Device arrays ===
	// Fields
	float* d_Hx;
	float* d_Bx;
	float* d_Hy;
	float* d_By;
	float* d_Ez;
	float* d_Dz;
	float* d_Dzx;
	float* d_Dzy;

	// Permeability, permittivity and conductance.
	float* d_urHx;
	float* d_urHy;
	float* d_erEz;

	float* d_ScmHx;
	float* d_ScmHy;
	float* d_Sc;
	// =============================

	bool cpu;		// Should simulation run on CPU or GPU?
	unsigned int flagHalf;

public:
	CFDTD2D ();
	int Initialize ();						// Initialize with default parameters.
	int initializeFDTD2DKernel ();
	int runFDTD2DKernels (bool=true);
	int RunSimulationCPU (bool=true);		// Run simulation on CPU single-threaded.

	void DisplaySimulationParameters ();

	// Timing.
	inline void StartClock () { tStart = clock(); }
	inline void StopClock () { tEnd = clock(); }
	inline float GetElapsedTime () { return (float)(tEnd-tStart)/CLOCKS_PER_SEC; }

	int Cleanup ();
	~CFDTD2D ();
};


#endif  /* #ifndef FDTD2D_H_ */
