#ifndef FDTD1DDNG_H_
#define FDTD1DDNG_H_

// Constants.
#define c0	299792458.
#define PI	3.14159265358979323846

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

	const double e0;
	const double u0;
	const double dt;
	const double dz;
	const double Sc;

	// Frequency, wavelength, wave number.
	const double l;
	const double f;
	const double fmax;
	const double w;
	const double k0;
	const double fp; // Ricker wavelet peak frequency.
	const double dr; // Ricker wavelet delay.

	// Data arrays.
	double *Ex_;
	double *Dx_;
	double *Hy_;
	double *By_;
	unsigned int frame;

	// Incident and transmitted fields.
	double *Exi;
	double *Ext;
	double *Extt;
	const unsigned int x1;

	// Refractive index.
	const unsigned int Z1;
	const double z1;
	const unsigned int Z2;
	const double z2;
	double *Exz1;
	double *Exz2;

	// Drude material parameters.
	double *einf;
	double *uinf;
	double *wpesq;
	double *wpmsq;
	double *ge;
	double *gm;

	// Auxiliary field scalars.
	double *ae0, *ae, *be, *ce, *de, *ee;
	double *am0, *am, *bm, *cm, *dm, *em;

	// Time indices.
	unsigned int nf, n0, np;

	// Timer variables.
	__int64 tStart;
	__int64 tEnd;

public:
	CFDTD1DDNG(unsigned int=4U*1024U, unsigned int=10U, unsigned int=16U, unsigned int=1U);

	// Space calculations.
	unsigned long SimSize();
	unsigned long HDDSpace();

	// Initialisation and memory allocation.
	void AllocateMemoryCPU();
	void InitialiseCPU();
	void InitialiseExHyCPU();

	// Simulations.
	int DryRunCPU();
	int RunSimulationCPU(bool=true);

	// Timing.
	void StartTimer();
	void StopTimer();
	double GetElapsedTime();

	~CFDTD1DDNG();
};
#endif // #ifndef FDTD1DDNG_H_
