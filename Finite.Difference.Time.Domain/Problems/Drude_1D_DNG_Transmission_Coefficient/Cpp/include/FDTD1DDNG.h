#include "Timer.h"
// Constants.
#define c0	299792458.
#define PI	3.14159265358979323846

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

	// Choice of source.
	// 1. Gaussian pulse 2. Sine wave 3. Ricker wavelet
	const unsigned int SourceChoice;
	const double fp; // Ricker wavelet peak frequency.
	const double dr; // Ricker wavelet delay.

	// Data arrays.
	double *Ex;
	double *Dx;
	double *Hy;
	double *By;
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
	unsigned int n0, n1;

	// Timer variables.
	__int64 tStart;
	__int64 tEnd;

public:
	CFDTD1DDNG();
	unsigned long SimSize();
	void AllocateMemory();
	void Initialise();
	int RunSimulationCPU();
	// Timing.
	void StartTimer();
	void StopTimer();
	~CFDTD1DDNG();
};
