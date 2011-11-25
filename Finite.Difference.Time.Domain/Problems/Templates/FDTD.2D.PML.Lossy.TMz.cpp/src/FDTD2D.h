#ifdef WIN32
#include <fstream>
#else
#include <fcntl.h>
#endif
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <ctime>

typedef unsigned int uint;

class CFDTD2D
{
private:

	// Generic simulation parameters.
	const uint I;				// Width.
	const uint J;				// Height.
	const double c;				// Speed of light.
	const double delta;			// dx and dy.
	const double dx;			// dx if being used.
	const double dy;			// dy if being used.
	const double dtscalar;		// dt scale factor. dt will be divided by this factor.
	const double dt;			// dt.
	const uint PMLw;			// Width of PML layer.

	const uint NMax;			// Maximum n

	// Constants.
	const double f;				// frequency
	const double pi;			// pi
	const double e0;			// epsilon nought
	const double u0;			// mu nought

	// miscellenaeous
	const double Two_pi_f_deltat;
	const double NHW;			// Half wave cycle.
	const uint Js;				// J position of plane wave front.
	const uint Is;				// I position of plane wave front.
	uint n0, n1, n2;			// past/present/future time indices.
	const uint tResolution;		// Snapshots will be saved after this much time.
	const uint xResolution;		// Resolution of plotted field will be divided by this factor.
	const uint yResolution;

	// TMz parameters.
	const uint IHx;
	const uint JHx;
	const uint IHy;
	const uint JHy;
	const uint IEz;
	const uint JEz;

	// Geometry parameters.
	const uint XCenter;
	const uint YCenter;
	const double ra;
	const double rb;
	const double imp0;			// Impedence of free space

	// =========== Data Arrays ==========
	double *Hx;					// Hx, magnetic field.
	double *Bx;
	double *Hy;					// Hy, magnetic field.
	double *By;

	double *Ez;					// Ez, electric field.
	double *Dz;

	double *Dzx;
	double *Dzy;
	double *EzSnapshots;		// For storing Ez snapshots.

	// ========= Field specific arrays =======

	// Permeability and permittivity.
	double *urHx;
	double *urHy;
	double *erEz;

	// Magnetic and electric conductances.
	double *smHx;
	double *smHy;
	double *sEz;

	// s*dt/(2*er) and sm*dt/(2*ur)
	double *Sc;
	double *ScmHx;
	double *ScmHy;

	// PML conductance arrays.
	double *sex;	// sigma ex
	double *sey;	// sigma ey
	double *smx;	// sigma mx
	double *smy;	// sigma my

	double *Scsx;
	double *Scsy;
	double *ScmsmxHy;
	double *ScmsmyHx;

	// Timing.
	clock_t tStart;
	clock_t tEnd;
	clock_t tElapsed;

public:
	CFDTD2D ();
	void Initialize ();			// Initialize with default parameters.
	void RunSimulation (bool=true);

	// Timing.
	inline void StartClock () { tStart = clock(); }
	inline void StopClock () { tEnd = clock(); }
	inline double GetElapsedTime () { return (double)(tEnd-tStart)/CLOCKS_PER_SEC; }

	~CFDTD2D ();
};
