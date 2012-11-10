#include "FDTD1DDNG.h"
#include <cmath>
#include <iostream>
using namespace std;
CFDTD1DDNG::CFDTD1DDNG():
							// Simulation parameters.
							Size(4*1024),
							MaxTime(4*1024),
							PulseWidth(Size/8),
							td(PulseWidth),
							SourceLocation(10),
							SlabLeft(Size/3),
							SlabRight(2*Size/3),
							SnapshotInterval(16),
							e0((1e-9)/(36*PI)),
							u0((1e-7)*4*PI),
							dt(0.5e-11),
							dz(3e-3),
							Sc(c0*dt/dz),
							// Frequency, wavelength, wave number.
							l(PulseWidth*dz),
							f(c0/l),
							fmax(1/(2*dt)),
							w(2*PI*f),
							k0(w/c0),
							// Source choice.
							SourceChoice(1),
							fp(f),
							dr(PulseWidth*dt*2),
							// Data arrays.
							Ex(NULL), Dx(NULL), Hy(NULL), By(NULL),
							frame(0),
							// Incident and transmitted fields.
							Exi(NULL), Ext(NULL), Extt(NULL),
							x1(SlabLeft+1),
							// Refractive index.
							Z1(SlabLeft+50),
							z1(Z1*dz),
							Z2(SlabLeft+60),
							z2(Z2*dz),
							Exz1(NULL),
							Exz2(NULL),
							// Drude material parameters.
							einf(NULL), uinf(NULL), wpesq(NULL), wpmsq(NULL), ge(NULL), gm(NULL),
							// Auxiliary field scalars.
							ae0(NULL), ae(NULL), be(NULL), ce(NULL), de(NULL), ee(NULL),
							am0(NULL), am(NULL), bm(NULL), cm(NULL), dm(NULL), em(NULL),
							// Time indices.
							nf(2), n0(1), np(0),
							// Timer variables.
							tStart(0LL), tEnd(0LL)
{
}
unsigned long CFDTD1DDNG::SimSize()
{
	return sizeof(*this)+8*(30*Size+5*MaxTime);
}
void CFDTD1DDNG::AllocateMemoryCPU()
{
	// Field arrays.
	Ex = new double[Size*3];
	Dx = new double[Size*3];
	Hy = new double[Size*3];
	By = new double[Size*3];
	// Incident and transmitted fields.
	Exi = new double[MaxTime];
	Ext = new double[MaxTime];
	Extt = new double[MaxTime];
	// Refractive index.
	Exz1 = new double[MaxTime];
	Exz2 = new double[MaxTime];
	// Drude parameter arrays.
	einf = new double[Size];
	uinf = new double[Size];
	wpesq = new double[Size];
	wpmsq = new double[Size];
	ge = new double[Size];
	gm = new double[Size];
	// Auxiliary field scalars.
	ae0 = new double[Size];
	ae = new double[Size];
	be = new double[Size];
	ce = new double[Size];
	de = new double[Size];
	ee = new double[Size];
	am0 = new double[Size];
	am = new double[Size];
	bm = new double[Size];
	cm = new double[Size];
	dm = new double[Size];
	em = new double[Size];
}
void CFDTD1DDNG::InitialiseCPU()
{
	for (unsigned int i=0; i<Size*3; i++)
	{
		Ex[i] = 0.;
		Dx[i] = 0.;
		Hy[i] = 0.;
		By[i] = 0.;

		// Parameteric and auxiliary arrays.
		if (i<Size)
		{
			// Outside Slab.
			if (i<SlabLeft || i>SlabRight)
			{
				// Drude parameters.
				einf[i] = 1.;
				uinf[i] = 1.;
				wpesq[i] = 0.;
				wpmsq[i] = 0.;
				ge[i] = 0.;
				gm[i] = 0.;
			}
			// Inside Slab.
			else
			{
				// Drude parameters.
				einf[i] = 1.;
				uinf[i] = 1.;
				wpesq[i] = 2*pow(w, 2);
				wpmsq[i] = 2*pow(w, 2);
				ge[i] = w/32.;
				gm[i] = w/32.;
			}
			// Auxiliary scalars.
			ae0[i] = (4.*pow(dt,2))/(e0*(4.*einf[i]+pow(dt,2)*wpesq[i]+2.*dt*einf[i]*ge[i]));
			ae[i] = (1./pow(dt,2))*ae0[i];
			be[i] = (1./(2.*dt))*ge[i]*ae0[i];
			ce[i] = (e0/pow(dt,2))*einf[i]*ae0[i];
			de[i] = (-1.*e0/4.)*wpesq[i]*ae0[i];
			ee[i] = (1./(2.*dt))*e0*einf[i]*ge[i]*ae0[i];
			am0[i] = (4.*pow(dt,2))/(u0*(4.*uinf[i]+pow(dt,2)*wpmsq[i]+2.*dt*uinf[i]*gm[i]));;
			am[i] = (1./pow(dt,2))*am0[i];
			bm[i] = (1./(2.*dt))*gm[i]*am0[i];
			cm[i] = (u0/pow(dt,2))*uinf[i]*am0[i];
			dm[i] = (-1.*u0/4.)*wpmsq[i]*am0[i];
			em[i] = (1./(2.*dt))*u0*uinf[i]*gm[i]*am0[i];
		}
	}
	for (unsigned int i=0; i<MaxTime; i++)
	{
		Exi[i] = 0.;
		Ext[i] = 0.;
		Extt[i] = 0.;
		Exz1[i] = 0.;
		Exz2[i] = 0.;
	}
}

// Useful for indexing arrays.
#define Ex(i,n) Ex[(i)+Size*(n)]
#define Dx(i,n) Dx[(i)+Size*(n)]
#define Hy(i,n) Hy[(i)+Size*(n)]
#define By(i,n) By[(i)+Size*(n)]

int CFDTD1DDNG::DryRunCPU()
{
	for (unsigned int n=0; n<MaxTime; n++)
	{
		// Calculation of Hy using update difference equation for Hy. This is time step n.
		for (unsigned int i=0; i<Size-1; i++)
		{
			Hy(i,nf) = Hy(i,n0) + (Ex(i,n0)-Ex(i+1,n0))*dt/(u0*dz);
		}
		// ABC for Hy at i=Size-1.
		Hy(Size-1,nf) = Hy(Size-2,n0) + (Sc-1)/(Sc+1)*(Hy(Size-2,nf)-Hy(Size-1,n0));

		// Calculation of Ex using update difference equation for Ex. This is time step n+1/2.
		for (unsigned int i=1; i<Size; i++)
		{
			Ex(i,nf) = Ex(i,n0) + (Hy(i-1,nf)-Hy(i,nf))*dt/(e0*dz);
		}
		// ABC for Ex at i=0;
		Ex(0,nf) = Ex(1,n0) + (Sc-1)/(Sc+1)*(Ex(1,nf)-Ex(1,n0));

		// Source.
		if (SourceChoice == 1)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + exp(-1.*pow((n-td)/(PulseWidth/4.),2)) * Sc;
		}
		else if (SourceChoice == 2)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + sin(2.*PI*f*n*dt) * Sc;
		}
		else if (SourceChoice == 3)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + (1.-2.*pow(PI*fp*(n*dt-dr),2))*exp(-1.*pow(PI*fp*(n*dt-dr),2)) * Sc;
		}
		// Recording incident field.
		Exi[n] = Ex(x1,nf);

		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}
	return 0;
}
int CFDTD1DDNG::RunSimulationCPU()
{
	return 0;
}
// Timing.
void CFDTD1DDNG::StartTimer()
{
	tStart = GetTimeus64();
}
void CFDTD1DDNG::StopTimer()
{
	tEnd = GetTimeus64();
	cout << "Time taken = " << ((double)(tEnd-tStart))/(1000000.) << " seconds." << endl;
}
CFDTD1DDNG::~CFDTD1DDNG()
{
	// Field arrays.
	if (Ex != NULL) delete Ex;
	if (Dx != NULL) delete Dx;
	if (Hy != NULL) delete Hy;
	if (By != NULL) delete By;
	// Incident and transmitted fields.
	if (Exi != NULL) delete Exi;
	if (Ext != NULL) delete Ext;
	if (Extt != NULL) delete Extt;
	// Refractive index.
	if (Exz1 != NULL) delete Exz1;
	if (Exz2 != NULL) delete Exz2;
	// Drude parameter arrays.
	if (einf != NULL) delete einf;
	if (uinf != NULL) delete uinf;
	if (wpesq != NULL) delete wpesq;
	if (wpmsq != NULL) delete wpmsq;
	if (ge != NULL) delete ge;
	if (gm != NULL) delete gm;
	// Auxiliary field scalars.
	if (ae0 != NULL) delete ae0;
	if (ae != NULL) delete ae;
	if (be != NULL) delete be;
	if (ce != NULL) delete ce;
	if (de != NULL) delete de;
	if (ee != NULL) delete ee;
	if (am0 != NULL) delete am0;
	if (am != NULL) delete am;
	if (bm != NULL) delete bm;
	if (cm != NULL) delete cm;
	if (dm != NULL) delete dm;
	if (em != NULL) delete em;
}
