#include "FDTD1DDNG.h"
#include <iostream>
using namespace std;
CFDTD1DDNG::CFDTD1DDNG():
							// Simulation parameters.
							Size(1024),
							MaxTime(1024),
							PulseWidth(Size/8),
							td(PulseWidth),
							SourceLocation(10),
							SlabLeft(Size/3),
							SlabRight(2*Size/3),
							SnapshotInterval(32),
							e0((1e-9)/(36*pi)),
							u0((1e-7)*4*pi),
							dt(0.5e-11),
							dz(3e-3),
							Sc(c0*dt/dz),
							// Frequency, wavelength, wave number.
							l(PulseWidth*dz),
							f(c0/l),
							fmax(1/(2*dt)),
							w(2*pi*f),
							k0(w/c0),
							// Source choice.
							SourceChoice(1),
							fp(f),
							dr(PulseWidth*dt*2),
							// Data arrays.
							Ex(NULL), Dx(NULL), Hy(NULL), By(NULL),
							ExSnapshots(NULL),
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
							n0(0), n1(1)
{
}
CFDTD1DDNG::~CFDTD1DDNG()
{
}
