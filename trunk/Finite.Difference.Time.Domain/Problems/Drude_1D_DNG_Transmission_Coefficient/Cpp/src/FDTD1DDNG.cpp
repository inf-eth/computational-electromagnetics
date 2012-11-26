// Useful for indexing arrays.
#define Ex(i,n) Ex_[(i)+Size*(n)]
#define Dx(i,n) Dx_[(i)+Size*(n)]
#define Hy(i,n) Hy_[(i)+Size*(n)]
#define By(i,n) By_[(i)+Size*(n)]

#include <FDTD1DDNG.h>
#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;

CFDTD1DDNG::CFDTD1DDNG(unsigned int pSize, unsigned int pSourceLocation, unsigned int pSnapshotInterval, unsigned int pSourceChoice):
							// Simulation parameters.
							Size(pSize),
							MaxTime(4*Size),
							PulseWidth(Size/8),
							td(PulseWidth),
							SourceLocation(pSourceLocation),
							SlabLeft(Size/3),
							SlabRight(2*Size/3),
							SnapshotInterval(pSnapshotInterval),
							SourceChoice(pSourceChoice),
							e0((1e-9)/(36.*PI)),
							u0((1e-7)*4.*PI),
							dt(0.5e-11),
							dz(3e-3),
							Sc(c0*dt/dz),
							// Frequency, wavelength, wave number.
							l(PulseWidth*dz),
							f(c0/l),
							fmax(1/(2*dt)),
							w(2*PI*f),
							k0(w/c0),
							fp(f),
							dr(PulseWidth*dt*2),
							// Data arrays.
							Ex_(NULL), Dx_(NULL), Hy_(NULL), By_(NULL),
							frame(0),
							// Incident and transmitted fields.
							Exi(NULL), Ext(NULL), Extt(NULL),
							x1(SlabLeft+1),
							// Refractive index.
							Z1(SlabLeft+50),
							z1((double)Z1*dz),
							Z2(SlabLeft+60),
							z2((double)Z2*dz),
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
	// Creating directory and writing simulation parameters.
#if defined __linux__ || defined __CYGWIN__
	system("if test -d FieldData; then rm -rf FieldData; fi");
	system("mkdir -p FieldData");
#else
	system("IF exist FieldData (rd FieldData /s /q)");
	system("mkdir FieldData");
#endif
	fstream parametersfile;
	parametersfile.open("FieldData/Parameters.smp", std::ios::out|std::ios::binary);
	parametersfile.write((char*)&dt, sizeof(double));
	parametersfile.write((char*)&k0, sizeof(double));
	parametersfile.write((char*)&z1, sizeof(double));
	parametersfile.write((char*)&z2, sizeof(double));
	parametersfile.write((char*)&Size, sizeof(unsigned int));
	parametersfile.write((char*)&MaxTime, sizeof(unsigned int));
	parametersfile.write((char*)&SnapshotInterval, sizeof(unsigned int));
	parametersfile.write((char*)&SlabLeft, sizeof(unsigned int));
	parametersfile.write((char*)&SlabRight, sizeof(unsigned int));
	parametersfile.close();

	// Printing simulation parameters.
	cout << "Size      = " << Size << endl;
	cout << "MaxTime   = " << MaxTime << endl;
	cout << "frequency = " << f << " Hz (" << f/1e9 << " GHz)" << endl;
	cout << "fmax      = " << fmax << " Hz (" << fmax/1e9 << " GHz)" << endl;
	cout << "Sc        = " << Sc << endl;
	cout << "Slab left = " << SlabLeft << endl;
	cout << "Slab right= " << SlabRight << endl;
}
unsigned long CFDTD1DDNG::SimSize()
{
	return (unsigned long)sizeof(*this)+8UL*(30UL*(unsigned long)Size+5UL*(unsigned long)MaxTime);
}
unsigned long CFDTD1DDNG::HDDSpace()
{
	return 8UL*(5UL*(unsigned long)MaxTime+(unsigned long)Size*((unsigned long)MaxTime/(unsigned long)SnapshotInterval+1UL));
}
void CFDTD1DDNG::AllocateMemoryCPU()
{
	// Field arrays.
	Ex_ = new double[Size*3];
	Dx_ = new double[Size*3];
	Hy_ = new double[Size*3];
	By_ = new double[Size*3];
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
	for (unsigned int i=0; i<3*Size; i++)
	{
		Ex_[i] = 0.;
		Dx_[i] = 0.;
		Hy_[i] = 0.;
		By_[i] = 0.;

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
				wpesq[i] = 2.*pow(w, 2);
				wpmsq[i] = 2.*pow(w, 2);
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
void CFDTD1DDNG::InitialiseExHyCPU()
{
	for (unsigned int i=0; i<3*Size; i++)
	{
		Ex_[i] = 0.;
		Hy_[i] = 0.;
	}
}
int CFDTD1DDNG::DryRunCPU()
{
	cout << "Dry run (CPU) started..." << endl;
	for (unsigned int n=0; n<MaxTime; n++)
	{
		if (n%SnapshotInterval == 0)
			cout << "\r" << setprecision(4) << (float)n*100/(MaxTime-1) << "%  " << flush;

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
		Ex(0,nf) = Ex(1,n0) + (Sc-1)/(Sc+1)*(Ex(1,nf)-Ex(0,n0));

		// Source.
		if (SourceChoice == 1)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + exp(-1.*pow(((double)n-(double)td)/((double)PulseWidth/4.),2)) * Sc;
		}
		else if (SourceChoice == 2)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + sin(2.*PI*f*(double)n*dt) * Sc;
		}
		else if (SourceChoice == 3)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + (1.-2.*pow(PI*fp*((double)n*dt-dr),2))*exp(-1.*pow(PI*fp*((double)n*dt-dr),2)) * Sc;
		}
		// Recording incident field.
		Exi[n] = Ex(x1,nf);

		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}
	cout << endl << "Dry run (CPU) completed!" << endl;
	return 0;
}
int CFDTD1DDNG::RunSimulationCPU(bool SaveFields)
{
	stringstream framestream;
	string basename = "FieldData/Ex";
	string filename;
	fstream snapshot;

	cout << "Simulation (CPU) started..." << endl;
	for (unsigned int n=0; n<MaxTime; n++)
	{
		if (n%SnapshotInterval == 0)
			cout << "\r" << setprecision(4) << (float)n*100/(MaxTime-1) << "%  " << flush;

		// Calculation of By using update difference equation for Hy. This is time step n.
		for (unsigned int i=0; i<Size-1; i++)
		{
			By(i,nf) = By(i,n0) + (Ex(i,n0)-Ex(i+1,n0))*dt/dz;
			Hy(i,nf) = am[i]*(By(i,nf)-2*By(i,n0)+By(i,np)) + bm[i]*(By(i,nf)-By(i,np)) + cm[i]*(2*Hy(i,n0)-Hy(i,np)) + dm[i]*(2*Hy(i,n0)+Hy(i,np)) + em[i]*(Hy(i,np));
		}
		// ABC for Hy at i=Size-1.
		Hy(Size-1,nf) = Hy(Size-2,n0) + (Sc-1)/(Sc+1)*(Hy(Size-2,nf)-Hy(Size-1,n0));
		By(Size-1,nf) = u0*Hy(Size-1,nf);

		// Calculation of Ex using update difference equation for Ex. This is time step n+1/2.
		for (unsigned int i=1; i<Size; i++)
		{
			Dx(i,nf) = Dx(i,n0) + (Hy(i-1,nf)-Hy(i,nf))*dt/dz;
			Ex(i,nf) = ae[i]*(Dx(i,nf)-2*Dx(i,n0)+Dx(i,np)) + be[i]*(Dx(i,nf)-Dx(i,np)) + ce[i]*(2*Ex(i,n0)-Ex(i,np)) + de[i]*(2*Ex(i,n0)+Ex(i,np)) + ee[i]*(Ex(i,np));
		}
		// ABC for Ex at i=0;
		Ex(0,nf) = Ex(1,n0) + (Sc-1)/(Sc+1)*(Ex(1,nf)-Ex(0,n0));
		Dx(0,nf) = e0*Ex(0,nf);

		// Source.
		if (SourceChoice == 1)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + exp(-1.*pow(((double)n-(double)td)/((double)PulseWidth/4.),2)) * Sc;
		}
		else if (SourceChoice == 2)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + sin(2.*PI*f*(double)n*dt) * Sc;
		}
		else if (SourceChoice == 3)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + (1.-2.*pow(PI*fp*((double)n*dt-dr),2))*exp(-1.*pow(PI*fp*((double)n*dt-dr),2)) * Sc;
		}
		Dx(SourceLocation,nf) = e0*Ex(SourceLocation,nf);
		// Recording transmitted fields.
		Ext[n] = Ex(x1,nf);
		Extt[n] = Ex(SlabRight+10,nf);
		// Fields for refractive index.
		Exz1[n] = Ex(Z1, nf);
		Exz2[n] = Ex(Z2, nf);

		// Saving electric field snapshot.
		if (n%SnapshotInterval == 0 && SaveFields == true)
		{
			// Write E-field to file.
			framestream.str(std::string());			// Clearing stringstream contents.
			framestream << ++frame;
			filename = basename + framestream.str() + ".fdt";
			snapshot.open(filename.c_str(), std::ios::out|std::ios::binary);
			snapshot.write((char*)&(Ex(0,nf)), sizeof(double)*Size);
			snapshot.close();
		}
		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}
	// Saving electric field snapshot.
	if (SaveFields == true)
	{
		fstream parametersfile;
		parametersfile.open("FieldData/Parameters.smp", std::ios::out|std::ios::binary|std::ios::app);
		parametersfile.write((char*)&(frame), sizeof(unsigned int));
		parametersfile.close();
		// Write saved fields to files.
		snapshot.open("FieldData/Exi.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Exi, sizeof(double)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ext.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Ext, sizeof(double)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Extt.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Extt, sizeof(double)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Exz1.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Exz1, sizeof(double)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Exz2.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Exz2, sizeof(double)*MaxTime);
		snapshot.close();
	}
	cout << endl << "Simulation (CPU) completed!" << endl;
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
}
double CFDTD1DDNG::GetElapsedTime()
{
	return ((double)(tEnd-tStart))/(1000000.);
}
CFDTD1DDNG::~CFDTD1DDNG()
{
	// Field arrays.
	if (Ex_ != NULL) delete[] Ex_;
	if (Dx_ != NULL) delete[] Dx_;
	if (Hy_ != NULL) delete[] Hy_;
	if (By_ != NULL) delete[] By_;
	// Incident and transmitted fields.
	if (Exi != NULL) delete[] Exi;
	if (Ext != NULL) delete[] Ext;
	if (Extt != NULL) delete[] Extt;
	// Refractive index.
	if (Exz1 != NULL) delete[] Exz1;
	if (Exz2 != NULL) delete[] Exz2;
	// Drude parameter arrays.
	if (einf != NULL) delete[] einf;
	if (uinf != NULL) delete[] uinf;
	if (wpesq != NULL) delete[] wpesq;
	if (wpmsq != NULL) delete[] wpmsq;
	if (ge != NULL) delete[] ge;
	if (gm != NULL) delete[] gm;
	// Auxiliary field scalars.
	if (ae0 != NULL) delete[] ae0;
	if (ae != NULL) delete[] ae;
	if (be != NULL) delete[] be;
	if (ce != NULL) delete[] ce;
	if (de != NULL) delete[] de;
	if (ee != NULL) delete[] ee;
	if (am0 != NULL) delete[] am0;
	if (am != NULL) delete[] am;
	if (bm != NULL) delete[] bm;
	if (cm != NULL) delete[] cm;
	if (dm != NULL) delete[] dm;
	if (em != NULL) delete[] em;
}
