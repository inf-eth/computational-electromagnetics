// Field arrays indexing macros.
#define Ez(i,j,n) Ez_[(i)+IEz*(j)+IEz*JEz*(n)]
#define Dz(i,j,n) Dz_[(i)+IEz*(j)+IEz*JEz*(n)]
#define Hx(i,j,n) Hx_[(i)+IHx*(j)+IHx*JHx*(n)]
#define Bx(i,j,n) Bx_[(i)+IHx*(j)+IHx*JHx*(n)]
#define Hy(i,j,n) Hy_[(i)+IHy*(j)+IHy*JHy*(n)]
#define By(i,j,n) By_[(i)+IHy*(j)+IHy*JHy*(n)]
// Drude parameter arrays indexing macros.
#define einf(i,j) einf_[(i)+IHx*(j)]
#define uinf(i,j) uinf_[(i)+IHx*(j)]
#define wpesq(i,j) wpesq_[(i)+IHx*(j)]
#define wpmsq(i,j) wpmsq_[(i)+IHx*(j)]
#define ge(i,j) ge_[(i)+IHx*(j)]
#define gm(i,j) gm_[(i)+IHx*(j)]
// Auxiliary scalar arrays indexing macros.
#define ae0(i,j) ae0_[(i)+IHx*(j)]
#define ae(i,j) ae_[(i)+IHx*(j)]
#define be(i,j) be_[(i)+IHx*(j)]
#define ce(i,j) ce_[(i)+IHx*(j)]
#define de(i,j) de_[(i)+IHx*(j)]
#define ee(i,j) ee_[(i)+IHx*(j)]
#define am0(i,j) am0_[(i)+IHx*(j)]
#define am(i,j) am_[(i)+IHx*(j)]
#define bm(i,j) bm_[(i)+IHx*(j)]
#define cm(i,j) cm_[(i)+IHx*(j)]
#define dm(i,j) dm_[(i)+IHx*(j)]
#define em(i,j) em_[(i)+IHx*(j)]
// PML arrays indexing macros.
#define PsiEzX(i,j) PsiEzX_[(i)+IEz*(j)]
#define PsiEzY(i,j) PsiEzY_[(i)+IEz*(j)]
#define PsiHyX(i,j) PsiHyX_[(i)+IHy*(j)]
#define PsiHxY(i,j) PsiHxY_[(i)+IHx*(j)]

#include <FDTD2DDNG.h>
#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream>
using namespace std;

CFDTD2DDNG::CFDTD2DDNG(
						unsigned int pI,
						unsigned int pJ,
						unsigned int pPMLw,
						unsigned int MaxTime,
						unsigned int pSnapshotResolution,
						unsigned int pSnapshotInterval,
						unsigned int pSourceChoice,
						unsigned int pSourcePlane,
						unsigned int pSourceLocationX,
						unsigned int pSourceLocationY):

							// Simulation parameters.
							I(pI),
							J(pJ),
							PMLw(pPMLw),
							SlabLeft(J/3+PMLw),
							SlabRight(2*J/3+PMLw),
							MaxTime(4*J),
							PulseWidth(J/8),
							td(PulseWidth),
							SnapshotResolution(pSnapshotResolution),
							SnapshotInterval(pSnapshotInterval),
							SourceChoice(pSourceChoice),
							SourcePlane(pSourcePlane),
							SourceLocationX(pSourceLocationX),
							SourceLocationY(pSourceLocationY),
							// Constants.
							c(299792458.),
							pi(3.14159265358979323846),
							e0((1e-9)/(36.*pi)),
							u0((1e-7)*4.*pi),
							dt(0.5e-11),
							delta(3e-3),
							Sc(c*dt/delta),
							// Frequency, wavelength, wave number.
							l(PulseWidth*delta),
							f(c/l),
							fmax(1/(2*dt)),
							w(2*pi*f),
							k0(w/c),
							fp(f),
							dr(PulseWidth*dt*2),
							// Data array sizes.
							IEz(I), JEz(J+2*PMLw),
							IHx(I), JHx(J+2*PMLw+1),
							IHy(I), JHy(J+2*PMLw),
							// Data arrays.
							Ez_(NULL), Dz_(NULL), Hx_(NULL), Bx_(NULL), Hy_(NULL), By_(NULL),
							// Incident and transmitted fields.
							Ezi(NULL), Ezt(NULL), Eztt(NULL),
							x1(SlabLeft+1),
							// Refractive index.
							Y1(SlabLeft+5),
							y1((PRECISION)Y1*delta),
							Y2(SlabLeft+6),
							y2((PRECISION)Y1*delta),
							Ezy1(NULL),
							Ezy2(NULL),
							// Drude material parameters.
							einf_(NULL), uinf_(NULL), wpesq_(NULL), wpmsq_(NULL), ge_(NULL), gm_(NULL),
							// Auxiliary field scalars.
							ae0_(NULL), ae_(NULL), be_(NULL), ce_(NULL), de_(NULL), ee_(NULL),
							am0_(NULL), am_(NULL), bm_(NULL), cm_(NULL), dm_(NULL), em_(NULL),
							// PML arrays.
							PsiEzX_(NULL), PsiEzY_(NULL),
							PsiHyX_(NULL), PsiHxY_(NULL),
							// PML parameters.
							kappex(1.), kappey(kappex), kappmx(kappex), kappmy(kappex),
							aex(0.0004), aey(aex), amx(aex), amy(aex),
							sigex(0.), sigey(0.045), sigmx(0.), sigmy(u0/e0*sigey),
							bex(exp(-1*(aex/e0+sigex/(kappex*e0))*dt)), bey(exp(-1*(aey/e0+sigey/(kappey*e0))*dt)), bmx(exp(-1*(amx/u0+sigmx/(kappmx*u0))*dt)), bmy(exp(-1*(amy/u0+sigmy/(kappmy*u0))*dt)),
							Cex((bex-1)*sigex/(sigex*kappex+pow(kappex, 2)*aex)), Cey((bey-1)*sigey/(sigey*kappey+pow(kappey, 2)*aey)), Cmx((bmx-1)*sigmx/(sigmx*kappmx+pow(kappmx, 2)*amx)), Cmy((bmy-1)*sigmy/(sigmy*kappmy+pow(kappmy, 2)*amy)),
							// Snapshot frame number.
							frame(0),
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
	parametersfile.open("FieldData/Parameters.smp", ios::out|ios::binary);
	parametersfile.write((char*)&dt, sizeof(PRECISION));
	parametersfile.write((char*)&k0, sizeof(PRECISION));
	parametersfile.write((char*)&y1, sizeof(PRECISION));
	parametersfile.write((char*)&y2, sizeof(PRECISION));
	parametersfile.write((char*)&I, sizeof(unsigned int));
	parametersfile.write((char*)&J, sizeof(unsigned int));
	parametersfile.write((char*)&PMLw, sizeof(unsigned int));
	parametersfile.write((char*)&MaxTime, sizeof(unsigned int));
	parametersfile.write((char*)&SnapshotResolution, sizeof(unsigned int));
	parametersfile.write((char*)&SnapshotInterval, sizeof(unsigned int));
	parametersfile.write((char*)&SlabLeft, sizeof(unsigned int));
	parametersfile.write((char*)&SlabRight, sizeof(unsigned int));
	parametersfile.close();

	// Printing simulation parameters.
	cout << "I      = " << I << endl;
	cout << "J      = " << J << endl;
	cout << "PMLw       = " << PMLw << endl;
	cout << "MaxTime    = " << MaxTime << endl;
	cout << "frequency  = " << f << " Hz (" << f/1e9 << " GHz)" << endl;
	cout << "fmax       = " << fmax << " Hz (" << fmax/1e9 << " GHz)" << endl;
	cout << "Sc         = " << Sc << endl;
	cout << "Slab left  = " << SlabLeft << endl;
	cout << "Slab right = " << SlabRight << endl;
}
unsigned long CFDTD2DDNG::SimSize()
{
	return (unsigned long)sizeof(*this)+(unsigned long)sizeof(PRECISION)*(8UL*(unsigned long)IEz*(unsigned long)JEz+7UL*(unsigned long)IHy*(unsigned long)JHy+25UL*(unsigned long)IHx*(unsigned long)JHx+5UL*(unsigned long)MaxTime);
}
unsigned long CFDTD2DDNG::HDDSpace()
{
	return (unsigned long)sizeof(PRECISION)*(5UL*(unsigned long)MaxTime+(unsigned long)IEz*(unsigned long)JEz*((unsigned long)MaxTime/(unsigned long)SnapshotInterval+1UL));
}
int CFDTD2DDNG::AllocateMemoryCPU()
{
	// Field arrays.
	Ez_ = new PRECISION[IEz*JEz*3];
	Dz_ = new PRECISION[IEz*JEz*3];
	Hx_ = new PRECISION[IHx*JHx*3];
	Bx_ = new PRECISION[IHx*JHx*3];
	Hy_ = new PRECISION[IHy*JHy*3];
	By_ = new PRECISION[IHy*JHy*3];
	// Incident and transmitted fields.
	Ezi = new PRECISION[MaxTime];
	Ezt = new PRECISION[MaxTime];
	Eztt = new PRECISION[MaxTime];
	// Refractive index.
	Ezy1 = new PRECISION[MaxTime];
	Ezy2 = new PRECISION[MaxTime];
	// Drude parameter arrays.
	einf_ = new PRECISION[IHx*JHx];
	uinf_ = new PRECISION[IHx*JHx];
	wpesq_ = new PRECISION[IHx*JHx];
	wpmsq_ = new PRECISION[IHx*JHx];
	ge_ = new PRECISION[IHx*JHx];
	gm_ = new PRECISION[IHx*JHx];
	// Auxiliary field scalars.
	ae0_ = new PRECISION[IHx*JHx];
	ae_ = new PRECISION[IHx*JHx];
	be_ = new PRECISION[IHx*JHx];
	ce_ = new PRECISION[IHx*JHx];
	de_ = new PRECISION[IHx*JHx];
	ee_ = new PRECISION[IHx*JHx];
	am0_ = new PRECISION[IHx*JHx];
	am_ = new PRECISION[IHx*JHx];
	bm_ = new PRECISION[IHx*JHx];
	cm_ = new PRECISION[IHx*JHx];
	dm_ = new PRECISION[IHx*JHx];
	em_ = new PRECISION[IHx*JHx];
	// PML arrays.
	PsiEzX_ = new PRECISION[IEz*JEz];
	PsiEzY_ = new PRECISION[IEz*JEz];
	PsiHyX_ = new PRECISION[IHy*JHy];
	PsiHxY_ = new PRECISION[IHx*JHx];

	return 0;
}
int CFDTD2DDNG::InitialiseCPU()
{
	for (unsigned int i=0; i<IHx; i++)
	{
		for (unsigned int j=0; j<JHx; j++)
		{
			for (unsigned int n=0; n<3; n++)
			{
				if (i<IEz && j<JEz)
				{
					Ez(i,j,n) = 0.;
					Dz(i,j,n) = 0.;
				}
				Hx(i,j,n) = 0.;
				Bx(i,j,n) = 0.;
				if (i<IHy && j<JHy)
				{
					Hy(i,j,n) = 0.;
					By(i,j,n) = 0.;
				}
			}
			// Parameteric and auxiliary arrays.
			if (j<SlabLeft || j>SlabRight)
			{
				// Outside slab.
				einf(i,j) = 1.;
				uinf(i,j) = 1.;
				wpesq(i,j) = 0.;
				wpmsq(i,j) = 0.;
				ge(i,j) = 0.;
				gm(i,j) = 0.;
			}
			else // Inside slab.
			{
				einf(i,j) = 1.;
				uinf(i,j) = 1.;
				wpesq(i,j) = 2.*pow(w, 2);;
				wpmsq(i,j) = 2.*pow(w, 2);;
				ge(i,j) = w/32.;
				gm(i,j) = w/32.;
			}
			// Auxiliary scalars.
			ae0(i,j) = (4.*pow(dt,2))/(e0*(4.*einf(i,j)+pow(dt,2)*wpesq(i,j)+2.*dt*einf(i,j)*ge(i,j)));
			ae(i,j) = (1./pow(dt,2))*ae0(i,j);
			be(i,j) = (1./(2.*dt))*ge(i,j)*ae0(i,j);
			ce(i,j) = (e0/pow(dt,2))*einf(i,j)*ae0(i,j);
			de(i,j) = (-1.*e0/4.)*wpesq(i,j)*ae0(i,j);
			ee(i,j) = (1./(2.*dt))*e0*einf(i,j)*ge(i,j)*ae0(i,j);
			am0(i,j) = (4.*pow(dt,2))/(u0*(4.*uinf(i,j)+pow(dt,2)*wpmsq(i,j)+2.*dt*uinf(i,j)*gm(i,j)));;
			am(i,j) = (1./pow(dt,2))*am0(i,j);
			bm(i,j) = (1./(2.*dt))*gm(i,j)*am0(i,j);
			cm(i,j) = (u0/pow(dt,2))*uinf(i,j)*am0(i,j);
			dm(i,j) = (-1.*u0/4.)*wpmsq(i,j)*am0(i,j);
			em(i,j) = (1./(2.*dt))*u0*uinf(i,j)*gm(i,j)*am0(i,j);
			// PML Psi arrays.
			if (i<IEz && j<JEz)
			{
				PsiEzX(i,j) = 0.;
				PsiEzY(i,j) = 0.;
			}
			PsiHxY(i,j) = 0.;
			if (i<IHy && j<JHy)
				PsiHyX(i,j) = 0.;
		}
	}
	for (unsigned int i=0; i<MaxTime; i++)
	{
		Ezi[i] = 0.;
		Ezt[i] = 0.;
		Eztt[i] = 0.;
		Ezy1[i] = 0.;
		Ezy2[i] = 0.;
	}

	return 0;
}
int CFDTD2DDNG::InitialiseForSimulationCPU()
{
	for (unsigned int i=0; i<IHx; i++)
	{
		for (unsigned int j=0; j<JHx; j++)
		{
			for (unsigned int n=0; n<3; n++)
			{
				if (i<IEz && j<JEz)
				{
					Ez(i,j,n) = 0.;
					Dz(i,j,n) = 0.;
				}
				Hx(i,j,n) = 0.;
				Bx(i,j,n) = 0.;
				if (i<IHy && j<JHy)
				{
					Hy(i,j,n) = 0.;
					By(i,j,n) = 0.;
				}
			}
			// PML Psi arrays.
			if (i<IEz && j<JEz)
			{
				PsiEzX(i,j) = 0.;
				PsiEzY(i,j) = 0.;
			}
			PsiHxY(i,j) = 0.;
			if (i<IHy && j<JHy)
				PsiHyX(i,j) = 0.;
		}
	}

	return 0;
}
int CFDTD2DDNG::DryRunCPU()
{
	cout << "Dry run (CPU) started..." << endl;/*
	for (unsigned int n=0; n<MaxTime; n++)
	{
		cout << "\r\t\t\r" << n*100/(MaxTime-1) << "%";
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
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
		}
		else if (SourceChoice == 2)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + sin(2.*PI*f*(PRECISION)n*dt) * Sc;
		}
		else if (SourceChoice == 3)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + (1.-2.*pow(PI*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(PI*fp*((PRECISION)n*dt-dr),2)) * Sc;
		}
		// Recording incident field.
		Exi[n] = Ex(x1,nf);

		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}*/
	cout << endl << "Dry run (CPU) completed!" << endl;
	return 0;
}
int CFDTD2DDNG::RunSimulationCPU(bool SaveFields)
{
	stringstream framestream;
	string basename = "FieldData/Ex";
	string filename;
	fstream snapshot;

	cout << "Simulation (CPU) started..." << endl;/*
	for (unsigned int n=0; n<MaxTime; n++)
	{
		cout << "\r\t\t\r" << n*100/(MaxTime-1) << "%";
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
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
		}
		else if (SourceChoice == 2)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + sin(2.*PI*f*(PRECISION)n*dt) * Sc;
		}
		else if (SourceChoice == 3)
		{
			Ex(SourceLocation,nf) = Ex(SourceLocation,nf) + (1.-2.*pow(PI*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(PI*fp*((PRECISION)n*dt-dr),2)) * Sc;
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
			snapshot.write((char*)&(Ex(0,nf)), sizeof(PRECISION)*Size);
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
		snapshot.write((char*)Exi, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ext.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Ext, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Extt.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Extt, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Exz1.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Exz1, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Exz2.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Exz2, sizeof(PRECISION)*MaxTime);
		snapshot.close();
	}*/
	cout << endl << "Simulation (CPU) completed!" << endl;
	return 0;
}
int CFDTD2DDNG::CompleteRunCPU(bool SaveFields)
{
	cout << "Memory required for simulation = " << SimSize() << " bytes (" << (double)SimSize()/1024UL << "kB/" << (double)SimSize()/1024UL/1024UL << "MB)." << endl;
	cout << "HDD space required for data storage = " << HDDSpace() << " bytes (" << (double)HDDSpace()/1024UL << "kB/" << (double)HDDSpace()/1024UL/1024UL << "MB)." << endl;
	SafeCall(AllocateMemoryCPU(), "Error: Allocating memory on CPU.");
	SafeCall(InitialiseCPU(), "Error: Initialising data on CPU.");
	SafeCall(DryRunCPU(), "Error: Dry run (CPU).");
	SafeCall(InitialiseForSimulationCPU(), "Error: Initalising for simulation (CPU).");
	SafeCall(RunSimulationCPU(SaveFields), "Error: Running simulation (CPU).");
	SafeCall(CleanupCPU(), "Error: Cleaning up CPU.");

	return 0;
}
// Timing.
void CFDTD2DDNG::StartTimer()
{
	tStart = GetTimeus64();
}
void CFDTD2DDNG::StopTimer()
{
	tEnd = GetTimeus64();
}
PRECISION CFDTD2DDNG::GetElapsedTime()
{
	return ((PRECISION)(tEnd-tStart))/(1000000.);
}
int CFDTD2DDNG::SafeCall(int Status, const char *Error)
{
	if (Status != 0)
	{
		if (Error!=NULL) cout << Error << endl;
		exit(Status);
	}
	return Status;
}
int CFDTD2DDNG::CleanupCPU()
{
	// Field arrays.
	DeleteArray(Ez_);
	DeleteArray(Dz_);
	DeleteArray(Hx_);
	DeleteArray(Bx_);
	DeleteArray(Hy_);
	DeleteArray(By_);
	// Incident and transmitted fields.
	DeleteArray(Ezi);
	DeleteArray(Ezt);
	DeleteArray(Eztt);
	// Refractive index.
	DeleteArray(Ezy1);
	DeleteArray(Ezy2);
	// Drude parameter arrays.
	DeleteArray(einf_);
	DeleteArray(uinf_);
	DeleteArray(wpesq_);
	DeleteArray(wpmsq_);
	DeleteArray(ge_);
	DeleteArray(gm_);
	// Auxiliary field scalars.
	DeleteArray(ae0_);
	DeleteArray(ae_);
	DeleteArray(be_);
	DeleteArray(ce_);
	DeleteArray(de_);
	DeleteArray(ee_);
	DeleteArray(am0_);
	DeleteArray(am_);
	DeleteArray(bm_);
	DeleteArray(cm_);
	DeleteArray(dm_);
	DeleteArray(em_);
	// PML arrays.
	DeleteArray(PsiEzX_);
	DeleteArray(PsiEzY_);
	DeleteArray(PsiHyX_);
	DeleteArray(PsiHxY_);

	return 0;
}
CFDTD2DDNG::~CFDTD2DDNG()
{
	// Field arrays.
	DeleteArray(Ez_);
	DeleteArray(Dz_);
	DeleteArray(Hx_);
	DeleteArray(Bx_);
	DeleteArray(Hy_);
	DeleteArray(By_);
	// Incident and transmitted fields.
	DeleteArray(Ezi);
	DeleteArray(Ezt);
	DeleteArray(Eztt);
	// Refractive index.
	DeleteArray(Ezy1);
	DeleteArray(Ezy2);
	// Drude parameter arrays.
	DeleteArray(einf_);
	DeleteArray(uinf_);
	DeleteArray(wpesq_);
	DeleteArray(wpmsq_);
	DeleteArray(ge_);
	DeleteArray(gm_);
	// Auxiliary field scalars.
	DeleteArray(ae0_);
	DeleteArray(ae_);
	DeleteArray(be_);
	DeleteArray(ce_);
	DeleteArray(de_);
	DeleteArray(ee_);
	DeleteArray(am0_);
	DeleteArray(am_);
	DeleteArray(bm_);
	DeleteArray(cm_);
	DeleteArray(dm_);
	DeleteArray(em_);
	// PML arrays.
	DeleteArray(PsiEzX_);
	DeleteArray(PsiEzY_);
	DeleteArray(PsiHyX_);
	DeleteArray(PsiHxY_);
}
