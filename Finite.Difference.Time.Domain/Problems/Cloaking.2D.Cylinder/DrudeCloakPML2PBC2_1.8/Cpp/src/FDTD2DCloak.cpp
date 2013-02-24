// Field arrays indexing macros.
#define Ez(i,j,n) Ez_[(i)+IEz*(j)+IEz*JEz*(n)]
#define Dz(i,j,n) Dz_[(i)+IEz*(j)+IEz*JEz*(n)]
#define EzMask(i,j) EzMask_[(i)+IEz*(j)]
#define Hx(i,j,n) Hx_[(i)+IHx*(j)+IHx*JHx*(n)]
#define Bx(i,j,n) Bx_[(i)+IHx*(j)+IHx*JHx*(n)]
#define BxAve(i,j,n) BxAve_[(i)+IHy*(j)+IHy*JHy*(n)]
#define Hy(i,j,n) Hy_[(i)+IHy*(j)+IHy*JHy*(n)]
#define By(i,j,n) By_[(i)+IHy*(j)+IHy*JHy*(n)]
#define ByAve(i,j,n) ByAve_[(i)+IHx*(j)+IHx*JHx*(n)]
// Auxiliary scalar arrays indexing macros.
#define ax(i,j,k) ax_[(i)+IHx*(j)+IHx*JHx*(k)]
#define ay(i,j,k) ay_[(i)+IHy*(j)+IHy*JHy*(k)]
#define az(i,j,k) az_[(i)+IEz*(j)+IEz*JEz*(k)]
// PML arrays indexing macros.
#define PsiEzX(i,j) PsiEzX_[(i)+IEz*(j)]
#define PsiEzY(i,j) PsiEzY_[(i)+IEz*(j)]
#define PsiHyX(i,j) PsiHyX_[(i)+IHy*(j)]
#define PsiHxY(i,j) PsiHxY_[(i)+IHx*(j)]

#include <FDTD2DCloak.h>
#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;

CFDTD2DCloak::CFDTD2DCloak(
						unsigned int pI,
						unsigned int pJ,
						unsigned int pPMLw,
						unsigned int pMaxTime,
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
							ra(0.1),
							rb(0.2),
							MaxTime(pMaxTime),
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
							dt(0.2e-11),
							delta(2.5e-3),
							Sc(c*dt/delta),
							x0(delta*(PRECISION)I/2.),
							y0(delta*((PRECISION)PMLw+(PRECISION)I/2.)),
							// Frequency, wavelength, wave number.
							l(PulseWidth*delta),
							f(2e9),
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
							Ez_(NULL), Dz_(NULL), EzMask_(NULL), Hx_(NULL), Bx_(NULL), BxAve_(NULL), Hy_(NULL), By_(NULL), ByAve_(NULL),
							// Incident and transmitted fields.
							Ezi(NULL), Ezt(NULL), Eztt(NULL),
							x1(J-PMLw-10),
							// Refractive index.
							Y1(J-PMLw-10),
							y1((PRECISION)Y1*delta),
							Y2(J-PMLw-9),
							y2((PRECISION)Y2*delta),
							Ezy1(NULL),
							Ezy2(NULL),
							// Auxiliary field scalars.
							ax_(NULL), ay_(NULL), az_(NULL),
							// PML arrays.
							PsiEzX_(NULL), PsiEzY_(NULL),
							PsiHyX_(NULL), PsiHxY_(NULL),
							// PML parameters.
							kappe(1.), kappm(1.),
							kappex(1.), kappey(kappex), kappmx(kappex), kappmy(kappex),
							aex(0.0004), aey(aex), amx(aex), amy(aex),
							sigex(0.), sigey(0.045), sigmx(0.), sigmy(u0/e0*sigey),
							bex(exp(-1.*(aex/e0+sigex/(kappex*e0))*dt)), bey(exp(-1.*(aey/e0+sigey/(kappey*e0))*dt)), bmx(exp(-1.*(amx/u0+sigmx/(kappmx*u0))*dt)), bmy(exp(-1.*(amy/u0+sigmy/(kappmy*u0))*dt)),
							Cex((bex-1.)*sigex/(sigex*kappex+pow(kappe, 2)*aex)), Cey((bey-1.)*sigey/(sigey*kappey+pow(kappe, 2)*aey)), Cmx((bmx-1.)*sigmx/(sigmx*kappmx+pow(kappm, 2)*amx)), Cmy((bmy-1.)*sigmy/(sigmy*kappmy+pow(kappm, 2)*amy)),
							// Snapshot frame number.
							frame(0),
							// Time indices.
							nf(2), n0(1), np(0),
							// Timer variables.
							tStart(0LL), tEnd(0LL),
							tDelta(0LL), tPaused(true)
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
	parametersfile.write((char*)&x1, sizeof(unsigned int));
	parametersfile.write((char*)&x1, sizeof(unsigned int));
	parametersfile.close();

	// Printing simulation parameters.
	cout << "I      = " << I << endl;
	cout << "J      = " << J << endl;
	cout << "PMLw       = " << PMLw << endl;
	cout << "MaxTime    = " << MaxTime << endl;
	cout << "frequency  = " << f << " Hz (" << f/1e9 << " GHz)" << endl;
	cout << "fmax       = " << fmax << " Hz (" << fmax/1e9 << " GHz)" << endl;
	cout << "Sc         = " << Sc << endl;
	cout << "Slab left  = " << x1 << endl;
	cout << "Slab right = " << x1 << endl;
}
unsigned long CFDTD2DCloak::SimSize()
{
	return (unsigned long)sizeof(*this)+(unsigned long)sizeof(PRECISION)*(14UL*(unsigned long)IEz*(unsigned long)JEz+19UL*(unsigned long)IHy*(unsigned long)JHy+19UL*(unsigned long)IHx*(unsigned long)JHx+5UL*(unsigned long)MaxTime);
}
unsigned long CFDTD2DCloak::HDDSpace()
{
	return (unsigned long)sizeof(PRECISION)*(5UL*(unsigned long)MaxTime+(unsigned long)IEz*(unsigned long)JEz*((unsigned long)MaxTime/(unsigned long)SnapshotInterval+1UL));
}
int CFDTD2DCloak::AllocateMemoryCPU()
{
	// Field arrays.
	Ez_ = new PRECISION[IEz*JEz*3];
	Dz_ = new PRECISION[IEz*JEz*3];
	EzMask_ = new PRECISION[IEz*JEz];
	Hx_ = new PRECISION[IHx*JHx*3];
	Bx_ = new PRECISION[IHx*JHx*3];
	BxAve_ = new PRECISION[IHy*JHy*3];
	Hy_ = new PRECISION[IHy*JHy*3];
	By_ = new PRECISION[IHy*JHy*3];
	ByAve_ = new PRECISION[IHx*JHx*3];
	// Incident and transmitted fields.
	Ezi = new PRECISION[MaxTime];
	Ezt = new PRECISION[MaxTime];
	Eztt = new PRECISION[MaxTime];
	// Refractive index.
	Ezy1 = new PRECISION[MaxTime];
	Ezy2 = new PRECISION[MaxTime];
	// Auxiliary field scalars.
	ax_ = new PRECISION[IHx*JHx*9];
	ay_ = new PRECISION[IHy*JHy*9];
	az_ = new PRECISION[IEz*JEz*5];
	// PML arrays.
	PsiEzX_ = new PRECISION[IEz*JEz];
	PsiEzY_ = new PRECISION[IEz*JEz];
	PsiHyX_ = new PRECISION[IHy*JHy];
	PsiHxY_ = new PRECISION[IHx*JHx];

	return 0;
}
int CFDTD2DCloak::InitialiseCPU()
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
					if (n<1)
						EzMask(i,j) = mask((PRECISION)i,(PRECISION)j-0.5);
				}
				Hx(i,j,n) = 0.;
				Bx(i,j,n) = 0.;
				ByAve(i,j,n) = 0.;
				if (i<IHy && j<JHy)
				{
					Hy(i,j,n) = 0.;
					By(i,j,n) = 0.;
					BxAve(i,j,n) = 0.;
				}
			}
			// Auxiliary scalars.
			ax(i,j,0) = pow(sinphi((PRECISION)i,(PRECISION)j),2)*(1./pow(dt,2)+wpmsq((PRECISION)i,(PRECISION)j)/4.)+uphi((PRECISION)i,(PRECISION)j)*pow(cosphi((PRECISION)i,(PRECISION)j),2)*(1./pow(dt,2));
			ax(i,j,1) = pow(sinphi((PRECISION)i,(PRECISION)j),2)*(-2./pow(dt,2)+wpmsq((PRECISION)i,(PRECISION)j)/2.)-uphi((PRECISION)i,(PRECISION)j)*pow(cosphi((PRECISION)i,(PRECISION)j),2)*(2./pow(dt,2));
			ax(i,j,2) = ax(i,j,0);
			ax(i,j,3) = (uphi((PRECISION)i,(PRECISION)j)*(1./pow(dt,2))-(1./pow(dt,2)+wpmsq((PRECISION)i,(PRECISION)j)/4.))*sinphi((PRECISION)i,(PRECISION)j)*cosphi((PRECISION)i,(PRECISION)j);
			ax(i,j,4) = (uphi((PRECISION)i,(PRECISION)j)*(-2./pow(dt,2))-(-2./pow(dt,2)+wpmsq((PRECISION)i,(PRECISION)j)/2.))*sinphi((PRECISION)i,(PRECISION)j)*cosphi((PRECISION)i,(PRECISION)j);
			ax(i,j,5) = ax(i,j,3);
			ax(i,j,6) = u0*uphi((PRECISION)i,(PRECISION)j)*(-2./pow(dt,2)+wpmsq((PRECISION)i,(PRECISION)j)/2.);
			ax(i,j,7) = u0*uphi((PRECISION)i,(PRECISION)j)*(1./pow(dt,2)+wpmsq((PRECISION)i,(PRECISION)j)/4.);
			ax(i,j,8) = ax(i,j,7);
			if (j < JHx-1)
			{
				ay(i,j,0) = pow(cosphi((PRECISION)i-0.5,(PRECISION)j-0.5),2)*(1./pow(dt,2)+wpmsq((PRECISION)i-0.5,(PRECISION)j-0.5)/4.)+uphi((PRECISION)i-0.5,(PRECISION)j-0.5)*pow(sinphi((PRECISION)i-0.5,(PRECISION)j-0.5),2)*(1./pow(dt,2));
				ay(i,j,1) = pow(cosphi((PRECISION)i-0.5,(PRECISION)j-0.5),2)*(-2./pow(dt,2)+wpmsq((PRECISION)i-0.5,(PRECISION)j-0.5)/2.)-uphi((PRECISION)i-0.5,(PRECISION)j-0.5)*pow(sinphi((PRECISION)i-0.5,(PRECISION)j-0.5),2)*(2./pow(dt,2));
				ay(i,j,2) = ay(i,j,0);
				ay(i,j,3) = (uphi((PRECISION)i-0.5,(PRECISION)j-0.5)*(1./pow(dt,2))-(1./pow(dt,2)+wpmsq((PRECISION)i-0.5,(PRECISION)j-0.5)/4.))*sinphi((PRECISION)i-0.5,(PRECISION)j-0.5)*cosphi((PRECISION)i-0.5,(PRECISION)j-0.5);
				ay(i,j,4) = (uphi((PRECISION)i-0.5,(PRECISION)j-0.5)*(-2./pow(dt,2))-(-2./pow(dt,2)+wpmsq((PRECISION)i-0.5,(PRECISION)j-0.5)/2.))*sinphi((PRECISION)i-0.5,(PRECISION)j-0.5)*cosphi((PRECISION)i-0.5,(PRECISION)j-0.5);
				ay(i,j,5) = ay(i,j,3);
				ay(i,j,6) = u0*uphi((PRECISION)i-0.5,(PRECISION)j-0.5)*(-2./pow(dt,2)+wpmsq((PRECISION)i-0.5,(PRECISION)j-0.5)/2.);
				ay(i,j,7) = u0*uphi((PRECISION)i-0.5,(PRECISION)j-0.5)*(1./pow(dt,2)+wpmsq((PRECISION)i-0.5,(PRECISION)j-0.5)/4.);
				ay(i,j,8) = ay(i,j,7);

				PRECISION a0 = (4.*pow(dt,2))/(e0*(4.*einf((PRECISION)i,(PRECISION)j-0.5)+pow(dt,2)*wpesq((PRECISION)i,(PRECISION)j-0.5)+2.*dt*einf((PRECISION)i,(PRECISION)j-0.5)*ge((PRECISION)i,(PRECISION)j-0.5)));
				az(i,j,0) = (1./pow(dt,2))*a0/A((PRECISION)i,(PRECISION)j-0.5);
				az(i,j,1) = (1./(2.*dt))*ge((PRECISION)i,(PRECISION)j-0.5)*a0/A((PRECISION)i,(PRECISION)j-0.5);
				az(i,j,2) = (e0/pow(dt,2))*einf((PRECISION)i,(PRECISION)j-0.5)*a0;
				az(i,j,3) = (-1.*e0/4.)*wpesq((PRECISION)i,(PRECISION)j-0.5)*a0;
				az(i,j,4) = (1./(2.*dt))*e0*einf((PRECISION)i,(PRECISION)j-0.5)*ge((PRECISION)i,(PRECISION)j-0.5)*a0;
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
int CFDTD2DCloak::InitialiseForSimulationCPU()
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
int CFDTD2DCloak::DryRunCPU()
{
	cout << "Dry run (CPU) started..." << endl;
	for (unsigned int n=0; n<MaxTime; n++)
	{
		if (n%SnapshotInterval == 0)
			cout << "\r" << setprecision(4) << (float)n*100/(MaxTime-1) << "%  " << flush;

		// ========================== Bx and Hx ==========================
		for (unsigned int i=0; i<IHx; i++)
		{
			// Calculation of PsiHxY.
			for (unsigned int j=1; j<JHx-1; j++)
				PsiHxY(i,j) = (Cmy/delta)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) + bmy*PsiHxY(i,j);

			// Bx in normal space.
			for (unsigned int j=1+PMLw; j<JHx-PMLw-1; j++)
			{
				Bx(i,j,nf) = Bx(i,j,n0) + (-Ez(i,j,n0) + Ez(i,j-1,n0)) * dt/delta;
				Hx(i,j,nf) = Bx(i,j,nf)/u0;
			}
			// Bx in lower PML.
			for (unsigned int j=1; j<PMLw+1; j++)
			{
				Bx(i,j,nf) = Bx(i,j,n0) + dt*((1./kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1./delta + PsiHxY(i,j));
				Hx(i,j,nf) = Bx(i,j,nf)/u0;
			}
			// Bx in upper PML.
			for (unsigned int j=JHx-PMLw-1; j<JHx-1; j++)
			{
				Bx(i,j,nf) = Bx(i,j,n0) + dt*((1./kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1./delta + PsiHxY(i,j));
				Hx(i,j,nf) = Bx(i,j,nf)/u0;
			}
		}
		
		// ========================== By and Hy ==========================
		for (unsigned int i=0; i<IHy-1; i++)
		{
			// PsiHyX arrays.
			for (unsigned int j=0; j<JHy; j++)
			{
				PsiHyX(i,j) = (Cmx/delta)*(Ez(i+1,j,n0)-Ez(i,j,n0)) + bmx*PsiHyX(i,j);
				if (i==0)
					PsiHyX(IHy-1,j) = (Cmx/delta)*(Ez(0,j,n0)-Ez(IHy-1,j,n0)) + bmx*PsiHyX(IHy-1,j); // PBC
			}
			// By in normal space.
			for (unsigned int j=PMLw; j<JHy-PMLw; j++)
			{
				By(i,j,nf) = By(i,j,n0) + (Ez(i+1,j,n0) - Ez(i,j,n0)) * dt/delta;
				Hy(i,j,nf) = By(i,j,nf)/u0;
				if (i==0)
				{
					By(IHy-1,j,nf) = By(IHy-1,j,n0) + (Ez(0,j,n0) - Ez(IHy-1,j,n0)) * dt/delta; // PBC
					Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/u0; // PBC
				}
			}
			// By in Lower PML.
			for (unsigned int j=0; j<PMLw; j++)
			{
				By(i,j,nf) = By(i,j,n0) + dt*((1./kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1./delta + PsiHyX(i,j));
				Hy(i,j,nf) = By(i,j,nf)/u0;
				if (i==0)
				{
					By(IHy-1,j,nf) = By(IHy-1,j,n0) + dt*((1./kappmx)*(Ez(0,j,n0) - Ez(IHy-1,j,n0)) * 1./delta + PsiHyX(IHy-1,j)); // PBC
					Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/u0; // PBC
				}
			}
			// By in upper PML.
			for (unsigned int j=JHy-PMLw; j<JHy; j++)
			{
				By(i,j,nf) = By(i,j,n0) + dt*((1./kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1./delta + PsiHyX(i,j));
				Hy(i,j,nf) = By(i,j,nf)/u0;
				if (i==0)
				{
					By(IHy-1,j,nf) = By(IHy-1,j,n0) + dt*((1./kappmx)*(Ez(0,j,n0) - Ez(IHy-1,j,n0)) * 1./delta + PsiHyX(IHy-1,j)); // PBC
					Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/u0; // PBC
				}
			}
		}

		// ========================== Dz and Ez ==========================
		for (unsigned int i=1; i<IEz; i++)
		{
			// Psi arrays.
			for (unsigned int j=0; j<JEz; j++)
			{
				PsiEzX(i,j) = (Cex/delta)*(Hy(i,j,nf)-Hy(i-1,j,nf)) + bex*PsiEzX(i,j);
				if (i==1)
					PsiEzX(0,j) = (Cex/delta)*(Hy(0,j,nf)-Hy(IEz-1,j,nf)) + bex*PsiEzX(0,j); // PBC
				PsiEzY(i,j) = (Cey/delta)*(-Hx(i,j+1,nf)+Hx(i,j,nf)) + bey*PsiEzY(i,j);
				if (i==1)
					PsiEzY(0,j) = (Cey/delta)*(-Hx(0,j+1,nf)+Hx(0,j,nf)) + bey*PsiEzY(0,j); // PBC
			}
			// Dz in normal space.
			for (unsigned int j=PMLw; j<JEz-PMLw; j++)
			{
				Dz(i,j,nf) = Dz(i,j,n0) + (Hy(i,j,nf)-Hy(i-1,j,nf)-Hx(i,j+1,nf)+Hx(i,j,nf)) * dt/delta;
				Ez(i,j,nf) = Dz(i,j,nf)/e0;
				if (i==1)
				{
					Dz(0,j,nf) = Dz(0,j,n0) + (Hy(0,j,nf)-Hy(IEz-1,j,nf)-Hx(0,j+1,nf)+Hx(0,j,nf)) * dt/delta; // PBC
					Ez(0,j,nf) = Dz(0,j,nf)/e0; // PBC
				}
			}
			// Dz in lower PML.
			for (unsigned int j=0; j<PMLw; j++)
			{
				Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1./kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1./kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1./delta + PsiEzX(i,j) + PsiEzY(i,j));
				Ez(i,j,nf) = Dz(i,j,nf)/e0;
				if (i==1)
				{
					Dz(0,j,nf) = Dz(0,j,n0) + dt*(((1./kappex)*(Hy(0,j,nf)-Hy(IEz-1,j,nf))+(1./kappey)*(-Hx(0,j+1,nf)+Hx(0,j,nf))) * 1./delta + PsiEzX(0,j) + PsiEzY(0,j)); // PBC
					Ez(0,j,nf) = Dz(0,j,nf)/e0; // PBC
				}
			}
			// Dz in upper PML.
			for (unsigned int j=JEz-PMLw; j<JEz; j++)
			{
				Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1./kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1./kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1./delta + PsiEzX(i,j) + PsiEzY(i,j));
				Ez(i,j,nf) = Dz(i,j,nf)/e0;
				if (i==1)
				{
					Dz(0,j,nf) = Dz(0,j,n0) + dt*(((1./kappex)*(Hy(0,j,nf)-Hy(IEz-1,j,nf))+(1./kappey)*(-Hx(0,j+1,nf)+Hx(0,j,nf))) * 1./delta + PsiEzX(0,j) + PsiEzY(0,j)); // PBC
					Ez(0,j,nf) = Dz(0,j,nf)/e0; // PBC
				}
			}
		}

		// ====================== Source ======================
		if (SourcePlane == 1)
		{
			unsigned int j=SourceLocationY;
			for (unsigned int i=0; i<IEz; i++)
			{
				if (SourceChoice == 1)
					Ez(i,j,nf) += exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
				else if (SourceChoice == 2)
					Ez(i,j,nf) += sin(2.*pi*f*(PRECISION)n*dt) * Sc;
				else if (SourceChoice == 3)
					Ez(i,j,nf) += (1.-2.*pow(pi*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(pi*fp*((PRECISION)n*dt-dr),2)) * Sc;

				Dz(i,j,nf) = e0*Ez(i,j,nf);
			}
		}
		else
		{
			unsigned int i=SourceLocationX;
			unsigned int j=SourceLocationY;
			if (SourceChoice == 1)
				Ez(i,j,nf) += exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
			else if (SourceChoice == 2)
				Ez(i,j,nf) += sin(2.*pi*f*(PRECISION)n*dt) * Sc;
			else if (SourceChoice == 3)
				Ez(i,j,nf) += (1.-2.*pow(pi*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(pi*fp*((PRECISION)n*dt-dr),2)) * Sc;

			Dz(i,j,nf) = e0*Ez(i,j,nf);
		}
		Ezi[n] = Ez(IEz/2,x1,nf); // Incident field.

		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}
	cout << endl << "Dry run (CPU) completed!" << endl;
	return 0;
}
int CFDTD2DCloak::RunSimulationCPU(bool SaveFields)
{
	stringstream framestream;
	string basename = "FieldData/Ez";
	string filename;
	fstream snapshot;
	frame = 0U;

	cout << "Simulation (CPU) started..." << endl;
	for (unsigned int n=0; n<MaxTime; n++)
	{
		if (n%SnapshotInterval == 0)
			cout << "\r" << setprecision(4) << (float)n*100/(MaxTime-1) << "%  " << flush;

		// ========================== Bx ==========================
		for (unsigned int i=0; i<IHx; i++)
		{
			// Calculation of PsiHxY.
			for (unsigned int j=1; j<JHx-1; j++)
				PsiHxY(i,j) = (Cmy/delta)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) + bmy*PsiHxY(i,j);

			// Bx in normal space.
			for (unsigned int j=1+PMLw; j<JHx-PMLw-1; j++)
				Bx(i,j,nf) = Bx(i,j,n0) + (-Ez(i,j,n0) + Ez(i,j-1,n0)) * dt/delta;

			// Bx in lower PML.
			for (unsigned int j=1; j<PMLw+1; j++)
				Bx(i,j,nf) = Bx(i,j,n0) + dt*((1./kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1./delta + PsiHxY(i,j));

			// Bx in upper PML.
			for (unsigned int j=JHx-PMLw-1; j<JHx-1; j++)
				Bx(i,j,nf) = Bx(i,j,n0) + dt*((1./kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1./delta + PsiHxY(i,j));
		}

		// ========================== By ==========================
		for (unsigned int i=0; i<IHy-1; i++)
		{
			// PsiHyX arrays.
			for (unsigned int j=0; j<JHy; j++)
			{
				PsiHyX(i,j) = (Cmx/delta)*(Ez(i+1,j,n0)-Ez(i,j,n0)) + bmx*PsiHyX(i,j);
				if (i==0)
					PsiHyX(IHy-1,j) = (Cmx/delta)*(Ez(0,j,n0)-Ez(IHy-1,j,n0)) + bmx*PsiHyX(IHy-1,j); // PBC
			}
			// By in normal space.
			for (unsigned int j=PMLw; j<JHy-PMLw; j++)
			{
				By(i,j,nf) = By(i,j,n0) + (Ez(i+1,j,n0) - Ez(i,j,n0)) * dt/delta;
				if (i==0)
					By(IHy-1,j,nf) = By(IHy-1,j,n0) + (Ez(0,j,n0) - Ez(IHy-1,j,n0)) * dt/delta; // PBC
			}
			// By in Lower PML.
			for (unsigned int j=0; j<PMLw; j++)
			{
				By(i,j,nf) = By(i,j,n0) + dt*((1./kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1./delta + PsiHyX(i,j));
				if (i==0)
					By(IHy-1,j,nf) = By(IHy-1,j,n0) + dt*((1./kappmx)*(Ez(0,j,n0) - Ez(IHy-1,j,n0)) * 1./delta + PsiHyX(IHy-1,j)); // PBC
			}
			// By in upper PML.
			for (unsigned int j=JHy-PMLw; j<JHy; j++)
			{
				By(i,j,nf) = By(i,j,n0) + dt*((1./kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1./delta + PsiHyX(i,j));
				if (i==0)
					By(IHy-1,j,nf) = By(IHy-1,j,n0) + dt*((1./kappmx)*(Ez(0,j,n0) - Ez(IHy-1,j,n0)) * 1./delta + PsiHyX(IHy-1,j)); // PBC
			}
		}

		// --- Sync ---

		for (unsigned int i=2; i<IHy-2; i++)
		{
			for (unsigned int j=3+PMLw; j<JHx-PMLw-3; j++)
			{
				BxAve(i,j,nf) = (Bx(i,j,nf)+Bx(i+1,j,nf)+Bx(i,j+1,nf)+Bx(i+1,j+1,nf))/4.;
				ByAve(i,j,nf) = (By(i,j,nf)+By(i-1,j,nf)+By(i,j-1,nf)+By(i-1,j-1,nf))/4.;
			}
		}

		// ========================== Hx ==========================
		for (unsigned int i=0; i<IHx; i++)
		{
			// Hx in normal space.
			for (unsigned int j=1+PMLw; j<JHx-PMLw-1; j++)
				Hx(i,j,nf) = (ax(i,j,0)*Bx(i,j,nf)+ax(i,j,1)*Bx(i,j,n0)+ax(i,j,2)*Bx(i,j,np)+ax(i,j,3)*ByAve(i,j,nf)+ax(i,j,4)*ByAve(i,j,n0)+ax(i,j,5)*ByAve(i,j,np)-ax(i,j,6)*Hx(i,j,n0)-ax(i,j,7)*Hx(i,j,np))/ax(i,j,8);

			// Hx in lower PML.
			for (unsigned int j=1; j<PMLw+1; j++)
				Hx(i,j,nf) = Bx(i,j,nf)/u0;

			// Hx in upper PML.
			for (unsigned int j=JHx-PMLw-1; j<JHx-1; j++)
				Hx(i,j,nf) = Bx(i,j,nf)/u0;
		}

		// ========================== Hy ==========================
		for (unsigned int i=0; i<IHy-1; i++)
		{
			// Hy in normal space.
			for (unsigned int j=PMLw; j<JHy-PMLw; j++)
				Hy(i,j,nf) = (ay(i,j,0)*By(i,j,nf)+ay(i,j,1)*By(i,j,n0)+ay(i,j,2)*By(i,j,np)+ay(i,j,3)*BxAve(i,j,nf)+ay(i,j,4)*BxAve(i,j,n0)+ay(i,j,5)*BxAve(i,j,np)-ay(i,j,6)*Hy(i,j,n0)-ay(i,j,7)*Hy(i,j,np))/ay(i,j,8);

			// Hy in Lower PML.
			for (unsigned int j=0; j<PMLw; j++)
				Hy(i,j,nf) = By(i,j,nf)/u0;

			// Hy in upper PML.
			for (unsigned int j=JHy-PMLw; j<JHy; j++)
				Hy(i,j,nf) = By(i,j,nf)/u0;
		}

		// ========================== Dz and Ez ==========================
		for (unsigned int i=1; i<IEz; i++)
		{
			// Psi arrays.
			for (unsigned int j=0; j<JEz; j++)
			{
				PsiEzX(i,j) = (Cex/delta)*(Hy(i,j,nf)-Hy(i-1,j,nf)) + bex*PsiEzX(i,j);
				PsiEzY(i,j) = (Cey/delta)*(-Hx(i,j+1,nf)+Hx(i,j,nf)) + bey*PsiEzY(i,j);
				if (i==1)
				{
					PsiEzX(0,j) = (Cex/delta)*(Hy(0,j,nf)-Hy(IEz-1,j,nf)) + bex*PsiEzX(0,j); // PBC
					PsiEzY(0,j) = (Cey/delta)*(-Hx(0,j+1,nf)+Hx(0,j,nf)) + bey*PsiEzY(0,j); // PBC
				}
			}
			// Dz in normal space.
			for (unsigned int j=PMLw; j<JEz-PMLw; j++)
			{
				Dz(i,j,nf) = Dz(i,j,n0) + (Hy(i,j,nf)-Hy(i-1,j,nf)-Hx(i,j+1,nf)+Hx(i,j,nf)) * dt/delta;
				Ez(i,j,nf) = EzMask(i,j) * (az(i,j,0)*(Dz(i,j,nf)-2.*Dz(i,j,n0)+Dz(i,j,np))+az(i,j,1)*(Dz(i,j,nf)-Dz(i,j,np))+az(i,j,2)*(2.*Ez(i,j,n0)-Ez(i,j,np))+az(i,j,3)*(2.*Ez(i,j,n0)+Ez(i,j,np))+az(i,j,4)*Ez(i,j,np));
				if (i==1)
				{
					Dz(0,j,nf) = Dz(0,j,n0) + (Hy(0,j,nf)-Hy(IEz-1,j,nf)-Hx(0,j+1,nf)+Hx(0,j,nf)) * dt/delta; // PBC
					Ez(0,j,nf) = EzMask(0,j) * (az(0,j,0)*(Dz(0,j,nf)-2.*Dz(0,j,n0)+Dz(0,j,np))+az(0,j,1)*(Dz(0,j,nf)-Dz(0,j,np))+az(0,j,2)*(2.*Ez(0,j,n0)-Ez(0,j,np))+az(0,j,4)*(2.*Ez(0,j,n0)+Ez(0,j,np))+az(0,j,4)*Ez(0,j,np)); // PBC
				}
			}
			// Dz in lower PML.
			for (unsigned int j=0; j<PMLw; j++)
			{
				Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1./kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1./kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1./delta + PsiEzX(i,j) + PsiEzY(i,j));
				Ez(i,j,nf) = EzMask(i,j) * Dz(i,j,nf)/e0;
				if (i==1)
				{
					Dz(0,j,nf) = Dz(0,j,n0) + dt*(((1./kappex)*(Hy(0,j,nf)-Hy(IEz-1,j,nf))+(1./kappey)*(-Hx(0,j+1,nf)+Hx(0,j,nf))) * 1./delta + PsiEzX(0,j) + PsiEzY(0,j)); // PBC
					Ez(0,j,nf) = EzMask(0,j) * Dz(0,j,nf)/e0; // PBC
				}
			}
			// Dz in upper PML.
			for (unsigned int j=JEz-PMLw; j<JEz; j++)
			{
				Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1./kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1./kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1./delta + PsiEzX(i,j) + PsiEzY(i,j));
				Ez(i,j,nf) = EzMask(0,j) * Dz(i,j,nf)/e0;
				if (i==1)
				{
					Dz(0,j,nf) = Dz(0,j,n0) + dt*(((1./kappex)*(Hy(0,j,nf)-Hy(IEz-1,j,nf))+(1./kappey)*(-Hx(0,j+1,nf)+Hx(0,j,nf))) * 1./delta + PsiEzX(0,j) + PsiEzY(0,j)); // PBC
					Ez(0,j,nf) = EzMask(0,j) * Dz(0,j,nf)/e0; // PBC
				}
			}
		}

		// ====================== Source ======================
		if (SourcePlane == 1)
		{
			unsigned int j=SourceLocationY;
			for (unsigned int i=0; i<IEz; i++)
			{
				if (SourceChoice == 1)
					Ez(i,j,nf) += exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
				else if (SourceChoice == 2)
					Ez(i,j,nf) += sin(2.*pi*f*(PRECISION)n*dt) * Sc;
				else if (SourceChoice == 3)
					Ez(i,j,nf) += (1.-2.*pow(pi*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(pi*fp*((PRECISION)n*dt-dr),2)) * Sc;

				Dz(i,j,nf) = e0*Ez(i,j,nf);
			}
		}
		else
		{
			unsigned int i=SourceLocationX;
			unsigned int j=SourceLocationY;
			if (SourceChoice == 1)
				Ez(i,j,nf) += exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
			else if (SourceChoice == 2)
				Ez(i,j,nf) += sin(2.*pi*f*(PRECISION)n*dt) * Sc;
			else if (SourceChoice == 3)
				Ez(i,j,nf) += (1.-2.*pow(pi*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(pi*fp*((PRECISION)n*dt-dr),2)) * Sc;

			Dz(i,j,nf) = e0*Ez(i,j,nf);
		}
		// Recording transmitted fields.
		Ezt[n] = Ez(IEz/2,x1,nf);
		Eztt[n] = Ez(IEz/2,x1+1,nf);
		// Fields for refractive index.
		Ezy1[n] = Ez(IEz/2,Y1,nf);
		Ezy2[n] = Ez(IEz/2,Y2,nf);

		// Saving electric field snapshot.
		if (n%SnapshotInterval == 0 && SaveFields == true)
		{
			// Write E-field to file.
			framestream.str(std::string());			// Clearing stringstream contents.
			framestream << ++frame;
			filename = basename + framestream.str() + ".fdt";
			snapshot.open(filename.c_str(), std::ios::out|std::ios::binary);
			snapshot.write((char*)&(Ez(0,0,nf)), sizeof(PRECISION)*IEz*JEz);
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
		snapshot.open("FieldData/Ezi.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Ezi, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ezt.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Ezt, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Eztt.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Eztt, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ezy1.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Ezy1, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ezy2.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Ezy2, sizeof(PRECISION)*MaxTime);
		snapshot.close();
	}
	cout << endl << "Simulation (CPU) completed!" << endl;
	return 0;
}
int CFDTD2DCloak::CompleteRunCPU(bool SaveFields)
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
// Cylinder related functions.
PRECISION CFDTD2DCloak::A(PRECISION i, PRECISION j)
{
	PRECISION x = (i*delta) - x0;
	PRECISION y = (j*delta) - y0;
	PRECISION r = sqrt(pow(x,2) + pow(y,2));
	PRECISION rA = rb/(rb-ra);

	if (r < rb && r > ra)
	{
		return rA;
	}
	else
		return 1.;
}
PRECISION CFDTD2DCloak::mask(PRECISION i, PRECISION j)
{
	PRECISION x = (i*delta) - x0;
	PRECISION y = (j*delta) - y0;
	PRECISION r = sqrt(pow(x,2) + pow(y,2));

	if (r <= ra)
	{
		return 0.;
	}
	else
		return 1.;
}
PRECISION CFDTD2DCloak::sinphi(PRECISION i, PRECISION j)
{
	PRECISION x = (i*delta) - x0;
	PRECISION y = (j*delta) - y0;
	PRECISION r = sqrt(pow(x,2) + pow(y,2));

	if (r == 0.)
	{
		return 0.;
	}
	else
		return y/r;
}
PRECISION CFDTD2DCloak::cosphi(PRECISION i, PRECISION j)
{
	PRECISION x = (i*delta) - x0;
	PRECISION y = (j*delta) - y0;
	PRECISION r = sqrt(pow(x,2) + pow(y,2));

	if (r == 0.)
	{
		return 0.;
	}
	else
		return x/r;
}
PRECISION CFDTD2DCloak::einf(PRECISION i, PRECISION j)
{
	PRECISION x = (i*delta) - x0;
	PRECISION y = (j*delta) - y0;
	PRECISION r = sqrt(pow(x,2) + pow(y,2));
	PRECISION reinf = 1.;

	if (r < rb && r > ra)
	{
		return reinf;
	}
	else
		return 1.;
}
PRECISION CFDTD2DCloak::uinf(PRECISION i, PRECISION j)
{
	PRECISION x = (i*delta) - x0;
	PRECISION y = (j*delta) - y0;
	PRECISION r = sqrt(pow(x,2) + pow(y,2));
	PRECISION ruinf = 1.;

	if (r < rb && r > ra)
	{
		return ruinf;
	}
	else
		return 1.;
}
PRECISION CFDTD2DCloak::ge(PRECISION i, PRECISION j)
{
	PRECISION x = (i*delta) - x0;
	PRECISION y = (j*delta) - y0;
	PRECISION r = sqrt(pow(x,2) + pow(y,2));
	PRECISION rge = 0.;

	if (r < rb && r > ra)
	{
		return rge;
	}
	else
		return 0.;
}
PRECISION CFDTD2DCloak::ez(PRECISION i, PRECISION j)
{
	PRECISION x = (i*delta) - x0;
	PRECISION y = (j*delta) - y0;
	PRECISION r = sqrt(pow(x,2) + pow(y,2));

	PRECISION rez = pow(rb/(rb-ra),2) * ((r-ra)/r); // Ideal.
	// PRECISION rez = pow(rb/(rb-ra),2.); // Pendry reduced
	// PRECISION rez = rb/(rb-ra); // Zhao reduced.
	// PRECISION rez = 1.; // Free space

	if (r < rb && r > ra)
	{
		return rez;
	}
	else
		return 1.;
}
PRECISION CFDTD2DCloak::ur(PRECISION i, PRECISION j)
{
	PRECISION x = (i*delta) - x0;
	PRECISION y = (j*delta) - y0;
	PRECISION r = sqrt(pow(x,2) + pow(y,2));

	PRECISION rur = (r-ra)/r; // Ideal.
	// PRECISION rur = rb/(rb-1.) * pow((r-ra)/r,2.); // Pendry reduced
	// PRECISION rur = pow(r/(r-ra),2.); // Zhao reduced.
	// PRECISION rur = 1.; // Free space

	if (r < rb && r > ra && rur == rur)
	{
		return rur;
	}
	else
		return 1.;
}
PRECISION CFDTD2DCloak::uphi(PRECISION i, PRECISION j)
{
	PRECISION x = (i*delta) - x0;
	PRECISION y = (j*delta) - y0;
	PRECISION r = sqrt(pow(x,2) + pow(y,2));

	PRECISION ruphi = r/(r-ra); // Ideal.
	// PRECISION ruphi = pow(r/(r-ra),2.); // Pendry reduced
	// PRECISION ruphi = 1.; // Zhao reduced.
	// PRECISION ruphi = 1.; // Free space

	if (r < rb && r > ra)
	{
		return ruphi;
	}
	else
		return 1.;
}
PRECISION CFDTD2DCloak::wpesq(PRECISION i, PRECISION j)
{
	return pow(w,2) * (einf(i,j)-ez(i,j)/A(i,j));
}
PRECISION CFDTD2DCloak::wpmsq(PRECISION i, PRECISION j) // Corrected.
{
	return (2.*sin(w*dt/2.)*(-2.*(ur(i,j)-1.)*sin(w*dt/2.)))/(pow(dt,2)*pow(cos(w*dt/2.),2));
}
// Timing.
void CFDTD2DCloak::StartTimer()
{
	if (tPaused == true)
	{
		tStart = GetTimeus64();
		tPaused = false;
	}
}
void CFDTD2DCloak::StopTimer()
{
	if (tPaused == false)
	{
		tEnd = GetTimeus64();
		tDelta += tEnd - tStart;
		tStart = tEnd;
		tPaused = true;
	}
}
void CFDTD2DCloak::ResetTimer()
{
	if (tPaused == true)
		tStart = tEnd;
	else
		tStart = GetTimeus64();

	tDelta = 0UL;
}
double CFDTD2DCloak::GetElapsedTime()
{
	if (tPaused == false)
		tEnd = GetTimeus64();

	return ((double)(tEnd-tStart+tDelta))/(1000000.);
}
int CFDTD2DCloak::SafeCall(int Status, const char *Error)
{
	if (Status != 0)
	{
		if (Error!=NULL) cout << Error << endl;
		exit(Status);
	}
	return Status;
}
int CFDTD2DCloak::CleanupCPU()
{
	// Field arrays.
	DeleteArray(Ez_);
	DeleteArray(Dz_);
	DeleteArray(EzMask_);
	DeleteArray(Hx_);
	DeleteArray(Bx_);
	DeleteArray(BxAve_);
	DeleteArray(Hy_);
	DeleteArray(By_);
	DeleteArray(ByAve_);
	// Incident and transmitted fields.
	DeleteArray(Ezi);
	DeleteArray(Ezt);
	DeleteArray(Eztt);
	// Refractive index.
	DeleteArray(Ezy1);
	DeleteArray(Ezy2);
	// Auxiliary field scalars.
	DeleteArray(ax_);
	DeleteArray(ay_);
	DeleteArray(az_);
	// PML arrays.
	DeleteArray(PsiEzX_);
	DeleteArray(PsiEzY_);
	DeleteArray(PsiHyX_);
	DeleteArray(PsiHxY_);

	return 0;
}
CFDTD2DCloak::~CFDTD2DCloak()
{
	// Field arrays.
	DeleteArray(Ez_);
	DeleteArray(Dz_);
	DeleteArray(EzMask_);
	DeleteArray(Hx_);
	DeleteArray(Bx_);
	DeleteArray(BxAve_);
	DeleteArray(Hy_);
	DeleteArray(By_);
	DeleteArray(ByAve_);
	// Incident and transmitted fields.
	DeleteArray(Ezi);
	DeleteArray(Ezt);
	DeleteArray(Eztt);
	// Refractive index.
	DeleteArray(Ezy1);
	DeleteArray(Ezy2);
	// Auxiliary field scalars.
	DeleteArray(ax_);
	DeleteArray(ay_);
	DeleteArray(az_);
	// PML arrays.
	DeleteArray(PsiEzX_);
	DeleteArray(PsiEzY_);
	DeleteArray(PsiHyX_);
	DeleteArray(PsiHxY_);
}
