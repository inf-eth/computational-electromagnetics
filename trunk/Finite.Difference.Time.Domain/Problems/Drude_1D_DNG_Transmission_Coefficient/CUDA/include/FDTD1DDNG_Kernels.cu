// Macros for indexing data arrays.
#define Ex(i,n) Ex_[(i)+Size*(n)]
#define Dx(i,n) Dx_[(i)+Size*(n)]
#define Hy(i,n) Hy_[(i)+Size*(n)]
#define By(i,n) By_[(i)+Size*(n)]

#include <FDTD1DDNG.hpp>

// Dry run kernel.
template <unsigned int BlockX, unsigned int BlockY> __global__ void FDTD1DDNGKernel_DryRun(
							const unsigned int Size,
							const unsigned int PulseWidth,
							const unsigned int td,
							const unsigned int SourceLocation,
							const unsigned int SourceChoice,
							const PRECISION e0,
							const PRECISION u0,
							const PRECISION dt,
							const PRECISION dz,
							const PRECISION Sc,
							// Frequency, wavelength, wave number.
							const PRECISION f,
							const PRECISION fp,
							const PRECISION dr,
							// Data arrays.
							PRECISION *Ex_, PRECISION *Hy_,
							// Incident field.
							PRECISION *Exi,
							const unsigned int x1,
							// Time indices.
							const unsigned int n,
							const unsigned int np,
							const unsigned int n0,
							const unsigned int nf)
{
	unsigned int i = BlockX*blockIdx.x+threadIdx.x;

	if (i != Size-1) // Normal update equation.
	{
		Hy(i,nf) = Hy(i,n0) + (Ex(i,n0)-Ex(i+1,n0))*dt/(u0*dz);
	}
	else // ABC
	{
		__syncthreads();
		Hy(i,nf) = Hy(i-1,n0) + (Sc-1)/(Sc+1)*(Hy(i-1,nf)-Hy(i,n0));
	}

	__syncthreads();

	if (i != 0)
	{
		Ex(i,nf) = Ex(i,n0) + (Hy(i-1,nf)-Hy(i,nf))*dt/(e0*dz);
	}
	else // ABC
	{
		__syncthreads();
		Ex(i,nf) = Ex(i+1,n0) + (Sc-1)/(Sc+1)*(Ex(i+1,nf)-Ex(i,n0));
	}

	// Source.
	if (i == SourceLocation)
	{
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
	}
	// Recording incident Field.
	if (i == x1)
		Exi[n] = Ex(i,nf);
}

#define ae_(i) (1./pow(dt,2))*ae0[i]
#define be_(i) (1./(2.*dt))*ge[i]*ae0[i]
#define ce_(i) (e0/pow(dt,2))*einf[i]*ae0[i]
#define de_(i) (-1.*e0/4.)*wpesq[i]*ae0[i]
#define ee_(i) (1./(2.*dt))*e0*einf[i]*ge[i]*ae0[i]
#define am_(i) (1./pow(dt,2))*am0[i]
#define bm_(i) (1./(2.*dt))*gm[i]*am0[i]
#define cm_(i) (u0/pow(dt,2))*uinf[i]*am0[i]
#define dm_(i) (-1.*u0/4.)*wpmsq[i]*am0[i]
#define em_(i) (1./(2.*dt))*u0*uinf[i]*gm[i]*am0[i]
			
// Simulation kernel.
template <unsigned int BlockX, unsigned int BlockY> __global__ void FDTD1DDNGKernel_Simulation(
							const unsigned int Size,
							const unsigned int PulseWidth,
							const unsigned int td,
							const unsigned int SourceLocation,
							const unsigned int SourceChoice,
							const PRECISION e0,
							const PRECISION u0,
							const PRECISION dt,
							const PRECISION dz,
							const PRECISION Sc,
							// Frequency, wavelength, wave number.
							const PRECISION f,
							const PRECISION fp,
							const PRECISION dr,
							// Data arrays.
							PRECISION *Ex_, PRECISION *Dx_, PRECISION *Hy_, PRECISION *By_,
							// Drude material parameters.
							PRECISION *einf, PRECISION *uinf, PRECISION *wpesq, PRECISION *wpmsq, PRECISION *ge, PRECISION *gm,
							// Drude scalars.
							PRECISION *ae0, //PRECISION *ae, PRECISION *be, PRECISION *ce, PRECISION *de, PRECISION *ee,
							PRECISION *am0, //PRECISION *am, PRECISION *bm, PRECISION *cm, PRECISION *dm, PRECISION *em,
							// Incident field.
							PRECISION *Ext,
							PRECISION *Extt,
							PRECISION *Exz1,
							PRECISION *Exz2,
							const unsigned int x1,
							const unsigned int Z1,
							const unsigned int Z2,
							// Time indices.
							const unsigned int n,
							const unsigned int np,
							const unsigned int n0,
							const unsigned int nf)
{
	unsigned int i = BlockX*blockIdx.x+threadIdx.x;

	if (i != Size-1) // Normal update equation.
	{
		By(i,nf) = By(i,n0) + (Ex(i,n0)-Ex(i+1,n0))*dt/dz;
		Hy(i,nf) = am_(i)*(By(i,nf)-2*By(i,n0)+By(i,np)) + bm_(i)*(By(i,nf)-By(i,np)) + cm_(i)*(2*Hy(i,n0)-Hy(i,np)) + dm_(i)*(2*Hy(i,n0)+Hy(i,np)) + em_(i)*(Hy(i,np));
	}
	else // ABC
	{
		__syncthreads();
		Hy(i,nf) = Hy(i-1,n0) + (Sc-1)/(Sc+1)*(Hy(i-1,nf)-Hy(i,n0));
		By(i,nf) = u0*Hy(i,nf);
	}

	__syncthreads();

	if (i != 0)
	{
		Dx(i,nf) = Dx(i,n0) + (Hy(i-1,nf)-Hy(i,nf))*dt/dz;
		Ex(i,nf) = ae_(i)*(Dx(i,nf)-2*Dx(i,n0)+Dx(i,np)) + be_(i)*(Dx(i,nf)-Dx(i,np)) + ce_(i)*(2*Ex(i,n0)-Ex(i,np)) + de_(i)*(2*Ex(i,n0)+Ex(i,np)) + ee_(i)*(Ex(i,np));
	}
	else // ABC
	{
		__syncthreads();
		Ex(i,nf) = Ex(i+1,n0) + (Sc-1)/(Sc+1)*(Ex(i+1,nf)-Ex(i,n0));
		Dx(i,nf) = e0*Ex(i,nf);
	}

	// Source.
	if (i == SourceLocation)
	{
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
	}
	// Recording transmitted Fields.
	if (i == x1)
		Ext[n] = Ex(i,nf);
	if (i == (Size-(2*Size/3))+10)
		Extt[n] = Ex(i,nf);
	if (i == Z1)
		Exz1[n] = Ex(i,nf);
	if (i == Z2)
		Exz2[n] = Ex(i,nf);
}
