// Macros for indexing data arrays.
#define Ex(i,n) Ex_[(i)+Size*(n)]
#define Dx(i,n) Dx_[(i)+Size*(n)]
#define Hy(i,n) Hy_[(i)+Size*(n)]
#define By(i,n) By_[(i)+Size*(n)]

#define PI	3.14159265358979323846
#define PRECISION double

#ifdef KHR_DP_EXTENSION
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#else
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

// Number of const uint or double arguments does NOT have any significant impact on performance.
__kernel void FDTD1DDNGKernel_DryRun(
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
							__global PRECISION *Ex_, __global PRECISION *Hy_,
							// Incident field.
							__global PRECISION *Exi,
							const unsigned int x1,
							// Time indices.
							const unsigned int n,
							const unsigned int np,
							const unsigned int n0,
							const unsigned int nf)
{
	unsigned int i = get_global_id(0);

	if (i != Size-1) // Normal update equation.
		Hy(i,nf) = Hy(i,n0) + (Ex(i,n0)-Ex(i+1,n0))*dt/(u0*dz);

	barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);

	// ABC
	if (i == Size-1)
		Hy(i,nf) = Hy(i-1,n0) + (Sc-1)/(Sc+1)*(Hy(i-1,nf)-Hy(i,n0));

	barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);

	if (i != 0)
		Ex(i,nf) = Ex(i,n0) + (Hy(i-1,nf)-Hy(i,nf))*dt/(e0*dz);

	barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);

	// ABC
	if (i == 0)
		Ex(i,nf) = Ex(i+1,n0) + (Sc-1)/(Sc+1)*(Ex(i+1,nf)-Ex(i,n0));

	barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);

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
// Simulation kernel
__kernel void FDTD1DDNGKernel_Simulation(
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
							__global PRECISION *Ex_, __global PRECISION *Dx_, __global PRECISION *Hy_, __global PRECISION *By_,
							// Drude material parameters.
							__global PRECISION *einf, __global PRECISION *uinf, __global PRECISION *wpesq, __global PRECISION *wpmsq, __global PRECISION *ge, __global PRECISION *gm,
							// Drude scalars.
							__global PRECISION *ae0, __global PRECISION *ae, __global PRECISION *be, __global PRECISION *ce, __global PRECISION *de, __global PRECISION *ee,
							__global PRECISION *am0, __global PRECISION *am, __global PRECISION *bm, __global PRECISION *cm, __global PRECISION *dm, __global PRECISION *em,
							// Incident field.
							__global PRECISION *Ext,
							__global PRECISION *Extt,
							__global PRECISION *Exz1,
							__global PRECISION *Exz2,
							const unsigned int x1,
							const unsigned int Z1,
							const unsigned int Z2,
							// Time indices.
							const unsigned int n,
							const unsigned int np,
							const unsigned int n0,
							const unsigned int nf)
{
	unsigned int i = get_global_id(0);

	if (i != Size-1) // Normal update equation.
	{
		By(i,nf) = By(i,n0) + (Ex(i,n0)-Ex(i+1,n0))*dt/dz;
		Hy(i,nf) = am[i]*(By(i,nf)-2*By(i,n0)+By(i,np)) + bm[i]*(By(i,nf)-By(i,np)) + cm[i]*(2*Hy(i,n0)-Hy(i,np)) + dm[i]*(2*Hy(i,n0)+Hy(i,np)) + em[i]*(Hy(i,np));
	}
	barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);

	// ABC
	if (i == Size-1)
	{
		Hy(i,nf) = Hy(i-1,n0) + (Sc-1)/(Sc+1)*(Hy(i-1,nf)-Hy(i,n0));
		By(i,nf) = u0*Hy(i,nf);
	}
	barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);

	if (i != 0)
	{
		Dx(i,nf) = Dx(i,n0) + (Hy(i-1,nf)-Hy(i,nf))*dt/dz;
		Ex(i,nf) = ae[i]*(Dx(i,nf)-2*Dx(i,n0)+Dx(i,np)) + be[i]*(Dx(i,nf)-Dx(i,np)) + ce[i]*(2*Ex(i,n0)-Ex(i,np)) + de[i]*(2*Ex(i,n0)+Ex(i,np)) + ee[i]*(Ex(i,np));
	}
	barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);

	// ABC
	if (i == 0)
	{
		Ex(i,nf) = Ex(i+1,n0) + (Sc-1)/(Sc+1)*(Ex(i+1,nf)-Ex(i,n0));
		Dx(i,nf) = e0*Ex(i,nf);
	}
	barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);

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
