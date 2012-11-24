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

#define PRECISION double

#ifdef KHR_DP_EXTENSION
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#else
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

// Number of const uint or double arguments does NOT have any significant impact on performance.
__kernel void FDTD2DDNGKernel_DryRun_M(
							const unsigned int I,
							const unsigned int J,
							const unsigned int PMLw,
							const unsigned int PulseWidth,
							const unsigned int td,
							const unsigned int SourceChoice,
							const unsigned int SourcePlane,
							const unsigned int SourceLocationX,
							const unsigned int SourceLocationY,
							const PRECISION c,
							const PRECISION pi,
							const PRECISION e0,
							const PRECISION u0,
							const PRECISION dt,
							const PRECISION delta,
							const PRECISION Sc,
							const PRECISION f,
							const PRECISION fp,
							const PRECISION dr,
							const unsigned int IEz, const unsigned int JEz,
							const unsigned int IHx, const unsigned int JHx,
							const unsigned int IHy, const unsigned int JHy,
							const __global PRECISION *Ez_, const __global PRECISION *Dz_, __global PRECISION *Hx_, __global PRECISION *Bx_, __global PRECISION *Hy_, __global PRECISION *By_,
							const __global PRECISION *PsiEzX_, const __global PRECISION *PsiEzY_, __global PRECISION *PsiHyX_, __global PRECISION *PsiHxY_,
							const __global PRECISION *einf_, const __global PRECISION *uinf_,
							const PRECISION kappex, const PRECISION kappey, const PRECISION kappmx, const PRECISION kappmy,
							const PRECISION bex, const PRECISION bey, const PRECISION bmx, const PRECISION bmy,
							const PRECISION Cex, const PRECISION Cey, const PRECISION Cmx, const PRECISION Cmy,
							__global PRECISION *Ezi,
							const unsigned int x1,
							const unsigned int n,
							const unsigned int np,
							const unsigned int n0,
							const unsigned int nf)
{
	const unsigned int i = get_global_id(0);
	const unsigned int j = get_global_id(1);

	// Calculation of PsiHxY.
	if (j>0 && j<JHx-1)
		PsiHxY(i,j) = (Cmy/delta)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) + bmy*PsiHxY(i,j);

	// Bx in normal space.
	if( j>PMLw && j<JHx-PMLw-1)
	{
		Bx(i,j,nf) = Bx(i,j,n0) + (-Ez(i,j,n0) + Ez(i,j-1,n0)) * dt/delta;
		Hx(i,j,nf) = Bx(i,j,nf)/(u0*uinf(i,j));
	}
	// Bx in lower PML.
	if (j>0 && j<PMLw+1)
	{
		Bx(i,j,nf) = Bx(i,j,n0) + dt*((1./kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1./delta + PsiHxY(i,j));
		Hx(i,j,nf) = Bx(i,j,nf)/(u0*uinf(i,j));
	}
	// Bx in upper PML.
	if (j>JHx-PMLw-2 && j<JHx-1)
	{
		Bx(i,j,nf) = Bx(i,j,n0) + dt*((1./kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1./delta + PsiHxY(i,j));
		Hx(i,j,nf) = Bx(i,j,nf)/(u0*uinf(i,j));
	}

	// ========================== By and Hy ==========================
	if (i<IHy-1)
	{
		// PsiHyX array.
		if (j<JHy)
			PsiHyX(i,j) = (Cmx/delta)*(Ez(i+1,j,n0)-Ez(i,j,n0)) + bmx*PsiHyX(i,j);

		// By in normal space.
		if (j>PMLw-1 && j<JHy-PMLw)
		{
			By(i,j,nf) = By(i,j,n0) + (Ez(i+1,j,n0) - Ez(i,j,n0)) * dt/delta;
			Hy(i,j,nf) = By(i,j,nf)/(u0*uinf(i,j));
		}
		// By in Lower PML.
		if (j<PMLw)
		{
			By(i,j,nf) = By(i,j,n0) + dt*((1./kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1./delta + PsiHyX(i,j));
			Hy(i,j,nf) = By(i,j,nf)/(u0*uinf(i,j));
		}
		// By in upper PML.
		if (j>JHy-PMLw-1 && j<JHy)
		{
			By(i,j,nf) = By(i,j,n0) + dt*((1./kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1./delta + PsiHyX(i,j));
			Hy(i,j,nf) = By(i,j,nf)/(u0*uinf(i,j));
		}
	}
	else
	{
		// PsiHyX array.
		if (j<JHy)
			PsiHyX(IHy-1,j) = (Cmx/delta)*(Ez(0,j,n0)-Ez(IHy-1,j,n0)) + bmx*PsiHyX(IHy-1,j); // PBC

		// By in normal space.
		if (j>PMLw-1 && j<JHy-PMLw)
		{
			By(IHy-1,j,nf) = By(IHy-1,j,n0) + (Ez(0,j,n0) - Ez(IHy-1,j,n0)) * dt/delta; // PBC
			Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/(u0*uinf(IHy-1,j)); // PBC
		}
		// By in Lower PML.
		if (j<PMLw)
		{
			By(IHy-1,j,nf) = By(IHy-1,j,n0) + dt*((1./kappmx)*(Ez(0,j,n0) - Ez(IHy-1,j,n0)) * 1./delta + PsiHyX(IHy-1,j)); // PBC
			Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/(u0*uinf(IHy-1,j)); // PBC
		}
		// By in upper PML.
		if (j>JHy-PMLw-1 && j<JHy)
		{
			By(IHy-1,j,nf) = By(IHy-1,j,n0) + dt*((1./kappmx)*(Ez(0,j,n0) - Ez(IHy-1,j,n0)) * 1./delta + PsiHyX(IHy-1,j)); // PBC
			Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/(u0*uinf(IHy-1,j)); // PBC
		}
	}
}
__kernel void FDTD2DDNGKernel_DryRun_E(
							const unsigned int I,
							const unsigned int J,
							const unsigned int PMLw,
							const unsigned int PulseWidth,
							const unsigned int td,
							const unsigned int SourceChoice,
							const unsigned int SourcePlane,
							const unsigned int SourceLocationX,
							const unsigned int SourceLocationY,
							const PRECISION c,
							const PRECISION pi,
							const PRECISION e0,
							const PRECISION u0,
							const PRECISION dt,
							const PRECISION delta,
							const PRECISION Sc,
							const PRECISION f,
							const PRECISION fp,
							const PRECISION dr,
							const unsigned int IEz, const unsigned int JEz,
							const unsigned int IHx, const unsigned int JHx,
							const unsigned int IHy, const unsigned int JHy,
							__global PRECISION *Ez_, __global PRECISION *Dz_, const __global PRECISION *Hx_, const __global PRECISION *Bx_, const __global PRECISION *Hy_, const __global PRECISION *By_,
							__global PRECISION *PsiEzX_, __global PRECISION *PsiEzY_, const __global PRECISION *PsiHyX_, const __global PRECISION *PsiHxY_,
							const __global PRECISION *einf_, const __global PRECISION *uinf_,
							const PRECISION kappex, const PRECISION kappey, const PRECISION kappmx, const PRECISION kappmy,
							const PRECISION bex, const PRECISION bey, const PRECISION bmx, const PRECISION bmy,
							const PRECISION Cex, const PRECISION Cey, const PRECISION Cmx, const PRECISION Cmy,
							__global PRECISION *Ezi,
							const unsigned int x1,
							const unsigned int n,
							const unsigned int np,
							const unsigned int n0,
							const unsigned int nf)
{
	const unsigned int i = get_global_id(0);
	const unsigned int j = get_global_id(1);

	// ========================== Dz and Ez ==========================
	if (i>0)
	{
		// Psi arrays.
		if (j<JEz)
		{
			PsiEzX(i,j) = (Cex/delta)*(Hy(i,j,nf)-Hy(i-1,j,nf)) + bex*PsiEzX(i,j);
			PsiEzY(i,j) = (Cey/delta)*(-Hx(i,j+1,nf)+Hx(i,j,nf)) + bey*PsiEzY(i,j);
		}
		// Dz in normal space.
		if (j>PMLw-1 && j<JEz-PMLw)
		{
			Dz(i,j,nf) = Dz(i,j,n0) + (Hy(i,j,nf)-Hy(i-1,j,nf)-Hx(i,j+1,nf)+Hx(i,j,nf)) * dt/delta;
			Ez(i,j,nf) = Dz(i,j,nf)/(e0*einf(i,j));
		}
		// Dz in lower PML.
		if (j<PMLw)
		{
			Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1./kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1./kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1./delta + PsiEzX(i,j) + PsiEzY(i,j));
			Ez(i,j,nf) = Dz(i,j,nf)/(e0*einf(i,j));
		}
		// Dz in upper PML.
		if (j>JEz-PMLw-1 && j<JEz)
		{
			Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1./kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1./kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1./delta + PsiEzX(i,j) + PsiEzY(i,j));
			Ez(i,j,nf) = Dz(i,j,nf)/(e0*einf(i,j));
		}
	}
	else
	{
		// Psi arrays.
		if (j<JEz)
		{
			PsiEzX(0,j) = (Cex/delta)*(Hy(0,j,nf)-Hy(IEz-1,j,nf)) + bex*PsiEzX(0,j); // PBC
			PsiEzY(0,j) = (Cey/delta)*(-Hx(0,j+1,nf)+Hx(0,j,nf)) + bey*PsiEzY(0,j); // PBC
		}
		// Dz in normal space.
		if (j>PMLw-1 && j<JEz-PMLw)
		{
			Dz(0,j,nf) = Dz(0,j,n0) + (Hy(0,j,nf)-Hy(IEz-1,j,nf)-Hx(0,j+1,nf)+Hx(0,j,nf)) * dt/delta; // PBC
			Ez(0,j,nf) = Dz(0,j,nf)/(e0*einf(0,j)); // PBC
		}
		// Dz in lower PML.
		if (j<PMLw)
		{
			Dz(0,j,nf) = Dz(0,j,n0) + dt*(((1./kappex)*(Hy(0,j,nf)-Hy(IEz-1,j,nf))+(1./kappey)*(-Hx(0,j+1,nf)+Hx(0,j,nf))) * 1./delta + PsiEzX(0,j) + PsiEzY(0,j)); // PBC
			Ez(0,j,nf) = Dz(0,j,nf)/(e0*einf(0,j)); // PBC
		}
		// Dz in upper PML.
		if (j>JEz-PMLw-1 && j<JEz)
		{
			Dz(0,j,nf) = Dz(0,j,n0) + dt*(((1./kappex)*(Hy(0,j,nf)-Hy(IEz-1,j,nf))+(1./kappey)*(-Hx(0,j+1,nf)+Hx(0,j,nf))) * 1./delta + PsiEzX(0,j) + PsiEzY(0,j)); // PBC
			Ez(0,j,nf) = Dz(0,j,nf)/(e0*einf(0,j)); // PBC
		}
	}

	// ====================== Source ======================
	if (SourcePlane == 1 && i<IEz && j == SourceLocationY)
	{
		if (SourceChoice == 1)
			Ez(i,j,nf) += exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
		else if (SourceChoice == 2)
			Ez(i,j,nf) += sin(2.*pi*f*(PRECISION)n*dt) * Sc;
		else if (SourceChoice == 3)
			Ez(i,j,nf) += (1.-2.*pow(pi*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(pi*fp*((PRECISION)n*dt-dr),2)) * Sc;

		Dz(i,j,nf) = e0*Ez(i,j,nf);
	}
	else if (i == SourceLocationX && j == SourceLocationY)
	{
		if (SourceChoice == 1)
			Ez(i,j,nf) += exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
		else if (SourceChoice == 2)
			Ez(i,j,nf) += sin(2.*pi*f*(PRECISION)n*dt) * Sc;
		else if (SourceChoice == 3)
			Ez(i,j,nf) += (1.-2.*pow(pi*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(pi*fp*((PRECISION)n*dt-dr),2)) * Sc;

		Dz(i,j,nf) = e0*Ez(i,j,nf);
	}
	if (j==x1)
		Ezi[n] = Ez(IEz/2,x1,nf); // Incident field.
}
// Simulation kernel
__kernel void FDTD2DDNGKernel_Simulation_M(
							const unsigned int I,
							const unsigned int J,
							const unsigned int PMLw,
							const unsigned int PulseWidth,
							const unsigned int td,
							const unsigned int SourceChoice,
							const unsigned int SourcePlane,
							const unsigned int SourceLocationX,
							const unsigned int SourceLocationY,
							const PRECISION c,
							const PRECISION pi,
							const PRECISION e0,
							const PRECISION u0,
							const PRECISION dt,
							const PRECISION delta,
							const PRECISION Sc,
							const PRECISION f,
							const PRECISION fp,
							const PRECISION dr,
							const unsigned int IEz, const unsigned int JEz,
							const unsigned int IHx, const unsigned int JHx,
							const unsigned int IHy, const unsigned int JHy,
							const __global PRECISION *Ez_, const __global PRECISION *Dz_, __global PRECISION *Hx_, __global PRECISION *Bx_, __global PRECISION *Hy_, __global PRECISION *By_,
							const __global PRECISION *einf_, const __global PRECISION *uinf_, const __global PRECISION *wpesq_, const __global PRECISION *wpmsq_, const __global PRECISION *ge_, const __global PRECISION *gm_,
							const __global PRECISION *ae0_, const __global PRECISION *ae_, const __global PRECISION *be_, const __global PRECISION *ce_, const __global PRECISION *de_, const __global PRECISION *ee_,
							const __global PRECISION *am0_, const __global PRECISION *am_, const __global PRECISION *bm_, const __global PRECISION *cm_, const __global PRECISION *dm_, const __global PRECISION *em_,
							const __global PRECISION *PsiEzX_, const __global PRECISION *PsiEzY_, __global PRECISION *PsiHyX_, __global PRECISION *PsiHxY_,
							const PRECISION kappex, const PRECISION kappey, const PRECISION kappmx, const PRECISION kappmy,
							const PRECISION bex, const PRECISION bey, const PRECISION bmx, const PRECISION bmy,
							const PRECISION Cex, const PRECISION Cey, const PRECISION Cmx, const PRECISION Cmy,
							__global PRECISION *Ezt, __global PRECISION *Eztt, __global PRECISION *Ezy1, __global PRECISION *Ezy2,
							const unsigned int x1, const unsigned int Y1, const unsigned int Y2,
							const unsigned int n,
							const unsigned int np,
							const unsigned int n0,
							const unsigned int nf)
{
	const unsigned int i = get_global_id(0);
	const unsigned int j = get_global_id(1);

	// Calculation of PsiHxY.
	if (j>0 && j<JHx-1)
		PsiHxY(i,j) = (Cmy/delta)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) + bmy*PsiHxY(i,j);

	// Bx in normal space.
	if (j>PMLw && j<JHx-PMLw-1)
	{
		Bx(i,j,nf) = Bx(i,j,n0) + (-Ez(i,j,n0) + Ez(i,j-1,n0)) * dt/delta;
		Hx(i,j,nf) = am(i,j)*(Bx(i,j,nf)-2.*Bx(i,j,n0)+Bx(i,j,np))+bm(i,j)*(Bx(i,j,nf)-Bx(i,j,np))+cm(i,j)*(2.*Hx(i,j,n0)-Hx(i,j,np))+dm(i,j)*(2.*Hx(i,j,n0)+Hx(i,j,np))+em(i,j)*Hx(i,j,np);
	}
	// Bx in lower PML.
	if (j>0 && j<PMLw+1)
	{
		Bx(i,j,nf) = Bx(i,j,n0) + dt*((1./kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1./delta + PsiHxY(i,j));
		Hx(i,j,nf) = Bx(i,j,nf)/(u0*uinf(i,j));
	}
	// Bx in upper PML.
	if (j>JHx-PMLw-2 && j<JHx-1)
	{
		Bx(i,j,nf) = Bx(i,j,n0) + dt*((1./kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1./delta + PsiHxY(i,j));
		Hx(i,j,nf) = Bx(i,j,nf)/(u0*uinf(i,j));
	}


	// ========================== By and Hy ==========================
	if (i<IHy-1)
	{
		// PsiHyX arrays.
		if (j<JHy)
			PsiHyX(i,j) = (Cmx/delta)*(Ez(i+1,j,n0)-Ez(i,j,n0)) + bmx*PsiHyX(i,j);

		// By in normal space.
		if (j>PMLw-1 && j<JHy-PMLw)
		{
			By(i,j,nf) = By(i,j,n0) + (Ez(i+1,j,n0) - Ez(i,j,n0)) * dt/delta;
			Hy(i,j,nf) = am(i,j)*(By(i,j,nf)-2.*By(i,j,n0)+By(i,j,np))+bm(i,j)*(By(i,j,nf)-By(i,j,np))+cm(i,j)*(2.*Hy(i,j,n0)-Hy(i,j,np))+dm(i,j)*(2.*Hy(i,j,n0)+Hy(i,j,np))+em(i,j)*Hy(i,j,np);
		}
		// By in Lower PML.
		if (j<PMLw)
		{
			By(i,j,nf) = By(i,j,n0) + dt*((1./kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1./delta + PsiHyX(i,j));
			Hy(i,j,nf) = By(i,j,nf)/(u0*uinf(i,j));
		}
		// By in upper PML.
		if (j>JHy-PMLw-1 && j<JHy)
		{
			By(i,j,nf) = By(i,j,n0) + dt*((1./kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1./delta + PsiHyX(i,j));
			Hy(i,j,nf) = By(i,j,nf)/(u0*uinf(i,j));
		}
	}
	else
	{
		// PsiHyX array.
		if (j<JHy)
			PsiHyX(IHy-1,j) = (Cmx/delta)*(Ez(0,j,n0)-Ez(IHy-1,j,n0)) + bmx*PsiHyX(IHy-1,j); // PBC

		// By in normal space.
		if (j>PMLw-1 && j<JHy-PMLw)
		{
			By(IHy-1,j,nf) = By(IHy-1,j,n0) + (Ez(0,j,n0) - Ez(IHy-1,j,n0)) * dt/delta; // PBC
			Hy(IHy-1,j,nf) = am(IHy-1,j)*(By(IHy-1,j,nf)-2.*By(IHy-1,j,n0)+By(IHy-1,j,np))+bm(IHy-1,j)*(By(IHy-1,j,nf)-By(IHy-1,j,np))+cm(IHy-1,j)*(2.*Hy(IHy-1,j,n0)-Hy(IHy-1,j,np))+dm(IHy-1,j)*(2.*Hy(IHy-1,j,n0)+Hy(IHy-1,j,np))+em(IHy-1,j)*Hy(IHy-1,j,np); // PBC
		}
		// By in Lower PML.
		if (j<PMLw)
		{
			By(IHy-1,j,nf) = By(IHy-1,j,n0) + dt*((1./kappmx)*(Ez(0,j,n0) - Ez(IHy-1,j,n0)) * 1./delta + PsiHyX(IHy-1,j)); // PBC
			Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/(u0*uinf(IHy-1,j)); // PBC
		}
		// By in upper PML.
		if (j>JHy-PMLw-1 && j<JHy)
		{
			By(IHy-1,j,nf) = By(IHy-1,j,n0) + dt*((1./kappmx)*(Ez(0,j,n0) - Ez(IHy-1,j,n0)) * 1./delta + PsiHyX(IHy-1,j)); // PBC
			Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/(u0*uinf(IHy-1,j)); // PBC
		}
	}
}
__kernel void FDTD2DDNGKernel_Simulation_E(
							const unsigned int I,
							const unsigned int J,
							const unsigned int PMLw,
							const unsigned int PulseWidth,
							const unsigned int td,
							const unsigned int SourceChoice,
							const unsigned int SourcePlane,
							const unsigned int SourceLocationX,
							const unsigned int SourceLocationY,
							const PRECISION c,
							const PRECISION pi,
							const PRECISION e0,
							const PRECISION u0,
							const PRECISION dt,
							const PRECISION delta,
							const PRECISION Sc,
							const PRECISION f,
							const PRECISION fp,
							const PRECISION dr,
							const unsigned int IEz, const unsigned int JEz,
							const unsigned int IHx, const unsigned int JHx,
							const unsigned int IHy, const unsigned int JHy,
							__global PRECISION *Ez_, __global PRECISION *Dz_, const __global PRECISION *Hx_, const __global PRECISION *Bx_, const __global PRECISION *Hy_, const __global PRECISION *By_,
							const __global PRECISION *einf_, const __global PRECISION *uinf_, const __global PRECISION *wpesq_, const __global PRECISION *wpmsq_, const __global PRECISION *ge_, const __global PRECISION *gm_,
							const __global PRECISION *ae0_, const __global PRECISION *ae_, const __global PRECISION *be_, const __global PRECISION *ce_, const __global PRECISION *de_, const __global PRECISION *ee_,
							const __global PRECISION *am0_, const __global PRECISION *am_, const __global PRECISION *bm_, const __global PRECISION *cm_, const __global PRECISION *dm_, const __global PRECISION *em_,
							__global PRECISION *PsiEzX_, __global PRECISION *PsiEzY_, const __global PRECISION *PsiHyX_, const __global PRECISION *PsiHxY_,
							const PRECISION kappex, const PRECISION kappey, const PRECISION kappmx, const PRECISION kappmy,
							const PRECISION bex, const PRECISION bey, const PRECISION bmx, const PRECISION bmy,
							const PRECISION Cex, const PRECISION Cey, const PRECISION Cmx, const PRECISION Cmy,
							__global PRECISION *Ezt, __global PRECISION *Eztt, __global PRECISION *Ezy1, __global PRECISION *Ezy2,
							const unsigned int x1, const unsigned int Y1, const unsigned int Y2,
							const unsigned int n,
							const unsigned int np,
							const unsigned int n0,
							const unsigned int nf)
{
	const unsigned int i = get_global_id(0);
	const unsigned int j = get_global_id(1);

	// ========================== Dz and Ez ==========================
	if (i>0)
	{
		// Psi arrays.
		if (j<JEz)
		{
			PsiEzX(i,j) = (Cex/delta)*(Hy(i,j,nf)-Hy(i-1,j,nf)) + bex*PsiEzX(i,j);
			PsiEzY(i,j) = (Cey/delta)*(-Hx(i,j+1,nf)+Hx(i,j,nf)) + bey*PsiEzY(i,j);
		}
		// Dz in normal space.
		if (j>PMLw-1 && j<JEz-PMLw)
		{
			Dz(i,j,nf) = Dz(i,j,n0) + (Hy(i,j,nf)-Hy(i-1,j,nf)-Hx(i,j+1,nf)+Hx(i,j,nf)) * dt/delta;
			Ez(i,j,nf) = ae(i,j)*(Dz(i,j,nf)-2.*Dz(i,j,n0)+Dz(i,j,np))+be(i,j)*(Dz(i,j,nf)-Dz(i,j,np))+ce(i,j)*(2.*Ez(i,j,n0)-Ez(i,j,np))+de(i,j)*(2.*Ez(i,j,n0)+Ez(i,j,np))+ee(i,j)*Ez(i,j,np);
		}
		// Dz in lower PML.
		if (j<PMLw)
		{
			Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1./kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1./kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1./delta + PsiEzX(i,j) + PsiEzY(i,j));
			Ez(i,j,nf) = Dz(i,j,nf)/(e0*einf(i,j));
		}
		// Dz in upper PML.
		if (j>JEz-PMLw-1 && j<JEz)
		{
			Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1./kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1./kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1./delta + PsiEzX(i,j) + PsiEzY(i,j));
			Ez(i,j,nf) = Dz(i,j,nf)/(e0*einf(i,j));
		}
	}
	else
	{
		// Psi arrays.
		if (j<JEz)
		{
			PsiEzX(0,j) = (Cex/delta)*(Hy(0,j,nf)-Hy(IEz-1,j,nf)) + bex*PsiEzX(0,j); // PBC
			PsiEzY(0,j) = (Cey/delta)*(-Hx(0,j+1,nf)+Hx(0,j,nf)) + bey*PsiEzY(0,j); // PBC
		}
		// Dz in normal space.
		if (j>PMLw-1 && j<JEz-PMLw)
		{
			Dz(0,j,nf) = Dz(0,j,n0) + (Hy(0,j,nf)-Hy(IEz-1,j,nf)-Hx(0,j+1,nf)+Hx(0,j,nf)) * dt/delta; // PBC
			Ez(0,j,nf) = ae(0,j)*(Dz(0,j,nf)-2.*Dz(0,j,n0)+Dz(0,j,np))+be(0,j)*(Dz(0,j,nf)-Dz(0,j,np))+ce(0,j)*(2.*Ez(0,j,n0)-Ez(0,j,np))+de(0,j)*(2.*Ez(0,j,n0)+Ez(0,j,np))+ee(0,j)*Ez(0,j,np); // PBC
		}
		// Dz in lower PML.
		if (j<PMLw)
		{
			Dz(0,j,nf) = Dz(0,j,n0) + dt*(((1./kappex)*(Hy(0,j,nf)-Hy(IEz-1,j,nf))+(1./kappey)*(-Hx(0,j+1,nf)+Hx(0,j,nf))) * 1./delta + PsiEzX(0,j) + PsiEzY(0,j)); // PBC
			Ez(0,j,nf) = Dz(0,j,nf)/(e0*einf(0,j)); // PBC
		}
		// Dz in upper PML.
		if (j>JEz-PMLw-1 && j<JEz)
		{
			Dz(0,j,nf) = Dz(0,j,n0) + dt*(((1./kappex)*(Hy(0,j,nf)-Hy(IEz-1,j,nf))+(1./kappey)*(-Hx(0,j+1,nf)+Hx(0,j,nf))) * 1./delta + PsiEzX(0,j) + PsiEzY(0,j)); // PBC
			Ez(0,j,nf) = Dz(0,j,nf)/(e0*einf(0,j)); // PBC
		}
	}

	// ====================== Source ======================
	if (SourcePlane == 1 && i<IEz && j == SourceLocationY)
	{
		if (SourceChoice == 1)
			Ez(i,j,nf) += exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
		else if (SourceChoice == 2)
			Ez(i,j,nf) += sin(2.*pi*f*(PRECISION)n*dt) * Sc;
		else if (SourceChoice == 3)
			Ez(i,j,nf) += (1.-2.*pow(pi*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(pi*fp*((PRECISION)n*dt-dr),2)) * Sc;

		Dz(i,j,nf) = e0*Ez(i,j,nf);
	}
	else if (i == SourceLocationX && j == SourceLocationY)
	{
		if (SourceChoice == 1)
			Ez(i,j,nf) += exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
		else if (SourceChoice == 2)
			Ez(i,j,nf) += sin(2.*pi*f*(PRECISION)n*dt) * Sc;
		else if (SourceChoice == 3)
			Ez(i,j,nf) += (1.-2.*pow(pi*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(pi*fp*((PRECISION)n*dt-dr),2)) * Sc;

		Dz(i,j,nf) = e0*Ez(i,j,nf);
	}
	// Recording transmitted fields.
	if (i==IEz/2 && j==x1)
		Ezt[n] = Ez(i,j,nf);
	if (i==IEz/2 && j==2*J/3+PMLw+1)
		Eztt[n] = Ez(i,j,nf);
	// Fields for refractive index.
	if (i==IEz/2 && j==Y1)
		Ezy1[n] = Ez(i,j,nf);
	if (i==IEz/2 && j==Y2)
		Ezy2[n] = Ez(i,j,nf);
}
