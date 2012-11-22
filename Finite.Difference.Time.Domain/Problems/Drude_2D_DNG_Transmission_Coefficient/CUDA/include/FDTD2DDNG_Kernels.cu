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

#include <FDTD2DDNG.hpp>

// Dry run kernel.
template <unsigned int BlockX, unsigned int BlockY> __global__ void FDTD2DDNGKernel_DryRun_M(
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
							const PRECISION *Ez_, const PRECISION *Dz_, PRECISION *Hx_, PRECISION *Bx_, PRECISION *Hy_, PRECISION *By_,
							const PRECISION *PsiEzX_, const PRECISION *PsiEzY_, PRECISION *PsiHyX_, PRECISION *PsiHxY_,
							const PRECISION *einf_, const PRECISION *uinf_,
							const PRECISION kappex, const PRECISION kappey, const PRECISION kappmx, const PRECISION kappmy,
							const PRECISION bex, const PRECISION bey, const PRECISION bmx, const PRECISION bmy,
							const PRECISION Cex, const PRECISION Cey, const PRECISION Cmx, const PRECISION Cmy,
							PRECISION *Ezi,
							const unsigned int x1,
							const unsigned int n,
							const unsigned int np,
							const unsigned int n0,
							const unsigned int nf)
{
	const unsigned int i = BlockX*blockIdx.x+threadIdx.x;
	const unsigned int j = BlockY*blockIdx.y+threadIdx.y;

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
	// PsiHyX arrays.
	if (i<IHy-1)
	{

		if (j<JHy)
		{
			PsiHyX(i,j) = (Cmx/delta)*(Ez(i+1,j,n0)-Ez(i,j,n0)) + bmx*PsiHyX(i,j);
			if (i==0)
				PsiHyX(IHy-1,j) = (Cmx/delta)*(Ez(0,j,n0)-Ez(IHy-1,j,n0)) + bmx*PsiHyX(IHy-1,j); // PBC
		}
		// By in normal space.
		if (j>PMLw-1 && j<JHy-PMLw)
		{
			By(i,j,nf) = By(i,j,n0) + (Ez(i+1,j,n0) - Ez(i,j,n0)) * dt/delta;
			Hy(i,j,nf) = By(i,j,nf)/(u0*uinf(i,j));
			if (i==0)
			{
				By(IHy-1,j,nf) = By(IHy-1,j,n0) + (Ez(0,j,n0) - Ez(IHy-1,j,n0)) * dt/delta; // PBC
				Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/(u0*uinf(IHy-1,j)); // PBC
			}
		}
		// By in Lower PML.
		if (j<PMLw)
		{
			By(i,j,nf) = By(i,j,n0) + dt*((1./kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1./delta + PsiHyX(i,j));
			Hy(i,j,nf) = By(i,j,nf)/(u0*uinf(i,j));
			if (i==0)
			{
				By(IHy-1,j,nf) = By(IHy-1,j,n0) + dt*((1./kappmx)*(Ez(0,j,n0) - Ez(IHy-1,j,n0)) * 1./delta + PsiHyX(IHy-1,j)); // PBC
				Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/(u0*uinf(IHy-1,j)); // PBC
			}
		}
		// By in upper PML.
		if (j>JHy-PMLw-1 && j<JHy)
		{
			By(i,j,nf) = By(i,j,n0) + dt*((1./kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1./delta + PsiHyX(i,j));
			Hy(i,j,nf) = By(i,j,nf)/(u0*uinf(i,j));
			if (i==0)
			{
				By(IHy-1,j,nf) = By(IHy-1,j,n0) + dt*((1./kappmx)*(Ez(0,j,n0) - Ez(IHy-1,j,n0)) * 1./delta + PsiHyX(IHy-1,j)); // PBC
				Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/(u0*uinf(IHy-1,j)); // PBC
			}
		}
	}
}
template <unsigned int BlockX, unsigned int BlockY> __global__ void FDTD2DDNGKernel_DryRun_E(
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
							PRECISION *Ez_, PRECISION *Dz_, const PRECISION *Hx_, const PRECISION *Bx_, const PRECISION *Hy_, const PRECISION *By_,
							PRECISION *PsiEzX_, PRECISION *PsiEzY_, const PRECISION *PsiHyX_, const PRECISION *PsiHxY_,
							const PRECISION *einf_, const PRECISION *uinf_,
							const PRECISION kappex, const PRECISION kappey, const PRECISION kappmx, const PRECISION kappmy,
							const PRECISION bex, const PRECISION bey, const PRECISION bmx, const PRECISION bmy,
							const PRECISION Cex, const PRECISION Cey, const PRECISION Cmx, const PRECISION Cmy,
							PRECISION *Ezi,
							const unsigned int x1,
							const unsigned int n,
							const unsigned int np,
							const unsigned int n0,
							const unsigned int nf)
{
	const unsigned int i = BlockX*blockIdx.x+threadIdx.x;
	const unsigned int j = BlockY*blockIdx.y+threadIdx.y;

	// ========================== Dz and Ez ==========================
	if (i>0 && i<IEz)
	{
		// Psi arrays.
		if (j<JEz)
		{
			PsiEzX(i,j) = (Cex/delta)*(Hy(i,j,nf)-Hy(i-1,j,nf)) + bex*PsiEzX(i,j);
			if (i==1)
				PsiEzX(0,j) = (Cex/delta)*(Hy(0,j,nf)-Hy(IEz-1,j,nf)) + bex*PsiEzX(0,j); // PBC
			PsiEzY(i,j) = (Cey/delta)*(-Hx(i,j+1,nf)+Hx(i,j,nf)) + bey*PsiEzY(i,j);
			if (i==1)
				PsiEzY(0,j) = (Cey/delta)*(-Hx(0,j+1,nf)+Hx(0,j,nf)) + bey*PsiEzY(0,j); // PBC
		}
		// Dz in normal space.
		if (j>PMLw-1 && j<JEz-PMLw)
		{
			Dz(i,j,nf) = Dz(i,j,n0) + (Hy(i,j,nf)-Hy(i-1,j,nf)-Hx(i,j+1,nf)+Hx(i,j,nf)) * dt/delta;
			Ez(i,j,nf) = Dz(i,j,nf)/(e0*einf(i,j));
			if (i==1)
			{
				Dz(0,j,nf) = Dz(0,j,n0) + (Hy(0,j,nf)-Hy(IEz-1,j,nf)-Hx(0,j+1,nf)+Hx(0,j,nf)) * dt/delta; // PBC
				Ez(0,j,nf) = Dz(0,j,nf)/(e0*einf(0,j)); // PBC
			}
		}
		// Dz in lower PML.
		if (j<PMLw)
		{
			Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1./kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1./kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1./delta + PsiEzX(i,j) + PsiEzY(i,j));
			Ez(i,j,nf) = Dz(i,j,nf)/(e0*einf(i,j));
			if (i==1)
			{
				Dz(0,j,nf) = Dz(0,j,n0) + dt*(((1./kappex)*(Hy(0,j,nf)-Hy(IEz-1,j,nf))+(1./kappey)*(-Hx(0,j+1,nf)+Hx(0,j,nf))) * 1./delta + PsiEzX(0,j) + PsiEzY(0,j)); // PBC
				Ez(0,j,nf) = Dz(0,j,nf)/(e0*einf(0,j)); // PBC
			}
		}
		// Dz in upper PML.
		if (j>JEz-PMLw-1 && j<JEz)
		{
			Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1./kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1./kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1./delta + PsiEzX(i,j) + PsiEzY(i,j));
			Ez(i,j,nf) = Dz(i,j,nf)/(e0*einf(i,j));
			if (i==1)
			{
				Dz(0,j,nf) = Dz(0,j,n0) + dt*(((1./kappex)*(Hy(0,j,nf)-Hy(IEz-1,j,nf))+(1./kappey)*(-Hx(0,j+1,nf)+Hx(0,j,nf))) * 1./delta + PsiEzX(0,j) + PsiEzY(0,j)); // PBC
				Ez(0,j,nf) = Dz(0,j,nf)/(e0*einf(0,j)); // PBC
			}
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
// Simulation kernel.
template <unsigned int BlockX, unsigned int BlockY> __global__ void FDTD2DDNGKernel_Simulation_M(
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
							const PRECISION *Ez_, const PRECISION *Dz_, PRECISION *Hx_, PRECISION *Bx_, PRECISION *Hy_, PRECISION *By_,
							const PRECISION *einf_, const PRECISION *uinf_, const PRECISION *wpesq_, const PRECISION *wpmsq_, const PRECISION *ge_, const PRECISION *gm_,
							const PRECISION *ae0_, const PRECISION *ae_, const PRECISION *be_, const PRECISION *ce_, const PRECISION *de_, const PRECISION *ee_,
							const PRECISION *am0_, const PRECISION *am_, const PRECISION *bm_, const PRECISION *cm_, const PRECISION *dm_, const PRECISION *em_,
							const PRECISION *PsiEzX_, const PRECISION *PsiEzY_, PRECISION *PsiHyX_, PRECISION *PsiHxY_,
							const PRECISION kappex, const PRECISION kappey, const PRECISION kappmx, const PRECISION kappmy,
							const PRECISION bex, const PRECISION bey, const PRECISION bmx, const PRECISION bmy,
							const PRECISION Cex, const PRECISION Cey, const PRECISION Cmx, const PRECISION Cmy,
							PRECISION *Ezt, PRECISION *Eztt, PRECISION *Ezy1, PRECISION *Ezy2,
							const unsigned int x1, const unsigned int Y1, const unsigned int Y2,
							const unsigned int n,
							const unsigned int np,
							const unsigned int n0,
							const unsigned int nf)
{
	const unsigned int i = BlockX*blockIdx.x+threadIdx.x;
	const unsigned int j = BlockY*blockIdx.y+threadIdx.y;

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
		{
			PsiHyX(i,j) = (Cmx/delta)*(Ez(i+1,j,n0)-Ez(i,j,n0)) + bmx*PsiHyX(i,j);
			if (i==0)
				PsiHyX(IHy-1,j) = (Cmx/delta)*(Ez(0,j,n0)-Ez(IHy-1,j,n0)) + bmx*PsiHyX(IHy-1,j); // PBC
		}
		// By in normal space.
		if (j>PMLw-1 && j<JHy-PMLw)
		{
			By(i,j,nf) = By(i,j,n0) + (Ez(i+1,j,n0) - Ez(i,j,n0)) * dt/delta;
			Hy(i,j,nf) = am(i,j)*(By(i,j,nf)-2.*By(i,j,n0)+By(i,j,np))+bm(i,j)*(By(i,j,nf)-By(i,j,np))+cm(i,j)*(2.*Hy(i,j,n0)-Hy(i,j,np))+dm(i,j)*(2.*Hy(i,j,n0)+Hy(i,j,np))+em(i,j)*Hy(i,j,np);
			if (i==0)
			{
				By(IHy-1,j,nf) = By(IHy-1,j,n0) + (Ez(0,j,n0) - Ez(IHy-1,j,n0)) * dt/delta; // PBC
				Hy(IHy-1,j,nf) = am(IHy-1,j)*(By(IHy-1,j,nf)-2.*By(IHy-1,j,n0)+By(IHy-1,j,np))+bm(IHy-1,j)*(By(IHy-1,j,nf)-By(IHy-1,j,np))+cm(IHy-1,j)*(2.*Hy(IHy-1,j,n0)-Hy(IHy-1,j,np))+dm(IHy-1,j)*(2.*Hy(IHy-1,j,n0)+Hy(IHy-1,j,np))+em(IHy-1,j)*Hy(IHy-1,j,np); // PBC
			}
		}
		// By in Lower PML.
		if (j<PMLw)
		{
			By(i,j,nf) = By(i,j,n0) + dt*((1./kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1./delta + PsiHyX(i,j));
			Hy(i,j,nf) = By(i,j,nf)/(u0*uinf(i,j));
			if (i==0)
			{
				By(IHy-1,j,nf) = By(IHy-1,j,n0) + dt*((1./kappmx)*(Ez(0,j,n0) - Ez(IHy-1,j,n0)) * 1./delta + PsiHyX(IHy-1,j)); // PBC
				Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/(u0*uinf(IHy-1,j)); // PBC
			}
		}
		// By in upper PML.
		if (j>JHy-PMLw-1 && j<JHy)
		{
			By(i,j,nf) = By(i,j,n0) + dt*((1./kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1./delta + PsiHyX(i,j));
			Hy(i,j,nf) = By(i,j,nf)/(u0*uinf(i,j));
			if (i==0)
			{
				By(IHy-1,j,nf) = By(IHy-1,j,n0) + dt*((1./kappmx)*(Ez(0,j,n0) - Ez(IHy-1,j,n0)) * 1./delta + PsiHyX(IHy-1,j)); // PBC
				Hy(IHy-1,j,nf) = By(IHy-1,j,nf)/(u0*uinf(IHy-1,j)); // PBC
			}
		}
	}
}
template <unsigned int BlockX, unsigned int BlockY> __global__ void FDTD2DDNGKernel_Simulation_E(
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
							PRECISION *Ez_, PRECISION *Dz_, const PRECISION *Hx_, const PRECISION *Bx_, const PRECISION *Hy_, const PRECISION *By_,
							const PRECISION *einf_, const PRECISION *uinf_, const PRECISION *wpesq_, const PRECISION *wpmsq_, const PRECISION *ge_, const PRECISION *gm_,
							const PRECISION *ae0_, const PRECISION *ae_, const PRECISION *be_, const PRECISION *ce_, const PRECISION *de_, const PRECISION *ee_,
							const PRECISION *am0_, const PRECISION *am_, const PRECISION *bm_, const PRECISION *cm_, const PRECISION *dm_, const PRECISION *em_,
							PRECISION *PsiEzX_, PRECISION *PsiEzY_, const PRECISION *PsiHyX_, const PRECISION *PsiHxY_,
							const PRECISION kappex, const PRECISION kappey, const PRECISION kappmx, const PRECISION kappmy,
							const PRECISION bex, const PRECISION bey, const PRECISION bmx, const PRECISION bmy,
							const PRECISION Cex, const PRECISION Cey, const PRECISION Cmx, const PRECISION Cmy,
							PRECISION *Ezt, PRECISION *Eztt, PRECISION *Ezy1, PRECISION *Ezy2,
							const unsigned int x1, const unsigned int Y1, const unsigned int Y2,
							const unsigned int n,
							const unsigned int np,
							const unsigned int n0,
							const unsigned int nf)
{
	const unsigned int i = BlockX*blockIdx.x+threadIdx.x;
	const unsigned int j = BlockY*blockIdx.y+threadIdx.y;

	// ========================== Dz and Ez ==========================
	if (i>0 && i<IEz)
	{
		// Psi arrays.
		if (j<JEz)
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
		if (j>PMLw-1 && j<JEz-PMLw)
		{
			Dz(i,j,nf) = Dz(i,j,n0) + (Hy(i,j,nf)-Hy(i-1,j,nf)-Hx(i,j+1,nf)+Hx(i,j,nf)) * dt/delta;
			Ez(i,j,nf) = ae(i,j)*(Dz(i,j,nf)-2.*Dz(i,j,n0)+Dz(i,j,np))+be(i,j)*(Dz(i,j,nf)-Dz(i,j,np))+ce(i,j)*(2.*Ez(i,j,n0)-Ez(i,j,np))+de(i,j)*(2.*Ez(i,j,n0)+Ez(i,j,np))+ee(i,j)*Ez(i,j,np);
			if (i==1)
			{
				Dz(0,j,nf) = Dz(0,j,n0) + (Hy(0,j,nf)-Hy(IEz-1,j,nf)-Hx(0,j+1,nf)+Hx(0,j,nf)) * dt/delta; // PBC
				Ez(0,j,nf) = ae(0,j)*(Dz(0,j,nf)-2.*Dz(0,j,n0)+Dz(0,j,np))+be(0,j)*(Dz(0,j,nf)-Dz(0,j,np))+ce(0,j)*(2.*Ez(0,j,n0)-Ez(0,j,np))+de(0,j)*(2.*Ez(0,j,n0)+Ez(0,j,np))+ee(0,j)*Ez(0,j,np); // PBC
			}
		}
		// Dz in lower PML.
		if (j<PMLw)
		{
			Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1./kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1./kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1./delta + PsiEzX(i,j) + PsiEzY(i,j));
			Ez(i,j,nf) = Dz(i,j,nf)/(e0*einf(i,j));
			if (i==1)
			{
				Dz(0,j,nf) = Dz(0,j,n0) + dt*(((1./kappex)*(Hy(0,j,nf)-Hy(IEz-1,j,nf))+(1./kappey)*(-Hx(0,j+1,nf)+Hx(0,j,nf))) * 1./delta + PsiEzX(0,j) + PsiEzY(0,j)); // PBC
				Ez(0,j,nf) = Dz(0,j,nf)/(e0*einf(0,j)); // PBC
			}
		}
		// Dz in upper PML.
		if (j>JEz-PMLw-1 && j<JEz)
		{
			Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1./kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1./kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1./delta + PsiEzX(i,j) + PsiEzY(i,j));
			Ez(i,j,nf) = Dz(i,j,nf)/(e0*einf(i,j));
			if (i==1)
			{
				Dz(0,j,nf) = Dz(0,j,n0) + dt*(((1./kappex)*(Hy(0,j,nf)-Hy(IEz-1,j,nf))+(1./kappey)*(-Hx(0,j+1,nf)+Hx(0,j,nf))) * 1./delta + PsiEzX(0,j) + PsiEzY(0,j)); // PBC
				Ez(0,j,nf) = Dz(0,j,nf)/(e0*einf(0,j)); // PBC
			}
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
