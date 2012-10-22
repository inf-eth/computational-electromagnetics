#define uint unsigned int
#define HX(i,j,n) Hx[i+IHx*(j)+IHx*JHx*n]
#define BX(i,j,n) Bx[i+IHx*(j)+IHx*JHx*n]
#define HY(i,j,n) Hy[i+IHy*(j)+IHy*JHy*n]
#define BY(i,j,n) By[i+IHy*(j)+IHy*JHy*n]
#define EZ(i,j,n) Ez[i+IEz*(j)+IEz*JEz*n]
#define DZ(i,j,n) Dz[i+IEz*(j)+IEz*JEz*n]
#define DZX(i,j,n) Dzx[i+IEz*(j)+IEz*JEz*n]
#define DZY(i,j,n) Dzy[i+IEz*(j)+IEz*JEz*n]

template <unsigned int BlockX, unsigned int BlockY> __global__ void FDTD2DKernel(
							float *Hx,
							float *Bx,
							float *Hy,
							float *By,
							float *Ez,
							float *Dz,
							float *Dzx,
							float *Dzy,
							float *urHx,
							float *urHy,
							float *erEz,
							float *ScmHx,
							float *ScmHy,
							float *Sc,
							float *Scsx,
							float *Scsy,
							float *ScmsmxHy,
							float *ScmsmyHx,
							const float delta,
							const float dtscalar,
							const float dt,
							const uint PMLw,
							const float e0,
							const float u0,
							const float Two_pi_f_deltat,
							const uint NHW,
							const uint Is,
							const uint Js,
							const uint IHx,
							const uint JHx,
							const uint IHy,
							const uint JHy,
							const uint IEz,
							const uint JEz,
							const uint n,
							const uint n0,
							const uint n1,
							const uint n2,
							const uint flag)
{
	uint i = BlockX*blockIdx.x+threadIdx.x;
	uint j = BlockY*blockIdx.y+threadIdx.y;

	// Half time step flag is either 0 or 1 indicating whether magnetic field or electric field is to be calculated, respectively.
	if (flag == 0)
	{
		if (i < IHx)
		{
			// Normal space.
			if (j >= PMLw && j < JHx-PMLw)
			{
				BX(i,j,n2) = (1-ScmHx[i+IHx*j])/(1+ScmHx[i+IHx*j]) * BX(i,j,n1) + ( (dt/delta)/(1+ScmHx[i+IHx*j]) * (EZ(i,j,n1)-EZ(i,j+1,n1)) );
				HX(i,j,n2) = BX(i,j,n2)/(u0*urHx[i+IHx*j]);

				BY(i+1,j+1,n2) = (1-ScmHy[(i+1)+IHy*(j+1)])/(1+ScmHy[(i+1)+IHy*(j+1)]) * BY(i+1,j+1,n1) + ( (dt/delta)/(1+ScmHy[(i+1)+IHy*(j+1)]) * (EZ(i+1,j+1,n1)-EZ(i,j+1,n1)) );
				HY(i+1,j+1,n2) = BY(i+1,j+1,n2)/(u0*urHy[(i+1)+IHy*(j+1)]);
			}
			// Lower PML region.
			if (j < PMLw)
			{
				BX(i,j,n2) = (1-ScmsmyHx[i+IHx*j])/(1+ScmsmyHx[i+IHx*j]) * BX(i,j,n1) + ( (dt/delta)/(1+ScmsmyHx[i+IHx*j]) * (EZ(i,j,n1)-EZ(i,j+1,n1)) );
				HX(i,j,n2) = BX(i,j,n2)/(u0*urHx[i+IHx*j]);

				BY(i+1,j+1,n2) = (1-ScmsmxHy[(i+1)+IHy*(j+1)])/(1+ScmsmxHy[(i+1)+IHy*(j+1)]) * BY(i+1,j+1,n1) + ( (dt/delta)/(1+ScmsmxHy[(i+1)+IHy*(j+1)]) * (EZ(i+1,j+1,n1)-EZ(i,j+1,n1)) );
				HY(i+1,j+1,n2) = BY(i+1,j+1,n2)/(u0*urHy[(i+1)+IHy*(j+1)]);
			}
			// Upper PML region.
			if (j >= JHx-PMLw && j < JHx)
			{
				BX(i,j,n2) = (1-ScmsmyHx[i+IHx*j])/(1+ScmsmyHx[i+IHx*j]) * BX(i,j,n1) + ( (dt/delta)/(1+ScmsmyHx[i+IHx*j]) * (EZ(i,j,n1)-EZ(i,j+1,n1)) );
				HX(i,j,n2) = BX(i,j,n2)/(u0*urHx[i+IHx*j]);

				BY(i+1,j+1,n2) = (1-ScmsmxHy[(i+1)+IHy*(j+1)])/(1+ScmsmxHy[(i+1)+IHy*(j+1)]) * BY(i+1,j+1,n1) + ( (dt/delta)/(1+ScmsmxHy[(i+1)+IHy*(j+1)]) * (EZ(i+1,j+1,n1)-EZ(i,j+1,n1)) );
				HY(i+1,j+1,n2) = BY(i+1,j+1,n2)/(u0*urHy[(i+1)+IHy*(j+1)]);
			}
		}
	}
	else

	{
		if (i < IEz)
		{
			if (j != 0 && j < JEz-1 )
			{
				DZ(i,j,n2) = (1-Sc[i+IEz*j])/(1+Sc[i+IEz*j]) * DZ(i,j,n1) + ( (dt/delta)/(1+Sc[i+IEz*j]) * ( HY(i+1,j,n2) - HY(i,j,n2) - HX(i,j,n2) + HX(i,j-1,n2)) );
				EZ(i,j,n2) = DZ(i,j,n2)/(e0*erEz[i+IEz*j]);
			}
			// Source.
			if (j == Js && n < NHW)
			{
				EZ(i,j,n2) = EZ(i,j,n2) + 1 * sin (Two_pi_f_deltat * n) / dtscalar;
				DZ(i,j,n2) = e0 * EZ(i,j,n2);
			}
			// Lower PML region.
			if (j > 0 && j < PMLw+1)
			{
				DZX(i,j,n2) = (1-Scsx[i+IEz*j])/(1+Scsx[i+IEz*j]) * DZX(i,j,n1) + ( (dt/delta)/(1+Scsx[i+IEz*j]) * ( HY(i+1,j,n2) - HY(i,j,n2)) );
				DZY(i,j,n2) = (1-Scsy[i+IEz*j])/(1+Scsy[i+IEz*j]) * DZY(i,j,n1) + ( (dt/delta)/(1+Scsy[i+IEz*j]) * (- HX(i,j,n2) + HX(i,j-1,n2)) );
				DZ(i,j,n2) = DZX(i,j,n2) + DZY(i,j,n2);
				EZ(i,j,n2) = DZ(i,j,n2)/(e0*erEz[i+IEz*j]);
			}
			// Upper PML region.
			if (j >= JEz-PMLw-1 && j < JEz-1)
			{
				DZX(i,j,n2) = (1-Scsx[i+IEz*j])/(1+Scsx[i+IEz*j]) * DZX(i,j,n1) + ( (dt/delta)/(1+Scsx[i+IEz*j]) * ( HY(i+1,j,n2) - HY(i,j,n2)) );
				DZY(i,j,n2) = (1-Scsy[i+IEz*j])/(1+Scsy[i+IEz*j]) * DZY(i,j,n1) + ( (dt/delta)/(1+Scsy[i+IEz*j]) * (- HX(i,j,n2) + HX(i,j-1,n2)) );
				DZ(i,j,n2) = DZX(i,j,n2) + DZY(i,j,n2);
				EZ(i,j,n2) = DZ(i,j,n2)/(e0*erEz[i+IEz*j]);
			}
		}
	}
}
