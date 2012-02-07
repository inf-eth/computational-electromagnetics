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
				Bx[i+IHx*j+IHx*JHx*n2] = (1-ScmHx[i+IHx*j])/(1+ScmHx[i+IHx*j]) * Bx[i+IHx*j+IHx*JHx*n1] + ( (dt/delta)/(1+ScmHx[i+IHx*j]) * (Ez[i+IEz*j+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hx[i+IHx*j+IHx*JHx*n2] = Bx[i+IHx*j+IHx*JHx*n2]/(u0*urHx[i+IHx*j]);

				By[(i+1)+IHy*(j+1)+IHy*JHy*n2] = (1-ScmHy[(i+1)+IHy*(j+1)])/(1+ScmHy[(i+1)+IHy*(j+1)]) * By[(i+1)+IHy*(j+1)+IHy*JHy*n1] + ( (dt/delta)/(1+ScmHy[(i+1)+IHy*(j+1)]) * (Ez[(i+1)+IEz*(j+1)+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hy[(i+1)+IHy*(j+1)+IHy*JHy*n2] = By[(i+1)+IHy*(j+1)+IHy*JHy*n2]/(u0*urHy[(i+1)+IHy*(j+1)]);
			}
			// Lower PML region.
			if (j < PMLw)
			{			
				Bx[i+IHx*j+IHx*JHx*n2] = (1-ScmsmyHx[i+IHx*j])/(1+ScmsmyHx[i+IHx*j]) * Bx[i+IHx*j+IHx*JHx*n1] + ( (dt/delta)/(1+ScmsmyHx[i+IHx*j]) * (Ez[i+IEz*j+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hx[i+IHx*j+IHx*JHx*n2] = Bx[i+IHx*j+IHx*JHx*n2]/(u0*urHx[i+IHx*j]);

				By[(i+1)+IHy*(j+1)+IHy*JHy*n2] = (1-ScmsmxHy[(i+1)+IHy*(j+1)])/(1+ScmsmxHy[(i+1)+IHy*(j+1)]) * By[(i+1)+IHy*(j+1)+IHy*JHy*n1] + ( (dt/delta)/(1+ScmsmxHy[(i+1)+IHy*(j+1)]) * (Ez[(i+1)+IEz*(j+1)+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hy[(i+1)+IHy*(j+1)+IHy*JHy*n2] = By[(i+1)+IHy*(j+1)+IHy*JHy*n2]/(u0*urHy[(i+1)+IHy*(j+1)]);			
			}
			// Upper PML region.
			if (j >= JHx-PMLw && j < JHx)
			{
				Bx[i+IHx*j+IHx*JHx*n2] = (1-ScmsmyHx[i+IHx*j])/(1+ScmsmyHx[i+IHx*j]) * Bx[i+IHx*j+IHx*JHx*n1] + ( (dt/delta)/(1+ScmsmyHx[i+IHx*j]) * (Ez[i+IEz*j+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hx[i+IHx*j+IHx*JHx*n2] = Bx[i+IHx*j+IHx*JHx*n2]/(u0*urHx[i+IHx*j]);

				By[(i+1)+IHy*(j+1)+IHy*JHy*n2] = (1-ScmsmxHy[(i+1)+IHy*(j+1)])/(1+ScmsmxHy[(i+1)+IHy*(j+1)]) * By[(i+1)+IHy*(j+1)+IHy*JHy*n1] + ( (dt/delta)/(1+ScmsmxHy[(i+1)+IHy*(j+1)]) * (Ez[(i+1)+IEz*(j+1)+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hy[(i+1)+IHy*(j+1)+IHy*JHy*n2] = By[(i+1)+IHy*(j+1)+IHy*JHy*n2]/(u0*urHy[(i+1)+IHy*(j+1)]);			
			}
		}
	}
	else

	{
		if (i < IEz)
		{
			if (j != 0 && j < JEz-1 )
			{
				Dz[i+IEz*j+IEz*JEz*n2] = (1-Sc[i+IEz*j])/(1+Sc[i+IEz*j]) * Dz[i+IEz*j+IEz*JEz*n1] + ( (dt/delta)/(1+Sc[i+IEz*j]) * ( Hy[(i+1)+IHy*j+IHy*JHy*n2] - Hy[i+IHy*j+IHy*JHy*n2] - Hx[i+IHx*j+IHx*JHx*n2] + Hx[i+IHx*(j-1)+IHx*JHx*n2]) );
				Ez[i+IEz*j+IEz*JEz*n2] = Dz[i+IEz*j+IEz*JEz*n2]/(e0*erEz[i+IEz*j]);
			}
			// Source.
			if (j == Js && n < NHW)
			{
				Ez[i+IEz*j+IEz*JEz*n2] = Ez[i+IEz*j+IEz*JEz*n2] + 1 * sin (Two_pi_f_deltat * n) / dtscalar;
				Dz[i+IEz*j+IEz*JEz*n2] = e0 * Ez[i+IEz*j+IEz*JEz*n2];
			}
			// Lower PML region.
			if (j > 0 && j < PMLw+1)
			{
				Dzx[i+IEz*j+IEz*JEz*n2] = (1-Scsx[i+IEz*j])/(1+Scsx[i+IEz*j]) * Dzx[i+IEz*j+IEz*JEz*n1] + ( (dt/delta)/(1+Scsx[i+IEz*j]) * ( Hy[(i+1)+IHy*j+IHy*JHy*n2] - Hy[i+IHy*j+IHy*JHy*n2]) );
				Dzy[i+IEz*j+IEz*JEz*n2] = (1-Scsy[i+IEz*j])/(1+Scsy[i+IEz*j]) * Dzy[i+IEz*j+IEz*JEz*n1] + ( (dt/delta)/(1+Scsy[i+IEz*j]) * (- Hx[i+IHx*j+IHx*JHx*n2] + Hx[i+IHx*(j-1)+IHx*JHx*n2]) );
				Dz[i+IEz*j+IEz*JEz*n2] = Dzx[i+IEz*j+IEz*JEz*n2] + Dzy[i+IEz*j+IEz*JEz*n2];
				Ez[i+IEz*j+IEz*JEz*n2] = Dz[i+IEz*j+IEz*JEz*n2]/(e0*erEz[i+IEz*j]);
			}
			// Upper PML region.
			if (j >= JEz-PMLw-1 && j < JEz-1)
			{
				Dzx[i+IEz*j+IEz*JEz*n2] = (1-Scsx[i+IEz*j])/(1+Scsx[i+IEz*j]) * Dzx[i+IEz*j+IEz*JEz*n1] + ( (dt/delta)/(1+Scsx[i+IEz*j]) * ( Hy[(i+1)+IHy*j+IHy*JHy*n2] - Hy[i+IHy*j+IHy*JHy*n2]) );
				Dzy[i+IEz*j+IEz*JEz*n2] = (1-Scsy[i+IEz*j])/(1+Scsy[i+IEz*j]) * Dzy[i+IEz*j+IEz*JEz*n1] + ( (dt/delta)/(1+Scsy[i+IEz*j]) * (- Hx[i+IHx*j+IHx*JHx*n2] + Hx[i+IHx*(j-1)+IHx*JHx*n2]) );
				Dz[i+IEz*j+IEz*JEz*n2] = Dzx[i+IEz*j+IEz*JEz*n2] + Dzy[i+IEz*j+IEz*JEz*n2];
				Ez[i+IEz*j+IEz*JEz*n2] = Dz[i+IEz*j+IEz*JEz*n2]/(e0*erEz[i+IEz*j]);
			}
					
		}
	}
}
