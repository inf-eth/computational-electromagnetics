#include <FDTD2D.hpp>

CFDTD2D::CFDTD2D () : 
						I(256),
						J(256),
						c(299792458.),
						delta(2.5e-3),
						dx(delta),
						dy(delta),
						dtscalar(2.),
						dt(delta/(sqrt(2.)*c) /dtscalar),
						PMLw(32),
						NMax(1000),
						f(2.e9),
						pi(4*atan(1.)),
						e0(1.e-9/(36.*pi)),
						u0(1.e-7*4.*pi),
						Two_pi_f_deltat(2*pi*f*dt),
						NHW((unsigned int)1./(2.*f*dt)),
						Js(20+PMLw), 
						Is(2),
						n(1),
						n0(0),
						n1(1),
						n2(2),
						tResolution(2),
						xResolution(2),
						yResolution(2),
						IHx(I),
						JHx(J+2*PMLw-1),
						IHy(I+1),
						JHy(J+2*PMLw),
						IEz(I),
						JEz(J+2*PMLw),
						XCenter((I+1)/2),
						YCenter((J-1)/2),
						ra(0.1),
						rb(0.2),
						imp0(377.0),

						// Data arrays
						Hx(NULL),
						Bx(NULL),
						Hy(NULL),
						By(NULL),
						Ez(NULL),
						Dz(NULL),
						Dzx(NULL),
						Dzy(NULL),
						EzSnapshots(NULL),
						urHx(NULL),
						urHy(NULL),
						erEz(NULL),
						smHx(NULL),
						smHy(NULL),
						sEz(NULL),
						Sc(NULL),
						ScmHx(NULL),
						ScmHy(NULL),
						sex(NULL),
						sey(NULL),
						smx(NULL),
						smy(NULL),
						Scsx(NULL),
						Scsy(NULL),
						ScmsmxHy(NULL),
						ScmsmyHx(NULL),
						cpu(false),
						flagHalf(0)
{
	uint tbytes = 0.5*I*(J+2*PMLw)*17*3*8+50*8;		// Dynamic arrays + predefined variables.
	std::cout << "Approximate memory required for simulation: " << (float)tbytes/1024 << " Kbytes (" << (float)tbytes/(1024*1024) << " MB)." << std::endl;

	// Writing simulation parameters for plotting.
	#ifdef WIN32
	std::fstream parametersfile;
	parametersfile.open ("./FieldData/Parameters.smp", std::ios::out|std::ios::binary);
	parametersfile.write ((char*)&I, sizeof(unsigned int));
	parametersfile.write ((char*)&J, sizeof(unsigned int));
	parametersfile.write ((char*)&tResolution, sizeof(unsigned int));
	parametersfile.write ((char*)&xResolution, sizeof(unsigned int));
	parametersfile.write ((char*)&yResolution, sizeof(unsigned int));
	parametersfile.write ((char*)&NMax, sizeof(unsigned int));
	parametersfile.write ((char*)&PMLw, sizeof(uint));
	parametersfile.close ();
	#else
	int fdp;
	fdp = open ( "./FieldData/Parameters.smp", O_CREAT|O_WRONLY|O_TRUNC, S_IRWXU );
	write (fdp, (void*)&I, sizeof(unsigned int));
	write (fdp, (void*)&J, sizeof(unsigned int));
	write (fdp, (void*)&tResolution, sizeof(unsigned int));
	write (fdp, (void*)&xResolution, sizeof(unsigned int));
	write (fdp, (void*)&yResolution, sizeof(unsigned int));
	write (fdp, (void*)&NMax, sizeof(unsigned int));
	write (fdp, (void*)&PMLw, sizeof(uint));
	close (fdp);
	#endif
}
// Initialize data arrays.
int CFDTD2D::Initialize ()
{
	Hx = new float[IHx*JHx*3];
	Bx = new float[IHx*JHx*3];
	Hy = new float[IHy*JHy*3];
	By = new float[IHy*JHy*3];
	Ez = new float[IEz*JEz*3];
	Dz = new float[IEz*JEz*3];
	Dzx = new float[IEz*JEz*3];
	Dzy = new float[IEz*JEz*3];
	//EzSnapshots = new float[(IEz/xResolution)*(JEz/yResolution)];

	urHx = new float[IHx*JHx];
	urHy = new float[IHy*JHy];
	erEz = new float[IEz*JEz];
	smHx = new float[IHx*JHx];
	smHy = new float[IHy*JHy];
	sEz = new float[IEz*JEz];
	Sc = new float[IEz*JEz];
	ScmHx = new float[IHx*JHx];
	ScmHy = new float[IHy*JHy];

	sex = new float[IEz*JEz];
	sey = new float[IEz*JEz];
	smx = new float[IHy*JHy];
	smy = new float[IHx*JHx];
	Scsx = new float[IEz*JEz];
	Scsy = new float[IEz*JEz];
	ScmsmxHy = new float[IHy*JHy];
	ScmsmyHx = new float[IHx*JHx];

	unsigned int i, j, t;
	for (t=0; t<3; t++)
	{
		for (i=0; i<IHy; i++)
		{
			for(j=0; j<JHy; j++)
			{
				// Field initialization.

				Hy[i+IHy*j+t*IHy*JHy] = 0.;
				By[i+IHy*j+t*IHy*JHy] = 0.;

				if (i<IHx && j<JHx)
				{
					Hx[i+IHx*j+t*IHx*JHx] = 0.;
					Bx[i+IHx*j+t*IHx*JHx] = 0.;
				}
				if (i<IEz && J<JEz)
				{
					Ez[i+IEz*j+t*IEz*JEz] = 0.;
					Dz[i+IEz*j+t*IEz*JEz] = 0.;
					Dzx[i+IEz*j+t*IEz*JEz] = 0.;
					Dzy[i+IEz*j+t*IEz*JEz] = 0.;
				}

				if (t==0)
				{
					// Hy related arrays.
					urHy[i+IHy*j] = 1.;
					smHy[i+IHy*j] = 0.;
					ScmHy[i+IHy*j] = (dt*smHy[i+IHy*j])/(2*urHy[i+IHy*j]);
					smx[i+IHy*j] = 0.;
					ScmsmxHy[i+IHy*j] = (dt*smx[i+IHy*j])/(2*urHy[i+IHy*j]);

					// Hx related arrays.
					if (i<IHx && j<JHx)
					{
						urHx[i+IHx*j] = 1.;
						smHx[i+IHx*j] = 0.;
						ScmHx[i+IHx*j] = (dt*smHx[i+IHx*j])/(2*urHx[i+IHx*j]);

						// Initializing PML conductances.
						if (j < PMLw+1 || j > JHx-PMLw-1)
						{
							smy[i+IHx*j] = 1.6e10;
						}
						else
							smy[i+IHx*j] = 0.;

						ScmsmyHx[i+IHx*j] = (dt*smy[i+IHx*j])/(2*urHx[i+IHx*j]);
					}

					// Ez related arrays.
					if (i<IEz && j<JEz)
					{
						erEz[i+IEz*j] = 1.;
						sEz[i+IEz*j] = 0.;
						Sc[i+IEz*j] = (dt*sEz[i+IEz*j])/(2*erEz[i+IEz*j]);
						sex[i+IEz*j] = 0.;

						// Initializing PML conductances.
						if (j < PMLw+1 || j > JEz-PMLw-1)
						{
							sey[i+IEz*j] = 1.6e10;
						}
						else
							sey[i+IEz*j] = 0.;

						Scsx[i+IEz*j] = (dt*sex[i+IEz*j])/(2*erEz[i+IEz*j]);
						Scsy[i+IEz*j] = (dt*sey[i+IEz*j])/(2*erEz[i+IEz*j]);
					}
				}
			}
		}
	}
	return 0;
}

int CFDTD2D::initializeFDTD2DKernel ()
{
	// Device memory allocation
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Hx, sizeof(float) * IHx*JHx*3) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Bx, sizeof(float) * IHx*JHx*3) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Hy, sizeof(float) * IHy*JHy*3) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_By, sizeof(float) * IHy*JHy*3) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Ez, sizeof(float) * IEz*JEz*3) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Dz, sizeof(float) * IEz*JEz*3) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Dzx, sizeof(float) * IEz*JEz*3) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Dzy, sizeof(float) * IEz*JEz*3) );
	
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_urHx, sizeof(float) * IHx*JHx) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_urHy, sizeof(float) * IHy*JHy) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_erEz, sizeof(float) * IEz*JEz) );
	
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_smHx, sizeof(float) * IHx*JHx) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_smHy, sizeof(float) * IHy*JHy) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_sEz, sizeof(float) * IEz*JEz) );

	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Sc, sizeof(float) * IEz*JEz) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_ScmHx, sizeof(float) * IHx*JHx) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_ScmHy, sizeof(float) * IHy*JHy) );
	
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Scsx, sizeof(float) * IEz*JEz) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Scsy, sizeof(float) * IEz*JEz) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_ScmsmxHy, sizeof(float) * IHx*JHx) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_ScmsmyHx, sizeof(float) * IHy*JHy) );
		
	// Device data initialization.
	CUDA_SAFE_CALL( cudaMemcpy(d_Hx, Hx, sizeof(float) * IHx*JHx*3, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_Bx, Bx, sizeof(float) * IHx*JHx*3, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_Hy, Hy, sizeof(float) * IHy*JHy*3, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_By, By, sizeof(float) * IHy*JHy*3, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_Ez, Ez, sizeof(float) * IEz*JEz*3, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_Dz, Dz, sizeof(float) * IEz*JEz*3, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_Dzx, Dzx, sizeof(float) * IEz*JEz*3, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_Dzy, Dzy, sizeof(float) * IEz*JEz*3, cudaMemcpyHostToDevice) );
	
	CUDA_SAFE_CALL( cudaMemcpy(d_urHx, urHx, sizeof(float) * IHx*JHx, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_urHy, urHy, sizeof(float) * IHy*JHy, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_erEz, erEz, sizeof(float) * IEz*JEz, cudaMemcpyHostToDevice) );
	
	CUDA_SAFE_CALL( cudaMemcpy(d_smHx, smHx, sizeof(float) * IHx*JHx, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_smHy, smHy, sizeof(float) * IHy*JHy, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_sEz, sEz, sizeof(float) * IEz*JEz, cudaMemcpyHostToDevice) );

	CUDA_SAFE_CALL( cudaMemcpy(d_ScmHx, ScmHx, sizeof(float) * IHx*JHx, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_ScmHy, ScmHy, sizeof(float) * IHy*JHy, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_Sc, Sc, sizeof(float) * IEz*JEz, cudaMemcpyHostToDevice) );

	CUDA_SAFE_CALL( cudaMemcpy(d_Scsx, Scsx, sizeof(float) * IEz*JEz, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_Scsy, Scsy, sizeof(float) * IEz*JEz, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_ScmsmxHy, ScmsmxHy, sizeof(float) * IHx*JHx, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_ScmsmyHx, ScmsmyHx, sizeof(float) * IHy*JHy, cudaMemcpyHostToDevice) );

	return 0;
}

int CFDTD2D::runFDTD2DKernels (bool SaveFields)
{
/*	cl_int status;
	unsigned int maxDims;
	cl_event events[2];
	size_t globalThreads[2];
	size_t localThreads[2];
	size_t maxWorkGroupSize;
	size_t maxWorkItemSizes[3];
*/
	// Total local threads in a block. Can be thought of as Block dimensions.
	unsigned int ThreadsX = 16;
	unsigned int ThreadsY = 16;
	
	// Total blocks in simulation grid. Can be thought of as no. of blocks in grid.
	// Obviously, I and J should be divisible by block dimensions.
	unsigned int BlocksX = I/ThreadsX;
	unsigned int BlocksY = (J+2*PMLw)/ThreadsY;

	// Kernel parameters.
	dim3 Blocks(BlocksX, BlocksY);
	dim3 Threads(ThreadsX, ThreadsY);

	std::cout << "Block dimensions: " << ThreadsX << "x" << ThreadsY << std::endl;
	std::cout << "Grid dimensions: " << BlocksX << "x" << BlocksY << std::endl;
/*	if(localThreads[0] > maxWorkGroupSize ||
		localThreads[0] > maxWorkItemSizes[0])
	{
		std::cout<<"Unsupported: Device does not support requested number of work items.";
		return 1;
	}
*/
	// ====== Set appropriate arguments to the kernel ======
	
	// Input Field arrays.

	// Permeability and permittivity.

	// ==========================================================
	
//	cl_ulong startTime, endTime;
//	cl_ulong kernelExecTimeNs;
//	cl_ulong kernelExecTimeNsT = 0;

	unsigned int hTimer;
	CUT_SAFE_CALL (cutCreateTimer(&hTimer));
	CUT_SAFE_CALL (cutResetTimer(hTimer));
		
	// File handling from chapter 3 of Understanding FDTD. J. B. Schneider
	std::stringstream framestream;
	std::string basename = "./FieldData/Ez";
	std::string filename;
	unsigned int frame = 1;

	#ifdef WIN32
	std::fstream snapshot;
	#else
	int fd;
	#endif

	unsigned int ProgressResolution;
	NMax > 3000 ? ProgressResolution = NMax/100 : ProgressResolution = 1;
	std::cout << "Simulation started..." << std::endl;

	for (n = 1; n < NMax; n++)
	{
		CUT_SAFE_CALL( cutStartTimer(hTimer) );
		for ( unsigned int step=0; step < 2; step++)
		{
			// Kernel call. 14 data pointers. 23 non-pointer arguments.
			FDTD2DKernel <16, 16> <<<Blocks, Threads>>>(
												d_Hx,
												d_Bx,
												d_Hy,
												d_By,
												d_Ez,
												d_Dz,
												d_Dzx,
												d_Dzy,
												d_urHx,
												d_urHy,
												d_erEz,
												d_ScmHx,
												d_ScmHy,
												d_Sc,
												d_Scsx,
												d_Scsy,
												d_ScmsmxHy,
												d_ScmsmyHx,
												delta,
												dtscalar,
												dt,
												PMLw,
												e0,
												u0,
												Two_pi_f_deltat,
												NHW,
												Is,
												Js,
												IHx,
												JHx,
												IHy,
												JHy,
												IEz,
												JEz,
												n,
												n0,
												n1,
												n2,
												flagHalf);


			CUT_CHECK_ERROR("FDTD2DKernel() execution failed\n");
			CUDA_SAFE_CALL( cudaThreadSynchronize() );
			flagHalf = !flagHalf;
		}
		CUT_SAFE_CALL( cutStopTimer(hTimer) );
		
		// Write field snapshot.
		if (n % tResolution == 0 && SaveFields == true)
		{
			// Copy the data back to the host
			CUDA_SAFE_CALL( cudaMemcpy(Ez, d_Ez, sizeof(float) * IEz*JEz*3, cudaMemcpyDeviceToHost) );
		
			framestream.str(std::string());			// Clearing stringstream contents.
			framestream << frame;
			filename = basename + framestream.str() + ".fdt";

			#ifdef WIN32
			snapshot.open (filename.c_str(), std::ios::out|std::ios::binary);
			snapshot.write ( (char*)&(Ez[IEz*JEz*n2]), sizeof(double)*IEz*JEz);
			snapshot.close();
			#else
			fd = open (filename.c_str(), O_CREAT|O_WRONLY|O_TRUNC, S_IRWXU);
			write (fd, (void*)&(Ez[IEz*JEz*n2]), sizeof(double)*IEz*JEz);
			close (fd);
			#endif

			frame++;
		}
		n0 = (n0+1)%3;
		n1 = (n1+1)%3;
		n2 = (n2+1)%3;

		if (n%ProgressResolution == 0)
			std::cout << std::setprecision(4) << std::setw(3) << "\r" << (float)n/(NMax-1)*100 << "%";
	}
	std::cout << "\r" << "Simulation complete!" << std::endl;
	std::cout << "kernel execution time = " << cutGetTimerValue(hTimer) << " ms." << std::endl;
	CUT_SAFE_CALL(cutDeleteTimer(hTimer));

	return 0;
}
int CFDTD2D::RunSimulationCPU (bool SaveFields)
{
	// Time loop.
	unsigned int n, i, j;

	// File Handling.
	std::stringstream framestream;
	std::string basename = "../../FieldData/Ez";
	std::string filename;
	unsigned int frame = 1;

	#ifdef WIN32
	std::fstream snapshot;
	#else
	int fd;
	#endif

	unsigned int ProgressResolution;
	NMax > 3000 ? ProgressResolution = NMax/100 : ProgressResolution = 1;
	std::cout << "Simulation started..." << std::endl;

	for (n=0; n < NMax-1; n++)
	{
		if (n%ProgressResolution == 0)
			std::cout << std::setprecision(4) << std::setw(3) << "\r" << (float)n/(NMax-1)*100 << "%";

		// t = 1/2.
		// Magnetic field. IHx and JHx are one less than IHy and JHy.
		for (i=0; i < IHx; i++)
		{
			for (j=0; j < JHx; j++)
			{
				Bx[i+IHx*j+IHx*JHx*n2] = (1-ScmHx[i+IHx*j])/(1+ScmHx[i+IHx*j]) * Bx[i+IHx*j+IHx*JHx*n1] + ( (dt/delta)/(1+ScmHx[i+IHx*j]) * (Ez[i+IEz*j+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hx[i+IHx*j+IHx*JHx*n2] = Bx[i+IHx*j+IHx*JHx*n2]/(u0*urHx[i+IHx*j]);

				By[(i+1)+IHy*(j+1)+IHy*JHy*n2] = (1-ScmHy[(i+1)+IHy*(j+1)])/(1+ScmHy[(i+1)+IHy*(j+1)]) * By[(i+1)+IHy*(j+1)+IHy*JHy*n1] + ( (dt/delta)/(1+ScmHy[(i+1)+IHy*(j+1)]) * (Ez[(i+1)+IEz*(j+1)+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hy[(i+1)+IHy*(j+1)+IHy*JHy*n2] = By[(i+1)+IHy*(j+1)+IHy*JHy*n2]/(u0*urHy[(i+1)+IHy*(j+1)]);
			}
		}
		// t = 1.
		// Electric field.
		for (i=0; i < IEz; i++)
		{
			for (j=1; j < JEz-1; j++)
			{
				Dz[i+IEz*j+IEz*JEz*n2] = (1-Sc[i+IEz*j])/(1+Sc[i+IEz*j]) * Dz[i+IEz*j+IEz*JEz*n1] + ( (dt/delta)/(1+Sc[i+IEz*j]) * ( Hy[(i+1)+IHy*j+IHy*JHy*n2] - Hy[i+IHy*j+IHy*JHy*n2] - Hx[i+IHx*j+IHx*JHx*n2] + Hx[i+IHx*(j-1)+IHx*JHx*n2]) );
				Ez[i+IEz*j+IEz*JEz*n2] = Dz[i+IEz*j+IEz*JEz*n2]/(e0*erEz[i+IEz*j]);

				// Source.
				if (j == Js && n < NHW)
				{
					Ez[i+IEz*j+IEz*JEz*n2] = Ez[i+IEz*j+IEz*JEz*n2] + 1 * sin (Two_pi_f_deltat * n) / dtscalar;
					Dz[i+IEz*j+IEz*JEz*n2] = e0 * Ez[i+IEz*j+IEz*JEz*n2];
				}
			}
		}

		// Write field snapshot.
		if (n % tResolution == 0 && SaveFields == true)
		{
			framestream.str(std::string());			// Clearing stringstream contents.
			framestream << frame;
			filename = basename + framestream.str() + ".fdt";

			#ifdef WIN32
			snapshot.open (filename.c_str(), std::ios::out|std::ios::binary);
			snapshot.write ( (char*)&(Ez[IEz*JEz*n2]), sizeof(float)*IEz*JEz);
			snapshot.close();
			#else
			fd = open (filename.c_str(), O_CREAT|O_WRONLY|O_TRUNC, S_IRWXU);
			write (fd, (void*)&(Ez[IEz*JEz*n2]), sizeof(float)*IEz*JEz);
			close (fd);
			#endif

			frame++;
		}
		n0 = (n0+1)%3;
		n1 = (n1+1)%3;
		n2 = (n2+1)%3;
	}
	std::cout << "\r" << "Simulation complete!" << std::endl;
	return 0;
}

void CFDTD2D::DisplaySimulationParameters ()
{
  	std::cout << "======= Simulation Parameters =======" << std::endl;
	std::cout << "I = " << I << std::endl;
	std::cout << "J = " << J << std::endl;
	std::cout << "NMax = " << NMax << std::endl;
	std::cout << "f = " << f << std::endl;
	std::cout << "delta = " << delta << std::endl;
	std::cout << "dt = " << dt << std::endl;
	std::cout << "t rez = " << tResolution << std::endl;
	std::cout << "x rez = " << xResolution << std::endl;
	std::cout << "y rez = " << yResolution << std::endl;
	std::cout << "=====================================" << std::endl;
}

int CFDTD2D::Cleanup ()
{
	// Device cleanup.
	CUDA_SAFE_CALL( cudaFree(d_Hx ) );
	CUDA_SAFE_CALL( cudaFree(d_Bx) );
	CUDA_SAFE_CALL( cudaFree(d_Hy) );
	CUDA_SAFE_CALL( cudaFree(d_By) );
	CUDA_SAFE_CALL( cudaFree(d_Ez) );
	CUDA_SAFE_CALL( cudaFree(d_Dz) );
	CUDA_SAFE_CALL( cudaFree(d_Dzx) );
	CUDA_SAFE_CALL( cudaFree(d_Dzy) );
	
	CUDA_SAFE_CALL( cudaFree(d_urHx) );
	CUDA_SAFE_CALL( cudaFree(d_urHy) );
	CUDA_SAFE_CALL( cudaFree(d_erEz) );
	
	CUDA_SAFE_CALL( cudaFree(d_ScmHx) );
	CUDA_SAFE_CALL( cudaFree(d_ScmHy) );
	CUDA_SAFE_CALL( cudaFree(d_Sc) );

	CUDA_SAFE_CALL( cudaFree(d_smHx) );
	CUDA_SAFE_CALL( cudaFree(d_smHy) );
	CUDA_SAFE_CALL( cudaFree(d_sEz) );
	
	CUDA_SAFE_CALL( cudaFree(d_Scsx) );
	CUDA_SAFE_CALL( cudaFree(d_Scsy) );
	CUDA_SAFE_CALL( cudaFree(d_ScmsmxHy) );
	CUDA_SAFE_CALL( cudaFree(d_ScmsmyHx) );
	
	return 0;
}

CFDTD2D::~CFDTD2D ()
{
	if (Hx != NULL) delete[] Hx;
	if (Bx != NULL) delete[] Bx;
	if (Hy != NULL) delete[] Hy;
	if (By != NULL) delete[] By;
	if (Ez != NULL) delete[] Ez;
	if (Dz != NULL) delete[] Dz;
	if (Dzx != NULL) delete[] Dzx;
	if (Dzy != NULL) delete[] Dzy;
	if (EzSnapshots != NULL) delete[] EzSnapshots;
	if (urHx != NULL) delete[] urHx;
	if (urHy != NULL) delete[] urHy;
	if (erEz != NULL) delete[] erEz;
	if (smHx != NULL) delete[] smHx;
	if (smHy != NULL) delete[] smHy;
	if (sEz != NULL) delete[] sEz;
	if (Sc != NULL) delete[] Sc;
	if (ScmHx != NULL) delete[] ScmHx;
	if (ScmHy != NULL) delete[] ScmHy;
	if (sex != NULL) delete[] sex;
	if (sey != NULL) delete[] sey;
	if (smx != NULL) delete[] smx;
	if (smy != NULL) delete[] smy;
	if (Scsx != NULL) delete[] Scsx;
	if (Scsy != NULL) delete[] Scsy;
	if (ScmsmxHy != NULL) delete[] ScmsmxHy;
	if (ScmsmyHx != NULL) delete[] ScmsmyHx;
}
