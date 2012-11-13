// Macros for indexing data arrays.
#define Ex(i,n) Ex_[(i)+Size*(n)]
#define Dx(i,n) Dx_[(i)+Size*(n)]
#define Hy(i,n) Hy_[(i)+Size*(n)]
#define By(i,n) By_[(i)+Size*(n)]

#include <FDTD1DDNG.hpp>
#include <Timer.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <string>
#include <FDTD1DDNG_Kernels.cu>
#include <cutil.h>
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
							z1((PRECISION)Z1*dz),
							Z2(SlabLeft+60),
							z2((PRECISION)Z2*dz),
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
	parametersfile.write((char*)&dt, sizeof(PRECISION));
	parametersfile.write((char*)&k0, sizeof(PRECISION));
	parametersfile.write((char*)&z1, sizeof(PRECISION));
	parametersfile.write((char*)&z2, sizeof(PRECISION));
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
	return (unsigned long)sizeof(*this)+(unsigned long)sizeof(PRECISION)*(30UL*(unsigned long)Size+5UL*(unsigned long)MaxTime);
}
unsigned long CFDTD1DDNG::HDDSpace()
{
	return (unsigned long)sizeof(PRECISION)*(5UL*(unsigned long)MaxTime+(unsigned long)Size*((unsigned long)MaxTime/(unsigned long)SnapshotInterval+1UL));
}
// Initialize data arrays.
void CFDTD1DDNG::AllocateMemoryCPU()
{
	// Field arrays.
	Ex_ = new PRECISION[Size*3];
	Dx_ = new PRECISION[Size*3];
	Hy_ = new PRECISION[Size*3];
	By_ = new PRECISION[Size*3];
	// Incident and transmitted fields.
	Exi = new PRECISION[MaxTime];
	Ext = new PRECISION[MaxTime];
	Extt = new PRECISION[MaxTime];
	// Refractive index.
	Exz1 = new PRECISION[MaxTime];
	Exz2 = new PRECISION[MaxTime];
	// Drude parameter arrays.
	einf = new PRECISION[Size];
	uinf = new PRECISION[Size];
	wpesq = new PRECISION[Size];
	wpmsq = new PRECISION[Size];
	ge = new PRECISION[Size];
	gm = new PRECISION[Size];
	// Auxiliary field scalars.
	ae0 = new PRECISION[Size];
	ae = new PRECISION[Size];
	be = new PRECISION[Size];
	ce = new PRECISION[Size];
	de = new PRECISION[Size];
	ee = new PRECISION[Size];
	am0 = new PRECISION[Size];
	am = new PRECISION[Size];
	bm = new PRECISION[Size];
	cm = new PRECISION[Size];
	dm = new PRECISION[Size];
	em = new PRECISION[Size];
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
int CFDTD1DDNG::AllocateMemoryGPU ()
{
	// Device memory allocation

	// Data arrays.
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_Ex_, sizeof(PRECISION)*Size*3));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_Dx_, sizeof(PRECISION)*Size*3));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_Hy_, sizeof(PRECISION)*Size*3));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_By_, sizeof(PRECISION)*Size*3));
	// Incident and transmitted fields.
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_Exi, sizeof(PRECISION)*MaxTime));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_Ext, sizeof(PRECISION)*MaxTime));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_Extt, sizeof(PRECISION)*MaxTime));
	// Refractive Index.
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_Exz1, sizeof(PRECISION)*MaxTime));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_Exz2, sizeof(PRECISION)*MaxTime));
	// Drude material parameters.
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_einf, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_uinf, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_wpesq, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_wpmsq, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_ge, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_gm, sizeof(PRECISION)*Size));
	// Auxiliary field scalars.
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_ae0, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_ae, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_be, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_ce, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_de, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_ee, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_am0, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_am, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_bm, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_cm, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_dm, sizeof(PRECISION)*Size));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_em, sizeof(PRECISION)*Size));

	return 0;
}
int CFDTD1DDNG::CopyDataCPUtoGPU()
{
	// Device data initialisation.

	// Data arrays.
	CUDA_SAFE_CALL(cudaMemcpy(d_Ex_, Ex_, sizeof(PRECISION)*Size*3, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_Dx_, Dx_, sizeof(PRECISION)*Size*3, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_Hy_, Hy_, sizeof(PRECISION)*Size*3, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_By_, By_, sizeof(PRECISION)*Size*3, cudaMemcpyHostToDevice));
	// Incident and transmitted fields.
	CUDA_SAFE_CALL(cudaMemcpy(d_Exi, Exi, sizeof(PRECISION)*MaxTime, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_Ext, Ext, sizeof(PRECISION)*MaxTime, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_Extt, Extt, sizeof(PRECISION)*MaxTime, cudaMemcpyHostToDevice));
	// Refractive Index.
	CUDA_SAFE_CALL(cudaMemcpy(d_Exz1, Exz1, sizeof(PRECISION)*MaxTime, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_Exz2, Exz2, sizeof(PRECISION)*MaxTime, cudaMemcpyHostToDevice));
	// Drude material parameters.
	CUDA_SAFE_CALL(cudaMemcpy(d_einf, einf, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_uinf, uinf, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_wpesq, wpesq, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_wpmsq, wpmsq, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_ge, ge, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_gm, gm, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	// Auxiliary field scalars.
	CUDA_SAFE_CALL(cudaMemcpy(d_ae0, ae0, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_ae, ae, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_be, be, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_ce, ce, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_de, de, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_ee, ee, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_am0, am0, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_am, am, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_bm, bm, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_cm, cm, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_dm, dm, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_em, em, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));

	return 0;
}
int CFDTD1DDNG::CopyExHyCPUtoGPU()
{
	CUDA_SAFE_CALL(cudaMemcpy(d_Ex_, Ex_, sizeof(PRECISION)*Size*3, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_Hy_, Hy_, sizeof(PRECISION)*Size*3, cudaMemcpyHostToDevice));

	return 0;
}
int CFDTD1DDNG::DryRunCPU()
{
	cout << "Dry run (CPU) started..." << endl;
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
		Ex(0,nf) = Ex(1,n0) + (Sc-1)/(Sc+1)*(Ex(1,nf)-Ex(1,n0));

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
	frame = 0U;

	cout << "Simulation (CPU) started..." << endl;
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
		Ex(0,nf) = Ex(1,n0) + (Sc-1)/(Sc+1)*(Ex(1,nf)-Ex(1,n0));
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
	// Saving electric fields.
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
	}
	cout << endl << "Simulation (CPU) completed!" << endl;
	return 0;
}
int CFDTD1DDNG::DryRunGPU()
{
	// Total local threads in a block. Can be thought of as Block dimensions.
	const unsigned int ThreadsX = 256;
	const unsigned int ThreadsY = 1;
	
	// Total blocks in simulation grid. Can be thought of as no. of blocks in grid.
	// Size should be divisible by 256.
	unsigned int BlocksX = Size/ThreadsX;
	unsigned int BlocksY = 1;

	// Kernel parameters.
	dim3 Blocks(BlocksX, BlocksY);
	dim3 Threads(ThreadsX, ThreadsY);

	cout << "Dry run (GPU) started..." << std::endl;

	cout << "Block dimensions: " << ThreadsX << "x" << ThreadsY << endl;
	cout << "Grid dimensions: " << BlocksX << "x" << BlocksY << endl;

	unsigned int hTimer;
	CUT_SAFE_CALL(cutCreateTimer(&hTimer));
	CUT_SAFE_CALL(cutResetTimer(hTimer));

	for (unsigned int n=0;n<MaxTime; n++)
	{
		if (n % (MaxTime/1024U) == 0)
			cout << "\r\t\t\r" << n*100/(MaxTime-1) << "%";
		CUT_SAFE_CALL(cutStartTimer(hTimer));

		// Kernel call.
		FDTD1DDNGKernel_DryRun <ThreadsX, ThreadsY> <<<Blocks, Threads>>>(
											Size,
											PulseWidth,
											td,
											SourceLocation,
											SourceChoice,
											e0,
											u0,
											dt,
											dz,
											Sc,
											f,
											fp,
											dr,
											d_Ex_,
											d_Hy_,
											d_Exi,
											x1,
											n,
											np,
											n0,
											nf);

		CUT_CHECK_ERROR("FDTDKernel() execution failed\n");
		CUDA_SAFE_CALL(cudaThreadSynchronize());
		CUT_SAFE_CALL( cutStopTimer(hTimer) );

		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}
	cout << "\r" << "Dry run kernel execution time = " << cutGetTimerValue(hTimer) << " ms." << endl;
	CUT_SAFE_CALL(cutDeleteTimer(hTimer));

	return 0;
}
int CFDTD1DDNG::RunSimulationGPU(bool SaveFields)
{
	// Total local threads in a block. Can be thought of as Block dimensions.
	const unsigned int ThreadsX = 256;
	const unsigned int ThreadsY = 1;
	
	// Total blocks in simulation grid. Can be thought of as no. of blocks in grid.
	// Size should be divisible by 256.
	unsigned int BlocksX = Size/ThreadsX;
	unsigned int BlocksY = 1;

	// Kernel parameters.
	dim3 Blocks(BlocksX, BlocksY);
	dim3 Threads(ThreadsX, ThreadsY);

	cout << "Simulation (GPU) started..." << endl;

	cout << "Block dimensions: " << ThreadsX << "x" << ThreadsY << endl;
	cout << "Grid dimensions: " << BlocksX << "x" << BlocksY << endl;

	unsigned int hTimer;
	CUT_SAFE_CALL(cutCreateTimer(&hTimer));
	CUT_SAFE_CALL(cutResetTimer(hTimer));

	stringstream framestream;
	string basename = "FieldData/Ex";
	string filename;
	fstream snapshot;
	frame = 0U;

	for (unsigned int n=0;n<MaxTime; n++)
	{
		if (n % (MaxTime/1024U) == 0)
			cout << "\r\t\t\r" << n*100/(MaxTime-1) << "%";
		CUT_SAFE_CALL(cutStartTimer(hTimer));

		// Kernel call.
		FDTD1DDNGKernel_Simulation <ThreadsX, ThreadsY> <<<Blocks, Threads>>>(
											Size,
											PulseWidth,
											td,
											SourceLocation,
											SourceChoice,
											e0,
											u0,
											dt,
											dz,
											Sc,
											f,
											fp,
											dr,
											d_Ex_, d_Dx_, d_Hy_, d_By_,
											d_einf, d_uinf, d_wpesq, d_wpmsq, d_ge, d_gm,
											d_ae0, d_ae, d_be, d_ce, d_de, d_ee,
											d_am0, d_am, d_bm, d_cm, d_dm, d_em,
											d_Ext, d_Extt, d_Exz1, d_Exz2,
											x1, Z1, Z2,
											n,
											np,
											n0,
											nf);

		CUT_CHECK_ERROR("FDTDKernel() execution failed\n");
		CUDA_SAFE_CALL(cudaThreadSynchronize());
		CUT_SAFE_CALL(cutStopTimer(hTimer));

		// Saving electric field snapshot.
		if (n%SnapshotInterval == 0 && SaveFields == true)
		{
			// Write E-field to file.
			framestream.str(std::string());			// Clearing stringstream contents.
			framestream << ++frame;
			filename = basename + framestream.str() + ".fdt";
			snapshot.open(filename.c_str(), std::ios::out|std::ios::binary);
			CUDA_SAFE_CALL(cudaMemcpy(Ex_, d_Ex_, sizeof(PRECISION)*Size*3, cudaMemcpyDeviceToHost));
			snapshot.write((char*)&(Ex(0,nf)), sizeof(PRECISION)*Size);
			snapshot.close();
		}

		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}
	cout << "\r" << "Simulation run kernel execution time = " << cutGetTimerValue(hTimer) << " ms." << endl;
	CUT_SAFE_CALL(cutDeleteTimer(hTimer));

	// Saving electric field data arrays.
	if (SaveFields == true)
	{
		fstream parametersfile;
		parametersfile.open("FieldData/Parameters.smp", std::ios::out|std::ios::binary|std::ios::app);
		parametersfile.write((char*)&(frame), sizeof(unsigned int));
		parametersfile.close();
		// Write saved fields to files.
		snapshot.open("FieldData/Exi.fdt", std::ios::out|std::ios::binary);
		CUDA_SAFE_CALL(cudaMemcpy(Exi, d_Exi, sizeof(PRECISION)*MaxTime, cudaMemcpyDeviceToHost));
		snapshot.write((char*)Exi, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ext.fdt", std::ios::out|std::ios::binary);
		CUDA_SAFE_CALL(cudaMemcpy(Ext, d_Ext, sizeof(PRECISION)*MaxTime, cudaMemcpyDeviceToHost));
		snapshot.write((char*)Ext, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Extt.fdt", std::ios::out|std::ios::binary);
		CUDA_SAFE_CALL(cudaMemcpy(Extt, d_Extt, sizeof(PRECISION)*MaxTime, cudaMemcpyDeviceToHost));
		snapshot.write((char*)Extt, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Exz1.fdt", std::ios::out|std::ios::binary);
		CUDA_SAFE_CALL(cudaMemcpy(Exz1, d_Exz1, sizeof(PRECISION)*MaxTime, cudaMemcpyDeviceToHost));
		snapshot.write((char*)Exz1, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Exz2.fdt", std::ios::out|std::ios::binary);
		CUDA_SAFE_CALL(cudaMemcpy(Exz2, d_Exz2, sizeof(PRECISION)*MaxTime, cudaMemcpyDeviceToHost));
		snapshot.write((char*)Exz2, sizeof(PRECISION)*MaxTime);
		snapshot.close();
	}

	return 0;
}
/*
int CFDTD1DDNG::runFDTD1DDNGKernels (bool SaveFields)
{
	// Total local threads in a block. Can be thought of as Block dimensions.
	unsigned int ThreadsX = 256;
	unsigned int ThreadsY = 1;
	
	// Total blocks in simulation grid. Can be thought of as no. of blocks in grid.
	// Obviously, I and J should be divisible by block dimensions.
	unsigned int BlocksX = I/ThreadsX;
	unsigned int BlocksY = (J+2*PMLw)/ThreadsY;

	// Kernel parameters.
	dim3 Blocks(BlocksX, BlocksY);
	dim3 Threads(ThreadsX, ThreadsY);

	cout << "Block dimensions: " << ThreadsX << "x" << ThreadsY << std::endl;
	cout << "Grid dimensions: " << BlocksX << "x" << BlocksY << std::endl;

	unsigned int hTimer;
	CUT_SAFE_CALL(cutCreateTimer(&hTimer));
	CUT_SAFE_CALL(cutResetTimer(hTimer));
		
	// File handling from chapter 3 of Understanding FDTD. J. B. Schneider
	fstream snapshot;
	stringstream framestream;
	string basename = "FieldData/Ex";
	string filename;
	unsigned int frame = 1;

	unsigned int ProgressResolution;
	NMax > 3000 ? ProgressResolution = NMax/100 : ProgressResolution = 1;
	cout << "Simulation (GPU) started..." << std::endl;

	for (unsigned int n=0;n<MaxTime; n++)
	{
		CUT_SAFE_CALL( cutStartTimer(hTimer) );
		for ( unsigned int step=0; step < 2; step++)
		{
			// Kernel call. 14 data pointers. 23 non-pointer arguments.
			FDTD1DDNGKernel <16, 16> <<<Blocks, Threads>>>(
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
			snapshot.write ( (char*)&(Ez[IEz*JEz*n2]), sizeof(PRECISION)*IEz*JEz);
			snapshot.close();
			#else
			fd = open (filename.c_str(), O_CREAT|O_WRONLY|O_TRUNC, S_IRWXU);
			write (fd, (void*)&(Ez[IEz*JEz*n2]), sizeof(PRECISION)*IEz*JEz);
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
int CFDTD1DDNG::RunSimulationCPU (bool SaveFields)
{
	// Time loop.
	unsigned int n, i, j;

	// File Handling.
	std::stringstream framestream;
	std::string basename = "../FieldData/Ez";
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

void CFDTD1DDNG::DisplaySimulationParameters ()
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
}*/
// Timing.
void CFDTD1DDNG::StartTimer()
{
	tStart = GetTimeus64();
}
void CFDTD1DDNG::StopTimer()
{
	tEnd = GetTimeus64();
}
PRECISION CFDTD1DDNG::GetElapsedTime()
{
	return ((PRECISION)(tEnd-tStart))/(1000000.);
}
int CFDTD1DDNG::CleanupGPU()
{
	// Device cleanup.
	CUDA_SAFE_CALL(cudaFree(d_Ex_));
	CUDA_SAFE_CALL(cudaFree(d_Dx_));
	CUDA_SAFE_CALL(cudaFree(d_Hy_));
	CUDA_SAFE_CALL(cudaFree(d_By_));

	CUDA_SAFE_CALL(cudaFree(d_Exi));
	CUDA_SAFE_CALL(cudaFree(d_Ext));
	CUDA_SAFE_CALL(cudaFree(d_Extt));
	CUDA_SAFE_CALL(cudaFree(d_Exz1));
	CUDA_SAFE_CALL(cudaFree(d_Exz2));

	CUDA_SAFE_CALL(cudaFree(d_einf));
	CUDA_SAFE_CALL(cudaFree(d_uinf));
	CUDA_SAFE_CALL(cudaFree(d_wpesq));
	CUDA_SAFE_CALL(cudaFree(d_wpmsq));
	CUDA_SAFE_CALL(cudaFree(d_ge));
	CUDA_SAFE_CALL(cudaFree(d_gm));

	CUDA_SAFE_CALL(cudaFree(d_ae0));
	CUDA_SAFE_CALL(cudaFree(d_ae));
	CUDA_SAFE_CALL(cudaFree(d_be));
	CUDA_SAFE_CALL(cudaFree(d_ce));
	CUDA_SAFE_CALL(cudaFree(d_de));
	CUDA_SAFE_CALL(cudaFree(d_ee));
	CUDA_SAFE_CALL(cudaFree(d_am0));
	CUDA_SAFE_CALL(cudaFree(d_am));
	CUDA_SAFE_CALL(cudaFree(d_bm));
	CUDA_SAFE_CALL(cudaFree(d_cm));
	CUDA_SAFE_CALL(cudaFree(d_dm));
	CUDA_SAFE_CALL(cudaFree(d_em));

	return 0;
}
CFDTD1DDNG::~CFDTD1DDNG ()
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
