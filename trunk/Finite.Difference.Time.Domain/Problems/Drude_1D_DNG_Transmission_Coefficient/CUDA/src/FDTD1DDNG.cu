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
#include <helper_cuda.h>
#include <helper_functions.h> // helper functions for SDK examples
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
// Allocate memory for data arrays.
int CFDTD1DDNG::AllocateMemoryCPU()
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

	return 0;
}
// Initialise CPU data.
int CFDTD1DDNG::InitialiseCPU()
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
	return 0;
}
int CFDTD1DDNG::InitialiseExHyCPU()
{
	for (unsigned int i=0; i<3*Size; i++)
	{
		Ex_[i] = 0.;
		Hy_[i] = 0.;
	}
	return 0;
}
int CFDTD1DDNG::AllocateMemoryGPU()
{
	// Device memory allocation

	// Data arrays.
	checkCudaErrors(cudaMalloc((void **)&d_Ex_, sizeof(PRECISION)*Size*3));
	checkCudaErrors(cudaMalloc((void **)&d_Dx_, sizeof(PRECISION)*Size*3));
	checkCudaErrors(cudaMalloc((void **)&d_Hy_, sizeof(PRECISION)*Size*3));
	checkCudaErrors(cudaMalloc((void **)&d_By_, sizeof(PRECISION)*Size*3));
	// Incident and transmitted fields.
	checkCudaErrors(cudaMalloc((void **)&d_Exi, sizeof(PRECISION)*MaxTime));
	checkCudaErrors(cudaMalloc((void **)&d_Ext, sizeof(PRECISION)*MaxTime));
	checkCudaErrors(cudaMalloc((void **)&d_Extt, sizeof(PRECISION)*MaxTime));
	// Refractive Index.
	checkCudaErrors(cudaMalloc((void **)&d_Exz1, sizeof(PRECISION)*MaxTime));
	checkCudaErrors(cudaMalloc((void **)&d_Exz2, sizeof(PRECISION)*MaxTime));
	// Drude material parameters.
	checkCudaErrors(cudaMalloc((void **)&d_einf, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_uinf, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_wpesq, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_wpmsq, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_ge, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_gm, sizeof(PRECISION)*Size));
	// Auxiliary field scalars.
	checkCudaErrors(cudaMalloc((void **)&d_ae0, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_ae, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_be, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_ce, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_de, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_ee, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_am0, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_am, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_bm, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_cm, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_dm, sizeof(PRECISION)*Size));
	checkCudaErrors(cudaMalloc((void **)&d_em, sizeof(PRECISION)*Size));

	return 0;
}
int CFDTD1DDNG::CopyDataCPUtoGPU()
{
	// Device data initialisation.

	// Data arrays.
	checkCudaErrors(cudaMemcpy(d_Ex_, Ex_, sizeof(PRECISION)*Size*3, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Dx_, Dx_, sizeof(PRECISION)*Size*3, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Hy_, Hy_, sizeof(PRECISION)*Size*3, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_By_, By_, sizeof(PRECISION)*Size*3, cudaMemcpyHostToDevice));
	// Incident and transmitted fields.
	checkCudaErrors(cudaMemcpy(d_Exi, Exi, sizeof(PRECISION)*MaxTime, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Ext, Ext, sizeof(PRECISION)*MaxTime, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Extt, Extt, sizeof(PRECISION)*MaxTime, cudaMemcpyHostToDevice));
	// Refractive Index.
	checkCudaErrors(cudaMemcpy(d_Exz1, Exz1, sizeof(PRECISION)*MaxTime, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Exz2, Exz2, sizeof(PRECISION)*MaxTime, cudaMemcpyHostToDevice));
	// Drude material parameters.
	checkCudaErrors(cudaMemcpy(d_einf, einf, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_uinf, uinf, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_wpesq, wpesq, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_wpmsq, wpmsq, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_ge, ge, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_gm, gm, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	// Auxiliary field scalars.
	checkCudaErrors(cudaMemcpy(d_ae0, ae0, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_ae, ae, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_be, be, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_ce, ce, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_de, de, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_ee, ee, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_am0, am0, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_am, am, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_bm, bm, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_cm, cm, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_dm, dm, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_em, em, sizeof(PRECISION)*Size, cudaMemcpyHostToDevice));

	return 0;
}
int CFDTD1DDNG::CopyExHyCPUtoGPU()
{
	checkCudaErrors(cudaMemcpy(d_Ex_, Ex_, sizeof(PRECISION)*Size*3, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Hy_, Hy_, sizeof(PRECISION)*Size*3, cudaMemcpyHostToDevice));

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
		if (n % (MaxTime/1024U) == 0)
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

	cout << "Dry run (GPU) started..." << endl;

	cout << "Block dimensions: " << ThreadsX << "x" << ThreadsY << endl;
	cout << "Grid dimensions: " << BlocksX << "x" << BlocksY << endl;

	StopWatchInterface *Timer = 0;
	sdkCreateTimer(&Timer);
	sdkResetTimer(&Timer);

	for (unsigned int n=0;n<MaxTime; n++)
	{
		if (n % (MaxTime/1024U) == 0)
			cout << "\r\t\t\r" << n*100/(MaxTime-1) << "%";
		sdkStartTimer(&Timer);

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

		getLastCudaError("Kernel execution failed");
		checkCudaErrors(cudaThreadSynchronize());
		sdkStopTimer(&Timer);

		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}
	cout << "\r" << "Dry run kernel execution time = " << sdkGetTimerValue(&Timer) << " ms." << endl;
	sdkDeleteTimer(&Timer);

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

	StopWatchInterface *Timer = 0;
	sdkCreateTimer(&Timer);
	sdkResetTimer(&Timer);

	stringstream framestream;
	string basename = "FieldData/Ex";
	string filename;
	fstream snapshot;
	frame = 0U;

	for (unsigned int n=0;n<MaxTime; n++)
	{
		if (n % (MaxTime/1024U) == 0)
			cout << "\r\t\t\r" << n*100/(MaxTime-1) << "%";
		sdkStartTimer(&Timer);

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

		getLastCudaError("Kernel execution failed");
		checkCudaErrors(cudaThreadSynchronize());
		sdkStopTimer(&Timer);

		// Saving electric field snapshot.
		if (n%SnapshotInterval == 0 && SaveFields == true)
		{
			// Write E-field to file.
			framestream.str(std::string());			// Clearing stringstream contents.
			framestream << ++frame;
			filename = basename + framestream.str() + ".fdt";
			snapshot.open(filename.c_str(), std::ios::out|std::ios::binary);
			checkCudaErrors(cudaMemcpy(Ex_, d_Ex_, sizeof(PRECISION)*Size*3, cudaMemcpyDeviceToHost));
			snapshot.write((char*)&(Ex(0,nf)), sizeof(PRECISION)*Size);
			snapshot.close();
		}

		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}
	cout << "\r" << "Simulation run kernel execution time = " << sdkGetTimerValue(&Timer) << " ms." << endl;
	sdkDeleteTimer(&Timer);

	// Saving electric field data arrays.
	if (SaveFields == true)
	{
		fstream parametersfile;
		parametersfile.open("FieldData/Parameters.smp", std::ios::out|std::ios::binary|std::ios::app);
		parametersfile.write((char*)&(frame), sizeof(unsigned int));
		parametersfile.close();
		// Write saved fields to files.
		snapshot.open("FieldData/Exi.fdt", std::ios::out|std::ios::binary);
		checkCudaErrors(cudaMemcpy(Exi, d_Exi, sizeof(PRECISION)*MaxTime, cudaMemcpyDeviceToHost));
		snapshot.write((char*)Exi, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ext.fdt", std::ios::out|std::ios::binary);
		checkCudaErrors(cudaMemcpy(Ext, d_Ext, sizeof(PRECISION)*MaxTime, cudaMemcpyDeviceToHost));
		snapshot.write((char*)Ext, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Extt.fdt", std::ios::out|std::ios::binary);
		checkCudaErrors(cudaMemcpy(Extt, d_Extt, sizeof(PRECISION)*MaxTime, cudaMemcpyDeviceToHost));
		snapshot.write((char*)Extt, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Exz1.fdt", std::ios::out|std::ios::binary);
		checkCudaErrors(cudaMemcpy(Exz1, d_Exz1, sizeof(PRECISION)*MaxTime, cudaMemcpyDeviceToHost));
		snapshot.write((char*)Exz1, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Exz2.fdt", std::ios::out|std::ios::binary);
		checkCudaErrors(cudaMemcpy(Exz2, d_Exz2, sizeof(PRECISION)*MaxTime, cudaMemcpyDeviceToHost));
		snapshot.write((char*)Exz2, sizeof(PRECISION)*MaxTime);
		snapshot.close();
	}

	return 0;
}
int CFDTD1DDNG::CompleteRunCPU(bool SaveFields)
{
	cout << "Memory required for simulation = " << SimSize() << " bytes (" << (double)SimSize()/1024UL << "kB/" << (double)SimSize()/1024UL/1024UL << "MB)." << endl;
	cout << "HDD space required for data storage = " << HDDSpace() << " bytes (" << (double)HDDSpace()/1024UL << "kB/" << (double)HDDSpace()/1024UL/1024UL << "MB)." << endl;
	SafeCall(AllocateMemoryCPU(), "Error: Allocating memory on CPU.");
	SafeCall(InitialiseCPU(), "Error: Initialising data on CPU.");
	SafeCall(DryRunCPU(), "Error: Dry run (CPU).");
	SafeCall(InitialiseExHyCPU(), "Error: Initalising Ex/Hy arrays (CPU).");
	SafeCall(RunSimulationCPU(SaveFields), "Error: Running simulation (CPU).");
	SafeCall(CleanupCPU(), "Error: Cleaning up CPU.");

	return 0;
}
int CFDTD1DDNG::CompleteRunGPU(bool SaveFields)
{
	cout << "Memory required for simulation = " << SimSize() << " bytes (" << (double)SimSize()/1024UL << "kB/" << (double)SimSize()/1024UL/1024UL << "MB)." << endl;
	cout << "HDD space required for data storage = " << HDDSpace() << " bytes (" << (double)HDDSpace()/1024UL << "kB/" << (double)HDDSpace()/1024UL/1024UL << "MB)." << endl;
	SafeCall(AllocateMemoryCPU(), "Error: Allocating memory on CPU.");
	SafeCall(InitialiseCPU(), "Error: Initialising data on CPU.");

	SafeCall(AllocateMemoryGPU(), "Error: Allocating memory on GPU.");
	SafeCall(CopyDataCPUtoGPU(), "Error: Copying data from CPU to GPU.");
	SafeCall(DryRunGPU(), "Error: Dry run (GPU).");
	SafeCall(InitialiseExHyCPU(), "Error: Initalising Ex/Hy arrays (CPU).");
	SafeCall(CopyExHyCPUtoGPU(), "Error: Copying Ex/Hy arrays from CPU to GPU.");
	SafeCall(RunSimulationGPU(SaveFields), "Error: Running simulation on GPU.");
	SafeCall(CleanupCPU(), "Error: Cleaning up CPU.");
	SafeCall(CleanupGPU(), "Error: Cleaning up GPU.");

	return 0;
}

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
int CFDTD1DDNG::SafeCall(int Status, const char *Error)
{
	if (Status != 0)
	{
		if (Error!=NULL) cout << Error << endl;
		exit(Status);
	}
	return Status;
}
template<typename T> void DeleteArray(T *&ptr)
{
	if (ptr != NULL)
	{
		delete[] ptr;
		ptr = NULL;
	}
}
int CFDTD1DDNG::CleanupCPU()
{
	// Field arrays.
	DeleteArray(Ex_);
	DeleteArray(Dx_);
	DeleteArray(Hy_);
	DeleteArray(By_);
	// Incident and transmitted fields.
	DeleteArray(Exi);
	DeleteArray(Ext);
	DeleteArray(Extt);
	// Refractive index.
	DeleteArray(Exz1);
	DeleteArray(Exz2);
	// Drude parameter arrays.
	DeleteArray(einf);
	DeleteArray(uinf);
	DeleteArray(wpesq);
	DeleteArray(wpmsq);
	DeleteArray(ge);
	DeleteArray(gm);
	// Auxiliary field scalars.
	DeleteArray(ae0);
	DeleteArray(ae);
	DeleteArray(be);
	DeleteArray(ce);
	DeleteArray(de);
	DeleteArray(ee);
	DeleteArray(am0);
	DeleteArray(am);
	DeleteArray(bm);
	DeleteArray(cm);
	DeleteArray(dm);
	DeleteArray(em);

	return 0;
}
int CFDTD1DDNG::CleanupGPU()
{
	// Device cleanup.
	checkCudaErrors(cudaFree(d_Ex_));
	checkCudaErrors(cudaFree(d_Dx_));
	checkCudaErrors(cudaFree(d_Hy_));
	checkCudaErrors(cudaFree(d_By_));

	checkCudaErrors(cudaFree(d_Exi));
	checkCudaErrors(cudaFree(d_Ext));
	checkCudaErrors(cudaFree(d_Extt));
	checkCudaErrors(cudaFree(d_Exz1));
	checkCudaErrors(cudaFree(d_Exz2));

	checkCudaErrors(cudaFree(d_einf));
	checkCudaErrors(cudaFree(d_uinf));
	checkCudaErrors(cudaFree(d_wpesq));
	checkCudaErrors(cudaFree(d_wpmsq));
	checkCudaErrors(cudaFree(d_ge));
	checkCudaErrors(cudaFree(d_gm));

	checkCudaErrors(cudaFree(d_ae0));
	checkCudaErrors(cudaFree(d_ae));
	checkCudaErrors(cudaFree(d_be));
	checkCudaErrors(cudaFree(d_ce));
	checkCudaErrors(cudaFree(d_de));
	checkCudaErrors(cudaFree(d_ee));
	checkCudaErrors(cudaFree(d_am0));
	checkCudaErrors(cudaFree(d_am));
	checkCudaErrors(cudaFree(d_bm));
	checkCudaErrors(cudaFree(d_cm));
	checkCudaErrors(cudaFree(d_dm));
	checkCudaErrors(cudaFree(d_em));

	cudaDeviceReset();

	return 0;
}
CFDTD1DDNG::~CFDTD1DDNG ()
{
	// Field arrays.
	DeleteArray(Ex_);
	DeleteArray(Dx_);
	DeleteArray(Hy_);
	DeleteArray(By_);
	// Incident and transmitted fields.
	DeleteArray(Exi);
	DeleteArray(Ext);
	DeleteArray(Extt);
	// Refractive index.
	DeleteArray(Exz1);
	DeleteArray(Exz2);
	// Drude parameter arrays.
	DeleteArray(einf);
	DeleteArray(uinf);
	DeleteArray(wpesq);
	DeleteArray(wpmsq);
	DeleteArray(ge);
	DeleteArray(gm);
	// Auxiliary field scalars.
	DeleteArray(ae0);
	DeleteArray(ae);
	DeleteArray(be);
	DeleteArray(ce);
	DeleteArray(de);
	DeleteArray(ee);
	DeleteArray(am0);
	DeleteArray(am);
	DeleteArray(bm);
	DeleteArray(cm);
	DeleteArray(dm);
	DeleteArray(em);
}
