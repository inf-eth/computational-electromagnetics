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
int CFDTD1DDNG::InitialiseCL()
{
	cl_int status = 0;
	size_t deviceListSize;

	/*
	* Have a look at the available platforms and pick either
	* the AMD one if available or a reasonable default.
	*/

	cl_uint numPlatforms;
	cl_platform_id platform = NULL;
	SafeCall(clGetPlatformIDs(0, NULL, &numPlatforms), "Error: Getting Platforms. (clGetPlatformsIDs)");

	char AMDPlatform[] = "Advanced Micro Devices, Inc.";
	char nVidiaPlatform[] = "NVIDIA Corporation";
	char *SelectedPlatform = NULL;

	char choice = '0';
	cout << "Choose a platform: " << endl;
	cout << "[1] Advanced Micro Devices, Inc. (default)" << endl;
	cout << "[2] NVIDIA Corporation" << endl;
	cout << ">>";
	StopTimer();
	cin >> choice;
	StartTimer();

	if (choice == '1')
		SelectedPlatform = AMDPlatform;
	else if (choice == '2')
		SelectedPlatform = nVidiaPlatform;
	else
	{
		cout << "Reverting to default platform..." << endl;
		SelectedPlatform = AMDPlatform;
	}

	cout << "Detecting platforms..." << endl;
	cout << "Available platforms are: " << endl;
	if(numPlatforms > 0)
	{
		cl_platform_id* platforms = new cl_platform_id[numPlatforms];
		SafeCall(clGetPlatformIDs(numPlatforms, platforms, NULL), "Error: Getting Platform Ids. (clGetPlatformsIDs)");

		for(unsigned int i=0; i < numPlatforms; ++i)
		{
			char pbuff[100];
			SafeCall(clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, sizeof(pbuff), pbuff, NULL), "Error: Getting Platform Info.(clGetPlatformInfo)");

			cout << "Platform " << i << " : " << pbuff << endl;
			if(!strcmp(pbuff, SelectedPlatform))
				platform = platforms[i];
		}
		delete platforms;
	}

	if(NULL == platform)
	{
		std::cout << "Selected platform not found so Exiting Application." << std::endl;
		return 1;
	}

	/*
	* If we could find our platform, use it. Otherwise use just available platform.
	*/
	cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0 };

	/////////////////////////////////////////////////////////////////
	// Create an OpenCL context
	/////////////////////////////////////////////////////////////////
	cl_device_type type;

	cout << "Emulate GPU run on CPU?" << endl;
	cout << "[1] Yes" << endl;
	cout << "[2] No (default)" << endl;
	cout << ">>";
	StopTimer();
	cin >> choice;
	StartTimer();

	if (choice == '1')
	{
		if(!strcmp(AMDPlatform, SelectedPlatform))
			cout << "Running on CPU with GPU emulation..." << endl;
		else
			cout << "Warning: Selected platform does not support GPU emulation on CPU." << endl;

		type = CL_DEVICE_TYPE_CPU;
	}
	else
	{
		cout << "Running on GPU..." << endl;
		type = CL_DEVICE_TYPE_GPU;
	}

	context = clCreateContextFromType(cps, type, NULL, NULL, &status);
	SafeCall(status, "Error: Creating Context. (clCreateContextFromType)");

	/* First, get the size of device list data */
	SafeCall(clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &deviceListSize), "Error: Getting Context Info (device list size, clGetContextInfo)");

	/////////////////////////////////////////////////////////////////
	// Detect OpenCL devices
	/////////////////////////////////////////////////////////////////
	devices = new cl_device_id[deviceListSize/sizeof(cl_device_id)];
	SafeCall(!devices, "Error: No devices found.");

	/* Now, get the device list data */
	SafeCall(clGetContextInfo(context, CL_CONTEXT_DEVICES, deviceListSize, devices, NULL), "Error: Getting Context Info (device list, clGetContextInfo)");

	char platformVendor[1024];
	SafeCall(clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, sizeof(platformVendor), platformVendor, NULL), "clGetPlatformInfo failed");
	cout << "Selected Platform Vendor : " << platformVendor << endl;

	// Get number of devices available 
	cl_uint deviceCount = 0;
	SafeCall(clGetDeviceIDs(platform, type, 0, NULL, &deviceCount), "clGetDeviceIDs failed");

	cl_device_id* deviceIds = (cl_device_id*)malloc(sizeof(cl_device_id) * deviceCount);
	SafeCall(!deviceIds, "Failed to allocate memory(deviceIds)");

	// Get device ids
	SafeCall(clGetDeviceIDs(platform, type, deviceCount, deviceIds, NULL), "clGetDeviceIDs failed");

	cout << "Available devices are: " << endl;
	// Print device index and device names
	for(cl_uint i = 0; i < deviceCount; ++i)
	{
		char deviceName[1024];
		SafeCall(clGetDeviceInfo(deviceIds[i], CL_DEVICE_NAME, sizeof(deviceName), deviceName, NULL), "clGetDeviceInfo failed");
		cout << "Device " << i << " : " << deviceName <<" Device ID is "<<deviceIds[i]<< endl;
	}
	free(deviceIds);
	/////////////////////////////////////////////////////////////////
	// Create an OpenCL command queue
	/////////////////////////////////////////////////////////////////
	cout << "Running on Device 0..." << endl;
	commandQueue = clCreateCommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE, &status);
	SafeCall(status, "Creating Command Queue. (clCreateCommandQueue)");

	return 0;
}

int CFDTD1DDNG::AllocateMemoryGPU()
{
	cl_int status;
	/////////////////////////////////////////////////////////////////
	// Create OpenCL memory buffers
	/////////////////////////////////////////////////////////////////
	// Fields.
	d_Ex_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size*3, Ex_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Ex_");
	d_Dx_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size*3, Dx_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Dx_");
	d_Hy_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size*3, Hy_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Hy_");
	d_By_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size*3, By_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for By_");
	// Incident and transmitted fields.
	d_Exi = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*MaxTime, Exi, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Exi");
	d_Ext = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*MaxTime, Ext, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Ext");
	d_Extt = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*MaxTime, Extt, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Extt");
	// Refractive Index.
	d_Exz1 = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*MaxTime, Exz1, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Exz1");
	d_Exz2 = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*MaxTime, Exz2, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Exz2");
	// Drude material parameters.
	d_einf = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, einf, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for einf");
	d_uinf = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, uinf, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for uinf");
	d_wpesq = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, wpesq, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for wpesq");
	d_wpmsq = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, wpmsq, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for wpmsq");
	d_ge = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, ge, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for ge");
	d_gm = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, gm, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for gm");
	// Auxiliary field scalars.
	d_ae0 = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, ae0, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for ae0");
	d_ae = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, ae, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for ae");
	d_be = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, be, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for be");
	d_ce = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, ce, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for ce");
	d_de = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, de, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for de");
	d_ee = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, ee, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for ee");

	d_am0 = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, am0, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for ae0");
	d_am = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, am, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for am");
	d_bm = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, bm, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for bm");
	d_cm = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, cm, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for cm");
	d_dm = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, dm, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for dm");
	d_em = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size, em, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for em");

	return 0;
}
int CFDTD1DDNG::InitialiseCLKernelsGPU()
{
	int status;
	/////////////////////////////////////////////////////////////////
	// Load CL file, build CL program object, create CL kernel object
	/////////////////////////////////////////////////////////////////
	const char *filename = "FDTD1DDNG_Kernels.cl";
	string sourceStr = convertToString(filename);
	const char *source = sourceStr.c_str();
	size_t sourceSize[] = {strlen(source)};

	program = clCreateProgramWithSource( context, 1, &source, sourceSize, &status);
	SafeCall(status, "Error: Loading Binary into cl_program (clCreateProgramWithBinary)\n");

	/* create a cl program executable for all the devices specified */
	status = clBuildProgram(program, 1, devices, NULL, NULL, NULL);
	if(status == CL_BUILD_PROGRAM_FAILURE)
	{
		cl_int logStatus;
		char *buildLog = NULL;
		size_t buildLogSize = 0;
		logStatus = clGetProgramBuildInfo (program, devices[0], CL_PROGRAM_BUILD_LOG, buildLogSize, buildLog, &buildLogSize);
		SafeCall(logStatus, "clGetProgramBuildInfo failed.");
		buildLog = new char[buildLogSize];
		SafeCall(!buildLog, "Failed to allocate host memory. (buildLog)");
		memset(buildLog, 0, buildLogSize);
		logStatus = clGetProgramBuildInfo (program, devices[0], CL_PROGRAM_BUILD_LOG, buildLogSize, buildLog, NULL);
		if (logStatus != CL_SUCCESS)
		{
			cout << "clGetProgramBuildInfo failed." << endl;
			free(buildLog);
			return -1;
		}
		cout << " \n\t\t\tBUILD LOG\n";
		cout << " ************************************************\n";
		cout << buildLog << std::endl;
		cout << " ************************************************\n";
		delete []buildLog;
	}

	// Attach kernel objects to respective kernel functions.
	DryRun_kernel_M = clCreateKernel(program, "FDTD1DDNGKernel_DryRun_M", &status);
	SafeCall(status, "Error: Creating dry run Kernel from program. (clCreateKernel)");
	DryRun_kernel_E = clCreateKernel(program, "FDTD1DDNGKernel_DryRun_E", &status);
	SafeCall(status, "Error: Creating dry run Kernel from program. (clCreateKernel)");
	Simulation_kernel_M = clCreateKernel(program, "FDTD1DDNGKernel_Simulation_M", &status);
	SafeCall(status, "Error: Creating simulation Kernel from program. (clCreateKernel)");
	Simulation_kernel_E = clCreateKernel(program, "FDTD1DDNGKernel_Simulation_E", &status);
	SafeCall(status, "Error: Creating simulation Kernel from program. (clCreateKernel)");

	// ====== Set appropriate arguments to the Dry Run kernel ======
	SafeCall(clSetKernelArg(DryRun_kernel_M, 0, sizeof(unsigned int), (void *)&Size), "Error: Setting kernel argument 'Size'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 1, sizeof(unsigned int), (void *)&PulseWidth), "Error: Setting kernel argument 'PulseWidth'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 2, sizeof(unsigned int), (void *)&td), "Error: Setting kernel argument 'td'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 3, sizeof(unsigned int), (void *)&SourceLocation), "Error: Setting kernel argument 'SourceLocation'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 4, sizeof(unsigned int), (void *)&SourceChoice), "Error: Setting kernel argument 'SourceChoice'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 5, sizeof(PRECISION), (void *)&e0), "Error: Setting kernel argument 'e0'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 6, sizeof(PRECISION), (void *)&u0), "Error: Setting kernel argument 'u0'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 7, sizeof(PRECISION), (void *)&dt), "Error: Setting kernel argument 'dt'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 8, sizeof(PRECISION), (void *)&dz), "Error: Setting kernel argument 'dz'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 9, sizeof(PRECISION), (void *)&Sc), "Error: Setting kernel argument 'Sc'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 10, sizeof(PRECISION), (void *)&f), "Error: Setting kernel argument 'f'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 11, sizeof(PRECISION), (void *)&fp), "Error: Setting kernel argument 'fp'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 12, sizeof(PRECISION), (void *)&dr), "Error: Setting kernel argument 'dr'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 13, sizeof(cl_mem), (void *)&d_Ex_), "Error: Setting kernel argument d_Ex_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 14, sizeof(cl_mem), (void *)&d_Hy_), "Error: Setting kernel argument d_Hy_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 15, sizeof(cl_mem), (void *)&d_Exi), "Error: Setting kernel argument d_Exi");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 16, sizeof(unsigned int), (void *)&td), "Error: Setting kernel argument 'x1'");

	SafeCall(clSetKernelArg(DryRun_kernel_E, 0, sizeof(unsigned int), (void *)&Size), "Error: Setting kernel argument 'Size'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 1, sizeof(unsigned int), (void *)&PulseWidth), "Error: Setting kernel argument 'PulseWidth'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 2, sizeof(unsigned int), (void *)&td), "Error: Setting kernel argument 'td'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 3, sizeof(unsigned int), (void *)&SourceLocation), "Error: Setting kernel argument 'SourceLocation'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 4, sizeof(unsigned int), (void *)&SourceChoice), "Error: Setting kernel argument 'SourceChoice'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 5, sizeof(PRECISION), (void *)&e0), "Error: Setting kernel argument 'e0'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 6, sizeof(PRECISION), (void *)&u0), "Error: Setting kernel argument 'u0'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 7, sizeof(PRECISION), (void *)&dt), "Error: Setting kernel argument 'dt'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 8, sizeof(PRECISION), (void *)&dz), "Error: Setting kernel argument 'dz'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 9, sizeof(PRECISION), (void *)&Sc), "Error: Setting kernel argument 'Sc'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 10, sizeof(PRECISION), (void *)&f), "Error: Setting kernel argument 'f'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 11, sizeof(PRECISION), (void *)&fp), "Error: Setting kernel argument 'fp'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 12, sizeof(PRECISION), (void *)&dr), "Error: Setting kernel argument 'dr'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 13, sizeof(cl_mem), (void *)&d_Ex_), "Error: Setting kernel argument d_Ex_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 14, sizeof(cl_mem), (void *)&d_Hy_), "Error: Setting kernel argument d_Hy_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 15, sizeof(cl_mem), (void *)&d_Exi), "Error: Setting kernel argument d_Exi");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 16, sizeof(unsigned int), (void *)&td), "Error: Setting kernel argument 'x1'");

	// ====== Set appropriate arguments to the Simulation kernel ======
	SafeCall(clSetKernelArg(Simulation_kernel_M, 0, sizeof(unsigned int), (void *)&Size), "Error: Setting kernel argument 'Size'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 1, sizeof(unsigned int), (void *)&PulseWidth), "Error: Setting kernel argument 'PulseWidth'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 2, sizeof(unsigned int), (void *)&td), "Error: Setting kernel argument 'td'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 3, sizeof(unsigned int), (void *)&SourceLocation), "Error: Setting kernel argument 'SourceLocation'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 4, sizeof(unsigned int), (void *)&SourceChoice), "Error: Setting kernel argument 'SourceChoice'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 5, sizeof(PRECISION), (void *)&e0), "Error: Setting kernel argument 'e0'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 6, sizeof(PRECISION), (void *)&u0), "Error: Setting kernel argument 'u0'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 7, sizeof(PRECISION), (void *)&dt), "Error: Setting kernel argument 'dt'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 8, sizeof(PRECISION), (void *)&dz), "Error: Setting kernel argument 'dz'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 9, sizeof(PRECISION), (void *)&Sc), "Error: Setting kernel argument 'Sc'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 10, sizeof(PRECISION), (void *)&f), "Error: Setting kernel argument 'f'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 11, sizeof(PRECISION), (void *)&fp), "Error: Setting kernel argument 'fp'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 12, sizeof(PRECISION), (void *)&dr), "Error: Setting kernel argument 'dr'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 13, sizeof(cl_mem), (void *)&d_Ex_), "Error: Setting kernel argument d_Ex_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 14, sizeof(cl_mem), (void *)&d_Dx_), "Error: Setting kernel argument d_Dx_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 15, sizeof(cl_mem), (void *)&d_Hy_), "Error: Setting kernel argument d_Hy_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 16, sizeof(cl_mem), (void *)&d_By_), "Error: Setting kernel argument d_By_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 17, sizeof(cl_mem), (void *)&d_einf), "Error: Setting kernel argument d_einf");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 18, sizeof(cl_mem), (void *)&d_uinf), "Error: Setting kernel argument d_uinf");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 19, sizeof(cl_mem), (void *)&d_wpesq), "Error: Setting kernel argument d_wpesq");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 20, sizeof(cl_mem), (void *)&d_wpmsq), "Error: Setting kernel argument d_wpmsq");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 21, sizeof(cl_mem), (void *)&d_ge), "Error: Setting kernel argument d_ge");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 22, sizeof(cl_mem), (void *)&d_gm), "Error: Setting kernel argument d_gm");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 23, sizeof(cl_mem), (void *)&d_ae0), "Error: Setting kernel argument d_ae0");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 24, sizeof(cl_mem), (void *)&d_ae), "Error: Setting kernel argument d_ae");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 25, sizeof(cl_mem), (void *)&d_be), "Error: Setting kernel argument d_be");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 26, sizeof(cl_mem), (void *)&d_ce), "Error: Setting kernel argument d_ce");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 27, sizeof(cl_mem), (void *)&d_de), "Error: Setting kernel argument d_de");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 28, sizeof(cl_mem), (void *)&d_ee), "Error: Setting kernel argument d_ee");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 29, sizeof(cl_mem), (void *)&d_am0), "Error: Setting kernel argument d_am0");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 30, sizeof(cl_mem), (void *)&d_am), "Error: Setting kernel argument d_am");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 31, sizeof(cl_mem), (void *)&d_bm), "Error: Setting kernel argument d_bm");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 32, sizeof(cl_mem), (void *)&d_cm), "Error: Setting kernel argument d_cm");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 33, sizeof(cl_mem), (void *)&d_dm), "Error: Setting kernel argument d_dm");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 34, sizeof(cl_mem), (void *)&d_em), "Error: Setting kernel argument d_em");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 35, sizeof(cl_mem), (void *)&d_Ext), "Error: Setting kernel argument d_Ext");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 36, sizeof(cl_mem), (void *)&d_Extt), "Error: Setting kernel argument d_Extt");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 37, sizeof(cl_mem), (void *)&d_Exz1), "Error: Setting kernel argument d_Exz1");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 38, sizeof(cl_mem), (void *)&d_Exz2), "Error: Setting kernel argument d_Exz2");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 39, sizeof(unsigned int), (void *)&x1), "Error: Setting kernel argument 'x1'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 40, sizeof(unsigned int), (void *)&Z1), "Error: Setting kernel argument 'Z1'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 41, sizeof(unsigned int), (void *)&Z2), "Error: Setting kernel argument 'Z2'");

	SafeCall(clSetKernelArg(Simulation_kernel_E, 0, sizeof(unsigned int), (void *)&Size), "Error: Setting kernel argument 'Size'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 1, sizeof(unsigned int), (void *)&PulseWidth), "Error: Setting kernel argument 'PulseWidth'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 2, sizeof(unsigned int), (void *)&td), "Error: Setting kernel argument 'td'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 3, sizeof(unsigned int), (void *)&SourceLocation), "Error: Setting kernel argument 'SourceLocation'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 4, sizeof(unsigned int), (void *)&SourceChoice), "Error: Setting kernel argument 'SourceChoice'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 5, sizeof(PRECISION), (void *)&e0), "Error: Setting kernel argument 'e0'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 6, sizeof(PRECISION), (void *)&u0), "Error: Setting kernel argument 'u0'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 7, sizeof(PRECISION), (void *)&dt), "Error: Setting kernel argument 'dt'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 8, sizeof(PRECISION), (void *)&dz), "Error: Setting kernel argument 'dz'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 9, sizeof(PRECISION), (void *)&Sc), "Error: Setting kernel argument 'Sc'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 10, sizeof(PRECISION), (void *)&f), "Error: Setting kernel argument 'f'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 11, sizeof(PRECISION), (void *)&fp), "Error: Setting kernel argument 'fp'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 12, sizeof(PRECISION), (void *)&dr), "Error: Setting kernel argument 'dr'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 13, sizeof(cl_mem), (void *)&d_Ex_), "Error: Setting kernel argument d_Ex_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 14, sizeof(cl_mem), (void *)&d_Dx_), "Error: Setting kernel argument d_Dx_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 15, sizeof(cl_mem), (void *)&d_Hy_), "Error: Setting kernel argument d_Hy_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 16, sizeof(cl_mem), (void *)&d_By_), "Error: Setting kernel argument d_By_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 17, sizeof(cl_mem), (void *)&d_einf), "Error: Setting kernel argument d_einf");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 18, sizeof(cl_mem), (void *)&d_uinf), "Error: Setting kernel argument d_uinf");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 19, sizeof(cl_mem), (void *)&d_wpesq), "Error: Setting kernel argument d_wpesq");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 20, sizeof(cl_mem), (void *)&d_wpmsq), "Error: Setting kernel argument d_wpmsq");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 21, sizeof(cl_mem), (void *)&d_ge), "Error: Setting kernel argument d_ge");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 22, sizeof(cl_mem), (void *)&d_gm), "Error: Setting kernel argument d_gm");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 23, sizeof(cl_mem), (void *)&d_ae0), "Error: Setting kernel argument d_ae0");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 24, sizeof(cl_mem), (void *)&d_ae), "Error: Setting kernel argument d_ae");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 25, sizeof(cl_mem), (void *)&d_be), "Error: Setting kernel argument d_be");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 26, sizeof(cl_mem), (void *)&d_ce), "Error: Setting kernel argument d_ce");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 27, sizeof(cl_mem), (void *)&d_de), "Error: Setting kernel argument d_de");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 28, sizeof(cl_mem), (void *)&d_ee), "Error: Setting kernel argument d_ee");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 29, sizeof(cl_mem), (void *)&d_am0), "Error: Setting kernel argument d_am0");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 30, sizeof(cl_mem), (void *)&d_am), "Error: Setting kernel argument d_am");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 31, sizeof(cl_mem), (void *)&d_bm), "Error: Setting kernel argument d_bm");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 32, sizeof(cl_mem), (void *)&d_cm), "Error: Setting kernel argument d_cm");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 33, sizeof(cl_mem), (void *)&d_dm), "Error: Setting kernel argument d_dm");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 34, sizeof(cl_mem), (void *)&d_em), "Error: Setting kernel argument d_em");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 35, sizeof(cl_mem), (void *)&d_Ext), "Error: Setting kernel argument d_Ext");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 36, sizeof(cl_mem), (void *)&d_Extt), "Error: Setting kernel argument d_Extt");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 37, sizeof(cl_mem), (void *)&d_Exz1), "Error: Setting kernel argument d_Exz1");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 38, sizeof(cl_mem), (void *)&d_Exz2), "Error: Setting kernel argument d_Exz2");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 39, sizeof(unsigned int), (void *)&x1), "Error: Setting kernel argument 'x1'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 40, sizeof(unsigned int), (void *)&Z1), "Error: Setting kernel argument 'Z1'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 41, sizeof(unsigned int), (void *)&Z2), "Error: Setting kernel argument 'Z2'");

	return 0;
}
int CFDTD1DDNG::ReinitialiseExHyGPU()
{
	int status;

	SafeCall(clReleaseMemObject(d_Ex_), "Error: clReleaseMemObject() cannot release memory input buffer for d_Ex_");
	SafeCall(clReleaseMemObject(d_Hy_), "Error: clReleaseMemObject() cannot release memory input buffer for d_Hy_");

	d_Ex_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size*3, Ex_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Ex_");
	d_Hy_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Size*3, Hy_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Hy_");

	SafeCall(clSetKernelArg(Simulation_kernel_M, 13, sizeof(cl_mem), (void *)&d_Ex_), "Error: Setting kernel argument d_Ex_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 15, sizeof(cl_mem), (void *)&d_Hy_), "Error: Setting kernel argument d_Hy_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 13, sizeof(cl_mem), (void *)&d_Ex_), "Error: Setting kernel argument d_Ex_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 15, sizeof(cl_mem), (void *)&d_Hy_), "Error: Setting kernel argument d_Hy_");

	return 0;
}
int CFDTD1DDNG::DryRunCPU()
{
	cout << "Dry run (CPU) started..." << endl;
	for (unsigned int n=0; n<MaxTime; n++)
	{
		if (n%SnapshotInterval == 0)
			cout << "\r" << setprecision(4) << (float)n*100/(MaxTime-1) << "%  " << flush;

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
		if (n%SnapshotInterval == 0)
			cout << "\r" << setprecision(4) << (float)n*100/(MaxTime-1) << "%  " << flush;

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
	cl_int status;
	cl_uint maxDims;
	cl_event events[2];
	size_t globalThreads[2];
	size_t localThreads[2];
	size_t maxWorkGroupSize;
	size_t maxWorkItemSizes[3];


	// Query device capabilities. Maximum work item dimensions and the maximmum work item sizes
	SafeCall(clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), (void*)&maxWorkGroupSize, NULL), "Error: Getting Device Info. (clGetDeviceInfo)");
	SafeCall(clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), (void*)&maxDims, NULL), "Error: Getting Device Info. (clGetDeviceInfo)");
	SafeCall(clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t)*maxDims, (void*)maxWorkItemSizes, NULL), "Error: Getting Device Info. (clGetDeviceInfo)");

	globalThreads[0] = Size;
	globalThreads[1] = 1;
	localThreads[0]  = 256;
	localThreads[1]  = 1;

	std::cout << "Max dimensions: " << maxDims << std::endl;
	std::cout << "Device maxWorkGroupSize = " << maxWorkGroupSize << std::endl;
	std::cout << "Device maxWorkItemSizes = " << maxWorkItemSizes[0] << std::endl;
	if(localThreads[0] > maxWorkGroupSize || localThreads[0] > maxWorkItemSizes[0])
	{
		cout<<"Unsupported: Device does not support requested number of work items." << endl;
		return 1;
	}

	cl_ulong startTime, endTime;
	cl_ulong kernelExecTimeNs;
	cl_ulong kernelExecTimeNsT = 0;

	np = 0;
	n0 = 1;
	nf = 2;

	cout << "Dry run (GPU) started..." << endl;
	cout << "Global threads: " << globalThreads[0] << "x" << globalThreads[1] << endl;
	cout << "Local threads: " << localThreads[0] << "x" << localThreads[1] << endl;

	for (unsigned int n=0;n<MaxTime; n++)
	{
		if (n%SnapshotInterval == 0)
			cout << "\r" << setprecision(4) << (float)n*100/(MaxTime-1) << "%  " << flush;

		SafeCall(clSetKernelArg(DryRun_kernel_M, 17, sizeof(unsigned int), (void *)&n), "Error: Setting kernel argument 'n'");
		SafeCall(clSetKernelArg(DryRun_kernel_M, 18, sizeof(unsigned int), (void *)&np), "Error: Setting kernel argument 'np'");
		SafeCall(clSetKernelArg(DryRun_kernel_M, 19, sizeof(unsigned int), (void *)&n0), "Error: Setting kernel argument 'n0'");
		SafeCall(clSetKernelArg(DryRun_kernel_M, 20, sizeof(unsigned int), (void *)&nf), "Error: Setting kernel argument 'nf'");

		// Enqueue a Dry Run kernel_M call.
		status = clEnqueueNDRangeKernel(commandQueue, DryRun_kernel_M, 2, NULL, globalThreads, localThreads, 0, NULL, &events[0]);
		if(status != CL_SUCCESS) 
		{ 
			cout << "Error: Enqueueing Dry Run kernel onto command queue (clEnqueueNDRangeKernel)" << endl;
			if ( status == CL_INVALID_COMMAND_QUEUE ) std::cout << "CL_INVALID_COMMAND_QUEUE." << endl;
			if ( status == CL_INVALID_PROGRAM_EXECUTABLE ) std::cout << "CL_INVALID_PROGRAM_EXECUTABLE." << endl;
			if ( status == CL_INVALID_KERNEL ) std::cout << "CL_INVALID_KERNEL." << endl;
			if ( status == CL_INVALID_WORK_DIMENSION ) std::cout << "CL_INVALID_WORK_DIMENSION." << endl;
			if ( status == CL_INVALID_CONTEXT ) std::cout << "CL_INVALID_CONTEXT." << endl;
			if ( status == CL_INVALID_KERNEL_ARGS ) std::cout << "CL_INVALID_KERNEL_ARGS." << endl;
			if ( status == CL_INVALID_WORK_GROUP_SIZE ) std::cout << "CL_INVALID_WORK_GROUP_SIZE." << endl;
			if ( status == CL_INVALID_WORK_ITEM_SIZE ) std::cout << "CL_INVALID_WORK_ITEM_SIZE." << endl;
			if ( status == CL_INVALID_GLOBAL_OFFSET ) std::cout << "CL_INVALID_GLOBAL_OFFSET." << endl;
			return 1;
		}

		// Wait for the Dry Run kernel_M call to finish execution.
		SafeCall(clWaitForEvents(1, &events[0]), "Error: Waiting for kernel run to finish. (clWaitForEvents)");

		clGetEventProfilingInfo(events[0], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, NULL);
		clGetEventProfilingInfo(events[0], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, NULL);
		kernelExecTimeNs = (cl_ulong)(1e-3*(endTime-startTime));
		kernelExecTimeNsT = kernelExecTimeNsT + kernelExecTimeNs;

		SafeCall(clSetKernelArg(DryRun_kernel_E, 17, sizeof(unsigned int), (void *)&n), "Error: Setting kernel argument 'n'");
		SafeCall(clSetKernelArg(DryRun_kernel_E, 18, sizeof(unsigned int), (void *)&np), "Error: Setting kernel argument 'np'");
		SafeCall(clSetKernelArg(DryRun_kernel_E, 19, sizeof(unsigned int), (void *)&n0), "Error: Setting kernel argument 'n0'");
		SafeCall(clSetKernelArg(DryRun_kernel_E, 20, sizeof(unsigned int), (void *)&nf), "Error: Setting kernel argument 'nf'");

		// Enqueue a Dry Run kernel_E call.
		status = clEnqueueNDRangeKernel(commandQueue, DryRun_kernel_E, 2, NULL, globalThreads, localThreads, 0, NULL, &events[0]);
		if(status != CL_SUCCESS) 
		{ 
			cout << "Error: Enqueueing Dry Run kernel onto command queue (clEnqueueNDRangeKernel)" << endl;
			if ( status == CL_INVALID_COMMAND_QUEUE ) std::cout << "CL_INVALID_COMMAND_QUEUE." << endl;
			if ( status == CL_INVALID_PROGRAM_EXECUTABLE ) std::cout << "CL_INVALID_PROGRAM_EXECUTABLE." << endl;
			if ( status == CL_INVALID_KERNEL ) std::cout << "CL_INVALID_KERNEL." << endl;
			if ( status == CL_INVALID_WORK_DIMENSION ) std::cout << "CL_INVALID_WORK_DIMENSION." << endl;
			if ( status == CL_INVALID_CONTEXT ) std::cout << "CL_INVALID_CONTEXT." << endl;
			if ( status == CL_INVALID_KERNEL_ARGS ) std::cout << "CL_INVALID_KERNEL_ARGS." << endl;
			if ( status == CL_INVALID_WORK_GROUP_SIZE ) std::cout << "CL_INVALID_WORK_GROUP_SIZE." << endl;
			if ( status == CL_INVALID_WORK_ITEM_SIZE ) std::cout << "CL_INVALID_WORK_ITEM_SIZE." << endl;
			if ( status == CL_INVALID_GLOBAL_OFFSET ) std::cout << "CL_INVALID_GLOBAL_OFFSET." << endl;
			return 1;
		}

		// Wait for the Dry Run kernel_E call to finish execution.
		SafeCall(clWaitForEvents(1, &events[0]), "Error: Waiting for kernel run to finish. (clWaitForEvents)");

		clGetEventProfilingInfo(events[0], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, NULL);
		clGetEventProfilingInfo(events[0], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, NULL);
		kernelExecTimeNs = (cl_ulong)(1e-3*(endTime-startTime));
		kernelExecTimeNsT = kernelExecTimeNsT + kernelExecTimeNs;

		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}
	std::cout << "\r" << "Dry run complete!" << std::endl;
	std::cout << "Dry Run kernel execution time = " << kernelExecTimeNsT/1e6 << "sec (" << kernelExecTimeNsT/1e3 << "ms or " << kernelExecTimeNsT << "us)" << std::endl;
	SafeCall(clReleaseEvent(events[0]), "Error: Release event object. (clReleaseEvent)\n");

	return 0;
}
int CFDTD1DDNG::RunSimulationGPU(bool SaveFields)
{
	cl_int status;
	cl_uint maxDims;
	cl_event events[2];
	size_t globalThreads[2];
	size_t localThreads[2];
	size_t maxWorkGroupSize;
	size_t maxWorkItemSizes[3];


	// Query device capabilities. Maximum work item dimensions and the maximmum work item sizes
	SafeCall(clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), (void*)&maxWorkGroupSize, NULL), "Error: Getting Device Info. (clGetDeviceInfo)");
	SafeCall(clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), (void*)&maxDims, NULL), "Error: Getting Device Info. (clGetDeviceInfo)");
	SafeCall(clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t)*maxDims, (void*)maxWorkItemSizes, NULL), "Error: Getting Device Info. (clGetDeviceInfo)");

	globalThreads[0] = Size;
	globalThreads[1] = 1;
	localThreads[0]  = 256;
	localThreads[1]  = 1;

	std::cout << "Max dimensions: " << maxDims << std::endl;
	std::cout << "Device maxWorkGroupSize = " << maxWorkGroupSize << std::endl;
	std::cout << "Device maxWorkItemSizes = " << maxWorkItemSizes[0] << std::endl;
	if(localThreads[0] > maxWorkGroupSize || localThreads[0] > maxWorkItemSizes[0])
	{
		cout<<"Unsupported: Device does not support requested number of work items." << endl;
		return 1;
	}

	cl_ulong startTime, endTime;
	cl_ulong kernelExecTimeNs;
	cl_ulong kernelExecTimeNsT = 0;

	np = 0;
	n0 = 1;
	nf = 2;

	stringstream framestream;
	string basename = "FieldData/Ex";
	string filename;
	fstream snapshot;
	frame = 0U;

	cout << "Simulation run (GPU) started..." << endl;
	cout << "Global threads: " << globalThreads[0] << "x" << globalThreads[1] << endl;
	cout << "Local threads: " << localThreads[0] << "x" << localThreads[1] << endl;

	for (unsigned int n=0;n<MaxTime; n++)
	{
		if (n%SnapshotInterval == 0)
			cout << "\r" << setprecision(4) << (float)n*100/(MaxTime-1) << "%  " << flush;

		SafeCall(clSetKernelArg(Simulation_kernel_M, 42, sizeof(unsigned int), (void *)&n), "Error: Setting kernel argument 'n'");
		SafeCall(clSetKernelArg(Simulation_kernel_M, 43, sizeof(unsigned int), (void *)&np), "Error: Setting kernel argument 'np'");
		SafeCall(clSetKernelArg(Simulation_kernel_M, 44, sizeof(unsigned int), (void *)&n0), "Error: Setting kernel argument 'n0'");
		SafeCall(clSetKernelArg(Simulation_kernel_M, 45, sizeof(unsigned int), (void *)&nf), "Error: Setting kernel argument 'nf'");

		// Enqueue a simulation kernel_M call.
		status = clEnqueueNDRangeKernel(commandQueue, Simulation_kernel_M, 2, NULL, globalThreads, localThreads, 0, NULL, &events[0]);
		if(status != CL_SUCCESS) 
		{ 
			cout << "Error: Enqueueing Dry Run kernel onto command queue (clEnqueueNDRangeKernel)" << endl;
			if ( status == CL_INVALID_COMMAND_QUEUE ) std::cout << "CL_INVALID_COMMAND_QUEUE." << endl;
			if ( status == CL_INVALID_PROGRAM_EXECUTABLE ) std::cout << "CL_INVALID_PROGRAM_EXECUTABLE." << endl;
			if ( status == CL_INVALID_KERNEL ) std::cout << "CL_INVALID_KERNEL." << endl;
			if ( status == CL_INVALID_WORK_DIMENSION ) std::cout << "CL_INVALID_WORK_DIMENSION." << endl;
			if ( status == CL_INVALID_CONTEXT ) std::cout << "CL_INVALID_CONTEXT." << endl;
			if ( status == CL_INVALID_KERNEL_ARGS ) std::cout << "CL_INVALID_KERNEL_ARGS." << endl;
			if ( status == CL_INVALID_WORK_GROUP_SIZE ) std::cout << "CL_INVALID_WORK_GROUP_SIZE." << endl;
			if ( status == CL_INVALID_WORK_ITEM_SIZE ) std::cout << "CL_INVALID_WORK_ITEM_SIZE." << endl;
			if ( status == CL_INVALID_GLOBAL_OFFSET ) std::cout << "CL_INVALID_GLOBAL_OFFSET." << endl;
			return 1;
		}

		// Wait for the simulation kernel_M call to finish execution.
		SafeCall(clWaitForEvents(1, &events[0]), "Error: Waiting for kernel run to finish. (clWaitForEvents)");

		clGetEventProfilingInfo(events[0], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, NULL);
		clGetEventProfilingInfo(events[0], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, NULL);
		kernelExecTimeNs = (cl_ulong)(1e-3*(endTime-startTime));
		kernelExecTimeNsT = kernelExecTimeNsT + kernelExecTimeNs;

		SafeCall(clSetKernelArg(Simulation_kernel_E, 42, sizeof(unsigned int), (void *)&n), "Error: Setting kernel argument 'n'");
		SafeCall(clSetKernelArg(Simulation_kernel_E, 43, sizeof(unsigned int), (void *)&np), "Error: Setting kernel argument 'np'");
		SafeCall(clSetKernelArg(Simulation_kernel_E, 44, sizeof(unsigned int), (void *)&n0), "Error: Setting kernel argument 'n0'");
		SafeCall(clSetKernelArg(Simulation_kernel_E, 45, sizeof(unsigned int), (void *)&nf), "Error: Setting kernel argument 'nf'");

		// Enqueue a simulation kernel_E call.
		status = clEnqueueNDRangeKernel(commandQueue, Simulation_kernel_E, 2, NULL, globalThreads, localThreads, 0, NULL, &events[0]);
		if(status != CL_SUCCESS) 
		{ 
			cout << "Error: Enqueueing Dry Run kernel onto command queue (clEnqueueNDRangeKernel)" << endl;
			if ( status == CL_INVALID_COMMAND_QUEUE ) std::cout << "CL_INVALID_COMMAND_QUEUE." << endl;
			if ( status == CL_INVALID_PROGRAM_EXECUTABLE ) std::cout << "CL_INVALID_PROGRAM_EXECUTABLE." << endl;
			if ( status == CL_INVALID_KERNEL ) std::cout << "CL_INVALID_KERNEL." << endl;
			if ( status == CL_INVALID_WORK_DIMENSION ) std::cout << "CL_INVALID_WORK_DIMENSION." << endl;
			if ( status == CL_INVALID_CONTEXT ) std::cout << "CL_INVALID_CONTEXT." << endl;
			if ( status == CL_INVALID_KERNEL_ARGS ) std::cout << "CL_INVALID_KERNEL_ARGS." << endl;
			if ( status == CL_INVALID_WORK_GROUP_SIZE ) std::cout << "CL_INVALID_WORK_GROUP_SIZE." << endl;
			if ( status == CL_INVALID_WORK_ITEM_SIZE ) std::cout << "CL_INVALID_WORK_ITEM_SIZE." << endl;
			if ( status == CL_INVALID_GLOBAL_OFFSET ) std::cout << "CL_INVALID_GLOBAL_OFFSET." << endl;
			return 1;
		}

		// Wait for the simulation kernel_E call to finish execution.
		SafeCall(clWaitForEvents(1, &events[0]), "Error: Waiting for kernel run to finish. (clWaitForEvents)");

		clGetEventProfilingInfo(events[0], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, NULL);
		clGetEventProfilingInfo(events[0], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, NULL);
		kernelExecTimeNs = (cl_ulong)(1e-3*(endTime-startTime));
		kernelExecTimeNsT = kernelExecTimeNsT + kernelExecTimeNs;

		// Saving electric field snapshot.
		if (n%SnapshotInterval == 0 && SaveFields == true)
		{
			// Write E-field to file.
			framestream.str(std::string());			// Clearing stringstream contents.
			framestream << ++frame;
			filename = basename + framestream.str() + ".fdt";
			snapshot.open(filename.c_str(), std::ios::out|std::ios::binary);

			// Enqueue read buffer.
			SafeCall(clEnqueueReadBuffer(commandQueue, d_Ex_, CL_TRUE, 0, sizeof(PRECISION)*Size*3, Ex_, 0, NULL, &events[1]), "Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)");
			// Wait for the read buffer to finish execution
			SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for read buffer call to finish. (clWaitForEvents)");

			snapshot.write((char*)&(Ex(0,nf)), sizeof(PRECISION)*Size);
			snapshot.close();
		}

		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}
	std::cout << "\r" << "Simulation complete!" << std::endl;
	std::cout << "Simulation kernel execution time = " << kernelExecTimeNsT/1e6 << "sec (" << kernelExecTimeNsT/1e3 << "ms or " << kernelExecTimeNsT << "us)" << std::endl;
	SafeCall(clReleaseEvent(events[0]), "Error: Release event object. (clReleaseEvent)\n");

	// Saving electric field data arrays.
	if (SaveFields == true)
	{
		fstream parametersfile;
		parametersfile.open("FieldData/Parameters.smp", std::ios::out|std::ios::binary|std::ios::app);
		parametersfile.write((char*)&(frame), sizeof(unsigned int));
		parametersfile.close();
		// Write saved fields to files.
		snapshot.open("FieldData/Exi.fdt", std::ios::out|std::ios::binary);
		SafeCall(clEnqueueReadBuffer(commandQueue, d_Exi, CL_TRUE, 0, sizeof(PRECISION)*Size, Exi, 0, NULL, &events[1]), "Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)");
		SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for read buffer call to finish. (clWaitForEvents)");
		snapshot.write((char*)Exi, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ext.fdt", std::ios::out|std::ios::binary);
		SafeCall(clEnqueueReadBuffer(commandQueue, d_Ext, CL_TRUE, 0, sizeof(PRECISION)*Size, Ext, 0, NULL, &events[1]), "Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)");
		SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for read buffer call to finish. (clWaitForEvents)");
		snapshot.write((char*)Ext, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Extt.fdt", std::ios::out|std::ios::binary);
		SafeCall(clEnqueueReadBuffer(commandQueue, d_Extt, CL_TRUE, 0, sizeof(PRECISION)*Size, Extt, 0, NULL, &events[1]), "Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)");
		SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for read buffer call to finish. (clWaitForEvents)");
		snapshot.write((char*)Extt, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Exz1.fdt", std::ios::out|std::ios::binary);
		SafeCall(clEnqueueReadBuffer(commandQueue, d_Exz1, CL_TRUE, 0, sizeof(PRECISION)*Size, Exz1, 0, NULL, &events[1]), "Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)");
		SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for read buffer call to finish. (clWaitForEvents)");
		snapshot.write((char*)Exz1, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Exz2.fdt", std::ios::out|std::ios::binary);
		SafeCall(clEnqueueReadBuffer(commandQueue, d_Exz2, CL_TRUE, 0, sizeof(PRECISION)*Size, Exz2, 0, NULL, &events[1]), "Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)");
		SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for read buffer call to finish. (clWaitForEvents)");
		snapshot.write((char*)Exz2, sizeof(PRECISION)*MaxTime);
		snapshot.close();

		SafeCall(clReleaseEvent(events[1]), "Error: Release event object. (clReleaseEvent)\n");
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

	SafeCall(InitialiseCL(), "Error: Initialiasing CL.");
	StartTimer();
	SafeCall(AllocateMemoryGPU(), "Error: Allocating memory on GPU.");
	SafeCall(InitialiseCLKernelsGPU(), "Error: Copying data from CPU to GPU.");
	SafeCall(DryRunGPU(), "Error: Dry run (GPU).");
	SafeCall(InitialiseExHyCPU(), "Error: Initalising Ex/Hy arrays (CPU).");
	SafeCall(ReinitialiseExHyGPU(), "Error: Reinitialising Ex/Hy arrays on GPU.");
	SafeCall(RunSimulationGPU(SaveFields), "Error: Running simulation on GPU.");
	SafeCall(CleanupCPU(), "Error: Cleaning up CPU.");
	SafeCall(CleanupCL(), "Error: Cleaning up CL.");
	SafeCall(CleanupGPU(), "Error: Cleaning up GPU.");

	return 0;
}
// Converts contents of a file into a string. From OPENCL examples.
string CFDTD1DDNG::convertToString(const char *filename)
{
	size_t size;
	char*  str;
	string s;
	fstream f(filename, (fstream::in | fstream::binary));

	if(f.is_open())
	{
		size_t fileSize;
		f.seekg(0, fstream::end);
		size = fileSize = (size_t)f.tellg();
		f.seekg(0, fstream::beg);

		str = new char[size+1];
		if(!str)
		{
			f.close();
			return NULL;
		}

		f.read(str, fileSize);
		f.close();
		str[size] = '\0';

		s = str;
		delete[] str;
		return s;
	}
	else
	{
		cout << "\nFile containg the kernel code(\".cl\") not found. Please copy the required file in the folder containg the executable." << endl;
		exit(1);
	}
	return NULL;
}
// Timing.
void CFDTD1DDNG::StartTimer()
{
	if (tPaused == true)
	{
		tStart = GetTimeus64();
		tPaused = false;
	}
}
void CFDTD1DDNG::StopTimer()
{
	if (tPaused == false)
	{
		tEnd = GetTimeus64();
		tDelta += tEnd - tStart;
		tStart = tEnd;
		tPaused = true;
	}
}
void CFDTD1DDNG::ResetTimer()
{
	if (tPaused == true)
		tStart = tEnd;
	else
		tStart = GetTimeus64();

	tDelta = 0UL;
}
double CFDTD1DDNG::GetElapsedTime()
{
	if (tPaused == false)
		tEnd = GetTimeus64();

	return ((double)(tEnd-tStart+tDelta))/(1000000.);
}
int CFDTD1DDNG::SafeCall(cl_int Status, const char *Error)
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
int CFDTD1DDNG::CleanupCL()
{
	SafeCall(clReleaseKernel(DryRun_kernel_M), "Error: In clReleaseKernel");
	SafeCall(clReleaseKernel(DryRun_kernel_E), "Error: In clReleaseKernel");
	SafeCall(clReleaseKernel(Simulation_kernel_M), "Error: In clReleaseKernel");
	SafeCall(clReleaseKernel(Simulation_kernel_E), "Error: In clReleaseKernel");
	SafeCall(clReleaseProgram(program), "Error: In clReleaseProgram");
	SafeCall(clReleaseCommandQueue(commandQueue), "Error: In clReleaseCommandQueue");
	SafeCall(clReleaseContext(context), "Error: In clReleaseContext");

	return 0;
}
int CFDTD1DDNG::CleanupGPU()
{
	SafeCall(clReleaseMemObject(d_Ex_), "Error: clReleaseMemObject() cannot release memory buffer for d_Ex_");
	SafeCall(clReleaseMemObject(d_Dx_), "Error: clReleaseMemObject() cannot release memory buffer for d_Dx_");
	SafeCall(clReleaseMemObject(d_Hy_), "Error: clReleaseMemObject() cannot release memory buffer for d_Hy_");
	SafeCall(clReleaseMemObject(d_By_), "Error: clReleaseMemObject() cannot release memory buffer for d_By_");

	SafeCall(clReleaseMemObject(d_Exi), "Error: clReleaseMemObject() cannot release memory buffer for d_Exi");
	SafeCall(clReleaseMemObject(d_Ext), "Error: clReleaseMemObject() cannot release memory buffer for d_Ext");
	SafeCall(clReleaseMemObject(d_Extt), "Error: clReleaseMemObject() cannot release memory buffer for d_Extt");
	SafeCall(clReleaseMemObject(d_Exz1), "Error: clReleaseMemObject() cannot release memory buffer for d_Exz1");
	SafeCall(clReleaseMemObject(d_Exz2), "Error: clReleaseMemObject() cannot release memory buffer for d_Exz2");

	SafeCall(clReleaseMemObject(d_einf), "Error: clReleaseMemObject() cannot release memory buffer for d_einf");
	SafeCall(clReleaseMemObject(d_uinf), "Error: clReleaseMemObject() cannot release memory buffer for d_uinf");
	SafeCall(clReleaseMemObject(d_wpesq), "Error: clReleaseMemObject() cannot release memory buffer for d_wpesq");
	SafeCall(clReleaseMemObject(d_wpmsq), "Error: clReleaseMemObject() cannot release memory buffer for d_wpmsq");
	SafeCall(clReleaseMemObject(d_ge), "Error: clReleaseMemObject() cannot release memory buffer for d_ge");
	SafeCall(clReleaseMemObject(d_gm), "Error: clReleaseMemObject() cannot release memory buffer for d_gm");

	SafeCall(clReleaseMemObject(d_ae0), "Error: clReleaseMemObject() cannot release memory buffer for d_ae0");
	SafeCall(clReleaseMemObject(d_ae), "Error: clReleaseMemObject() cannot release memory buffer for d_ae");
	SafeCall(clReleaseMemObject(d_be), "Error: clReleaseMemObject() cannot release memory buffer for d_be");
	SafeCall(clReleaseMemObject(d_ce), "Error: clReleaseMemObject() cannot release memory buffer for d_ce");
	SafeCall(clReleaseMemObject(d_de), "Error: clReleaseMemObject() cannot release memory buffer for d_de");
	SafeCall(clReleaseMemObject(d_ee), "Error: clReleaseMemObject() cannot release memory buffer for d_ee");
	SafeCall(clReleaseMemObject(d_am0), "Error: clReleaseMemObject() cannot release memory buffer for d_am0");
	SafeCall(clReleaseMemObject(d_am), "Error: clReleaseMemObject() cannot release memory buffer for d_am");
	SafeCall(clReleaseMemObject(d_bm), "Error: clReleaseMemObject() cannot release memory buffer for d_bm");
	SafeCall(clReleaseMemObject(d_cm), "Error: clReleaseMemObject() cannot release memory buffer for d_cm");
	SafeCall(clReleaseMemObject(d_dm), "Error: clReleaseMemObject() cannot release memory buffer for d_dm");
	SafeCall(clReleaseMemObject(d_em), "Error: clReleaseMemObject() cannot release memory buffer for d_em");

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
