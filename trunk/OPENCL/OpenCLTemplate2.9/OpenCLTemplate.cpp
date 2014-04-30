#include <OpenCLTemplate.hpp>
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

COpenCLTemplate::COpenCLTemplate(unsigned int pWidth, PRECISION pMultiplier): Width(pWidth), Multiplier(pMultiplier), input(NULL), output(NULL), tStart(0LL), tEnd(0LL), tDelta(0LL), tPaused(true)
{
	// Printing simulation size.
	cout << "Width      = " << Width << endl;
	cout << "Multiplier = " << Multiplier << endl;
}
// Allocate memory for data arrays.
int COpenCLTemplate::AllocateMemoryCPU()
{
	input = new PRECISION[Width];
	output = new PRECISION[Width];

	return 0;
}
// Initialise CPU data.
int COpenCLTemplate::InitialiseCPU()
{
	for (unsigned int i=0; i<Width; i++)
		input[i] = (PRECISION)(rand()%100);

	return 0;
}
int COpenCLTemplate::InitialiseCL()
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
		cout << "Selected platform not found so Exiting Application." << endl;
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
int COpenCLTemplate::AllocateMemoryGPU()
{
	cl_int status;

	d_input = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Width, input, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create input buffer");

	d_output = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*Width, output, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create output buffer");

	return 0;
}
int COpenCLTemplate::InitialiseCLKernelsGPU()
{
	int status;
	/////////////////////////////////////////////////////////////////
	// Load CL file, build CL program object, create CL kernel object
	/////////////////////////////////////////////////////////////////
	const char *filename = "OpenCLTemplate_Kernels.cl";
	string sourceStr = convertToString(filename);
	const char *source = sourceStr.c_str();
	size_t sourceSize[] = {strlen(source)};

	program = clCreateProgramWithSource( context, 1, &source, sourceSize, &status);
	SafeCall(status, "Error: Loading Binary into cl_program (clCreateProgramWithBinary)\n");

	// Create a cl program executable for all the devices specified.
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
		// Displaying build log in case of errors.
		cout << " \n\t\t\tBUILD LOG\n";
		cout << " ************************************************\n";
		cout << buildLog << endl;
		cout << " ************************************************\n";
		delete []buildLog;
	}

	// Attach kernel objects to respective kernel functions.
	kernel = clCreateKernel(program, "OpenCLTemplateKernel", &status);
	SafeCall(status, "Error: Creating Kernel from program. (clCreateKernel)");

	// ====== Set appropriate arguments to the kernel ======
	SafeCall(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&d_input), "Error: Setting kernel argument 'input'");
	SafeCall(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&d_output), "Error: Setting kernel argument 'output'");
	SafeCall(clSetKernelArg(kernel, 2, sizeof(PRECISION), (void*)&Multiplier), "Error: Setting kernel argument 'Multiplier'");

	return 0;
}
int COpenCLTemplate::RunCLKernels()
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

	globalThreads[0] = Width;
	globalThreads[1] = 1;
	localThreads[0]  = 256;
	localThreads[1]  = 1;

	cout << "Max dimensions: " << maxDims << endl;
	cout << "Device maxWorkGroupSize = " << maxWorkGroupSize << endl;
	cout << "Device maxWorkItemSizes = " << maxWorkItemSizes[0] << endl;
	if(localThreads[0] > maxWorkGroupSize || localThreads[0] > maxWorkItemSizes[0])
	{
		cout<<"Unsupported: Device does not support requested number of work items." << endl;
		return 1;
	}

	// Kernel timing variables.
	cl_ulong startTime, endTime;
	cl_ulong kernelExecTimeNs;
	cl_ulong kernelExecTimeNsT = 0;

	cout << "Launching CL Kernel..." << endl;
	cout << "Global threads: " << globalThreads[0] << "x" << globalThreads[1] << endl;
	cout << "Local threads: " << localThreads[0] << "x" << localThreads[1] << endl;

	// Enqueue a kernel call.
	status = clEnqueueNDRangeKernel(commandQueue, kernel, 2, NULL, globalThreads, localThreads, 0, NULL, &events[0]);
	if(status != CL_SUCCESS) 
	{ 
		cout << "Error: Enqueueing kernel onto command queue (clEnqueueNDRangeKernel)" << endl;
		if ( status == CL_INVALID_COMMAND_QUEUE ) cout << "CL_INVALID_COMMAND_QUEUE." << endl;
		if ( status == CL_INVALID_PROGRAM_EXECUTABLE ) cout << "CL_INVALID_PROGRAM_EXECUTABLE." << endl;
		if ( status == CL_INVALID_KERNEL ) cout << "CL_INVALID_KERNEL." << endl;
		if ( status == CL_INVALID_WORK_DIMENSION ) cout << "CL_INVALID_WORK_DIMENSION." << endl;
		if ( status == CL_INVALID_CONTEXT ) cout << "CL_INVALID_CONTEXT." << endl;
		if ( status == CL_INVALID_KERNEL_ARGS ) cout << "CL_INVALID_KERNEL_ARGS." << endl;
		if ( status == CL_INVALID_WORK_GROUP_SIZE ) cout << "CL_INVALID_WORK_GROUP_SIZE." << endl;
		if ( status == CL_INVALID_WORK_ITEM_SIZE ) cout << "CL_INVALID_WORK_ITEM_SIZE." << endl;
		if ( status == CL_INVALID_GLOBAL_OFFSET ) cout << "CL_INVALID_GLOBAL_OFFSET." << endl;
		return 1;
	}

	// Wait for the kernel call to finish execution.
	SafeCall(clWaitForEvents(1, &events[0]), "Error: Waiting for kernel run to finish. (clWaitForEvents)");

	clGetEventProfilingInfo(events[0], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, NULL);
	clGetEventProfilingInfo(events[0], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, NULL);
	kernelExecTimeNs = (cl_ulong)(1e-3*(endTime-startTime));
	kernelExecTimeNsT = kernelExecTimeNsT + kernelExecTimeNs;

	cout << "Kernel run complete!" << endl;
	cout << "Kernel execution time = " << kernelExecTimeNsT/1e6 << "sec (" << kernelExecTimeNsT/1e3 << "ms or " << kernelExecTimeNsT << "us)" << endl;
	SafeCall(clReleaseEvent(events[0]), "Error: Release event object. (clReleaseEvent)\n");

	// Copying output array data back to CPU.
	SafeCall(status = clEnqueueReadBuffer(commandQueue, /*Source: GPU*/d_output, CL_TRUE, 0,  /*Size of data in bytes*/sizeof(PRECISION)*Width, /*Destination: CPU*/output, 0, NULL, &events[1]),"Error: clEnqueueReadBuffer Failed for output array.");
	SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for output array to be copied back to CPU. (clWaitForEvents)");
	SafeCall(clReleaseEvent(events[1]), "Error: Release event object. (clReleaseEvent)\n");

	// Display results.
	cout << "Multiplier is " << Multiplier << endl;
	cout << "Input array is: " << endl;
	SafeCall(DisplayArray(Width, input), "Error: Displaying input array.");
	cout << "Output array is: " << endl;
	SafeCall(DisplayArray(Width, output), "Error: Displaying output array.");

	return 0;
}
int COpenCLTemplate::CompleteRun()
{
	SafeCall(AllocateMemoryCPU(), "Error: Allocating memory on CPU.");
	SafeCall(InitialiseCPU(), "Error: Initialising data on CPU.");
	SafeCall(InitialiseCL(), "Error: Initialiasing CL.");
	SafeCall(AllocateMemoryGPU(), "Error: Allocating memory on GPU.");
	SafeCall(InitialiseCLKernelsGPU(), "Error: Copying data from CPU to GPU.");
	SafeCall(RunCLKernels(), "Error: Running kernels (GPU).");
	SafeCall(CleanupCPU(), "Error: Cleaning up CPU.");
	SafeCall(CleanupCL(), "Error: Cleaning up CL.");
	SafeCall(CleanupGPU(), "Error: Cleaning up GPU.");

	return 0;
}
// Converts contents of a file into a string. From OPENCL examples.
string COpenCLTemplate::convertToString(const char *filename)
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

// Display array.
int COpenCLTemplate::DisplayArray(const unsigned int Size, PRECISION* Array)
{
	for (unsigned int i=0; i<Size; i++)
		cout << Array[i] << " ";
	cout << endl;

	return 0;
}

// Timing.
void COpenCLTemplate::StartTimer()
{
	if (tPaused == true)
	{
		tStart = GetTimeus64();
		tPaused = false;
	}
}
void COpenCLTemplate::StopTimer()
{
	if (tPaused == false)
	{
		tEnd = GetTimeus64();
		tDelta += tEnd - tStart;
		tStart = tEnd;
		tPaused = true;
	}
}
void COpenCLTemplate::ResetTimer()
{
	if (tPaused == true)
		tStart = tEnd;
	else
		tStart = GetTimeus64();

	tDelta = 0UL;
}
double COpenCLTemplate::GetElapsedTime()
{
	if (tPaused == false)
		tEnd = GetTimeus64();

	return ((double)(tEnd-tStart+tDelta))/(1000000.);
}
int COpenCLTemplate::SafeCall(cl_int Status, const char *Error)
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
int COpenCLTemplate::CleanupCPU()
{
	DeleteArray(input);
	DeleteArray(output);

	return 0;
}
int COpenCLTemplate::CleanupCL()
{
	SafeCall(clReleaseKernel(kernel), "Error: In clReleaseKernel");
	SafeCall(clReleaseProgram(program), "Error: In clReleaseProgram");
	SafeCall(clReleaseCommandQueue(commandQueue), "Error: In clReleaseCommandQueue");
	SafeCall(clReleaseContext(context), "Error: In clReleaseContext");

	return 0;
}
int COpenCLTemplate::CleanupGPU()
{
	SafeCall(clReleaseMemObject(d_input), "Error: clReleaseMemObject() cannot release input memory buffer");
	SafeCall(clReleaseMemObject(d_output), "Error: clReleaseMemObject() cannot release output memory buffer");

	return 0;
}
COpenCLTemplate::~COpenCLTemplate ()
{
	// Field arrays.
	DeleteArray(input);
	DeleteArray(output);
}
