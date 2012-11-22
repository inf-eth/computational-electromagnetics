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
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
using namespace std;

CFDTD2DDNG::CFDTD2DDNG(
						unsigned int pI,
						unsigned int pJ,
						unsigned int pPMLw,
						unsigned int MaxTime,
						unsigned int pSnapshotResolution,
						unsigned int pSnapshotInterval,
						unsigned int pSourceChoice,
						unsigned int pSourcePlane,
						unsigned int pSourceLocationX,
						unsigned int pSourceLocationY):

							// Simulation parameters.
							I(pI),
							J(pJ),
							PMLw(pPMLw),
							SlabLeft(J/3+PMLw),
							SlabRight(2*J/3+PMLw),
							MaxTime(4*J),
							PulseWidth(J/8),
							td(PulseWidth),
							SnapshotResolution(pSnapshotResolution),
							SnapshotInterval(pSnapshotInterval),
							SourceChoice(pSourceChoice),
							SourcePlane(pSourcePlane),
							SourceLocationX(pSourceLocationX),
							SourceLocationY(pSourceLocationY),
							// Constants.
							c(299792458.),
							pi(3.14159265358979323846),
							e0((1e-9)/(36.*pi)),
							u0((1e-7)*4.*pi),
							dt(0.5e-11),
							delta(3e-3),
							Sc(c*dt/delta),
							// Frequency, wavelength, wave number.
							l(PulseWidth*delta),
							f(c/l),
							fmax(1/(2*dt)),
							w(2*pi*f),
							k0(w/c),
							fp(f),
							dr(PulseWidth*dt*2),
							// Data array sizes.
							IEz(I), JEz(J+2*PMLw),
							IHx(I), JHx(J+2*PMLw+1),
							IHy(I), JHy(J+2*PMLw),
							// Data arrays.
							Ez_(NULL), Dz_(NULL), Hx_(NULL), Bx_(NULL), Hy_(NULL), By_(NULL),
							// Incident and transmitted fields.
							Ezi(NULL), Ezt(NULL), Eztt(NULL),
							x1(SlabLeft),
							// Refractive index.
							Y1(SlabLeft+5),
							y1((PRECISION)Y1*delta),
							Y2(SlabLeft+6),
							y2((PRECISION)Y2*delta),
							Ezy1(NULL),
							Ezy2(NULL),
							// Drude material parameters.
							einf_(NULL), uinf_(NULL), wpesq_(NULL), wpmsq_(NULL), ge_(NULL), gm_(NULL),
							// Auxiliary field scalars.
							ae0_(NULL), ae_(NULL), be_(NULL), ce_(NULL), de_(NULL), ee_(NULL),
							am0_(NULL), am_(NULL), bm_(NULL), cm_(NULL), dm_(NULL), em_(NULL),
							// PML arrays.
							PsiEzX_(NULL), PsiEzY_(NULL),
							PsiHyX_(NULL), PsiHxY_(NULL),
							// PML parameters.
							kappe(1.), kappm(1.),
							kappex(1.), kappey(kappex), kappmx(kappex), kappmy(kappex),
							aex(0.0004), aey(aex), amx(aex), amy(aex),
							sigex(0.), sigey(0.045), sigmx(0.), sigmy(u0/e0*sigey),
							bex(exp(-1.*(aex/e0+sigex/(kappex*e0))*dt)), bey(exp(-1.*(aey/e0+sigey/(kappey*e0))*dt)), bmx(exp(-1.*(amx/u0+sigmx/(kappmx*u0))*dt)), bmy(exp(-1.*(amy/u0+sigmy/(kappmy*u0))*dt)),
							Cex((bex-1.)*sigex/(sigex*kappex+pow(kappe, 2)*aex)), Cey((bey-1.)*sigey/(sigey*kappey+pow(kappe, 2)*aey)), Cmx((bmx-1.)*sigmx/(sigmx*kappmx+pow(kappm, 2)*amx)), Cmy((bmy-1.)*sigmy/(sigmy*kappmy+pow(kappm, 2)*amy)),
							// Snapshot frame number.
							frame(0),
							// Time indices.
							nf(2), n0(1), np(0),
							// Timer variables.
							tStart(0LL), tEnd(0LL),
							cpu(true)
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
	parametersfile.open("FieldData/Parameters.smp", ios::out|ios::binary);
	parametersfile.write((char*)&dt, sizeof(PRECISION));
	parametersfile.write((char*)&k0, sizeof(PRECISION));
	parametersfile.write((char*)&y1, sizeof(PRECISION));
	parametersfile.write((char*)&y2, sizeof(PRECISION));
	parametersfile.write((char*)&I, sizeof(unsigned int));
	parametersfile.write((char*)&J, sizeof(unsigned int));
	parametersfile.write((char*)&PMLw, sizeof(unsigned int));
	parametersfile.write((char*)&MaxTime, sizeof(unsigned int));
	parametersfile.write((char*)&SnapshotResolution, sizeof(unsigned int));
	parametersfile.write((char*)&SnapshotInterval, sizeof(unsigned int));
	parametersfile.write((char*)&SlabLeft, sizeof(unsigned int));
	parametersfile.write((char*)&SlabRight, sizeof(unsigned int));
	parametersfile.close();

	// Printing simulation parameters.
	cout << "I      = " << I << endl;
	cout << "J      = " << J << endl;
	cout << "PMLw       = " << PMLw << endl;
	cout << "MaxTime    = " << MaxTime << endl;
	cout << "frequency  = " << f << " Hz (" << f/1e9 << " GHz)" << endl;
	cout << "fmax       = " << fmax << " Hz (" << fmax/1e9 << " GHz)" << endl;
	cout << "Sc         = " << Sc << endl;
	cout << "Slab left  = " << SlabLeft << endl;
	cout << "Slab right = " << SlabRight << endl;
}
unsigned long CFDTD2DDNG::SimSize()
{
	return (unsigned long)sizeof(*this)+(unsigned long)sizeof(PRECISION)*(8UL*(unsigned long)IEz*(unsigned long)JEz+7UL*(unsigned long)IHy*(unsigned long)JHy+25UL*(unsigned long)IHx*(unsigned long)JHx+5UL*(unsigned long)MaxTime);
}
unsigned long CFDTD2DDNG::HDDSpace()
{
	return (unsigned long)sizeof(PRECISION)*(5UL*(unsigned long)MaxTime+(unsigned long)IEz*(unsigned long)JEz*((unsigned long)MaxTime/(unsigned long)SnapshotInterval+1UL));
}
// Allocate memory for data arrays.
int CFDTD2DDNG::AllocateMemoryCPU()
{
	// Field arrays.
	Ez_ = new PRECISION[IEz*JEz*3];
	Dz_ = new PRECISION[IEz*JEz*3];
	Hx_ = new PRECISION[IHx*JHx*3];
	Bx_ = new PRECISION[IHx*JHx*3];
	Hy_ = new PRECISION[IHy*JHy*3];
	By_ = new PRECISION[IHy*JHy*3];
	// Incident and transmitted fields.
	Ezi = new PRECISION[MaxTime];
	Ezt = new PRECISION[MaxTime];
	Eztt = new PRECISION[MaxTime];
	// Refractive index.
	Ezy1 = new PRECISION[MaxTime];
	Ezy2 = new PRECISION[MaxTime];
	// Drude parameter arrays.
	einf_ = new PRECISION[IHx*JHx];
	uinf_ = new PRECISION[IHx*JHx];
	wpesq_ = new PRECISION[IHx*JHx];
	wpmsq_ = new PRECISION[IHx*JHx];
	ge_ = new PRECISION[IHx*JHx];
	gm_ = new PRECISION[IHx*JHx];
	// Auxiliary field scalars.
	ae0_ = new PRECISION[IHx*JHx];
	ae_ = new PRECISION[IHx*JHx];
	be_ = new PRECISION[IHx*JHx];
	ce_ = new PRECISION[IHx*JHx];
	de_ = new PRECISION[IHx*JHx];
	ee_ = new PRECISION[IHx*JHx];
	am0_ = new PRECISION[IHx*JHx];
	am_ = new PRECISION[IHx*JHx];
	bm_ = new PRECISION[IHx*JHx];
	cm_ = new PRECISION[IHx*JHx];
	dm_ = new PRECISION[IHx*JHx];
	em_ = new PRECISION[IHx*JHx];
	// PML arrays.
	PsiEzX_ = new PRECISION[IEz*JEz];
	PsiEzY_ = new PRECISION[IEz*JEz];
	PsiHyX_ = new PRECISION[IHy*JHy];
	PsiHxY_ = new PRECISION[IHx*JHx];

	return 0;
}
// Initialise CPU data.
int CFDTD2DDNG::InitialiseCPU()
{
	for (unsigned int i=0; i<IHx; i++)
	{
		for (unsigned int j=0; j<JHx; j++)
		{
			for (unsigned int n=0; n<3; n++)
			{
				if (i<IEz && j<JEz)
				{
					Ez(i,j,n) = 0.;
					Dz(i,j,n) = 0.;
				}
				Hx(i,j,n) = 0.;
				Bx(i,j,n) = 0.;
				if (i<IHy && j<JHy)
				{
					Hy(i,j,n) = 0.;
					By(i,j,n) = 0.;
				}
			}
			// Parameteric and auxiliary arrays.
			if (j<SlabLeft || j>SlabRight)
			{
				// Outside slab.
				einf(i,j) = 1.;
				uinf(i,j) = 1.;
				wpesq(i,j) = 0.;
				wpmsq(i,j) = 0.;
				ge(i,j) = 0.;
				gm(i,j) = 0.;
			}
			else // Inside slab.
			{
				einf(i,j) = 1.;
				uinf(i,j) = 1.;
				wpesq(i,j) = 2.*pow(w, 2);;
				wpmsq(i,j) = 2.*pow(w, 2);;
				ge(i,j) = w/32.;
				gm(i,j) = w/32.;
			}
			// Auxiliary scalars.
			ae0(i,j) = (4.*pow(dt,2))/(e0*(4.*einf(i,j)+pow(dt,2)*wpesq(i,j)+2.*dt*einf(i,j)*ge(i,j)));
			ae(i,j) = (1./pow(dt,2))*ae0(i,j);
			be(i,j) = (1./(2.*dt))*ge(i,j)*ae0(i,j);
			ce(i,j) = (e0/pow(dt,2))*einf(i,j)*ae0(i,j);
			de(i,j) = (-1.*e0/4.)*wpesq(i,j)*ae0(i,j);
			ee(i,j) = (1./(2.*dt))*e0*einf(i,j)*ge(i,j)*ae0(i,j);
			am0(i,j) = (4.*pow(dt,2))/(u0*(4.*uinf(i,j)+pow(dt,2)*wpmsq(i,j)+2.*dt*uinf(i,j)*gm(i,j)));;
			am(i,j) = (1./pow(dt,2))*am0(i,j);
			bm(i,j) = (1./(2.*dt))*gm(i,j)*am0(i,j);
			cm(i,j) = (u0/pow(dt,2))*uinf(i,j)*am0(i,j);
			dm(i,j) = (-1.*u0/4.)*wpmsq(i,j)*am0(i,j);
			em(i,j) = (1./(2.*dt))*u0*uinf(i,j)*gm(i,j)*am0(i,j);
			// PML Psi arrays.
			if (i<IEz && j<JEz)
			{
				PsiEzX(i,j) = 0.;
				PsiEzY(i,j) = 0.;
			}
			PsiHxY(i,j) = 0.;
			if (i<IHy && j<JHy)
				PsiHyX(i,j) = 0.;
		}
	}
	for (unsigned int i=0; i<MaxTime; i++)
	{
		Ezi[i] = 0.;
		Ezt[i] = 0.;
		Eztt[i] = 0.;
		Ezy1[i] = 0.;
		Ezy2[i] = 0.;
	}

	return 0;
}
int CFDTD2DDNG::InitialiseForSimulationCPU()
{
	for (unsigned int i=0; i<IHx; i++)
	{
		for (unsigned int j=0; j<JHx; j++)
		{
			for (unsigned int n=0; n<3; n++)
			{
				if (i<IEz && j<JEz)
				{
					Ez(i,j,n) = 0.;
					Dz(i,j,n) = 0.;
				}
				Hx(i,j,n) = 0.;
				Bx(i,j,n) = 0.;
				if (i<IHy && j<JHy)
				{
					Hy(i,j,n) = 0.;
					By(i,j,n) = 0.;
				}
			}
			// PML Psi arrays.
			if (i<IEz && j<JEz)
			{
				PsiEzX(i,j) = 0.;
				PsiEzY(i,j) = 0.;
			}
			PsiHxY(i,j) = 0.;
			if (i<IHy && j<JHy)
				PsiHyX(i,j) = 0.;
		}
	}

	return 0;
}
int CFDTD2DDNG::InitialiseCL()
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

	if(numPlatforms > 0)
	{
		cl_platform_id* platforms = new cl_platform_id[numPlatforms];
		SafeCall(clGetPlatformIDs(numPlatforms, platforms, NULL), "Error: Getting Platform Ids. (clGetPlatformsIDs)");

		for(unsigned int i=0; i < numPlatforms; ++i)
		{
			char pbuff[100];
			SafeCall(clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, sizeof(pbuff), pbuff, NULL), "Error: Getting Platform Info.(clGetPlatformInfo)");

			platform = platforms[i];
			std::cout << "Device" << i << " = " << pbuff << std::endl;
			if(!strcmp(pbuff, "Advanced Micro Devices, Inc."))
			{
				break;
			}
		}
		delete platforms;
	}

	if(NULL == platform)
	{
		std::cout << "NULL platform found so Exiting Application." << std::endl;
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

	if (cpu == true)
	{
		std::cout << "Running on CPU (GPU emulation)..." << std::endl;
		type = CL_DEVICE_TYPE_CPU;
	}
	else
	{
		std::cout << "Running on GPU..." << std::endl;
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

	/////////////////////////////////////////////////////////////////
	// Create an OpenCL command queue
	/////////////////////////////////////////////////////////////////
	commandQueue = clCreateCommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE, &status);
	SafeCall(status, "Creating Command Queue. (clCreateCommandQueue)");

	return 0;
}

int CFDTD2DDNG::AllocateMemoryGPU()
{
	cl_int status;
	/////////////////////////////////////////////////////////////////
	// Create OpenCL memory buffers
	/////////////////////////////////////////////////////////////////
	// Fields.
	d_Ez_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IEz*JEz*3, Ez_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Ex_");
	d_Dz_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IEz*JEz*3, Dz_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Dx_");
	d_Hx_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx*3, Hx_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Hx_");
	d_Bx_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx*3, Bx_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Bx_");
	d_Hy_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHy*JHy*3, Hy_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Hy_");
	d_By_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHy*JHy*3, By_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for By_");
	// Incident and transmitted fields.
	d_Ezi = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*MaxTime, Ezi, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Ezi");
	d_Ezt = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*MaxTime, Ezt, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Ezt");
	d_Eztt = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*MaxTime, Eztt, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Eztt");
	// Refractive Index.
	d_Ezy1 = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*MaxTime, Ezy1, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Ezy1");
	d_Ezy2 = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*MaxTime, Ezy2, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for Ezy2");
	// Drude material parameters.
	d_einf_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, einf_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for einf");
	d_uinf_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, uinf_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for uinf");
	d_wpesq_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, wpesq_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for wpesq");
	d_wpmsq_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, wpmsq_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for wpmsq");
	d_ge_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, ge_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for ge");
	d_gm_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, gm_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for gm");
	// Auxiliary field scalars.
	d_ae0_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, ae0_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for ae0");
	d_ae_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, ae_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for ae");
	d_be_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, be_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for be");
	d_ce_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, ce_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for ce");
	d_de_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, de_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for de");
	d_ee_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, ee_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for ee");

	d_am0_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, am0_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for ae0");
	d_am_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, am_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for am");
	d_bm_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, bm_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for bm");
	d_cm_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, cm_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for cm");
	d_dm_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, dm_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for dm");
	d_em_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, em_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for em");

	d_PsiEzX_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IEz*JEz, PsiEzX_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for PsiEzX");
	d_PsiEzY_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IEz*JEz, PsiEzY_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for PsiEzY");
	d_PsiHyX_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHy*JHy, PsiHyX_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for PsiHyX");
	d_PsiHxY_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, PsiHxY_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot creae input buffer for PsiHxY");

	return 0;
}
int CFDTD2DDNG::InitialiseCLKernelsGPU()
{
	int status;
	/////////////////////////////////////////////////////////////////
	// Load CL file, build CL program object, create CL kernel object
	/////////////////////////////////////////////////////////////////
	const char *filename = "FDTD2DDNG_Kernels.cl";
	string sourceStr = convertToString(filename);
	const char *source = sourceStr.c_str();
	size_t sourceSize[] = {strlen(source)};

	program = clCreateProgramWithSource( context, 1, &source, sourceSize, &status);
	SafeCall(status, "Error: Loading Binary into cl_program (clCreateProgramWithBinary)\n");

	/* create a cl program executable for all the devices specified */
	SafeCall(clBuildProgram(program, 1, devices, NULL, NULL, NULL), "Error: Building Program (clBuildProgram)\n");

	// Attach kernel objects to respective kernel functions.
	DryRun_kernel_M = clCreateKernel(program, "FDTD2DDNGKernel_DryRun_M", &status);
	SafeCall(status, "Error: Creating dry run Kernel from program. (clCreateKernel)");
	DryRun_kernel_E = clCreateKernel(program, "FDTD2DDNGKernel_DryRun_E", &status);
	SafeCall(status, "Error: Creating dry run Kernel from program. (clCreateKernel)");
	Simulation_kernel_M = clCreateKernel(program, "FDTD2DDNGKernel_Simulation_M", &status);
	SafeCall(status, "Error: Creating simulation Kernel from program. (clCreateKernel)");
	Simulation_kernel_E = clCreateKernel(program, "FDTD2DDNGKernel_Simulation_E", &status);
	SafeCall(status, "Error: Creating simulation Kernel from program. (clCreateKernel)");

	// ====== Set appropriate arguments to the Dry Run kernel ======
	SafeCall(clSetKernelArg(DryRun_kernel_M, 0, sizeof(unsigned int), (void *)&I), "Error: Setting kernel argument 'I'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 1, sizeof(unsigned int), (void *)&J), "Error: Setting kernel argument 'J'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 2, sizeof(unsigned int), (void *)&PMLw), "Error: Setting kernel argument 'PMLw'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 3, sizeof(unsigned int), (void *)&PulseWidth), "Error: Setting kernel argument 'PulseWidth'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 4, sizeof(unsigned int), (void *)&td), "Error: Setting kernel argument 'td'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 5, sizeof(unsigned int), (void *)&SourceChoice), "Error: Setting kernel argument 'SourceChoice'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 6, sizeof(unsigned int), (void *)&SourcePlane), "Error: Setting kernel argument 'SourcePlane'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 7, sizeof(unsigned int), (void *)&SourceLocationX), "Error: Setting kernel argument 'SourceLocationX'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 8, sizeof(unsigned int), (void *)&SourceLocationY), "Error: Setting kernel argument 'SourceLocationY'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 9, sizeof(PRECISION), (void *)&c), "Error: Setting kernel argument 'c'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 10, sizeof(PRECISION), (void *)&pi), "Error: Setting kernel argument 'pi'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 11, sizeof(PRECISION), (void *)&e0), "Error: Setting kernel argument 'e0'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 12, sizeof(PRECISION), (void *)&u0), "Error: Setting kernel argument 'u0'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 13, sizeof(PRECISION), (void *)&dt), "Error: Setting kernel argument 'dt'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 14, sizeof(PRECISION), (void *)&delta), "Error: Setting kernel argument 'delta'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 15, sizeof(PRECISION), (void *)&Sc), "Error: Setting kernel argument 'Sc'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 16, sizeof(PRECISION), (void *)&f), "Error: Setting kernel argument 'f'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 17, sizeof(PRECISION), (void *)&fp), "Error: Setting kernel argument 'fp'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 18, sizeof(PRECISION), (void *)&dr), "Error: Setting kernel argument 'dr'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 19, sizeof(unsigned int), (void *)&IEz), "Error: Setting kernel argument 'IEz'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 20, sizeof(unsigned int), (void *)&JEz), "Error: Setting kernel argument 'JEz'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 21, sizeof(unsigned int), (void *)&IHx), "Error: Setting kernel argument 'IHx'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 22, sizeof(unsigned int), (void *)&JHx), "Error: Setting kernel argument 'JHx'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 23, sizeof(unsigned int), (void *)&IHy), "Error: Setting kernel argument 'IHy'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 24, sizeof(unsigned int), (void *)&JHy), "Error: Setting kernel argument 'JHy'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 25, sizeof(cl_mem), (void *)&d_Ez_), "Error: Setting kernel argument d_Ez_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 26, sizeof(cl_mem), (void *)&d_Dz_), "Error: Setting kernel argument d_Dz_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 27, sizeof(cl_mem), (void *)&d_Hx_), "Error: Setting kernel argument d_Hx_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 28, sizeof(cl_mem), (void *)&d_Bx_), "Error: Setting kernel argument d_Bx_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 29, sizeof(cl_mem), (void *)&d_Hy_), "Error: Setting kernel argument d_Hy_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 30, sizeof(cl_mem), (void *)&d_By_), "Error: Setting kernel argument d_By_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 31, sizeof(cl_mem), (void *)&d_PsiEzX_), "Error: Setting kernel argument d_PsiEzX_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 32, sizeof(cl_mem), (void *)&d_PsiEzY_), "Error: Setting kernel argument d_PsiEzY_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 33, sizeof(cl_mem), (void *)&d_PsiHyX_), "Error: Setting kernel argument d_PsiHyX_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 34, sizeof(cl_mem), (void *)&d_PsiHxY_), "Error: Setting kernel argument d_PsiHxY_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 35, sizeof(cl_mem), (void *)&d_einf_), "Error: Setting kernel argument d_einf_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 36, sizeof(cl_mem), (void *)&d_uinf_), "Error: Setting kernel argument d_uinf_");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 37, sizeof(PRECISION), (void *)&kappex), "Error: Setting kernel argument 'kappex'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 38, sizeof(PRECISION), (void *)&kappey), "Error: Setting kernel argument 'kappey'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 39, sizeof(PRECISION), (void *)&kappmx), "Error: Setting kernel argument 'kappmx'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 40, sizeof(PRECISION), (void *)&kappmy), "Error: Setting kernel argument 'kappmy'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 41, sizeof(PRECISION), (void *)&bex), "Error: Setting kernel argument 'bex'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 42, sizeof(PRECISION), (void *)&bey), "Error: Setting kernel argument 'bey'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 43, sizeof(PRECISION), (void *)&bmx), "Error: Setting kernel argument 'bmx'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 44, sizeof(PRECISION), (void *)&bmy), "Error: Setting kernel argument 'bmy'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 45, sizeof(PRECISION), (void *)&Cex), "Error: Setting kernel argument 'Cex'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 46, sizeof(PRECISION), (void *)&Cey), "Error: Setting kernel argument 'Cey'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 47, sizeof(PRECISION), (void *)&Cmx), "Error: Setting kernel argument 'Cmx'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 48, sizeof(PRECISION), (void *)&Cmy), "Error: Setting kernel argument 'Cmy'");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 49, sizeof(cl_mem), (void *)&d_Ezi), "Error: Setting kernel argument d_Ezi");
	SafeCall(clSetKernelArg(DryRun_kernel_M, 50, sizeof(unsigned int), (void *)&x1), "Error: Setting kernel argument 'x1'");

	SafeCall(clSetKernelArg(DryRun_kernel_E, 0, sizeof(unsigned int), (void *)&I), "Error: Setting kernel argument 'I'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 1, sizeof(unsigned int), (void *)&J), "Error: Setting kernel argument 'J'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 2, sizeof(unsigned int), (void *)&PMLw), "Error: Setting kernel argument 'PMLw'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 3, sizeof(unsigned int), (void *)&PulseWidth), "Error: Setting kernel argument 'PulseWidth'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 4, sizeof(unsigned int), (void *)&td), "Error: Setting kernel argument 'td'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 5, sizeof(unsigned int), (void *)&SourceChoice), "Error: Setting kernel argument 'SourceChoice'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 6, sizeof(unsigned int), (void *)&SourcePlane), "Error: Setting kernel argument 'SourcePlane'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 7, sizeof(unsigned int), (void *)&SourceLocationX), "Error: Setting kernel argument 'SourceLocationX'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 8, sizeof(unsigned int), (void *)&SourceLocationY), "Error: Setting kernel argument 'SourceLocationY'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 9, sizeof(PRECISION), (void *)&c), "Error: Setting kernel argument 'c'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 10, sizeof(PRECISION), (void *)&pi), "Error: Setting kernel argument 'pi'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 11, sizeof(PRECISION), (void *)&e0), "Error: Setting kernel argument 'e0'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 12, sizeof(PRECISION), (void *)&u0), "Error: Setting kernel argument 'u0'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 13, sizeof(PRECISION), (void *)&dt), "Error: Setting kernel argument 'dt'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 14, sizeof(PRECISION), (void *)&delta), "Error: Setting kernel argument 'delta'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 15, sizeof(PRECISION), (void *)&Sc), "Error: Setting kernel argument 'Sc'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 16, sizeof(PRECISION), (void *)&f), "Error: Setting kernel argument 'f'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 17, sizeof(PRECISION), (void *)&fp), "Error: Setting kernel argument 'fp'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 18, sizeof(PRECISION), (void *)&dr), "Error: Setting kernel argument 'dr'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 19, sizeof(unsigned int), (void *)&IEz), "Error: Setting kernel argument 'IEz'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 20, sizeof(unsigned int), (void *)&JEz), "Error: Setting kernel argument 'JEz'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 21, sizeof(unsigned int), (void *)&IHx), "Error: Setting kernel argument 'IHx'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 22, sizeof(unsigned int), (void *)&JHx), "Error: Setting kernel argument 'JHx'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 23, sizeof(unsigned int), (void *)&IHy), "Error: Setting kernel argument 'IHy'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 24, sizeof(unsigned int), (void *)&JHy), "Error: Setting kernel argument 'JHy'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 25, sizeof(cl_mem), (void *)&d_Ez_), "Error: Setting kernel argument d_Ez_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 26, sizeof(cl_mem), (void *)&d_Dz_), "Error: Setting kernel argument d_Dz_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 27, sizeof(cl_mem), (void *)&d_Hx_), "Error: Setting kernel argument d_Hx_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 28, sizeof(cl_mem), (void *)&d_Bx_), "Error: Setting kernel argument d_Bx_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 29, sizeof(cl_mem), (void *)&d_Hy_), "Error: Setting kernel argument d_Hy_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 30, sizeof(cl_mem), (void *)&d_By_), "Error: Setting kernel argument d_By_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 31, sizeof(cl_mem), (void *)&d_PsiEzX_), "Error: Setting kernel argument d_PsiEzX_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 32, sizeof(cl_mem), (void *)&d_PsiEzY_), "Error: Setting kernel argument d_PsiEzY_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 33, sizeof(cl_mem), (void *)&d_PsiHyX_), "Error: Setting kernel argument d_PsiHyX_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 34, sizeof(cl_mem), (void *)&d_PsiHxY_), "Error: Setting kernel argument d_PsiHxY_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 35, sizeof(cl_mem), (void *)&d_einf_), "Error: Setting kernel argument d_einf_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 36, sizeof(cl_mem), (void *)&d_uinf_), "Error: Setting kernel argument d_uinf_");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 37, sizeof(PRECISION), (void *)&kappex), "Error: Setting kernel argument 'kappex'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 38, sizeof(PRECISION), (void *)&kappey), "Error: Setting kernel argument 'kappey'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 39, sizeof(PRECISION), (void *)&kappmx), "Error: Setting kernel argument 'kappmx'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 40, sizeof(PRECISION), (void *)&kappmy), "Error: Setting kernel argument 'kappmy'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 41, sizeof(PRECISION), (void *)&bex), "Error: Setting kernel argument 'bex'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 42, sizeof(PRECISION), (void *)&bey), "Error: Setting kernel argument 'bey'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 43, sizeof(PRECISION), (void *)&bmx), "Error: Setting kernel argument 'bmx'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 44, sizeof(PRECISION), (void *)&bmy), "Error: Setting kernel argument 'bmy'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 45, sizeof(PRECISION), (void *)&Cex), "Error: Setting kernel argument 'Cex'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 46, sizeof(PRECISION), (void *)&Cey), "Error: Setting kernel argument 'Cey'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 47, sizeof(PRECISION), (void *)&Cmx), "Error: Setting kernel argument 'Cmx'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 48, sizeof(PRECISION), (void *)&Cmy), "Error: Setting kernel argument 'Cmy'");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 49, sizeof(cl_mem), (void *)&d_Ezi), "Error: Setting kernel argument d_Ezi");
	SafeCall(clSetKernelArg(DryRun_kernel_E, 50, sizeof(unsigned int), (void *)&x1), "Error: Setting kernel argument 'x1'");

	// ====== Set appropriate arguments to the Simulation kernel ======
	SafeCall(clSetKernelArg(Simulation_kernel_M, 0, sizeof(unsigned int), (void *)&I), "Error: Setting kernel argument 'I'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 1, sizeof(unsigned int), (void *)&J), "Error: Setting kernel argument 'J'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 2, sizeof(unsigned int), (void *)&PMLw), "Error: Setting kernel argument 'PMLw'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 3, sizeof(unsigned int), (void *)&PulseWidth), "Error: Setting kernel argument 'PulseWidth'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 4, sizeof(unsigned int), (void *)&td), "Error: Setting kernel argument 'td'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 5, sizeof(unsigned int), (void *)&SourceChoice), "Error: Setting kernel argument 'SourceChoice'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 6, sizeof(unsigned int), (void *)&SourcePlane), "Error: Setting kernel argument 'SourcePlane'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 7, sizeof(unsigned int), (void *)&SourceLocationX), "Error: Setting kernel argument 'SourceLocationX'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 8, sizeof(unsigned int), (void *)&SourceLocationY), "Error: Setting kernel argument 'SourceLocationY'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 9, sizeof(PRECISION), (void *)&c), "Error: Setting kernel argument 'c'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 10, sizeof(PRECISION), (void *)&pi), "Error: Setting kernel argument 'pi'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 11, sizeof(PRECISION), (void *)&e0), "Error: Setting kernel argument 'e0'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 12, sizeof(PRECISION), (void *)&u0), "Error: Setting kernel argument 'u0'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 13, sizeof(PRECISION), (void *)&dt), "Error: Setting kernel argument 'dt'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 14, sizeof(PRECISION), (void *)&delta), "Error: Setting kernel argument 'delta'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 15, sizeof(PRECISION), (void *)&Sc), "Error: Setting kernel argument 'Sc'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 16, sizeof(PRECISION), (void *)&f), "Error: Setting kernel argument 'f'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 17, sizeof(PRECISION), (void *)&fp), "Error: Setting kernel argument 'fp'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 18, sizeof(PRECISION), (void *)&dr), "Error: Setting kernel argument 'dr'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 19, sizeof(unsigned int), (void *)&IEz), "Error: Setting kernel argument 'IEz'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 20, sizeof(unsigned int), (void *)&JEz), "Error: Setting kernel argument 'JEz'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 21, sizeof(unsigned int), (void *)&IHx), "Error: Setting kernel argument 'IHx'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 22, sizeof(unsigned int), (void *)&JHx), "Error: Setting kernel argument 'JHx'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 23, sizeof(unsigned int), (void *)&IHy), "Error: Setting kernel argument 'IHy'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 24, sizeof(unsigned int), (void *)&JHy), "Error: Setting kernel argument 'JHy'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 25, sizeof(cl_mem), (void *)&d_Ez_), "Error: Setting kernel argument d_Ez_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 26, sizeof(cl_mem), (void *)&d_Dz_), "Error: Setting kernel argument d_Dz_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 27, sizeof(cl_mem), (void *)&d_Hx_), "Error: Setting kernel argument d_Hx_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 28, sizeof(cl_mem), (void *)&d_Bx_), "Error: Setting kernel argument d_Bx_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 29, sizeof(cl_mem), (void *)&d_Hy_), "Error: Setting kernel argument d_Hy_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 30, sizeof(cl_mem), (void *)&d_By_), "Error: Setting kernel argument d_By_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 31, sizeof(cl_mem), (void *)&d_einf_), "Error: Setting kernel argument d_einf");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 32, sizeof(cl_mem), (void *)&d_uinf_), "Error: Setting kernel argument d_uinf");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 33, sizeof(cl_mem), (void *)&d_wpesq_), "Error: Setting kernel argument d_wpesq");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 34, sizeof(cl_mem), (void *)&d_wpmsq_), "Error: Setting kernel argument d_wpmsq");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 35, sizeof(cl_mem), (void *)&d_ge_), "Error: Setting kernel argument d_ge");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 36, sizeof(cl_mem), (void *)&d_gm_), "Error: Setting kernel argument d_gm");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 37, sizeof(cl_mem), (void *)&d_ae0_), "Error: Setting kernel argument d_ae0");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 38, sizeof(cl_mem), (void *)&d_ae_), "Error: Setting kernel argument d_ae");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 39, sizeof(cl_mem), (void *)&d_be_), "Error: Setting kernel argument d_be");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 40, sizeof(cl_mem), (void *)&d_ce_), "Error: Setting kernel argument d_ce");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 41, sizeof(cl_mem), (void *)&d_de_), "Error: Setting kernel argument d_de");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 42, sizeof(cl_mem), (void *)&d_ee_), "Error: Setting kernel argument d_ee");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 43, sizeof(cl_mem), (void *)&d_am0_), "Error: Setting kernel argument d_am0");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 44, sizeof(cl_mem), (void *)&d_am_), "Error: Setting kernel argument d_am");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 45, sizeof(cl_mem), (void *)&d_bm_), "Error: Setting kernel argument d_bm");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 46, sizeof(cl_mem), (void *)&d_cm_), "Error: Setting kernel argument d_cm");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 47, sizeof(cl_mem), (void *)&d_dm_), "Error: Setting kernel argument d_dm");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 48, sizeof(cl_mem), (void *)&d_em_), "Error: Setting kernel argument d_em");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 49, sizeof(cl_mem), (void *)&d_PsiEzX_), "Error: Setting kernel argument d_PsiEzX_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 50, sizeof(cl_mem), (void *)&d_PsiEzY_), "Error: Setting kernel argument d_PsiEzY_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 51, sizeof(cl_mem), (void *)&d_PsiHyX_), "Error: Setting kernel argument d_PsiHyX_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 52, sizeof(cl_mem), (void *)&d_PsiHxY_), "Error: Setting kernel argument d_PsiHxY_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 53, sizeof(PRECISION), (void *)&kappex), "Error: Setting kernel argument 'kappex'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 54, sizeof(PRECISION), (void *)&kappey), "Error: Setting kernel argument 'kappey'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 55, sizeof(PRECISION), (void *)&kappmx), "Error: Setting kernel argument 'kappmx'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 56, sizeof(PRECISION), (void *)&kappmy), "Error: Setting kernel argument 'kappmy'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 57, sizeof(PRECISION), (void *)&bex), "Error: Setting kernel argument 'bex'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 58, sizeof(PRECISION), (void *)&bey), "Error: Setting kernel argument 'bey'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 59, sizeof(PRECISION), (void *)&bmx), "Error: Setting kernel argument 'bmx'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 60, sizeof(PRECISION), (void *)&bmy), "Error: Setting kernel argument 'bmy'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 61, sizeof(PRECISION), (void *)&Cex), "Error: Setting kernel argument 'Cex'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 62, sizeof(PRECISION), (void *)&Cey), "Error: Setting kernel argument 'Cey'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 63, sizeof(PRECISION), (void *)&Cmx), "Error: Setting kernel argument 'Cmx'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 64, sizeof(PRECISION), (void *)&Cmy), "Error: Setting kernel argument 'Cmy'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 65, sizeof(cl_mem), (void *)&d_Ezt), "Error: Setting kernel argument d_Ezt");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 66, sizeof(cl_mem), (void *)&d_Eztt), "Error: Setting kernel argument d_Eztt");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 67, sizeof(cl_mem), (void *)&d_Ezy1), "Error: Setting kernel argument d_Ezy1");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 68, sizeof(cl_mem), (void *)&d_Ezy2), "Error: Setting kernel argument d_Ezy2");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 69, sizeof(unsigned int), (void *)&x1), "Error: Setting kernel argument 'x1'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 70, sizeof(unsigned int), (void *)&Y1), "Error: Setting kernel argument 'Y1'");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 71, sizeof(unsigned int), (void *)&Y2), "Error: Setting kernel argument 'Y2'");

	SafeCall(clSetKernelArg(Simulation_kernel_E, 0, sizeof(unsigned int), (void *)&I), "Error: Setting kernel argument 'I'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 1, sizeof(unsigned int), (void *)&J), "Error: Setting kernel argument 'J'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 2, sizeof(unsigned int), (void *)&PMLw), "Error: Setting kernel argument 'PMLw'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 3, sizeof(unsigned int), (void *)&PulseWidth), "Error: Setting kernel argument 'PulseWidth'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 4, sizeof(unsigned int), (void *)&td), "Error: Setting kernel argument 'td'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 5, sizeof(unsigned int), (void *)&SourceChoice), "Error: Setting kernel argument 'SourceChoice'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 6, sizeof(unsigned int), (void *)&SourcePlane), "Error: Setting kernel argument 'SourcePlane'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 7, sizeof(unsigned int), (void *)&SourceLocationX), "Error: Setting kernel argument 'SourceLocationX'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 8, sizeof(unsigned int), (void *)&SourceLocationY), "Error: Setting kernel argument 'SourceLocationY'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 9, sizeof(PRECISION), (void *)&c), "Error: Setting kernel argument 'c'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 10, sizeof(PRECISION), (void *)&pi), "Error: Setting kernel argument 'pi'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 11, sizeof(PRECISION), (void *)&e0), "Error: Setting kernel argument 'e0'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 12, sizeof(PRECISION), (void *)&u0), "Error: Setting kernel argument 'u0'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 13, sizeof(PRECISION), (void *)&dt), "Error: Setting kernel argument 'dt'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 14, sizeof(PRECISION), (void *)&delta), "Error: Setting kernel argument 'delta'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 15, sizeof(PRECISION), (void *)&Sc), "Error: Setting kernel argument 'Sc'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 16, sizeof(PRECISION), (void *)&f), "Error: Setting kernel argument 'f'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 17, sizeof(PRECISION), (void *)&fp), "Error: Setting kernel argument 'fp'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 18, sizeof(PRECISION), (void *)&dr), "Error: Setting kernel argument 'dr'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 19, sizeof(unsigned int), (void *)&IEz), "Error: Setting kernel argument 'IEz'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 20, sizeof(unsigned int), (void *)&JEz), "Error: Setting kernel argument 'JEz'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 21, sizeof(unsigned int), (void *)&IHx), "Error: Setting kernel argument 'IHx'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 22, sizeof(unsigned int), (void *)&JHx), "Error: Setting kernel argument 'JHx'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 23, sizeof(unsigned int), (void *)&IHy), "Error: Setting kernel argument 'IHy'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 24, sizeof(unsigned int), (void *)&JHy), "Error: Setting kernel argument 'JHy'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 25, sizeof(cl_mem), (void *)&d_Ez_), "Error: Setting kernel argument d_Ez_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 26, sizeof(cl_mem), (void *)&d_Dz_), "Error: Setting kernel argument d_Dz_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 27, sizeof(cl_mem), (void *)&d_Hx_), "Error: Setting kernel argument d_Hx_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 28, sizeof(cl_mem), (void *)&d_Bx_), "Error: Setting kernel argument d_Bx_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 29, sizeof(cl_mem), (void *)&d_Hy_), "Error: Setting kernel argument d_Hy_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 30, sizeof(cl_mem), (void *)&d_By_), "Error: Setting kernel argument d_By_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 31, sizeof(cl_mem), (void *)&d_einf_), "Error: Setting kernel argument d_einf");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 32, sizeof(cl_mem), (void *)&d_uinf_), "Error: Setting kernel argument d_uinf");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 33, sizeof(cl_mem), (void *)&d_wpesq_), "Error: Setting kernel argument d_wpesq");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 34, sizeof(cl_mem), (void *)&d_wpmsq_), "Error: Setting kernel argument d_wpmsq");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 35, sizeof(cl_mem), (void *)&d_ge_), "Error: Setting kernel argument d_ge");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 36, sizeof(cl_mem), (void *)&d_gm_), "Error: Setting kernel argument d_gm");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 37, sizeof(cl_mem), (void *)&d_ae0_), "Error: Setting kernel argument d_ae0");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 38, sizeof(cl_mem), (void *)&d_ae_), "Error: Setting kernel argument d_ae");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 39, sizeof(cl_mem), (void *)&d_be_), "Error: Setting kernel argument d_be");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 40, sizeof(cl_mem), (void *)&d_ce_), "Error: Setting kernel argument d_ce");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 41, sizeof(cl_mem), (void *)&d_de_), "Error: Setting kernel argument d_de");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 42, sizeof(cl_mem), (void *)&d_ee_), "Error: Setting kernel argument d_ee");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 43, sizeof(cl_mem), (void *)&d_am0_), "Error: Setting kernel argument d_am0");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 44, sizeof(cl_mem), (void *)&d_am_), "Error: Setting kernel argument d_am");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 45, sizeof(cl_mem), (void *)&d_bm_), "Error: Setting kernel argument d_bm");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 46, sizeof(cl_mem), (void *)&d_cm_), "Error: Setting kernel argument d_cm");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 47, sizeof(cl_mem), (void *)&d_dm_), "Error: Setting kernel argument d_dm");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 48, sizeof(cl_mem), (void *)&d_em_), "Error: Setting kernel argument d_em");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 49, sizeof(cl_mem), (void *)&d_PsiEzX_), "Error: Setting kernel argument d_PsiEzX_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 50, sizeof(cl_mem), (void *)&d_PsiEzY_), "Error: Setting kernel argument d_PsiEzY_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 51, sizeof(cl_mem), (void *)&d_PsiHyX_), "Error: Setting kernel argument d_PsiHyX_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 52, sizeof(cl_mem), (void *)&d_PsiHxY_), "Error: Setting kernel argument d_PsiHxY_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 53, sizeof(PRECISION), (void *)&kappex), "Error: Setting kernel argument 'kappex'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 54, sizeof(PRECISION), (void *)&kappey), "Error: Setting kernel argument 'kappey'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 55, sizeof(PRECISION), (void *)&kappmx), "Error: Setting kernel argument 'kappmx'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 56, sizeof(PRECISION), (void *)&kappmy), "Error: Setting kernel argument 'kappmy'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 57, sizeof(PRECISION), (void *)&bex), "Error: Setting kernel argument 'bex'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 58, sizeof(PRECISION), (void *)&bey), "Error: Setting kernel argument 'bey'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 59, sizeof(PRECISION), (void *)&bmx), "Error: Setting kernel argument 'bmx'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 60, sizeof(PRECISION), (void *)&bmy), "Error: Setting kernel argument 'bmy'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 61, sizeof(PRECISION), (void *)&Cex), "Error: Setting kernel argument 'Cex'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 62, sizeof(PRECISION), (void *)&Cey), "Error: Setting kernel argument 'Cey'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 63, sizeof(PRECISION), (void *)&Cmx), "Error: Setting kernel argument 'Cmx'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 64, sizeof(PRECISION), (void *)&Cmy), "Error: Setting kernel argument 'Cmy'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 65, sizeof(cl_mem), (void *)&d_Ezt), "Error: Setting kernel argument d_Ezt");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 66, sizeof(cl_mem), (void *)&d_Eztt), "Error: Setting kernel argument d_Eztt");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 67, sizeof(cl_mem), (void *)&d_Ezy1), "Error: Setting kernel argument d_Ezy1");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 68, sizeof(cl_mem), (void *)&d_Ezy2), "Error: Setting kernel argument d_Ezy2");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 69, sizeof(unsigned int), (void *)&x1), "Error: Setting kernel argument 'x1'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 70, sizeof(unsigned int), (void *)&Y1), "Error: Setting kernel argument 'Y1'");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 71, sizeof(unsigned int), (void *)&Y2), "Error: Setting kernel argument 'Y2'");

	return 0;
}
int CFDTD2DDNG::InitialiseForSimulationGPU()
{
	int status;

	SafeCall(clReleaseMemObject(d_Ez_), "Error: clReleaseMemObject() cannot release memory input buffer for d_Ez_");
	SafeCall(clReleaseMemObject(d_Dz_), "Error: clReleaseMemObject() cannot release memory input buffer for d_Dz_");
	SafeCall(clReleaseMemObject(d_Hx_), "Error: clReleaseMemObject() cannot release memory input buffer for d_Hx_");
	SafeCall(clReleaseMemObject(d_Bx_), "Error: clReleaseMemObject() cannot release memory input buffer for d_Bx_");
	SafeCall(clReleaseMemObject(d_Hy_), "Error: clReleaseMemObject() cannot release memory input buffer for d_Hy_");
	SafeCall(clReleaseMemObject(d_By_), "Error: clReleaseMemObject() cannot release memory input buffer for d_By_");

	SafeCall(clReleaseMemObject(d_PsiEzX_), "Error: clReleaseMemObject() cannot release memory input buffer for d_PsiEzX_");
	SafeCall(clReleaseMemObject(d_PsiEzY_), "Error: clReleaseMemObject() cannot release memory input buffer for d_PsiEzY_");
	SafeCall(clReleaseMemObject(d_PsiHyX_), "Error: clReleaseMemObject() cannot release memory input buffer for d_PsiHyX_");
	SafeCall(clReleaseMemObject(d_PsiHxY_), "Error: clReleaseMemObject() cannot release memory input buffer for d_PsiHxY_");

	d_Ez_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IEz*JEz*3, Ez_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create input buffer for Ex_");
	d_Dz_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IEz*JEz*3, Dz_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create input buffer for Dx_");
	d_Hx_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx*3, Hx_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create input buffer for Hx_");
	d_Bx_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx*3, Bx_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create input buffer for Bx_");
	d_Hy_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHy*JHy*3, Hy_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create input buffer for Hy_");
	d_By_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHy*JHy*3, By_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create input buffer for By_");

	d_PsiEzX_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IEz*JEz, PsiEzX_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create input buffer for PsiEzX");
	d_PsiEzY_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IEz*JEz, PsiEzY_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create input buffer for PsiEzY");
	d_PsiHyX_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHy*JHy, PsiHyX_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create input buffer for PsiHyX");
	d_PsiHxY_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*IHx*JHx, PsiHxY_, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create input buffer for PsiHxY");

	SafeCall(clSetKernelArg(Simulation_kernel_M, 25, sizeof(cl_mem), (void *)&d_Ez_), "Error: Setting kernel argument d_Ez_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 26, sizeof(cl_mem), (void *)&d_Dz_), "Error: Setting kernel argument d_Dz_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 27, sizeof(cl_mem), (void *)&d_Hx_), "Error: Setting kernel argument d_Hx_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 28, sizeof(cl_mem), (void *)&d_Bx_), "Error: Setting kernel argument d_Bx_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 29, sizeof(cl_mem), (void *)&d_Hy_), "Error: Setting kernel argument d_Hy_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 30, sizeof(cl_mem), (void *)&d_By_), "Error: Setting kernel argument d_By_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 49, sizeof(cl_mem), (void *)&d_PsiEzX_), "Error: Setting kernel argument d_PsiEzX_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 50, sizeof(cl_mem), (void *)&d_PsiEzY_), "Error: Setting kernel argument d_PsiEzY_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 51, sizeof(cl_mem), (void *)&d_PsiHyX_), "Error: Setting kernel argument d_PsiHyX_");
	SafeCall(clSetKernelArg(Simulation_kernel_M, 52, sizeof(cl_mem), (void *)&d_PsiHxY_), "Error: Setting kernel argument d_PsiHxY_");

	SafeCall(clSetKernelArg(Simulation_kernel_E, 25, sizeof(cl_mem), (void *)&d_Ez_), "Error: Setting kernel argument d_Ez_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 26, sizeof(cl_mem), (void *)&d_Dz_), "Error: Setting kernel argument d_Dz_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 27, sizeof(cl_mem), (void *)&d_Hx_), "Error: Setting kernel argument d_Hx_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 28, sizeof(cl_mem), (void *)&d_Bx_), "Error: Setting kernel argument d_Bx_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 29, sizeof(cl_mem), (void *)&d_Hy_), "Error: Setting kernel argument d_Hy_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 30, sizeof(cl_mem), (void *)&d_By_), "Error: Setting kernel argument d_By_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 49, sizeof(cl_mem), (void *)&d_PsiEzX_), "Error: Setting kernel argument d_PsiEzX_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 50, sizeof(cl_mem), (void *)&d_PsiEzY_), "Error: Setting kernel argument d_PsiEzY_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 51, sizeof(cl_mem), (void *)&d_PsiHyX_), "Error: Setting kernel argument d_PsiHyX_");
	SafeCall(clSetKernelArg(Simulation_kernel_E, 52, sizeof(cl_mem), (void *)&d_PsiHxY_), "Error: Setting kernel argument d_PsiHxY_");

	return 0;
}
int CFDTD2DDNG::DryRunCPU()
{
	cout << "Dry run (CPU) started..." << endl;
	for (unsigned int n=0; n<MaxTime; n++)
	{
		if (n % (MaxTime/256U) == 0)
			cout << "\r\t\t\r" << n*100/(MaxTime-1) << "%";
		// ========================== Bx and Hx ==========================
		for (unsigned int i=0; i<IHx; i++)
		{
			// Calculation of PsiHxY.
			for (unsigned int j=1; j<JHx-1; j++)
				PsiHxY(i,j) = (Cmy/delta)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) + bmy*PsiHxY(i,j);

			// Bx in normal space.
			for (unsigned int j=1+PMLw; j<JHx-PMLw-1; j++)
			{
				Bx(i,j,nf) = Bx(i,j,n0) + (-Ez(i,j,n0) + Ez(i,j-1,n0)) * dt/delta;
				Hx(i,j,nf) = Bx(i,j,nf)/(u0*uinf(i,j));
			}
			// Bx in lower PML.
			for (unsigned int j=1; j<PMLw+1; j++)
			{
				Bx(i,j,nf) = Bx(i,j,n0) + dt*((1./kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1./delta + PsiHxY(i,j));
				Hx(i,j,nf) = Bx(i,j,nf)/(u0*uinf(i,j));
			}
			// Bx in upper PML.
			for (unsigned int j=JHx-PMLw-1; j<JHx-1; j++)
			{
				Bx(i,j,nf) = Bx(i,j,n0) + dt*((1./kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1./delta + PsiHxY(i,j));
				Hx(i,j,nf) = Bx(i,j,nf)/(u0*uinf(i,j));
			}
		}
		
		// ========================== By and Hy ==========================
		for (unsigned int i=0; i<IHy-1; i++)
		{
			// PsiHyX arrays.
			for (unsigned int j=0; j<JHy; j++)
			{
				PsiHyX(i,j) = (Cmx/delta)*(Ez(i+1,j,n0)-Ez(i,j,n0)) + bmx*PsiHyX(i,j);
				if (i==0)
					PsiHyX(IHy-1,j) = (Cmx/delta)*(Ez(0,j,n0)-Ez(IHy-1,j,n0)) + bmx*PsiHyX(IHy-1,j); // PBC
			}
			// By in normal space.
			for (unsigned int j=PMLw; j<JHy-PMLw; j++)
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
			for (unsigned int j=0; j<PMLw; j++)
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
			for (unsigned int j=JHy-PMLw; j<JHy; j++)
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

		// ========================== Dz and Ez ==========================
		for (unsigned int i=1; i<IEz; i++)
		{
			// Psi arrays.
			for (unsigned int j=0; j<JEz; j++)
			{
				PsiEzX(i,j) = (Cex/delta)*(Hy(i,j,nf)-Hy(i-1,j,nf)) + bex*PsiEzX(i,j);
				if (i==1)
					PsiEzX(0,j) = (Cex/delta)*(Hy(0,j,nf)-Hy(IEz-1,j,nf)) + bex*PsiEzX(0,j); // PBC
				PsiEzY(i,j) = (Cey/delta)*(-Hx(i,j+1,nf)+Hx(i,j,nf)) + bey*PsiEzY(i,j);
				if (i==1)
					PsiEzY(0,j) = (Cey/delta)*(-Hx(0,j+1,nf)+Hx(0,j,nf)) + bey*PsiEzY(0,j); // PBC
			}
			// Dz in normal space.
			for (unsigned int j=PMLw; j<JEz-PMLw; j++)
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
			for (unsigned int j=0; j<PMLw; j++)
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
			for (unsigned int j=JEz-PMLw; j<JEz; j++)
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
		if (SourcePlane == 1)
		{
			unsigned int j=SourceLocationY;
			for (unsigned int i=0; i<IEz; i++)
			{
				if (SourceChoice == 1)
					Ez(i,j,nf) += exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
				else if (SourceChoice == 2)
					Ez(i,j,nf) += sin(2.*pi*f*(PRECISION)n*dt) * Sc;
				else if (SourceChoice == 3)
					Ez(i,j,nf) += (1.-2.*pow(pi*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(pi*fp*((PRECISION)n*dt-dr),2)) * Sc;

				Dz(i,j,nf) = e0*Ez(i,j,nf);
			}
		}
		else
		{
			unsigned int i=SourceLocationX;
			unsigned int j=SourceLocationY;
			if (SourceChoice == 1)
				Ez(i,j,nf) += exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
			else if (SourceChoice == 2)
				Ez(i,j,nf) += sin(2.*pi*f*(PRECISION)n*dt) * Sc;
			else if (SourceChoice == 3)
				Ez(i,j,nf) += (1.-2.*pow(pi*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(pi*fp*((PRECISION)n*dt-dr),2)) * Sc;

			Dz(i,j,nf) = e0*Ez(i,j,nf);
		}
		Ezi[n] = Ez(IEz/2,x1,nf); // Incident field.

		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}
	cout << endl << "Dry run (CPU) completed!" << endl;
	return 0;
}
int CFDTD2DDNG::RunSimulationCPU(bool SaveFields)
{
	stringstream framestream;
	string basename = "FieldData/Ez";
	string filename;
	fstream snapshot;
	frame = 0U;

	cout << "Simulation (CPU) started..." << endl;
	for (unsigned int n=0; n<MaxTime; n++)
	{
		if (n % (MaxTime/256U) == 0)
			cout << "\r\t\t\r" << n*100/(MaxTime-1) << "%";
		// ========================== Bx and Hx ==========================
		for (unsigned int i=0; i<IHx; i++)
		{
			// Calculation of PsiHxY.
			for (unsigned int j=1; j<JHx-1; j++)
				PsiHxY(i,j) = (Cmy/delta)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) + bmy*PsiHxY(i,j);

			// Bx in normal space.
			for (unsigned int j=1+PMLw; j<JHx-PMLw-1; j++)
			{
				Bx(i,j,nf) = Bx(i,j,n0) + (-Ez(i,j,n0) + Ez(i,j-1,n0)) * dt/delta;
				Hx(i,j,nf) = am(i,j)*(Bx(i,j,nf)-2.*Bx(i,j,n0)+Bx(i,j,np))+bm(i,j)*(Bx(i,j,nf)-Bx(i,j,np))+cm(i,j)*(2.*Hx(i,j,n0)-Hx(i,j,np))+dm(i,j)*(2.*Hx(i,j,n0)+Hx(i,j,np))+em(i,j)*Hx(i,j,np);
			}
			// Bx in lower PML.
			for (unsigned int j=1; j<PMLw+1; j++)
			{
				Bx(i,j,nf) = Bx(i,j,n0) + dt*((1./kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1./delta + PsiHxY(i,j));
				Hx(i,j,nf) = Bx(i,j,nf)/(u0*uinf(i,j));
			}
			// Bx in upper PML.
			for (unsigned int j=JHx-PMLw-1; j<JHx-1; j++)
			{
				Bx(i,j,nf) = Bx(i,j,n0) + dt*((1./kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1./delta + PsiHxY(i,j));
				Hx(i,j,nf) = Bx(i,j,nf)/(u0*uinf(i,j));
			}
		}
		
		// ========================== By and Hy ==========================
		for (unsigned int i=0; i<IHy-1; i++)
		{
			// PsiHyX arrays.
			for (unsigned int j=0; j<JHy; j++)
			{
				PsiHyX(i,j) = (Cmx/delta)*(Ez(i+1,j,n0)-Ez(i,j,n0)) + bmx*PsiHyX(i,j);
				if (i==0)
					PsiHyX(IHy-1,j) = (Cmx/delta)*(Ez(0,j,n0)-Ez(IHy-1,j,n0)) + bmx*PsiHyX(IHy-1,j); // PBC
			}
			// By in normal space.
			for (unsigned int j=PMLw; j<JHy-PMLw; j++)
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
			for (unsigned int j=0; j<PMLw; j++)
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
			for (unsigned int j=JHy-PMLw; j<JHy; j++)
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

		// ========================== Dz and Ez ==========================
		for (unsigned int i=1; i<IEz; i++)
		{
			// Psi arrays.
			for (unsigned int j=0; j<JEz; j++)
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
			for (unsigned int j=PMLw; j<JEz-PMLw; j++)
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
			for (unsigned int j=0; j<PMLw; j++)
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
			for (unsigned int j=JEz-PMLw; j<JEz; j++)
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
		if (SourcePlane == 1)
		{
			unsigned int j=SourceLocationY;
			for (unsigned int i=0; i<IEz; i++)
			{
				if (SourceChoice == 1)
					Ez(i,j,nf) += exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
				else if (SourceChoice == 2)
					Ez(i,j,nf) += sin(2.*pi*f*(PRECISION)n*dt) * Sc;
				else if (SourceChoice == 3)
					Ez(i,j,nf) += (1.-2.*pow(pi*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(pi*fp*((PRECISION)n*dt-dr),2)) * Sc;

				Dz(i,j,nf) = e0*Ez(i,j,nf);
			}
		}
		else
		{
			unsigned int i=SourceLocationX;
			unsigned int j=SourceLocationY;
			if (SourceChoice == 1)
				Ez(i,j,nf) += exp(-1.*pow(((PRECISION)n-(PRECISION)td)/((PRECISION)PulseWidth/4.),2)) * Sc;
			else if (SourceChoice == 2)
				Ez(i,j,nf) += sin(2.*pi*f*(PRECISION)n*dt) * Sc;
			else if (SourceChoice == 3)
				Ez(i,j,nf) += (1.-2.*pow(pi*fp*((PRECISION)n*dt-dr),2))*exp(-1.*pow(pi*fp*((PRECISION)n*dt-dr),2)) * Sc;

			Dz(i,j,nf) = e0*Ez(i,j,nf);
		}
		// Recording transmitted fields.
		Ezt[n] = Ez(IEz/2,x1,nf);
		Eztt[n] = Ez(IEz/2,SlabRight+10,nf);
		// Fields for refractive index.
		Ezy1[n] = Ez(IEz/2,Y1,nf);
		Ezy2[n] = Ez(IEz/2,Y2,nf);

		// Saving electric field snapshot.
		if (n%SnapshotInterval == 0 && SaveFields == true)
		{
			// Write E-field to file.
			framestream.str(std::string());			// Clearing stringstream contents.
			framestream << ++frame;
			filename = basename + framestream.str() + ".fdt";
			snapshot.open(filename.c_str(), std::ios::out|std::ios::binary);
			snapshot.write((char*)&(Ez(0,0,nf)), sizeof(PRECISION)*IEz*JEz);
			snapshot.close();
		}

		np = (np+1)%3;
		n0 = (n0+1)%3;
		nf = (nf+1)%3;
	}
	// Saving electric field snapshot.
	if (SaveFields == true)
	{
		fstream parametersfile;
		parametersfile.open("FieldData/Parameters.smp", std::ios::out|std::ios::binary|std::ios::app);
		parametersfile.write((char*)&(frame), sizeof(unsigned int));
		parametersfile.close();
		// Write saved fields to files.
		snapshot.open("FieldData/Ezi.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Ezi, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ezt.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Ezt, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Eztt.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Eztt, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ezy1.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Ezy1, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ezy2.fdt", std::ios::out|std::ios::binary);
		snapshot.write((char*)Ezy2, sizeof(PRECISION)*MaxTime);
		snapshot.close();
	}
	cout << endl << "Simulation (CPU) completed!" << endl;

	return 0;
}
int CFDTD2DDNG::DryRunGPU()
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

	globalThreads[0] = I;
	globalThreads[1] = J+2*PMLw;
	localThreads[0]  = 16;
	localThreads[1]  = 16;

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
		if (n % (MaxTime/1024U) == 0)
			cout << "\r\t\t\r" << n*100/(MaxTime-1) << "%";

		SafeCall(clSetKernelArg(DryRun_kernel_M, 51, sizeof(unsigned int), (void *)&n), "Error: Setting kernel argument 'n'");
		SafeCall(clSetKernelArg(DryRun_kernel_M, 52, sizeof(unsigned int), (void *)&np), "Error: Setting kernel argument 'np'");
		SafeCall(clSetKernelArg(DryRun_kernel_M, 53, sizeof(unsigned int), (void *)&n0), "Error: Setting kernel argument 'n0'");
		SafeCall(clSetKernelArg(DryRun_kernel_M, 54, sizeof(unsigned int), (void *)&nf), "Error: Setting kernel argument 'nf'");

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
		kernelExecTimeNs = 1e-3*(endTime-startTime);
		kernelExecTimeNsT = kernelExecTimeNsT + kernelExecTimeNs;

		SafeCall(clSetKernelArg(DryRun_kernel_E, 51, sizeof(unsigned int), (void *)&n), "Error: Setting kernel argument 'n'");
		SafeCall(clSetKernelArg(DryRun_kernel_E, 52, sizeof(unsigned int), (void *)&np), "Error: Setting kernel argument 'np'");
		SafeCall(clSetKernelArg(DryRun_kernel_E, 53, sizeof(unsigned int), (void *)&n0), "Error: Setting kernel argument 'n0'");
		SafeCall(clSetKernelArg(DryRun_kernel_E, 54, sizeof(unsigned int), (void *)&nf), "Error: Setting kernel argument 'nf'");

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
		kernelExecTimeNs = 1e-3*(endTime-startTime);
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
int CFDTD2DDNG::RunSimulationGPU(bool SaveFields)
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

	globalThreads[0] = I;
	globalThreads[1] = J+2*PMLw;
	localThreads[0]  = 16;
	localThreads[1]  = 16;

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
	string basename = "FieldData/Ez";
	string filename;
	fstream snapshot;
	frame = 0U;

	cout << "Simulation run (GPU) started..." << endl;
	cout << "Global threads: " << globalThreads[0] << "x" << globalThreads[1] << endl;
	cout << "Local threads: " << localThreads[0] << "x" << localThreads[1] << endl;

	for (unsigned int n=0;n<MaxTime; n++)
	{
		if (n % (MaxTime/1024U) == 0)
			cout << "\r\t\t\r" << n*100/(MaxTime-1) << "%";

		SafeCall(clSetKernelArg(Simulation_kernel_M, 72, sizeof(unsigned int), (void *)&n), "Error: Setting kernel argument 'n'");
		SafeCall(clSetKernelArg(Simulation_kernel_M, 73, sizeof(unsigned int), (void *)&np), "Error: Setting kernel argument 'np'");
		SafeCall(clSetKernelArg(Simulation_kernel_M, 74, sizeof(unsigned int), (void *)&n0), "Error: Setting kernel argument 'n0'");
		SafeCall(clSetKernelArg(Simulation_kernel_M, 75, sizeof(unsigned int), (void *)&nf), "Error: Setting kernel argument 'nf'");

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
		kernelExecTimeNs = 1e-3*(endTime-startTime);
		kernelExecTimeNsT = kernelExecTimeNsT + kernelExecTimeNs;

		SafeCall(clSetKernelArg(Simulation_kernel_E, 72, sizeof(unsigned int), (void *)&n), "Error: Setting kernel argument 'n'");
		SafeCall(clSetKernelArg(Simulation_kernel_E, 73, sizeof(unsigned int), (void *)&np), "Error: Setting kernel argument 'np'");
		SafeCall(clSetKernelArg(Simulation_kernel_E, 74, sizeof(unsigned int), (void *)&n0), "Error: Setting kernel argument 'n0'");
		SafeCall(clSetKernelArg(Simulation_kernel_E, 75, sizeof(unsigned int), (void *)&nf), "Error: Setting kernel argument 'nf'");

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
		kernelExecTimeNs = 1e-3*(endTime-startTime);
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
			SafeCall(clEnqueueReadBuffer(commandQueue, d_Ez_, CL_TRUE, 0, sizeof(PRECISION)*IEz*JEz*3, Ez_, 0, NULL, &events[1]), "Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)");
			// Wait for the read buffer to finish execution
			SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for read buffer call to finish. (clWaitForEvents)");

			snapshot.write((char*)&(Ez(0,0,nf)), sizeof(PRECISION)*IEz*JEz);
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
		snapshot.open("FieldData/Ezi.fdt", std::ios::out|std::ios::binary);
		SafeCall(clEnqueueReadBuffer(commandQueue, d_Ezi, CL_TRUE, 0, sizeof(PRECISION)*MaxTime, Ezi, 0, NULL, &events[1]), "Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)");
		SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for read buffer call to finish. (clWaitForEvents)");
		snapshot.write((char*)Ezi, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ezt.fdt", std::ios::out|std::ios::binary);
		SafeCall(clEnqueueReadBuffer(commandQueue, d_Ezt, CL_TRUE, 0, sizeof(PRECISION)*MaxTime, Ezt, 0, NULL, &events[1]), "Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)");
		SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for read buffer call to finish. (clWaitForEvents)");
		snapshot.write((char*)Ezt, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Eztt.fdt", std::ios::out|std::ios::binary);
		SafeCall(clEnqueueReadBuffer(commandQueue, d_Eztt, CL_TRUE, 0, sizeof(PRECISION)*MaxTime, Eztt, 0, NULL, &events[1]), "Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)");
		SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for read buffer call to finish. (clWaitForEvents)");
		snapshot.write((char*)Eztt, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ezy1.fdt", std::ios::out|std::ios::binary);
		SafeCall(clEnqueueReadBuffer(commandQueue, d_Ezy1, CL_TRUE, 0, sizeof(PRECISION)*MaxTime, Ezy1, 0, NULL, &events[1]), "Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)");
		SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for read buffer call to finish. (clWaitForEvents)");
		snapshot.write((char*)Ezy1, sizeof(PRECISION)*MaxTime);
		snapshot.close();
		snapshot.open("FieldData/Ezy2.fdt", std::ios::out|std::ios::binary);
		SafeCall(clEnqueueReadBuffer(commandQueue, d_Ezy2, CL_TRUE, 0, sizeof(PRECISION)*MaxTime, Ezy2, 0, NULL, &events[1]), "Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)");
		SafeCall(clWaitForEvents(1, &events[1]), "Error: Waiting for read buffer call to finish. (clWaitForEvents)");
		snapshot.write((char*)Ezy2, sizeof(PRECISION)*MaxTime);
		snapshot.close();

		SafeCall(clReleaseEvent(events[1]), "Error: Release event object. (clReleaseEvent)\n");
	}

	return 0;
}
int CFDTD2DDNG::CompleteRunCPU(bool SaveFields)
{
	cout << "Memory required for simulation = " << SimSize() << " bytes (" << (double)SimSize()/1024UL << "kB/" << (double)SimSize()/1024UL/1024UL << "MB)." << endl;
	cout << "HDD space required for data storage = " << HDDSpace() << " bytes (" << (double)HDDSpace()/1024UL << "kB/" << (double)HDDSpace()/1024UL/1024UL << "MB)." << endl;
	SafeCall(AllocateMemoryCPU(), "Error: Allocating memory on CPU.");
	SafeCall(InitialiseCPU(), "Error: Initialising data on CPU.");
	SafeCall(DryRunCPU(), "Error: Dry run (CPU).");
	SafeCall(InitialiseForSimulationCPU(), "Error: Initalising for simulation (CPU).");
	SafeCall(RunSimulationCPU(SaveFields), "Error: Running simulation (CPU).");
	SafeCall(CleanupCPU(), "Error: Cleaning up CPU.");

	return 0;
}
int CFDTD2DDNG::CompleteRunGPU(bool SaveFields)
{
	cout << "Memory required for simulation = " << SimSize() << " bytes (" << (double)SimSize()/1024UL << "kB/" << (double)SimSize()/1024UL/1024UL << "MB)." << endl;
	cout << "HDD space required for data storage = " << HDDSpace() << " bytes (" << (double)HDDSpace()/1024UL << "kB/" << (double)HDDSpace()/1024UL/1024UL << "MB)." << endl;
	SafeCall(AllocateMemoryCPU(), "Error: Allocating memory on CPU.");
	SafeCall(InitialiseCPU(), "Error: Initialising data on CPU.");

	SafeCall(InitialiseCL(), "Error: Initialiasing CL.");
	SafeCall(AllocateMemoryGPU(), "Error: Allocating memory on GPU.");
	SafeCall(InitialiseCLKernelsGPU(), "Error: Copying data from CPU to GPU.");
	SafeCall(DryRunGPU(), "Error: Dry run (GPU).");
	SafeCall(InitialiseForSimulationCPU(), "Error: Initalising for simulation (CPU).");
	SafeCall(InitialiseForSimulationGPU(), "Error: Initalising for simulation (GPU).");
	SafeCall(RunSimulationGPU(SaveFields), "Error: Running simulation on GPU.");
	SafeCall(CleanupCPU(), "Error: Cleaning up CPU.");
	SafeCall(CleanupCL(), "Error: Cleaning up CL.");
	SafeCall(CleanupGPU(), "Error: Cleaning up GPU.");

	return 0;
}
// Converts contents of a file into a string. From OPENCL examples.
string CFDTD2DDNG::convertToString(const char *filename)
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
void CFDTD2DDNG::StartTimer()
{
	tStart = GetTimeus64();
}
void CFDTD2DDNG::StopTimer()
{
	tEnd = GetTimeus64();
}
PRECISION CFDTD2DDNG::GetElapsedTime()
{
	return ((PRECISION)(tEnd-tStart))/(1000000.);
}
int CFDTD2DDNG::SafeCall(int Status, const char *Error)
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
int CFDTD2DDNG::CleanupCPU()
{
	// Field arrays.
	DeleteArray(Ez_);
	DeleteArray(Dz_);
	DeleteArray(Hx_);
	DeleteArray(Bx_);
	DeleteArray(Hy_);
	DeleteArray(By_);
	// Incident and transmitted fields.
	DeleteArray(Ezi);
	DeleteArray(Ezt);
	DeleteArray(Eztt);
	// Refractive index.
	DeleteArray(Ezy1);
	DeleteArray(Ezy2);
	// Drude parameter arrays.
	DeleteArray(einf_);
	DeleteArray(uinf_);
	DeleteArray(wpesq_);
	DeleteArray(wpmsq_);
	DeleteArray(ge_);
	DeleteArray(gm_);
	// Auxiliary field scalars.
	DeleteArray(ae0_);
	DeleteArray(ae_);
	DeleteArray(be_);
	DeleteArray(ce_);
	DeleteArray(de_);
	DeleteArray(ee_);
	DeleteArray(am0_);
	DeleteArray(am_);
	DeleteArray(bm_);
	DeleteArray(cm_);
	DeleteArray(dm_);
	DeleteArray(em_);
	// PML arrays.
	DeleteArray(PsiEzX_);
	DeleteArray(PsiEzY_);
	DeleteArray(PsiHyX_);
	DeleteArray(PsiHxY_);

	return 0;
}
int CFDTD2DDNG::CleanupCL()
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
int CFDTD2DDNG::CleanupGPU()
{
	SafeCall(clReleaseMemObject(d_Ez_), "Error: clReleaseMemObject() cannot release memory buffer for d_Ez_");
	SafeCall(clReleaseMemObject(d_Dz_), "Error: clReleaseMemObject() cannot release memory buffer for d_Dz_");
	SafeCall(clReleaseMemObject(d_Hx_), "Error: clReleaseMemObject() cannot release memory buffer for d_Hx_");
	SafeCall(clReleaseMemObject(d_Bx_), "Error: clReleaseMemObject() cannot release memory buffer for d_Bx_");
	SafeCall(clReleaseMemObject(d_Hy_), "Error: clReleaseMemObject() cannot release memory buffer for d_Hy_");
	SafeCall(clReleaseMemObject(d_By_), "Error: clReleaseMemObject() cannot release memory buffer for d_By_");

	SafeCall(clReleaseMemObject(d_Ezi), "Error: clReleaseMemObject() cannot release memory buffer for d_Ezi");
	SafeCall(clReleaseMemObject(d_Ezt), "Error: clReleaseMemObject() cannot release memory buffer for d_Ezt");
	SafeCall(clReleaseMemObject(d_Eztt), "Error: clReleaseMemObject() cannot release memory buffer for d_Eztt");
	SafeCall(clReleaseMemObject(d_Ezy1), "Error: clReleaseMemObject() cannot release memory buffer for d_Ezy1");
	SafeCall(clReleaseMemObject(d_Ezy2), "Error: clReleaseMemObject() cannot release memory buffer for d_Ezy2");

	SafeCall(clReleaseMemObject(d_einf_), "Error: clReleaseMemObject() cannot release memory buffer for d_einf");
	SafeCall(clReleaseMemObject(d_uinf_), "Error: clReleaseMemObject() cannot release memory buffer for d_uinf");
	SafeCall(clReleaseMemObject(d_wpesq_), "Error: clReleaseMemObject() cannot release memory buffer for d_wpesq");
	SafeCall(clReleaseMemObject(d_wpmsq_), "Error: clReleaseMemObject() cannot release memory buffer for d_wpmsq");
	SafeCall(clReleaseMemObject(d_ge_), "Error: clReleaseMemObject() cannot release memory buffer for d_ge");
	SafeCall(clReleaseMemObject(d_gm_), "Error: clReleaseMemObject() cannot release memory buffer for d_gm");

	SafeCall(clReleaseMemObject(d_ae0_), "Error: clReleaseMemObject() cannot release memory buffer for d_ae0");
	SafeCall(clReleaseMemObject(d_ae_), "Error: clReleaseMemObject() cannot release memory buffer for d_ae");
	SafeCall(clReleaseMemObject(d_be_), "Error: clReleaseMemObject() cannot release memory buffer for d_be");
	SafeCall(clReleaseMemObject(d_ce_), "Error: clReleaseMemObject() cannot release memory buffer for d_ce");
	SafeCall(clReleaseMemObject(d_de_), "Error: clReleaseMemObject() cannot release memory buffer for d_de");
	SafeCall(clReleaseMemObject(d_ee_), "Error: clReleaseMemObject() cannot release memory buffer for d_ee");
	SafeCall(clReleaseMemObject(d_am0_), "Error: clReleaseMemObject() cannot release memory buffer for d_am0");
	SafeCall(clReleaseMemObject(d_am_), "Error: clReleaseMemObject() cannot release memory buffer for d_am");
	SafeCall(clReleaseMemObject(d_bm_), "Error: clReleaseMemObject() cannot release memory buffer for d_bm");
	SafeCall(clReleaseMemObject(d_cm_), "Error: clReleaseMemObject() cannot release memory buffer for d_cm");
	SafeCall(clReleaseMemObject(d_dm_), "Error: clReleaseMemObject() cannot release memory buffer for d_dm");
	SafeCall(clReleaseMemObject(d_em_), "Error: clReleaseMemObject() cannot release memory buffer for d_em");

	SafeCall(clReleaseMemObject(d_PsiEzX_), "Error: clReleaseMemObject() cannot release memory input buffer for d_PsiEzX_");
	SafeCall(clReleaseMemObject(d_PsiEzY_), "Error: clReleaseMemObject() cannot release memory input buffer for d_PsiEzY_");
	SafeCall(clReleaseMemObject(d_PsiHyX_), "Error: clReleaseMemObject() cannot release memory input buffer for d_PsiHyX_");
	SafeCall(clReleaseMemObject(d_PsiHxY_), "Error: clReleaseMemObject() cannot release memory input buffer for d_PsiHxY_");

	return 0;
}
CFDTD2DDNG::~CFDTD2DDNG ()
{
	// Field arrays.
	DeleteArray(Ez_);
	DeleteArray(Dz_);
	DeleteArray(Hx_);
	DeleteArray(Bx_);
	DeleteArray(Hy_);
	DeleteArray(By_);
	// Incident and transmitted fields.
	DeleteArray(Ezi);
	DeleteArray(Ezt);
	DeleteArray(Eztt);
	// Refractive index.
	DeleteArray(Ezy1);
	DeleteArray(Ezy2);
	// Drude parameter arrays.
	DeleteArray(einf_);
	DeleteArray(uinf_);
	DeleteArray(wpesq_);
	DeleteArray(wpmsq_);
	DeleteArray(ge_);
	DeleteArray(gm_);
	// Auxiliary field scalars.
	DeleteArray(ae0_);
	DeleteArray(ae_);
	DeleteArray(be_);
	DeleteArray(ce_);
	DeleteArray(de_);
	DeleteArray(ee_);
	DeleteArray(am0_);
	DeleteArray(am_);
	DeleteArray(bm_);
	DeleteArray(cm_);
	DeleteArray(dm_);
	DeleteArray(em_);
	// PML arrays.
	DeleteArray(PsiEzX_);
	DeleteArray(PsiEzY_);
	DeleteArray(PsiHyX_);
	DeleteArray(PsiHxY_);
}
