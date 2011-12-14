/* ============================================================

Copyright (c) 2009 Advanced Micro Devices, Inc.  All rights reserved.
 
Redistribution and use of this material is permitted under the following 
conditions:
 
Redistributions must retain the above copyright notice and all terms of this 
license.
 
In no event shall anyone redistributing or accessing or using this material 
commence or participate in any arbitration or legal action relating to this 
material against Advanced Micro Devices, Inc. or any copyright holders or 
contributors. The foregoing shall survive any expiration or termination of 
this license or any agreement or access or use related to this material. 

ANY BREACH OF ANY TERM OF THIS LICENSE SHALL RESULT IN THE IMMEDIATE REVOCATION 
OF ALL RIGHTS TO REDISTRIBUTE, ACCESS OR USE THIS MATERIAL.

THIS MATERIAL IS PROVIDED BY ADVANCED MICRO DEVICES, INC. AND ANY COPYRIGHT 
HOLDERS AND CONTRIBUTORS "AS IS" IN ITS CURRENT CONDITION AND WITHOUT ANY 
REPRESENTATIONS, GUARANTEE, OR WARRANTY OF ANY KIND OR IN ANY WAY RELATED TO 
SUPPORT, INDEMNITY, ERROR FREE OR UNINTERRUPTED OPERA TION, OR THAT IT IS FREE 
FROM DEFECTS OR VIRUSES.  ALL OBLIGATIONS ARE HEREBY DISCLAIMED - WHETHER 
EXPRESS, IMPLIED, OR STATUTORY - INCLUDING, BUT NOT LIMITED TO, ANY IMPLIED 
WARRANTIES OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, 
ACCURACY, COMPLETENESS, OPERABILITY, QUALITY OF SERVICE, OR NON-INFRINGEMENT. 
IN NO EVENT SHALL ADVANCED MICRO DEVICES, INC. OR ANY COPYRIGHT HOLDERS OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, PUNITIVE,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, REVENUE, DATA, OR PROFITS; OR 
BUSINESS INTERRUPTION) HOWEVER CAUSED OR BASED ON ANY THEORY OF LIABILITY 
ARISING IN ANY WAY RELATED TO THIS MATERIAL, EVEN IF ADVISED OF THE POSSIBILITY 
OF SUCH DAMAGE. THE ENTIRE AND AGGREGATE LIABILITY OF ADVANCED MICRO DEVICES, 
INC. AND ANY COPYRIGHT HOLDERS AND CONTRIBUTORS SHALL NOT EXCEED TEN DOLLARS 
(US $10.00). ANYONE REDISTRIBUTING OR ACCESSING OR USING THIS MATERIAL ACCEPTS 
THIS ALLOCATION OF RISK AND AGREES TO RELEASE ADVANCED MICRO DEVICES, INC. AND 
ANY COPYRIGHT HOLDERS AND CONTRIBUTORS FROM ANY AND ALL LIABILITIES, 
OBLIGATIONS, CLAIMS, OR DEMANDS IN EXCESS OF TEN DOLLARS (US $10.00). THE 
FOREGOING ARE ESSENTIAL TERMS OF THIS LICENSE AND, IF ANY OF THESE TERMS ARE 
CONSTRUED AS UNENFORCEABLE, FAIL IN ESSENTIAL PURPOSE, OR BECOME VOID OR 
DETRIMENTAL TO ADVANCED MICRO DEVICES, INC. OR ANY COPYRIGHT HOLDERS OR 
CONTRIBUTORS FOR ANY REASON, THEN ALL RIGHTS TO REDISTRIBUTE, ACCESS OR USE 
THIS MATERIAL SHALL TERMINATE IMMEDIATELY. MOREOVER, THE FOREGOING SHALL 
SURVIVE ANY EXPIRATION OR TERMINATION OF THIS LICENSE OR ANY AGREEMENT OR 
ACCESS OR USE RELATED TO THIS MATERIAL.

NOTICE IS HEREBY PROVIDED, AND BY REDISTRIBUTING OR ACCESSING OR USING THIS 
MATERIAL SUCH NOTICE IS ACKNOWLEDGED, THAT THIS MATERIAL MAY BE SUBJECT TO 
RESTRICTIONS UNDER THE LAWS AND REGULATIONS OF THE UNITED STATES OR OTHER 
COUNTRIES, WHICH INCLUDE BUT ARE NOT LIMITED TO, U.S. EXPORT CONTROL LAWS SUCH 
AS THE EXPORT ADMINISTRATION REGULATIONS AND NATIONAL SECURITY CONTROLS AS 
DEFINED THEREUNDER, AS WELL AS STATE DEPARTMENT CONTROLS UNDER THE U.S. 
MUNITIONS LIST. THIS MATERIAL MAY NOT BE USED, RELEASED, TRANSFERRED, IMPORTED,
EXPORTED AND/OR RE-EXPORTED IN ANY MANNER PROHIBITED UNDER ANY APPLICABLE LAWS, 
INCLUDING U.S. EXPORT CONTROL LAWS REGARDING SPECIFICALLY DESIGNATED PERSONS, 
COUNTRIES AND NATIONALS OF COUNTRIES SUBJECT TO NATIONAL SECURITY CONTROLS. 
MOREOVER, THE FOREGOING SHALL SURVIVE ANY EXPIRATION OR TERMINATION OF ANY 
LICENSE OR AGREEMENT OR ACCESS OR USE RELATED TO THIS MATERIAL.

NOTICE REGARDING THE U.S. GOVERNMENT AND DOD AGENCIES: This material is 
provided with "RESTRICTED RIGHTS" and/or "LIMITED RIGHTS" as applicable to 
computer software and technical data, respectively. Use, duplication, 
distribution or disclosure by the U.S. Government and/or DOD agencies is 
subject to the full extent of restrictions in all applicable regulations, 
including those found at FAR52.227 and DFARS252.227 et seq. and any successor 
regulations thereof. Use of this material by the U.S. Government and/or DOD 
agencies is acknowledgment of the proprietary rights of any copyright holders 
and contributors, including those of Advanced Micro Devices, Inc., as well as 
the provisions of FAR52.227-14 through 23 regarding privately developed and/or 
commercial computer software.

This license forms the entire agreement regarding the subject matter hereof and 
supersedes all proposals and prior discussions and writings between the parties 
with respect thereto. This license does not affect any ownership, rights, title,
or interest in, or relating to, this material. No terms of this license can be 
modified or waived, and no breach of this license can be excused, unless done 
so in a writing signed by all affected parties. Each term of this license is 
separately enforceable. If any term of this license is determined to be or 
becomes unenforceable or illegal, such term shall be reformed to the minimum 
extent necessary in order for this license to remain in effect in accordance 
with its terms as modified by such reformation. This license shall be governed 
by and construed in accordance with the laws of the State of Texas without 
regard to rules on conflicts of law of any state or jurisdiction or the United 
Nations Convention on the International Sale of Goods. All disputes arising out 
of this license shall be subject to the jurisdiction of the federal and state 
courts in Austin, Texas, and all defenses are hereby waived concerning personal 
jurisdiction and venue of these courts.

============================================================ */
#ifndef FDTD2D_H_
#define FDTD2D_H_

#ifdef WIN32
#include <fstream>
#else
#include <fcntl.h>
#endif

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <string>
#include <CL/cl.h>
/*
#include <SDKCommon.hpp>
#include <SDKApplication.hpp>
#include <SDKCommandArgs.hpp>
#include <SDKFile.hpp>
*/
class CFDTD2D
{
private:

	// Generic simulation parameters.
	const cl_uint I;				// Width.
	const cl_uint J;				// Height.
	const cl_double c;				// Speed of light.
	const cl_double delta;			// dx and dy.
	const cl_double dx;			// dx if being used.
	const cl_double dy;			// dy if being used.
	const cl_double dtscalar;		// dt scale factor. dt will be divided by this factor.
	const cl_double dt;			// dt.
	const cl_uint PMLw;			// Width of PML layer.

	const cl_uint NMax;			// Maximum n

	// Constants.
	const cl_double f;				// frequency
	const cl_double pi;			// pi
	const cl_double e0;			// epsilon nought
	const cl_double u0;			// mu nought

	// miscellenaeous
	const cl_double Two_pi_f_deltat;
	const cl_uint NHW;			// Half wave cycle.
	const cl_uint Js;				// J position of plane wave front.
	const cl_uint Is;				// I position of plane wave front.
	cl_uint n, n0, n1, n2;			// past/present/future time indices.
	const cl_uint tResolution;		// Snapshots will be saved after this much time.
	const cl_uint xResolution;		// Resolution of plotted field will be divided by this factor.
	const cl_uint yResolution;

	// TMz parameters.
	const cl_uint IHx;
	const cl_uint JHx;
	const cl_uint IHy;
	const cl_uint JHy;
	const cl_uint IEz;
	const cl_uint JEz;

	// Geometry parameters.
	const cl_uint XCenter;
	const cl_uint YCenter;
	const cl_double ra;
	const cl_double rb;
	const cl_double imp0;			// Impedence of free space

	// =========== Data Arrays ==========
	cl_double *Hx;					// Hx, magnetic field.
	cl_double *Bx;
	cl_double *Hy;					// Hy, magnetic field.
	cl_double *By;

	cl_double *Ez;					// Ez, electric field.
	cl_double *Dz;

	cl_double *Dzx;
	cl_double *Dzy;
	cl_double *EzSnapshots;		// For storing Ez snapshots.

	// ========= Field specific arrays =======

	// Permeability and permittivity.
	cl_double *urHx;
	cl_double *urHy;
	cl_double *erEz;

	// Magnetic and electric conductances.
	cl_double *smHx;
	cl_double *smHy;
	cl_double *sEz;

	// s*dt/(2*er) and sm*dt/(2*ur)
	cl_double *Sc;
	cl_double *ScmHx;
	cl_double *ScmHy;

	// PML conductance arrays.
	cl_double *sex;	// sigma ex
	cl_double *sey;	// sigma ey
	cl_double *smx;	// sigma mx
	cl_double *smy;	// sigma my

	cl_double *Scsx;
	cl_double *Scsy;
	cl_double *ScmsmxHy;
	cl_double *ScmsmyHx;

	// Timing.
	clock_t tStart;
	clock_t tEnd;
	clock_t tElapsed;

	// === OPENCL memory buffers ===
	// Fields
	cl_mem inputBufferHx;
	cl_mem inputBufferBx;
	cl_mem inputBufferHy;
	cl_mem inputBufferBy;
	cl_mem inputBufferEz;
	cl_mem inputBufferDz;
	cl_mem inputBufferDzx;
	cl_mem inputBufferDzy;

	// Permeability, permittivity and conductance.
	cl_mem inputBufferurHx;
	cl_mem inputBufferurHy;
	cl_mem inputBuffererEz;

	cl_mem inputBufferScmHx;
	cl_mem inputBufferScmHy;
	cl_mem inputBufferSc;
	// =============================

	// OPENCL related parameters.
	bool cpu;		// Is OPENCL using CPU or GPU?
	cl_uint flagHalf;

	// OPENCL context/device/program
	cl_context context;
	cl_device_id *devices;
	cl_command_queue commandQueue;
	cl_program program;

	/* This program uses only one kernel and this serves as a handle to it */
	cl_kernel  kernel;

public:
	CFDTD2D ();
	int Initialize ();						// Initialize with default parameters.
	int initializeCL ();
	int initializeFDTD2DKernel ();
	int runCLFDTD2DKernels (bool=true);
	int RunSimulationCPU (bool=true);		// Run simulation on CPU single-threaded.

	void DisplaySimulationParameters ();
	std::string convertToString(const char * filename);
	// Timing.
	inline void StartClock () { tStart = clock(); }
	inline void StopClock () { tEnd = clock(); }
	inline cl_double GetElapsedTime () { return (cl_double)(tEnd-tStart)/CLOCKS_PER_SEC; }

	int CleanupCL ();
	~CFDTD2D ();
};


#endif  /* #ifndef FDTD2D_H_ */
