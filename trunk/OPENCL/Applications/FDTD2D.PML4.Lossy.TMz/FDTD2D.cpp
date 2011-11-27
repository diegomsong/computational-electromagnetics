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


#include "FDTD2D.hpp"

CFDTD2D::CFDTD2D () : 
						I(200),
						J(200),
						c(299792458.),
						delta(2.5e-3),
						dx(delta),
						dy(delta),
						dtscalar(2.),
						dt(delta/(sqrt(2.)*c) /dtscalar),
						PMLw(0),
						NMax(256),
						f(2.e9),
						pi(4*atan(1.)),
						e0(1.e-9/(36.*pi)),
						u0(1.e-7*4.*pi),
						Two_pi_f_deltat(2*pi*f*dt),
						NHW(1./(2.*f*dt)),
						Js(20+PMLw),
						Is(2),
						n0(0),
						n1(1),
						n2(2),
						tResolution(4),
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
						ScmsmyHx(NULL)
{
	cl_double bytes = I*J*17*3*8+50*8;		// Dynamic arrays + predefined variables.
	std:: cout << "Approximate memory required for simulation: " << bytes/1024 << " Kbytes (" << bytes/(1024*1024) << " MB)." << std::endl;

	// Writing simulation parameters for plotting.
	#ifdef WIN32
	std::fstream parametersfile;
	parametersfile.open ("../../FieldData/Parameters.smp", std::ios::out|std::ios::binary);
	parametersfile.write ((char*)&I, sizeof(cl_uint));
	parametersfile.write ((char*)&J, sizeof(cl_uint));
	parametersfile.write ((char*)&tResolution, sizeof(cl_uint));
	parametersfile.write ((char*)&xResolution, sizeof(cl_uint));
	parametersfile.write ((char*)&yResolution, sizeof(cl_uint));
	parametersfile.write ((char*)&NMax, sizeof(cl_uint));
	parametersfile.close ();
	#else
	int fdp;
	fdp = open ( "../../FieldData/Parameters.smp", O_CREAT|O_WRONLY|O_TRUNC, S_IRWXU );
	write (fdp, (void*)&I, sizeof(cl_uint));
	write (fdp, (void*)&J, sizeof(cl_uint));
	write (fdp, (void*)&tResolution, sizeof(cl_uint));
	write (fdp, (void*)&xResolution, sizeof(cl_uint));
	write (fdp, (void*)&yResolution, sizeof(cl_uint));
	write (fdp, (void*)&NMax, sizeof(cl_uint));
	close (fdp);
	#endif
}
// Initialize data arrays.
int CFDTD2D::Initialize ()
{
	Hx = new cl_double[IHx*JHx*3];
	Bx = new cl_double[IHx*JHx*3];
	Hy = new cl_double[IHy*JHy*3];
	By = new cl_double[IHy*JHy*3];
	Ez = new cl_double[IEz*JEz*3];
	Dz = new cl_double[IEz*JEz*3];
	Dzx = new cl_double[IEz*JEz*3];
	Dzy = new cl_double[IEz*JEz*3];
	//EzSnapshots = new cl_double[(IEz/xResolution)*(JEz/yResolution)];

	urHx = new cl_double[IHx*JHx];
	urHy = new cl_double[IHy*JHy];
	erEz = new cl_double[IEz*JEz];
	smHx = new cl_double[IHx*JHx];
	smHy = new cl_double[IHy*JHy];
	sEz = new cl_double[IEz*JEz];
	Sc = new cl_double[IEz*JEz];
	ScmHx = new cl_double[IHx*JHx];
	ScmHy = new cl_double[IHy*JHy];

	sex = new cl_double[IEz*JEz];
	sey = new cl_double[IEz*JEz];
	smx = new cl_double[IHy*JHy];
	smy = new cl_double[IHx*JHx];
	Scsx = new cl_double[IEz*JEz];
	Scsy = new cl_double[IEz*JEz];
	ScmsmxHy = new cl_double[IHy*JHy];
	ScmsmyHx = new cl_double[IHx*JHx];

	cl_uint i, j, t;
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
						if (j < PMLw+1)
						{
							smy[i+IHx*j] = 1.7e10;
							if (J < PMLw)
								smy[i+IHx*(JHx-PMLw+j-1)] = 1.7e10;
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
						if (j < PMLw+1)
						{
							sey[i+IEz*j] = 1.7e10;
							if (J < PMLw)
								sey[i+IEz*(JEz-PMLw+j-1)] = 1.7e10;
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

int CFDTD2D::initializeCL ()
{

	return 0;
}

int CFDTD2D::initializeFDTD2DKernel ()
{
	return 0;
}

int CFDTD2D::RunSimulationCPU (bool SaveFields)
{
	// Time loop.
	cl_uint n, i, j;

	// File Handling.
	std::stringstream framestream;
	std::string basename = "../../FieldData/Ez";
	std::string filename;
	cl_uint frame = 1;

	#ifdef WIN32
	std::fstream snapshot;
	#else
	int fd;
	#endif

	cl_uint ProgressResolution;
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
			snapshot.write ( (char*)&(Ez[IEz*JEz*n2]), sizeof(cl_double)*IEz*JEz);
			snapshot.close();
			#else
			fd = open (filename.c_str(), O_CREAT|O_WRONLY|O_TRUNC, S_IRWXU);
			write (fd, (void*)&(Ez[IEz*JEz*n2]), sizeof(cl_double)*IEz*JEz);
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

// Converts contents of a file into a string. From OPENCL examples.
std::string CFDTD2D::convertToString(const char *filename)
{
	size_t size;
	char*  str;
	std::string s;
	std::fstream f(filename, (std::fstream::in | std::fstream::binary));

	if(f.is_open())
	{
		size_t fileSize;
		f.seekg(0, std::fstream::end);
		size = fileSize = (size_t)f.tellg();
		f.seekg(0, std::fstream::beg);

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
		std::cout << "\nFile containg the kernel code(\".cl\") not found. Please copy the required file in the folder containg the executable.\n";
		exit(1);
	}
	return NULL;
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

/*
 * \brief Host Initialization 
 *        Allocate and initialize memory 
 *        on the host. Print input array. 
 */
int
initializeHost(void)
{
	cpu			= false;	// Are we running on CPU or GPU?
	flagHalf	= 0;		// Flag to determine half time step.

	w			= 1024*1024;		// Width of 1D array.
	timeN		= 256;		// Total number of time steps.

    width		= 2*w;		// Size of Ez and Hy arrays.
    Ez			= NULL;
	Hy			= NULL;

	n0			= 0;
	n1			= 1;

	SaveFields	= false;	// Save fields?
	Binary		= true;	// Save fields as binary or text format?

    /////////////////////////////////////////////////////////////////
    // Allocate and initialize memory used by host 
    /////////////////////////////////////////////////////////////////
    cl_uint sizeInBytes = width * sizeof(cl_double);
    Ez = (cl_double *) malloc(sizeInBytes);
	Hy = (cl_double *) malloc(sizeInBytes);
    if(Ez == NULL || Hy == NULL)
    {
        std::cout<<"Error: Failed to allocate input memory on host\n";
        return 1; 
    }

    for(cl_uint i = 0; i < width; i++)
	{
		Ez[i] = 0;
		Hy[i] = 0;
	}

    // print input arrays
    //print1DArray(std::string("Input A").c_str(), Ez, width);
	//print1DArray(std::string("Input B").c_str(), Hy, width);
    return 0;
}

/*
 * Converts the contents of a file into a string
 */
std::string
convertToString(const char *filename)
{
    size_t size;
    char*  str;
    std::string s;

    std::fstream f(filename, (std::fstream::in | std::fstream::binary));

    if(f.is_open())
    {
        size_t fileSize;
        f.seekg(0, std::fstream::end);
        size = fileSize = (size_t)f.tellg();
        f.seekg(0, std::fstream::beg);

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
        std::cout << "\nFile containg the kernel code(\".cl\") not found. Please copy the required file in the folder containg the executable.\n";
        exit(1);
    }
    return NULL;
}

/*
 * \brief OpenCL related initialization 
 *        Create Context, Device list, Command Queue
 *        Create OpenCL memory buffer objects
 *        Load CL file, compile, link CL source 
 *		  Build program and kernel objects
 */
int
initializeCL(void)
{
    cl_int status = 0;
    size_t deviceListSize;
//	sampleCommon = new streamsdk::SDKCommon();
    /*
     * Have a look at the available platforms and pick either
     * the AMD one if available or a reasonable default.
     */

    cl_uint numPlatforms;
    cl_platform_id platform = NULL;
    status = clGetPlatformIDs(0, NULL, &numPlatforms);
    if(status != CL_SUCCESS)
    {
        std::cout << "Error: Getting Platforms. (clGetPlatformsIDs)\n";
        return 1;
    }
    
    if(numPlatforms > 0)
    {
        cl_platform_id* platforms = new cl_platform_id[numPlatforms];
        status = clGetPlatformIDs(numPlatforms, platforms, NULL);
        if(status != CL_SUCCESS)
        {
            std::cout << "Error: Getting Platform Ids. (clGetPlatformsIDs)\n";
            return 1;
        }
        for(unsigned int i=0; i < numPlatforms; ++i)
        {
            char pbuff[100];
            status = clGetPlatformInfo(
                        platforms[i],
                        CL_PLATFORM_VENDOR,
                        sizeof(pbuff),
                        pbuff,
                        NULL);
            if(status != CL_SUCCESS)
            {
                std::cout << "Error: Getting Platform Info.(clGetPlatformInfo)\n";
                return 1;
            }
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
		std::cout << "Running on CPU..." << std::endl;
		type = CL_DEVICE_TYPE_CPU;
	}
	else
	{
		std::cout << "Running on GPU..." << std::endl;
		type = CL_DEVICE_TYPE_GPU;
	}

    context = clCreateContextFromType(cps, 
                                      type, 
                                      NULL, 
                                      NULL, 
                                      &status);
    if(status != CL_SUCCESS) 
    {  
        std::cout<<"Error: Creating Context. (clCreateContextFromType)\n";
        return 1; 
    }

    /* First, get the size of device list data */
    status = clGetContextInfo(context, 
                              CL_CONTEXT_DEVICES, 
                              0, 
                              NULL, 
                              &deviceListSize);
    if(status != CL_SUCCESS) 
    {  
        std::cout<<
            "Error: Getting Context Info \
            (device list size, clGetContextInfo)\n";
        return 1;
    }

    /////////////////////////////////////////////////////////////////
    // Detect OpenCL devices
    /////////////////////////////////////////////////////////////////
    devices = (cl_device_id *)malloc(deviceListSize);
    if(devices == 0)
    {
        std::cout<<"Error: No devices found.\n";
        return 1;
    }

    /* Now, get the device list data */
    status = clGetContextInfo(
                 context, 
                 CL_CONTEXT_DEVICES, 
                 deviceListSize, 
                 devices, 
                 NULL);
    if(status != CL_SUCCESS) 
    { 
        std::cout<<
            "Error: Getting Context Info \
            (device list, clGetContextInfo)\n";
        return 1;
    }

    /////////////////////////////////////////////////////////////////
    // Create an OpenCL command queue
    /////////////////////////////////////////////////////////////////
    commandQueue = clCreateCommandQueue(
                       context, 
                       devices[0], 
                       CL_QUEUE_PROFILING_ENABLE, 
                       &status);
    if(status != CL_SUCCESS) 
    { 
        std::cout<<"Creating Command Queue. (clCreateCommandQueue)\n";
        return 1;
    }
	return 0;
}

int initializeFDTD1DKernel (void)
{
	cl_int status;
    /////////////////////////////////////////////////////////////////
    // Create OpenCL memory buffers
    /////////////////////////////////////////////////////////////////
    inputBufferEz = clCreateBuffer(
                      context, 
                      CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
                      sizeof(cl_double) * width,
                      Ez, 
                      &status);
    if(status != CL_SUCCESS) 
    { 
        std::cout<<"Error: clCreateBuffer (inputBufferEz)\n";
        return 1;
    }

	inputBufferHy = clCreateBuffer(
                      context, 
                      CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
                      sizeof(cl_double) * width,
                      Hy, 
                      &status);
    if(status != CL_SUCCESS) 
    { 
        std::cout<<"Error: clCreateBuffer (inputBufferHy)\n";
        return 1;
    }

    /////////////////////////////////////////////////////////////////
    // Load CL file, build CL program object, create CL kernel object
    /////////////////////////////////////////////////////////////////
    const char * filename  = "FDTD2D_Kernels.cl";
    std::string  sourceStr = convertToString(filename);
    const char * source    = sourceStr.c_str();
    size_t sourceSize[]    = { strlen(source) };

    program = clCreateProgramWithSource(
                  context, 
                  1, 
                  &source,
                  sourceSize,
                  &status);
    if(status != CL_SUCCESS) 
    { 
      std::cout<<
               "Error: Loading Binary into cl_program \
               (clCreateProgramWithBinary)\n";
      return 1;
    }

    /* create a cl program executable for all the devices specified */
    status = clBuildProgram(program, 1, devices, NULL, NULL, NULL);
    if(status != CL_SUCCESS) 
    { 
        std::cout<<"Error: Building Program (clBuildProgram)\n";
        return 1; 
    }

    /* get a kernel object handle for a kernel with the given name */
    kernel = clCreateKernel(program, "FDTD2DKernel", &status);
    if(status != CL_SUCCESS) 
    {  
        std::cout<<"Error: Creating Kernel from program. (clCreateKernel)\n";
        return 1;
    }

    return 0;
}


/*
 * \brief Run OpenCL program 
 *		  
 *        Bind host variables to kernel arguments 
 *		  Run the CL kernel
 */
int 
runCLKernels(void)
{
    cl_int   status;
    cl_uint maxDims;
    cl_event events[2];
    size_t globalThreads[1];
    size_t localThreads[1];
    size_t maxWorkGroupSize;
    size_t maxWorkItemSizes[3];

    /**
    * Query device capabilities. Maximum 
    * work item dimensions and the maximmum
    * work item sizes
    */ 
    status = clGetDeviceInfo(
        devices[0], 
        CL_DEVICE_MAX_WORK_GROUP_SIZE, 
        sizeof(size_t), 
        (void*)&maxWorkGroupSize, 
        NULL);
    if(status != CL_SUCCESS) 
    {  
        std::cout<<"Error: Getting Device Info. (clGetDeviceInfo)\n";
        return 1;
    }
    
    status = clGetDeviceInfo(
        devices[0], 
        CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, 
        sizeof(cl_uint), 
        (void*)&maxDims, 
        NULL);
    if(status != CL_SUCCESS) 
    {  
        std::cout<<"Error: Getting Device Info. (clGetDeviceInfo)\n";
        return 1;
    }

    status = clGetDeviceInfo(
        devices[0], 
        CL_DEVICE_MAX_WORK_ITEM_SIZES, 
        sizeof(size_t)*maxDims,
        (void*)maxWorkItemSizes,
        NULL);
    if(status != CL_SUCCESS) 
    {  
        std::cout<<"Error: Getting Device Info. (clGetDeviceInfo)\n";
        return 1;
    }
    
    globalThreads[0] = w;
    localThreads[0]  = 256;

	std::cout << "Max dimensions: " << maxDims << std::endl;
	std::cout << "Device maxWorkGroupSize = " << maxWorkGroupSize << std::endl;
	std::cout << "Device maxWorkItemSizes = " << maxWorkItemSizes[0] << std::endl;
    if(localThreads[0] > maxWorkGroupSize ||
        localThreads[0] > maxWorkItemSizes[0])
    {
        std::cout<<"Unsupported: Device does not support requested number of work items.";
        return 1;
    }

    /*** Set appropriate arguments to the kernel ***/

    /* the input array A to the kernel */
    status = clSetKernelArg(
                    kernel, 
                    0, 
                    sizeof(cl_mem), 
                    (void *)&inputBufferEz);
    if(status != CL_SUCCESS) 
    { 
        std::cout<<"Error: Setting kernel argument. (inputBufferEz)\n";
        return 1;
    }

    /* the input array B to the kernel*/
    status = clSetKernelArg(
                    kernel, 
                    1, 
                    sizeof(cl_mem), 
                    (void *)&inputBufferHy);
    if(status != CL_SUCCESS) 
    { 
        std::cout<< "Error: Setting kernel argument. (inputBufferHy)\n";
        return 1;
    }
	status = clSetKernelArg(
                    kernel, 
                    2, 
                    sizeof(cl_uint), 
                    (void *)&t);
    if(status != CL_SUCCESS) 
    { 
        std::cout<< "Error: Setting kernel argument. (t)\n";
        return 1;
    }
	status = clSetKernelArg(
                    kernel, 
                    3, 
                    sizeof(cl_uint), 
                    (void *)&w);
    if(status != CL_SUCCESS) 
    { 
        std::cout<< "Error: Setting kernel argument. (w)\n";
        return 1;
    }

	cl_ulong startTime, endTime;
	cl_ulong kernelExecTimeNs;
	cl_ulong kernelExecTimeNsT = 0;

	// File handling from chapter 3 of Understanding FDTD. J. B. Schneider
	char basename[20] = "./FieldData/Ez";
	char filename[30];

	#ifdef WIN32
	std::fstream snapshot;
	#else
	int fd;
	#endif

	unsigned int frame = 1;
	unsigned int SnapshotResolutiont = 1;	// Fields snapshots will be saved after this much interval.

	for (t = 1; t < timeN; t++)
	{
		status = clSetKernelArg(
                    kernel, 
                    2, 
                    sizeof(cl_uint), 
                    (void *)&t);
		if(status != CL_SUCCESS) 
		{ 
			std::cout<< "Error: Setting kernel argument. (t)\n";
			return 1;
		}

		for ( cl_uint step=0; step < 2; step++)
		{
			status = clSetKernelArg(
						kernel, 
						4, 
						sizeof(cl_uint), 
						(void *)&flagHalf);
			if(status != CL_SUCCESS) 
			{ 
				std::cout<< "Error: Setting kernel argument. (flag)\n";
				return 1;
			}
			status = clSetKernelArg(
						kernel, 
						5, 
						sizeof(cl_uint), 
						(void *)&n0);
			if(status != CL_SUCCESS) 
			{ 
				std::cout<< "Error: Setting kernel argument. (n0)\n";
				return 1;
			}
			status = clSetKernelArg(
						kernel, 
						6, 
						sizeof(cl_uint), 
						(void *)&n1);
			if(status != CL_SUCCESS) 
			{ 
				std::cout<< "Error: Setting kernel argument. (n1)\n";
				return 1;
			}
			/* 
			 * Enqueue a kernel run call.
			 */
			status = clEnqueueNDRangeKernel(
						 commandQueue,
						 kernel,
						 1,
						 NULL,
						 globalThreads,
						 localThreads,
						 0,
						 NULL,
						 &events[0]);
			if(status != CL_SUCCESS) 
			{ 
				std::cout<<
					"Error: Enqueueing kernel onto command queue. \
					(clEnqueueNDRangeKernel)\n";
				return 1;
			}


			/* wait for the kernel call to finish execution */
			status = clWaitForEvents(1, &events[0]);
			if(status != CL_SUCCESS) 
			{ 
				std::cout<<
					"Error: Waiting for kernel run to finish. \
					(clWaitForEvents)\n";
				return 1;
			}

			clGetEventProfilingInfo(events[0], CL_PROFILING_COMMAND_START,
			sizeof(cl_ulong), &startTime, NULL);
			clGetEventProfilingInfo(events[0], CL_PROFILING_COMMAND_END,
			sizeof(cl_ulong), &endTime, NULL);
			kernelExecTimeNs = 1e-3*(endTime-startTime);
			kernelExecTimeNsT = kernelExecTimeNsT + kernelExecTimeNs;

			flagHalf = !flagHalf;
		}
		/* Enqueue readBuffer*/
		status = clEnqueueReadBuffer(
					commandQueue,
					inputBufferEz,
					CL_TRUE,
					0,
					width * sizeof(cl_double),
					Ez,
					0,
					NULL,
					&events[1]);
    
		if(status != CL_SUCCESS) 
		{ 
			std::cout << 
				"Error: clEnqueueReadBuffer failed. \
				 (clEnqueueReadBuffer)\n";

			return 1;
		}

		// Write Ez snapshot.
		if (t%SnapshotResolutiont == 0 && SaveFields == true)
		{
			#ifdef WIN32
			sprintf_s (filename, "%s%d.fdt", basename, frame);
			snapshot.open (filename, std::ios::out|std::ios::binary);
			#else
			sprintf (filename, "%s%d.fdt", basename, frame);
			fd = open ( filename, O_CREAT|O_WRONLY|O_TRUNC, S_IRWXU );
			#endif

			if (Binary == true)
			{
				#ifdef WIN32
				snapshot.write ( (char*)(Ez+(n1*w)), sizeof(cl_double)*w);
				#else
				write (fd, (void*)(Ez+(n1*w)), sizeof(cl_double)*w);
				#endif
			}
			else
			{
				#ifdef WIN32
				for (unsigned int i=0; i<w; i++)
				{
					snapshot << Ez[i+n1*w] << std::endl;
				}
				#endif
			}

			#ifdef WIN32
			snapshot.close();
			#else
			close (fd);
			#endif

			frame++;
		}

		n0 = (n0+1)%2;
		n1 = (n1+1)%2;
	}
	std::cout << "kernel execution time = " << kernelExecTimeNsT << "us." << std::endl;
    status = clReleaseEvent(events[0]);
    if(status != CL_SUCCESS) 
    { 
        std::cout<<
            "Error: Release event object. \
            (clReleaseEvent)\n";
        return 1;
    }
    /* Wait for the read buffer to finish execution */
    status = clWaitForEvents(1, &events[1]);
    if(status != CL_SUCCESS) 
    { 
        std::cout<<
            "Error: Waiting for read buffer call to finish. \
            (clWaitForEvents)\n";
        return 1;
    }
    
	status = clEnqueueReadBuffer(
                commandQueue,
                inputBufferHy,
                CL_TRUE,
                0,
                width * sizeof(cl_double),
                Hy,
                0,
                NULL,
                &events[1]);
    
    if(status != CL_SUCCESS) 
    { 
        std::cout << 
            "Error: clEnqueueReadBuffer failed. \
             (clEnqueueReadBuffer)\n";

        return 1;
    }
    
    /* Wait for the read buffer to finish execution */
    status = clWaitForEvents(1, &events[1]);
    if(status != CL_SUCCESS) 
    { 
        std::cout<<
            "Error: Waiting for read buffer call to finish. \
            (clWaitForEvents)\n";
        return 1;
    }

    status = clReleaseEvent(events[1]);
    if(status != CL_SUCCESS) 
    { 
        std::cout<<
            "Error: Release event object. \
            (clReleaseEvent)\n";
        return 1;
    }
    return 0;
}

/*
 * \brief Release OpenCL resources (Context, Memory etc.) 
 */
int  
cleanupCL(void)
{
    cl_int status;

    status = clReleaseKernel(kernel);
    if(status != CL_SUCCESS)
    {
        std::cout<<"Error: In clReleaseKernel \n";
        return 1; 
    }
    status = clReleaseProgram(program);
    if(status != CL_SUCCESS)
    {
        std::cout<<"Error: In clReleaseProgram\n";
        return 1; 
    }
    status = clReleaseMemObject(inputBufferEz);
    if(status != CL_SUCCESS)
    {
        std::cout<<"Error: In clReleaseMemObject (inputBufferEz)\n";
        return 1; 
    }
	status = clReleaseMemObject(inputBufferHy);
    if(status != CL_SUCCESS)
    {
        std::cout<<"Error: In clReleaseMemObject (inputBufferHy)\n";
        return 1; 
    }
    status = clReleaseCommandQueue(commandQueue);
    if(status != CL_SUCCESS)
    {
        std::cout<<"Error: In clReleaseCommandQueue\n";
        return 1;
    }
    status = clReleaseContext(context);
    if(status != CL_SUCCESS)
    {
        std::cout<<"Error: In clReleaseContext\n";
        return 1;
    }

    return 0;
}


/* 
 * \brief Releases program's resources 
 */
void
cleanupHost(void)
{
    if(Ez != NULL)
    {
        free(Ez);
        Ez = NULL;
    }
	if(Hy != NULL)
    {
        free(Hy);
        Hy = NULL;
    }
    if(devices != NULL)
    {
        free(devices);
        devices = NULL;
    }
}


/*
 * \brief Print no more than 256 elements of the given array.
 *
 *        Print Array name followed by elements.
 */
void print1DArray(
         const std::string arrayName, 
         const cl_double * arrayData, 
         const unsigned int length)
{
    cl_uint i, j;

	std::cout << arrayName << ":" << std::endl;
	std::cout << std::endl;
    for(i = 0; i < timeN; ++i)
    {
		std:: cout << "t = " << i << " :" << std::endl;
		for (j = 0; j < w; j++ )
		{
			std::cout << arrayData[j + i*w] << " ";
		}
		std::cout << std::endl;
    }
    std::cout << std::endl;

}

int 
main(int argc, char * argv[])
{
    // Initialize Host application 
    if(initializeHost()==1)
        return 1;

    // Initialize OpenCL resources
    if(initializeCL()==1)
        return 1;

	if (initializeFDTD1DKernel () == 1)
		return 1;

    // Run the CL program
    if(runCLKernels()==1)
        return 1;

	//WriteSnapshots (1);
    // Print output array
    //print1DArray(std::string("Output"), outputC, width);
	//print1DArray(std::string("Ez"), Ez, width);
	//print1DArray(std::string("Hy"), Hy, width);

    // Verify output
    //verify();

    // Releases OpenCL resources 
    if(cleanupCL()==1)
        return 1;

    // Release host resources
    cleanupHost();

    return 0;
}