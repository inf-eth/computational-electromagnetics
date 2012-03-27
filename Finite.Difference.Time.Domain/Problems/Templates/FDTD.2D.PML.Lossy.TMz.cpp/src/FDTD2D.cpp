#include "FDTD2D.h"
CFDTD2D::CFDTD2D () : 
						I(256),
						J(256),
						c(299792458.),
						delta(2.5e-3),
						dx(delta),
						dy(delta),
						dtscalar(2.),
						dt(delta/(sqrt(2.)*c) /dtscalar),
						PMLw(20),
						NMax(1000),
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
						tResolution(1),
						xResolution(1),
						yResolution(1),
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
	double bytes = I*J*17*3*8+50*8;		// Dynamic arrays + predefined variables.
	std:: cout << "Approximate memory required for simulation: " << bytes/1024 << " Kbytes (" << bytes/(1024*1024) << " MB)." << std::endl;

	// Writing simulation parameters for plotting.
	#ifdef WIN32
	std::fstream parametersfile;
	parametersfile.open ("../../FieldData/Parameters.smp", std::ios::out|std::ios::binary);
	parametersfile.write ((char*)&I, sizeof(uint));
	parametersfile.write ((char*)&J, sizeof(uint));
	parametersfile.write ((char*)&tResolution, sizeof(uint));
	parametersfile.write ((char*)&xResolution, sizeof(uint));
	parametersfile.write ((char*)&yResolution, sizeof(uint));
	parametersfile.write ((char*)&NMax, sizeof(uint));
	parametersfile.write ((char*)&PMLw, sizeof(uint));
	parametersfile.close ();
	#else
	int fdp;
	fdp = open ( "../../FieldData/Parameters.smp", O_CREAT|O_WRONLY|O_TRUNC, S_IRWXU );
	write (fdp, (void*)&I, sizeof(uint));
	write (fdp, (void*)&J, sizeof(uint));
	write (fdp, (void*)&tResolution, sizeof(uint));
	write (fdp, (void*)&xResolution, sizeof(uint));
	write (fdp, (void*)&yResolution, sizeof(uint));
	write (fdp, (void*)&NMax, sizeof(uint));
	write (fdp, (void*)&PMLw, sizeof(uint));
	close (fdp);
	#endif
}
// Initialize data arrays.
void CFDTD2D::Initialize ()
{
	Hx = new double[IHx*JHx*3];
	Bx = new double[IHx*JHx*3];
	Hy = new double[IHy*JHy*3];
	By = new double[IHy*JHy*3];
	Ez = new double[IEz*JEz*3];
	Dz = new double[IEz*JEz*3];
	Dzx = new double[IEz*JEz*3];
	Dzy = new double[IEz*JEz*3];
	//EzSnapshots = new double[(IEz/xResolution)*(JEz/yResolution)];

	urHx = new double[IHx*JHx];
	urHy = new double[IHy*JHy];
	erEz = new double[IEz*JEz];
	smHx = new double[IHx*JHx];
	smHy = new double[IHy*JHy];
	sEz = new double[IEz*JEz];
	Sc = new double[IEz*JEz];
	ScmHx = new double[IHx*JHx];
	ScmHy = new double[IHy*JHy];

	sex = new double[IEz*JEz];
	sey = new double[IEz*JEz];
	smx = new double[IHy*JHy];
	smy = new double[IHx*JHx];
	Scsx = new double[IEz*JEz];
	Scsy = new double[IEz*JEz];
	ScmsmxHy = new double[IHy*JHy];
	ScmsmyHx = new double[IHx*JHx];

	uint i, j, t;
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
}

void CFDTD2D::RunSimulation (bool SaveFields)
{
	// Time loop.
	uint n, i, j;

	// File Handling.
	std::stringstream framestream;
	std::string basename = "../../FieldData/Ez";
	std::string filename;
	uint frame = 1;

	#ifdef WIN32
	std::fstream snapshot;
	#else
	int fd;
	#endif

	uint ProgressResolution;
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
			// Normal Space.
			for (j=PMLw; j < JHx-PMLw; j++)
			{
				Bx[i+IHx*j+IHx*JHx*n2] = (1-ScmHx[i+IHx*j])/(1+ScmHx[i+IHx*j]) * Bx[i+IHx*j+IHx*JHx*n1] + ( (dt/delta)/(1+ScmHx[i+IHx*j]) * (Ez[i+IEz*j+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hx[i+IHx*j+IHx*JHx*n2] = Bx[i+IHx*j+IHx*JHx*n2]/(u0*urHx[i+IHx*j]);

				By[(i+1)+IHy*(j+1)+IHy*JHy*n2] = (1-ScmHy[(i+1)+IHy*(j+1)])/(1+ScmHy[(i+1)+IHy*(j+1)]) * By[(i+1)+IHy*(j+1)+IHy*JHy*n1] + ( (dt/delta)/(1+ScmHy[(i+1)+IHy*(j+1)]) * (Ez[(i+1)+IEz*(j+1)+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hy[(i+1)+IHy*(j+1)+IHy*JHy*n2] = By[(i+1)+IHy*(j+1)+IHy*JHy*n2]/(u0*urHy[(i+1)+IHy*(j+1)]);
			}
			// Lower PML region.
			for (j=0; j < PMLw; j++)
			{
				Bx[i+IHx*j+IHx*JHx*n2] = (1-ScmsmyHx[i+IHx*j])/(1+ScmsmyHx[i+IHx*j]) * Bx[i+IHx*j+IHx*JHx*n1] + ( (dt/delta)/(1+ScmsmyHx[i+IHx*j]) * (Ez[i+IEz*j+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hx[i+IHx*j+IHx*JHx*n2] = Bx[i+IHx*j+IHx*JHx*n2]/(u0*urHx[i+IHx*j]);

				By[(i+1)+IHy*(j+1)+IHy*JHy*n2] = (1-ScmsmxHy[(i+1)+IHy*(j+1)])/(1+ScmsmxHy[(i+1)+IHy*(j+1)]) * By[(i+1)+IHy*(j+1)+IHy*JHy*n1] + ( (dt/delta)/(1+ScmsmxHy[(i+1)+IHy*(j+1)]) * (Ez[(i+1)+IEz*(j+1)+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hy[(i+1)+IHy*(j+1)+IHy*JHy*n2] = By[(i+1)+IHy*(j+1)+IHy*JHy*n2]/(u0*urHy[(i+1)+IHy*(j+1)]);
			}
			// Upper PML region.
			for (j=JHx-PMLw; j < JHx; j++)
			{
				Bx[i+IHx*j+IHx*JHx*n2] = (1-ScmsmyHx[i+IHx*j])/(1+ScmsmyHx[i+IHx*j]) * Bx[i+IHx*j+IHx*JHx*n1] + ( (dt/delta)/(1+ScmsmyHx[i+IHx*j]) * (Ez[i+IEz*j+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hx[i+IHx*j+IHx*JHx*n2] = Bx[i+IHx*j+IHx*JHx*n2]/(u0*urHx[i+IHx*j]);

				By[(i+1)+IHy*(j+1)+IHy*JHy*n2] = (1-ScmsmxHy[(i+1)+IHy*(j+1)])/(1+ScmsmxHy[(i+1)+IHy*(j+1)]) * By[(i+1)+IHy*(j+1)+IHy*JHy*n1] + ( (dt/delta)/(1+ScmsmxHy[(i+1)+IHy*(j+1)]) * (Ez[(i+1)+IEz*(j+1)+IEz*JEz*n1]-Ez[i+IEz*(j+1)+IEz*JEz*n1]) );
				Hy[(i+1)+IHy*(j+1)+IHy*JHy*n2] = By[(i+1)+IHy*(j+1)+IHy*JHy*n2]/(u0*urHy[(i+1)+IHy*(j+1)]);
			}
		}
		// t = 1.
		// Electric field.
		for (i=0; i < IEz; i++)
		{
			// Normal space.
			for (j=PMLw+1; j < JEz-PMLw-1; j++)
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
			// Lower PML.
			for (j=1; j < PMLw+1; j++)
			{
				Dzx[i+IEz*j+IEz*JEz*n2] = (1-Scsx[i+IEz*j])/(1+Scsx[i+IEz*j]) * Dzx[i+IEz*j+IEz*JEz*n1] + ( (dt/delta)/(1+Scsx[i+IEz*j]) * ( Hy[(i+1)+IHy*j+IHy*JHy*n2] - Hy[i+IHy*j+IHy*JHy*n2]) );
				Dzy[i+IEz*j+IEz*JEz*n2] = (1-Scsy[i+IEz*j])/(1+Scsy[i+IEz*j]) * Dzy[i+IEz*j+IEz*JEz*n1] + ( (dt/delta)/(1+Scsy[i+IEz*j]) * (- Hx[i+IHx*j+IHx*JHx*n2] + Hx[i+IHx*(j-1)+IHx*JHx*n2]) );
				Dz[i+IEz*j+IEz*JEz*n2] = Dzx[i+IEz*j+IEz*JEz*n2] + Dzy[i+IEz*j+IEz*JEz*n2];
				Ez[i+IEz*j+IEz*JEz*n2] = Dz[i+IEz*j+IEz*JEz*n2]/(e0*erEz[i+IEz*j]);
			}
			// Upper PML.
			for (j=JEz-PMLw-1; j < JEz-1; j++)
			{
				Dzx[i+IEz*j+IEz*JEz*n2] = (1-Scsx[i+IEz*j])/(1+Scsx[i+IEz*j]) * Dzx[i+IEz*j+IEz*JEz*n1] + ( (dt/delta)/(1+Scsx[i+IEz*j]) * ( Hy[(i+1)+IHy*j+IHy*JHy*n2] - Hy[i+IHy*j+IHy*JHy*n2]) );
				Dzy[i+IEz*j+IEz*JEz*n2] = (1-Scsy[i+IEz*j])/(1+Scsy[i+IEz*j]) * Dzy[i+IEz*j+IEz*JEz*n1] + ( (dt/delta)/(1+Scsy[i+IEz*j]) * (- Hx[i+IHx*j+IHx*JHx*n2] + Hx[i+IHx*(j-1)+IHx*JHx*n2]) );
				Dz[i+IEz*j+IEz*JEz*n2] = Dzx[i+IEz*j+IEz*JEz*n2] + Dzy[i+IEz*j+IEz*JEz*n2];
				Ez[i+IEz*j+IEz*JEz*n2] = Dz[i+IEz*j+IEz*JEz*n2]/(e0*erEz[i+IEz*j]);
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
	}
	std::cout << "\r" << "Simulation complete!" << std::endl;
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
