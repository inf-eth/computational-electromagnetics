#include "FDTD2D.h"

int main (int argc, char **argv)
{
	CFDTD2D FDTDSim;
	FDTDSim.StartClock ();
	FDTDSim.Initialize ();
	FDTDSim.StopClock ();
	std::cout << "Initialization Elapsed Time (sec): " << FDTDSim.GetElapsedTime () << std::endl;

	FDTDSim.StartClock ();
	FDTDSim.RunSimulation ();
	FDTDSim.StopClock ();
	std::cout << "Simulation Elapsed Time (sec): " << FDTDSim.GetElapsedTime () << std::endl;
	
	/*
	// File handling from chapter 3 of Understanding FDTD. J. B. Schneider
	char basename[20] = "../../FieldData/Ez";
	char filename[30];

	#ifdef WIN32
	std::fstream snapshot;
	#else
	int fd;
	#endif

	uint frame = 1;
	uint SnapshotResolutiont = 1;	// Fields snapshots will be saved after this much interval.
	bool SaveFields = true;		// Save field snapshots?
	bool Binary = true;			// Save fields in binary format?

	// Initialization.
	for (uint z = 0; z < (w*2); z++)
	{
		Ez[z] = 0;
		Hy[z] = 0;
	}

	uint i, t;

	// Present/future indices.
	uint n0 = 0;
	uint n1 = 1;

	clock_t start, end;
	start = clock();

	for (t = 1; t < TimeN; t++)
	{
		// Hy
		for (i=0; i < (w-1); i++)
		{
			Hy[i+n1*w] = Hy[i+n0*w] + ( (Ez[i+1+n0*w] - Ez[i+n0*w])/imp0 );
		}
		// Ez
		for (i=1; i<w; i++)
		{
			Ez[i+n1*w] = Ez[i+n0*w] + ( (Hy[i+n1*w] - Hy[i-1+n1*w])*imp0 );
		}
		Ez[0+n1*w] = exp ( -1 * pow((t-31.), 2)/100 );

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
				snapshot.write ( (char*)(Ez+(n1*w)), sizeof(double)*w);
				#else
				write (fd, (void*)(Ez+(n1*w)), sizeof(double)*w);
				#endif
			}
			else
			{
				#ifdef WIN32
				for (i=0; i<w; i++)
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
	
	end = clock();
	std::cout << "Time elapsed: " << (double)(end-start)/CLOCKS_PER_SEC << " seconds." << std::endl;*/
	return 0;
}
