#ifndef WIN32
#include <fcntl.h>
#endif
#include <iostream>
#include <cmath>
#include <ctime>

typedef unsigned int uint;

// Simulation parameters.
const uint w = 1024*1024;			// No. of spatial steps
const uint TimeN = 256;		// No. of time steps
const double imp0 = 377.0;		// Impedence of free space

// Data Arrays.
double *Ez = new double[w*2]; // z-component of E-field
double *Hy = new double[w*2]; // y-component of H-field

int main (int argc, char **argv)
{
	// File handling from chapter 3 of Understanding FDTD. J. B. Schneider
	char basename[20] = "../../FieldData/Ez";
	char filename[30];

	#ifdef WIN32
	FILE *snapshot;
	#else
	int fd;
	#endif

	uint frame = 1;
	uint SnapshotResolutiont = 1;	// Fields snapshots will be saved after this much interval.
	bool SaveFields = false;		// Save field snapshots?
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
			fopen_s (&snapshot, filename, "w");
			#else
			sprintf (filename, "%s%d.fdt", basename, frame);
			fd = open ( filename, O_CREAT|O_WRONLY|O_TRUNC, S_IRWXU );
			#endif

			if (Binary == true)
			{
				#ifdef WIN32
				fwrite ( (void*)(Ez+(n1*w)), sizeof(double), w, snapshot);
				#else
				write (fd, (void*)(Ez+(n1*w)), sizeof(double)*w);
				#endif
			}
			else
			{
				#ifdef WIN32
				for (i=0; i<w; i++)
				{
					fprintf_s (snapshot, "%g\n", Ez[i+n1*w]);
				}
				#endif
			}

			#ifdef WIN32
			fclose (snapshot);
			#else
			close (fd);
			#endif

			frame++;
		}
		n0 = (n0+1)%2;
		n1 = (n1+1)%2;
	}

	end = clock();
	std::cout << "Time elapsed: " << (double)(end-start)/CLOCKS_PER_SEC << " seconds." << std::endl;
	return 0;
}
