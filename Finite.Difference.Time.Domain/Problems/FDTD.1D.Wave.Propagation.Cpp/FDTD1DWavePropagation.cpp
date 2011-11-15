#include <iostream>
#include <cmath>
#include <ctime>

typedef unsigned int uint;

// Simulation parameters.
const uint w = 1024;			// No. of spatial steps
const uint TimeN = 30;		// No. of time steps
const double imp0 = 377.0;		// Impedence of free space

// Data Arrays.
double *Ez = new double[w*TimeN]; // z-component of E-field
double *Hy = new double[w*TimeN]; // y-component of H-field

int main (int argc, char **argv)
{
	// File handling from chapter 3 of Understanding FDTD. J. B. Schneider
	char basename[20] = "../../FieldData/Ez";
	char filename[30];
	uint frame = 1;
	FILE *snapshot;
	uint SnapshotResolutiont = 1;	// Fields snapshots will be saved after this much interval.
	bool SaveFields = true;		// Save field snapshots?
	bool Binary = true;			// Save fields in binary format?

	// Initialization.
	for (uint z = 0; z < (w*TimeN); z++)
	{
		Ez[z] = 0;
		Hy[z] = 0;
	}

	uint i, t;
	clock_t start, end;
	start = clock();
	for (t = 1; t < TimeN; t++)
	{
		// Hy
		for (i=0; i < (w-1); i++)
		{
			Hy[i+t*w] = Hy[i+(t-1)*w] + ( (Ez[i+1+(t-1)*w] - Ez[i+(t-1)*w])/imp0 );
		}
		// Ez
		for (i=1; i<w; i++)
		{
			Ez[i+t*w] = Ez[i+(t-1)*w] + ( (Hy[i+t*w] - Hy[i-1+t*w])*imp0 );
		}
		Ez[0+t*w] = Ez[0+t*w] + exp ( -1 * pow((t-30.), 2)/100 );

		// Write Ez snapshot.
		if (t%SnapshotResolutiont == 0 && SaveFields == true)
		{
			sprintf_s (filename, "%s%d.fdt", basename, frame);
			fopen_s (&snapshot, filename, "w");

			if (Binary == true)
			{
				fwrite ( (void*)&Ez[t*w], sizeof(double), w, snapshot);
			}
			else
			{
				for (i=0; i<w; i++)
				{
					fprintf_s (snapshot, "%g\n", Ez[i+t*w]);
				}
			}
			fclose (snapshot);
			frame++;
		}
	}
	end = clock();
	std::cout << "Time elapsed: " << (double)(end-start)/CLOCKS_PER_SEC << " seconds." << std::endl;
	return 0;
}