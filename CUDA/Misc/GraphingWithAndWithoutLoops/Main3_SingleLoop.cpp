#include <cstdio>
#include <cmath>

#define HMIN	11
#define HMAX	75
#define VMIN	11
#define VMAX	39
#define PI		3.14159

int main (int argc, char* argv[])
{
	int H, V;
	float n;
	
	int i;
	float r;
	int nf;	

	printf("Enter horizontal width (10-75 chars): ");
	scanf("%d", &H);
	H < 0 ? H = H*-1 : H = H;
	(H < HMIN || H > HMAX) ? (H < HMIN ? H = HMIN : H = HMAX) : (H%2 == 0 ? H++ : H = H);

	printf("Enter vertical width (11-39 chars): ");
	scanf("%d", &V);
	V < 0 ? V = V*-1 : V = V;
	(V < VMIN || V > VMAX) ? (V < VMIN ? V = VMIN : V = VMAX) : (V%2 == 0 ? V++ : V = V);

	printf("Enter number of cycles: ");
	scanf("%f", &n);
	n < 0.f ? n = n*-1.f : n == 0.f ? n = 1.f : n = n;

	printf("H = %d, V = %d, n = %.2f\n", H, V, n);

	for (i=0; i<V*H; i++)
	{
		(i != 0 && i%H == 0) ? printf("\n") : i=i;
		r = ((float)(sin(2.f*PI*((float)(i%H)/(H-1)*n)))*(V/2));
		nf = (int) ((r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5));		// Rounding off float to nearest int.
		((V/2)-nf == i/H) ? printf("*") : i%H == 0 ? printf("|") : (i/H == V/2 ? printf("-") : printf(" "));		
	}
	printf("\n");

	return 0;
}
