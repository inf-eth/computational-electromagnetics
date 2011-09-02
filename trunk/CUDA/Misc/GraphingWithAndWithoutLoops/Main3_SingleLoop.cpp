#include <iostream>
#include <cmath>

using std::cin;
using std::cout;
using std::endl;

#define HMIN	11
#define HMAX	75
#define VMIN	11
#define VMAX	39
#define PI		3.14159

int main (int argc, char* argv[])
{
	int H;
	int V;
	int n;
	int i;

	float t;
	float f;
	int nf;
	float r;

	cout << "Enter horizontal width (10-75 chars): ";
	cin >> H;
	H < 0 ? H = H*-1 : H = H;
	(H < HMIN || H > HMAX) ? (H < HMIN ? H = HMIN : H = HMAX) : (H%2 == 0 ? H++ : H = H);

	cout << "Enter vertical width (11-39 chars): ";
	cin >> V;
	V < 0 ? V = V*-1 : V = V;
	(V < VMIN || V > VMAX) ? (V < VMIN ? V = VMIN : V = VMAX) : (V%2 == 0 ? V++ : V = V);

	cout << "Enter number of cycles: ";
	cin >> n;
	n < 0 ? n = n*-1 : n == 0 ? n = 1 : n = n;

	cout << "H = " << H << ", V = " << V << ", n = " << n << endl;

	for (i=0; i<V*H; i++)
	{
		(i != 0 && i%H == 0) ? printf("\n") : i=i;
		t = ((float)(i%H)/(H-1)*n);
		f = sin(2.f*PI*t);
		r = ((float)f*(V/2));
		nf = (int) ((r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5));
		((V/2)-nf == i/H) ? printf("*") : i%H == 0 ? printf("|") : (i/H == V/2 ? printf("-") : printf(" "));		
	}
	printf("\n");

	return 0;
}
