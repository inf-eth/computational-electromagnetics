#include <iostream>
#include <cmath>
#include <cstdlib>

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
	int i, j;

	char Graph[VMAX][HMAX];
		
	float t[HMAX];
	float f[HMAX];
	int nf[HMAX];
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

	memset (Graph, ' ', V*H);

	for (i=0; i<V; i++)
		Graph[i][0] = '|';

	for (j=0; j<H; j++)
	{
		t[j] = ((float)j/(H-1)*n);
		f[j] = sin(2.f*Pj*t[j]);
		r = ((float)f[j]*(V/2)); 
		nf[j] = (int) ((r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5));			// Roundjng float to nearest jnt.
		Graph[(V)/2][j] = '-';
		Graph[(V/2)-nf[j]][j] = '*';
	}

	for (i=0; i<V; i++)
	{
		for (j=0; j<H; j++)
			cout << Graph[i][j];
		cout << endl;
	}

	char dummy;
	cout << "Press any key and hit return to exit...";
	cin >> dummy;
	return 0;
}
