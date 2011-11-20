#include "FDTD2D.h"
CFDTD2D::CFDTD2D () : 
						I(1024),
						J(1024),
						c(299792458),
						delta(5e-3),
						dx(delta),
						dy(delta),
						dt(delta/(sqrt(2.)*c)),
						PMLw(50),
						NMax(256),
						f(2e9),
						pi(4*atan(1.)),
						e0(1e-9/(36*pi)),
						u0(1e-7*4*pi),
						Two_pi_f_deltat(2*pi*f*dt),
						NHW(1/(2*f*dt)),
						Js(2+PMLw),
						Is(2),
						n0(0),
						n1(1),
						tResolution(1),
						xResolution(1),
						yResolution(1),
						n2(2),
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

}
CFDTD2D::~CFDTD2D ()
{

}
