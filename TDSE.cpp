#include <iostream>
#include <fstream>
//#include < string>
#include <sstream>
#include <time.h>
#include <stdio.h>
#include <arrayfire.h>
#define PI 3.141592653589793

using namespace af;
using namespace std;

//numerical solution of 3D TDSE following Chelkowski et.al(PRA. 46, R5342)
//Coulomb potential Vc(rou,z)=-sqrt(rou^2+(z-R/2)^2)-sqrt(rou^2+(z+R/2)^2)
//for H2+ in cylindrical coordinate system.
//Intgrnd=@(x) (sqrt(2)/L/besselj(1,xn(n))*besselj(0,xn(n)*x/L).*(sqrt(2)*exp(-sqrt(x.^2+z(m)^2)))).*x;

int main() {

	//af::info(); af::array myarray;

beginning:
	af::array imaginarynumber = af::array(0, 1);

	std::cout << "Input Total Longitudinal Grid Size (au): ";
	int lgp;
	cin >> lgp;
	int NumZ = 8 * lgp; //Number of grid points
	int NumZ1 = NumZ - 1;
	int Zmin = -lgp / 2; //Z + grid size
	int Zmax = -Zmin; //Z - grid size
	//double L=8.0;
	double L;
	double radius = 2.0;
	// std::cout << "Input radius: " ;
	// cin >> radius;
	std::cout << "Input Transverse Grid Size (max radial value)(au): ";
	int LL;
	cin >> LL;
	L = LL + 0.0;
	//cin >> L;
	//int gp=64*L;
	//std::cout << "Grid points per a.u. (transverse): ";
	int rgp = 8;
	//cin >> rgp;
	int gp = LL * rgp;
	//int gp=2048*(LL/8.0);
	double dgp = gp + 0.0;

	//Zeros of Bessel functions
	double xn[256] = { 2.404825557695773, 5.520078110286311, 8.653727912911011, 11.79153443901428, 14.93091770848779, 18.07106396791092, 21.21163662987926, 24.35247153074930, 27.49347913204025, 30.63460646843198, 33.77582021357357, 36.91709835366404, 40.05842576462824, 43.19979171317673, 46.34118837166181, 49.48260989739782, 52.62405184111500, 55.76551075501998, 58.90698392608094, 62.04846919022717, 65.18996480020686, 68.33146932985680, 71.47298160359373, 74.61450064370184, 77.75602563038806, 80.89755587113763, 84.03909077693819, 87.18062984364115, 90.32217263721048, 93.46371878194477, 96.60526795099627, 99.74681985868060, 102.8883742541948, 106.0299309164516, 109.1714896498054, 112.3130502804949, 115.4546126536669, 118.5961766308725, 121.7377420879510, 124.8793089132329, 128.0208770060083, 131.1624462752139, 134.3040166383055, 137.4455880202843, 140.5871603528543, 143.7287335736897, 146.8703076257966, 150.0118824569548, 153.1534580192279, 156.2950342685335, 159.4366111642631, 162.5781886689467, 165.7197667479550, 168.8613453692358, 172.0029245030782, 175.1445041219027, 178.2860842000738, 181.4276647137311, 184.5692456406387, 187.7108269600494, 190.8524086525815, 193.9939907001091, 197.1355730856614, 200.2771557933324, 203.4187388081986, 206.5603221162445, 209.7019057042941, 212.8434895599495, 215.9850736715340, 219.1266580280406, 222.2682426190843, 225.4098274348593, 228.5514124660988, 231.6929977040385, 234.8345831403832, 237.9761687672757, 241.1177545772680, 244.2593405632957, 247.4009267186528, 250.5425130369700, 253.6840995121931, 256.8256861385644, 259.9672729106045, 263.1088598230955, 266.2504468710659, 269.3920340497761, 272.5336213547049, 275.6752087815375, 278.8167963261531, 281.9583839846149, 285.0999717531596, 288.2415596281877, 291.3831476062552, 294.5247356840650, 297.6663238584589, 300.8079121264111, 303.9495004850206, 307.0910889315050, 310.2326774631950, 313.3742660775278, 316.5158547720429, 319.6574435443762, 322.7990323922556, 325.9406213134967, 329.0822103059986, 332.2237993677396, 335.3653884967741, 338.5069776912285, 341.6485669492981, 344.7901562692440, 347.9317456493902, 351.0733350881206, 354.2149245838764, 357.3565141351537, 360.4981037405011, 363.6396933985170, 366.7812831078483, 369.9228728671875, 373.0644626752712, 376.2060525308784, 379.3476424328284, 382.4892323799793, 385.6308223712263, 388.7724124055006, 391.9140024817673, 395.0555925990248, 398.1971827563028, 401.3387729526616, 404.4803631871904, 407.6219534590068, 410.7635437672553, 413.9051341111063, 417.0467244897553, 420.1883149024216, 423.3299053483481, 426.4714958267996, 429.6130863370627, 432.7546768784446, 435.8962674502723, 439.0378580518925, 442.1794486826700, 445.3210393419877, 448.4626300292460, 451.6042207438617, 454.7458114852677, 457.8874022529128, 461.0289930462604, 464.1705838647888, 467.3121747079900, 470.4537655753698, 473.5953564664471, 476.7369473807533, 479.8785383178323, 483.0201292772397, 486.1617202585426, 489.3033112613194, 492.4449022851590, 495.5864933296609, 498.7280843944346, 501.8696754790994, 505.0112665832841, 508.1528577066267, 511.2944488487740, 514.4360400093816, 517.5776311881133, 520.7192223846410, 523.8608135986446, 527.0024048298115, 530.1439960778367, 533.2855873424221, 536.4271786232769, 539.5687699201169, 542.7103612326644, 545.8519525606483, 548.9935439038036, 552.1351352618712, 555.2767266345981, 558.4183180217369, 561.5599094230457, 564.7015008382880, 567.8430922672325, 570.9846837096531, 574.1262751653285, 577.2678666340424, 580.4094581155830, 583.5510496097432, 586.6926411163202, 589.8342326351157, 592.9758241659354, 596.1174157085893, 599.2590072628912, 602.4005988286588, 605.5421904057138, 608.6837819938814, 611.8253735929903, 614.9669652028729, 618.1085568233650, 621.2501484543055, 624.3917400955367, 627.5333317469042, 630.6749234082564, 633.8165150794449, 636.9581067603241, 640.0996984507513, 643.2412901505867, 646.3828818596930, 649.5244735779358, 652.6660653051830, 655.8076570413054, 658.9492487861759, 662.0908405396701, 665.2324323016657, 668.3740240720429, 671.5156158506840, 674.6572076374736, 677.7987994322984, 680.9403912350472, 684.0819830456108, 687.2235748638821, 690.3651666897557, 693.5067585231285, 696.6483503638989, 699.7899422119673, 702.9315340672359, 706.0731259296085, 709.2147177989908, 712.3563096752900, 715.4979015584149, 718.6394934482762, 721.7810853447858, 724.9226772478572, 728.0642691574057, 731.2058610733475, 734.3474529956008, 737.4890449240848, 740.6306368587203, 743.7722287994293, 746.9138207461351, 750.0554126987625, 753.1970046572373, 756.3385966214867, 759.4801885914389, 762.6217805670236, 765.7633725481714, 768.9049645348141, 772.0465565268846, 775.1881485243170, 778.3297405270463, 781.4713325350086, 784.6129245481411, 787.7545165663819, 790.8961085896701, 794.0377006179459, 797.1792926511503, 800.3208846892252, 803.4624767321134 };
	int bz1;
	std::cout << "Insert Number of Bessel Zeros: ";
	cin >> bz1;
	const int bz = bz1;

	double rou0;
	std::cout << "Calculate flux at (radial dist a.u.): ";
	cin >> rou0;
	if (rou0 > LL) goto beginning;
	rou0 *= 8;
	double z0;
	std::cout << "Calculate flux at (longitudinal dist a.u.): ";
	cin >> z0;
	if (z0 > lgp) goto beginning;
	z0 *= 8;

	double dNumZ = NumZ + 0.0;
	std::cout << "Zmin = " << Zmin << endl;
	std::cout << "Zmax = " << Zmax << endl;
	std::cout << "Grid Points = " << NumZ << " x " << gp << endl;
	std::cout << "Flux at g.p. " << z0 << endl;
	double dZmin = Zmin + 0.0;
	double dZ = (Zmax - Zmin) / (dNumZ - 1);
	std::cout << "dZ = " << dZ << endl;

	//strings for coulomb matrices
	string cd1, cv1, lgth, b, gs;
	const char* cdc, * cvc;
	stringstream sL, sB, gsI;
	gsI << lgp;
	sL << L;
	sB << bz;
	gs = gsI.str();
	lgth = sL.str();
	b = sB.str();
	cd1 = "h_cd_"; cd1 += gs; cd1 += "au_"; cd1 += lgth; cd1 += "au_"; cd1 += b; cd1 += ".txt";
	cv1 = "h_cv_"; cv1 += gs; cv1 += "au_"; cv1 += lgth; cv1 += "au_"; cv1 += b; cv1 += ".txt";
	cdc = cd1.c_str();
	cvc = cv1.c_str();

	//check Coulomb Matrices
	ifstream dmat;  ifstream vmat;
	dmat.open(cdc); vmat.open(cvc);
	if (vmat.is_open() != 1 || dmat.is_open() != 1) { goto beginning; }

	//laser input
	double wavelen = 1064.0;//laser wavelength in nm
	std::cout << "Enter Wavelength (*100nm): ";
	cin >> wavelen;
	wavelen *= 100;
	std::cout << "Input Laser Intensity (*1E14): ";
	double intense;
	cin >> intense;
	int lp;
	std::cout << "Enter # laser periods : ";
	cin >> lp;
	//laser parameters
	double w = 2 * PI * 3E17 / 4.134E16 / wavelen; //laser frequency in atomic units
	double t0 = 729.0;  //pulse switch on time
	double inten; //laser peak intensity
	inten = intense * (1E14);
	double E0 = sqrt(inten / (3.5E16));//laser electric field

	//time
	double dT = 0.01;
	//std::cout << "Timestep: ";
	//cin >> dT;
	double tau = 2 * PI / w;
	const double TimeMax = tau * lp;
	const int Tint = ceil(tau * lp / dT); //==TimeMax/dT;
	//double Ttime=floor(tau/dT);
	//double Halfway=floor(Tint/2.0);
	double Ttime = Tint - 4;
	double Halfway = 0;
	std::cout << "Simulation time = " << TimeMax << " a.u." << endl;
	std::cout << "Flux start = " << Halfway / 100 << " a.u., end = " << (Halfway + Ttime) / 100 << " a.u." << endl;

	double* z1;  //z coordinate
	z1 = new double[NumZ];
	for (int n = 0; n != NumZ; n++) {
		z1[n] = n * dZ + dZmin;
	}
	//array z;
	//try{
	af::array z = af::constant(0, NumZ, f64);
	//}catch(af::exception& e){printf("%s\n",e.what());}
	//z=array(NumZ,z1);
	double dR = L / (dgp - 1);
	double* rou1;
	double* rou2;
	rou1 = new double[gp];
	rou2 = new double[gp * NumZ];

	for (int n = 0; n != gp; n++) {//radial coordinate
		rou1[n] = n * dR;
		for (int m = 0; m != NumZ; m++) {
			rou2[n + gp * m] = rou1[n] * pow(sqrt(rou1[n] * rou1[n] + z1[m] * z1[m]), -3);
		}
	}
	af::array rou = af::constant(0, gp, 1, f64);
	af::array zr3 = af::constant(0, gp, NumZ, f64);
	rou = af::array(gp, 1, rou1);
	zr3 = af::array(gp, NumZ, rou2);
	zr3(0, NumZ / 2) = 0;
	delete[] rou2;

	//Create output filename
	string fn, fnp, fnhg, dpn, pL, iS, wlen, srou0, sz0;
	stringstream slp, sI, Wavl, srou, sz;
	// std::cout << "Output Filename: " ;
	// cin >> fn;
	srou << rou0 / 8;
	sz << z0 / 8;
	slp << lp;
	Wavl << wavelen;
	sI << intense;
	srou0 = srou.str();
	sz0 = sz.str();
	pL = slp.str();
	wlen = Wavl.str();
	iS = sI.str();
	fn = "h_"; fn += wlen; fn += "nm_"; fn += gs; fn += "au_"; fn += lgth; fn += "au_"; fn += iS; fn += "I_"; fn += pL; fn += "t_"; fn += b; fn += "bz_z2p_"; fn += sz0; fn += "au"; fn += "_r_"; fn += srou0; fn += "au";
	fnp = fn;
	fnhg = fn;
	dpn = fn;
	fn += ".txt"; fnp += "_wf.txt"; fnhg += "_mom.txt"; dpn += "_dp.txt";
	const char* fn1, * fnp1, * hspec, * dpm;
	fn1 = fn.c_str();
	fnp1 = fnp.c_str();
	hspec = fnhg.c_str();
	dpm = dpn.c_str();

	//track calculation time
	clock_t start = clock();

	//boundary mask
	double* Mask1;
	Mask1 = new double[NumZ];
	af::array Mask = af::constant(0, NumZ, 1, f64);
	for (int n = 0; n != NumZ; n++) {
		if (abs(z1[n]) > (Zmax - 5.0)) {
			Mask1[n] = pow((sin((Zmax - abs(z1[n])) / 5.0 * PI / 2.0)), (0.125));
		}
		else {
			Mask1[n] = 1.0;
		}
	}
	Mask = af::array(NumZ, 1, Mask1);
	delete[] Mask1;

	double* rMask1;
	rMask1 = new double[gp];
	af::array rMask = af::constant(0, gp, 1, f64);
	for (int n = 0; n != gp; n++) {
		if (rou1[n] > (L - 5.0)) {
			rMask1[n] = pow((sin((L - rou1[n]) / 5.0 * PI / 2.0)), (0.125));
		}
		else {
			rMask1[n] = 1.0;
		}
	}
	rMask = af::array(gp, 1, rMask1);
	delete[] rMask1;

	//Set Initial Wave Functions
	af::array psi = af::constant(0, gp, NumZ, c64);
	//array psi0=zeros(gp,NumZ,c64);
	af::array psi2 = af::constant(0, bz, gp, c64);
	double* psi2d;
	psi2d = new double[bz * gp];

	for (int m = 0; m != bz; m++) {
		for (int n = 0; n != gp; n++) {
			psi2d[m + n * bz] = sqrt(2.0) / L / _j1(xn[m]) * _j0(xn[m] * rou1[n] / L);//bessel basis fn's
		}
	}
	psi2 = af::complex(af::array(bz, gp, psi2d));


	double* psid;
	psid = new double[bz * NumZ];
	for (int i = 0; i != bz; i++) {
		for (int j = 0; j != NumZ; j++) {
			psid[i + j * bz] = 0;
			for (int k = 0; k != gp; k++) {
				psid[i + j * bz] += psi2d[i + k * bz] * pow(PI, -0.5) * exp(-sqrt(z1[j] * z1[j] + rou1[k] * rou1[k])) * rou1[k] * dR;
			}
		}
	}
	af::array psi1f = af::constant(0, bz, NumZ, c64);
	af::array psiT = af::constant(0, NumZ, bz, c64);
	af::array psi_initial = af::constant(0, bz, NumZ, c64);
	psi1f = af::complex(af::array(bz, NumZ, psid));
	psi_initial = af::complex(af::array(bz, NumZ, psid));
	delete[] psi2d;
	delete[] psid;
	delete[] rou1;
	delete[] z1;

	double Norm = 0;
	double pint = 0;
	af::array Norm1 = af::constant(0, NumZ, c64);
	af::array Norm2 = af::constant(0, NumZ, c64);

	gfor(seq m, NumZ) {
		af::array cero = matmul(psi1f(span, m), af::conjg(psi1f(span, m)));
		af::array cero2 = matmul(psi1f(span, m), af::conjg(psi_initial(span, m)));

		Norm1(m) = sum(cero(span, 0));
		Norm2(m) = sum(cero2(span, 0));
	}
	Norm = sum<double>(real(Norm1(span)));
	pint = sum<double>(real(Norm2(span)));

	std::cout << "Norm = " << Norm << ' ' << pint << endl;
	psi1f = psi1f / sqrt(Norm);
	psi_initial = psi_initial / sqrt(Norm);

	//momentum grid
	af::array p = af::constant(0, NumZ, f64);
	double dP = 2 * PI / (Zmax - Zmin);
	double* p1 = new double[NumZ];
	double* p2 = new double[NumZ];
	for (int n = 0; n != NumZ; n++) {
		p1[n] = (n - NumZ / 2) * dP;
	}
	//FFTSHIFT (no fourier transform necessary)
	for (int n = 0; n != NumZ / 2; n++) {
		p2[n] = p1[NumZ / 2 + n];
		p2[NumZ / 2 + n] = p1[n];
	}
	p = af::array(NumZ, p2);
	delete[] p1;
	delete[] p2;

	double* creal, * cim, * dm;
	dm = new double[bz * bz * NumZ];

	std::cout << "cvmat" << endl;
	for (int n = 0; n != NumZ; n++) {
		for (int m = 0; m != bz; m++) {
			for (int l = 0; l != bz; l++) { //(n,m,l)==(NumZ,bz,bz), m=row, l=column, of 16x16 matrices.
				dmat >> dm[m + l * bz + n * bz * bz];
			}
		}
		if (n % 512 == 0) { std::cout << n / dNumZ << '\t'; }
	}
	std::cout << endl;
	dmat.close();
	af::array cd = af::constant(0, bz, bz, NumZ, c64);
	cd = af::complex(af::array(bz, bz, NumZ, dm));
	//std::cout << sum<double>(real(cd(bz/2,bz/2,NumZ/2))) << ' ' << dm[bz/2+bz*bz/2+bz*bz*NumZ/2] << endl;
	delete[] dm;

	creal = new double[bz * NumZ];
	cim = new double[bz * NumZ];
	//double dump;
	for (int n = 0; n != NumZ; n++) {
		for (int m = 0; m != bz; m++) {
			vmat >> creal[m + n * bz] >> cim[m + n * bz];
		}
		if (n % 512 == 0) { std::cout << n / dNumZ << '\t'; }
	}
	std::cout << endl;
	vmat.close();

	af::array cv0 = af::constant(0, bz, NumZ, c64);
	cv0 = af::complex(af::array(bz, NumZ, creal), af::array(bz, NumZ, cim));
	//std::cout << sum<double>(real(cv0(bz/2,NumZ/2))) << ' ' << sum<double>(imag(cv0(bz/2,NumZ/2))) << ' ' << creal[bz/2+bz*NumZ/2] << ' ' << cim[bz/2+bz*NumZ/2] << endl;
	delete[] creal;
	delete[] cim;

	//Time
	std::cout << "Completed Coulomb matrices: " << (clock() - start) / CLOCKS_PER_SEC << "sec" << endl;

	af::array mom = af::constant(0, NumZ, 1, c64);
	af::array las = af::constant(0, NumZ, 1, c64);

	mom(span, 0) = af::complex(cos(pow(p(span), 2) * 0.25 * dT), -sin(pow(p(span), 2) * 0.25 * dT));// def momentum operator


	af::array dsum = af::constant(0, NumZ, f64);
	af::array gNorm = af::constant(0, Tint + 4, f64);
	af::array gp0 = af::constant(0, Tint + 4, f64);
	af::array d0 = af::constant(0, 3, f64);
	af::array dp0 = af::constant(0, Tint + 4, f64);
	af::array d = af::constant(0, 3, f64);
	af::array dp = af::constant(0, Tint + 4, f64);
	af::array dA = af::constant(0, Tint + 4, f64);

	af::array Jrout = af::constant(0, Ttime, z0 * 2, c64);
	af::array Jzt = af::constant(0, Ttime, rou0, c64);
	af::array dblm = af::constant(0, rou0, c64);
	af::array dblm2 = af::constant(0, z0 * 2, c64);
	af::array Jz1 = af::constant(0, rou0, c64);
	af::array Jz2 = af::constant(0, rou0, c64);
	af::array Jrou1 = af::constant(0, z0 * 2, c64);
	af::array Jrou2 = af::constant(0, z0 * 2, c64);


	double t = 0.0;
	int time = 0;
	double v = 0.0;

	//propagation loop.
	psiT = psi1f.T();
	//psi0=psi2.T()*psi1f;
	af::array field = af::constant(0, 1, f64);
	//ofstream mout("mom.txt");

	//ofstream wf(fnp1);
	start = clock();
	af::array out = af::constant(0, 1, f64);
	while (t < TimeMax) {

		gfor(seq n, NumZ) {
			af::array cero = matmul(psi1f(span, n), af::conjg(psi1f(span, n)));
			Norm1(n) = sum(cero(span, 0));
			af::array cero2 = matmul(psi1f(span, n), af::conjg(psi_initial(span, n)));
			af::array cero3 = matmul(cero2, af::conjg(cero2));
			Norm2(n) = sum(sqrt(cero3(span, 0)));
		}
		gNorm(time) = sum(real(Norm1(span)));
		gp0(time) = sum(real(Norm2(span)));

		if (time % 100 == 0) {
			std::cout << t << '\t' << sum<double>(gNorm(time)) << '\t' << sum<double>(gp0(time)) << '\t';//Norm
			std::cout << (clock() - start) / CLOCKS_PER_SEC << " sec" << endl;//Time
		}

		gfor(seq n, bz) {
			af::array temp = matmul(fft(psiT(span, n)), mom(span, 0));
			psiT(span, n) = ifft(temp(span, 0)) / dNumZ;//no longer scaled by NumZ
		}

		//calculating coulomb coupling with diagonalizing matrix conjugate transpose
		psi1f = psiT.T();
		gfor(seq n, NumZ) {
			psi1f(span, n) = matmul(cd(span, span, n).H(), psi1f(span, n));
		}
		psi1f = matmul(cv0, psi1f);
		gfor(seq n, NumZ) {
			psi1f(span, n) = matmul(cd(span, span, n), psi1f(span, n));
		}
		psiT = psi1f.T();

		v = E0 * sin(w * (t + dT / 2)) * exp(-pow((t + dT / 2 - TimeMax / 2), 2) / (2 * pow(TimeMax / 6, 2)));// Gaussian Envelope == e^(-t^2/d^2), t=[-TimeMax/2,TimeMax/2], d~T
		field(0) = v * gNorm(time);

		//laser interaction
		las(span, 0) = af::complex(cos(v * z(span) * dT), -sin(v * z(span) * dT));
		gfor(seq n, bz) {
			psiT(span, n) = matmul(psiT(span, n), las(span, 0));
		}

		//transform to momentum space
		gfor(seq n, bz) {
			af::array temp = matmul(fft(psiT(span, n)), mom(span, 0));
			psiT(span, n) = matmul(ifft(temp(span, 0)), Mask(span, 0)) / dNumZ;//no longer scaled by NumZ

		}
		psi1f = psiT.T();

		psi = matmul(psi2.T(), psi1f); //convert to Cartesian coord's

		if (time < Tint) {
			if (time >= Halfway && time < (Halfway + Ttime)) {
				dblm2(span) = psi(rou0 + 1, seq(NumZ / 2 - z0, NumZ / 2 + z0 - 1)) - psi(rou0, seq(NumZ / 2 - z0, NumZ / 2 + z0 - 1));
				Jrou1(span) = matmul(psi(rou0, seq(NumZ / 2 - z0, NumZ / 2 + z0 - 1)).T(), af::conjg(dblm2(span)));
				Jrou2(span) = matmul(af::conjg(psi(rou0, seq(NumZ / 2 - z0, NumZ / 2 + z0 - 1))).T(), dblm2(span)) * rou0 / 8 / sqrt(z(seq(NumZ / 2 - z0, NumZ / 2 + z0 - 1)) * z(seq(NumZ / 2 - z0, NumZ / 2 + z0 - 1)) + rou0 * rou0 / 64);
				//}
				gfor(af::array i, rou0) {
					dblm(i) = psi(i, NumZ / 2 - z0 - 1) - psi(i, NumZ / 2 - z0);
					Jz1(i) = matmul(psi(i, NumZ / 2 - z0), af::conjg(dblm(i)));
					Jz2(i) = matmul(af::conjg(psi(i, NumZ / 2 - z0)), dblm(i)) * z0 / 8 / sqrt(rou(i) * rou(i) + z0 * z0 / 64);
				}
				Jrout(time - Halfway, span) = (Jrou1 - Jrou2) * af::conjg(imaginarynumber) * dR / 2.0;
				Jzt(time - Halfway, span) = (Jz1 - Jz2) * af::conjg(imaginarynumber) * dZ / 2.0;
				gfor(seq i, rou0) {
					dblm(i) = psi(i, NumZ / 2 + z0 + 1) - psi(i, NumZ / 2 - z0);
					Jz1(i) = matmul(psi(i, NumZ / 2 + z0), af::conjg(dblm(i)));
					Jz2(i) = matmul(af::conjg(psi(i, NumZ / 2 + z0)), dblm(i)) * z0 / 8 / sqrt(rou(i) * rou(i) + z0 * z0 / 64);
				}
				Jzt(time - Halfway, span) += (Jz1 - Jz2).T() * af::conjg(imaginarynumber) * dZ / 2.0;
			}
			else if (time > (Halfway + Ttime)) goto end;
		}

		gfor(af::array n, NumZ) {  //std. dipole
			af::array d3 = matmul(af::conjg(psi1f(span, n)), psi1f(span, n));
			af::array d4 = -2 * af::Pi * z(n) * dZ * d3;
			dsum(n) = sum(real(d4(span, 0)));
		}
		d(0) = sum(dsum(span));

		gfor(af::array n, NumZ) {//alt. dipole
			af::array d1 = matmul(af::conjg(psi_initial(span, n)), psi1f(span, n));
			af::array d2 = -2 * af::Pi * z(n) * dZ * sqrt(matmul(af::conjg(d1(span, 0)), d1(span, 0)));
			dsum(n) = sum(real(d2(span, 0)));
		}
		d0(0) = sum(dsum(span)); // dipole with initial wf

		gfor(af::array n, NumZ) {//acc. dipole
			af::array d5 = matmul(zr3(span, n), af::conjg(psi(span, n)));
			af::array d6 = matmul(d5(span, 0), psi(span, n));
			dsum(n) = -2 * af::Pi * dZ * dR * z(n) * sum(real(d6(span, 0)));
		}
		//psi1f=matmul(psi2R,psi*dR); //convert back to Bessel coord's

		if (time >= 4 && time <= Tint) {
			dp(time - 4) = d(0);

			dp0(time - 4) = d0(0);

			dA(time - 4) = sum(dsum(span));
			dA(time - 4) = dA(time - 4) - field(0);
		}
		t += dT;
		time++;

		if (time % 25000 == 0) std::cout << lgp << ' ' << L << ' ' << bz << ' ' << inten << endl;
	}

end:

	af::array JEz = af::constant(0, Ttime, f64);
	af::array JEr = af::constant(0, Ttime, f64);
	gfor(af::array n, z0 * 2) {
		Jrout(span, n) = fft(Jrout(span, n));
	}
	gfor(af::array n, Ttime) {
		JEr(n) = sum(real(matmul(Jrout(n, span), af::conjg(Jrout(n, span))))) * dZ;
	}
	gfor(af::array n, rou0) {
		Jzt(span, n) = fft(Jzt(span, n));
	}
	gfor(af::array n, Ttime) {
		JEz(n) = sum(real(matmul(Jzt(n, span), af::conjg(Jzt(n, span))))) * dR;
	}

	ofstream psinorm(fn1);
	for (int n = 0; n != Ttime; n++) {
		psinorm << sum<double>(gNorm(n)) << '\t' << sum<double>(gp0(n)) << '\t' << sum<double>(dp(n)) << '\t' << sum<double>(JEr(n)) << '\t' << sum<double>(JEz(n)) << '\t' << sum<double>(dp0(n)) << '\t' << "0" << '\t' << "0" << '\t' << sum<double>(dA(n)) << endl;
	}
	psinorm << endl << "Completed propagation\t" << (clock() - start) / CLOCKS_PER_SEC << endl;
	psinorm.close();

	std::cout << Ttime << endl;
	std::cout << "fin" << endl;

	return 0;

}
