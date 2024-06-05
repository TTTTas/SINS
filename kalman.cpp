#include "kalman.h"
#include"Atomsphere.h"

double KalmanFilter(OnePoch* op, KalRes* kalres)
{
	Kalman kal;
	Matrix<double, 3, 1>SP;
	Matrix<double, 3, 1>SV;
	//Ìî³ä¿¨¶ûÂüÂË²¨¾ØÕó
	Matrix<double, 3, 3>eye;
	Matrix<double, 8, 8>eye1;
	Matrix<double, Dynamic, 8>H;
	Matrix<double, Dynamic, Dynamic>R;
	Matrix<double, 3, 1>dr;
	Matrix<double, 3, 3>I;
	Matrix<double, Dynamic, Dynamic>K;
	Matrix<double, Dynamic, 1>dz;
	H.resize(2 * op->b, 8);
	R.resize(2 * op->b, 2 * op->b);
	K.resize(2 * op->b, 2 * op->b);
	dz.resize(2 * op->b, 1);
	eye.setIdentity();
	eye1.setIdentity();
	kal.OME << 0, -bomega, 0,
		bomega, 0, 0,
		0, 0, 0;
	kal.P = kalres->P;
	kal.X << kalres->X, kalres->Y, kalres->Z, kalres->Vx, kalres->Vy, kalres->Vz, kalres->b, kalres->db;
	kal.phi.setIdentity();
	kal.phi(0, 3) = 1;
	kal.phi(1, 4) = 1;
	kal.phi(2, 5) = 1;
	kal.phi(6, 7) = 1;
	kal.Q.setZero();
	kal.Q.block<3, 3>(0, 0) = eye * 10 / 3;
	kal.Q.block<3, 3>(0, 3) = eye * 10 / 2;
	kal.Q.block<3, 3>(3, 0) = eye * 10 / 2;
	kal.Q.block<3, 3>(3, 3) = eye * 10;
	kal.Q(6, 6) = 1 / 3 + 1;
	kal.Q(6, 7) = 0.5;
	kal.Q(7, 6) = 0.5;
	kal.Q(7, 7) = 1;
	H.setZero();
	R.setZero();


	kalres->Gtime = op->GpsT;
	//Kalman filter
	kal.X = kal.phi * kal.X;
	kal.P = kal.phi * kal.P * kal.phi.transpose() + kal.Q;

	int j = 0;
	for (int i = 0; i < op->bnum; i++)
	{
		if (op->bObs[i].Valid != -1)
		{
			int PRN = op->bObs[i].PRN;
			GNSSInfo B1 = op->bInfo[PRN];
			if (op->bInfo[PRN].flag == true)
			{
				if (j >= op->b)
				{
					return 0;
				}
				SP << B1.X, B1.Y, B1.Z;
				SV << B1.Vx, B1.Vy, B1.Vz;
				dr = SP - kal.X.block<3, 1>(0, 0);
				double r = sqrt((dr.transpose() * dr)(0, 0));
				XYZ xyz1(B1.X, B1.Y, B1.Z);
				NEU neu;
				XYZ xyz0(kalres->X, kalres->Y, kalres->Z);
				XYZ2NEU(WGS84, &xyz1, &xyz0, &neu);
				double Tr = Hopefield(StA, neu);
				double Pred = (op->bObs[i].IF - FB1 * FB1 * Cv * op->Bdse[PRN].tgd1 / (FB1 * FB1 - FB3 * FB3)) - kal.X(6, 0) - r - Tr + Cv * op->bInfo[PRN].dt;
				H.block<1, 3>(j, 0) = dr.transpose() / r;
				double rdot = ((H.block<1, 3>(j, 0)) * (eye * (SV + kal.OME * SP) - (kal.X.block<3, 1>(3, 0) + kal.OME * (kal.X.block<3, 1>(0, 0)))))(0, 0);
				double Pred1 = rdot + kal.X(7, 0);
				H.block<1, 3>(j, 0) = -H.block<1, 3>(j, 0);
				H.block<1, 3>(j + op->b, 3) = H.block<1, 3>(j, 0);
				H(j, 6) = 1;
				H(j + op->b, 7) = 1;
				R(j, j) = op->bObs[i].raw[0].psrsd * op->bObs[i].raw[0].psrsd;
				R(j + op->b, j + op->b) = 0.05 * 0.05;
				dz(j, 0) = Pred;
				dz(j + op->b, 0) = -op->bObs[i].raw[0].dopp * WLB1 - Pred1;
				j++;
			}
		}
	}
	K = kal.P * H.transpose() * (H * kal.P * H.transpose() + R).inverse();
	kal.X = kal.X + K * dz;
	kal.P = (eye1 - K * H) * kal.P;

	kalres->P = kal.P;
	kalres->X = kal.X(0, 0);
	kalres->Y = kal.X(1, 0);
	kalres->Z = kal.X(2, 0);
	kalres->Vx = kal.X(3, 0);
	kalres->Vy = kal.X(4, 0);
	kalres->Vz = kal.X(5, 0);
	kalres->b = kal.X(6, 0);
	kalres->db = kal.X(7, 0);
	BLH blh0;
	XYZ xyz1(kalres->X, kalres->Y, kalres->Z);
	XYZ2BLH(WGS84, &blh0, &xyz1);
	kalres->blh = blh0;
	NEU neu0;
	XYZ xyz(-2267804.526, 5009342.372, 3220991.863);
	XYZ2NEU(WGS84, &xyz1, &xyz, &neu0);
	kalres->neu = neu0;

	return 1;
}

int SaveTofileKal(ofstream& fout, KalRes r)
{
	fout << setw(5) << setfill(' ') << fixed << setprecision(0) << setiosflags(ios::right) << r.Gtime.Week;
	fout << setw(12) << setfill(' ') << fixed << setprecision(3) << r.Gtime.WeekSec;
	fout << setw(14) << setfill(' ') << fixed << setprecision(4) << r.X;
	fout << setw(14) << setfill(' ') << fixed << setprecision(4) << r.Y;
	fout << setw(14) << setfill(' ') << fixed << setprecision(4) << r.Z;
	fout << setw(14) << setfill(' ') << fixed << setprecision(3) << r.neu.e;
	fout << setw(14) << setfill(' ') << fixed << setprecision(3) << r.neu.n;
	fout << setw(14) << setfill(' ') << fixed << setprecision(3) << r.neu.u;
	fout << setw(13) << setfill(' ') << fixed << setprecision(8) << r.blh.b / PI * 180.0;
	fout << setw(13) << setfill(' ') << fixed << setprecision(8) << r.blh.l / PI * 180.0;
	fout << setw(13) << setfill(' ') << fixed << setprecision(8) << r.blh.h;
	fout << setw(10) << setfill(' ') << fixed << setprecision(3) << r.Vx;
	fout << setw(10) << setfill(' ') << fixed << setprecision(3) << r.Vy;
	fout << setw(10) << setfill(' ') << fixed << setprecision(3) << r.Vz << endl;

	return 1;
}