#include"transform.h"
#include"cal.h"

MJD UTC2MJD(UTC utc)//常规时间转化至简化儒略日
{
	MJD mjd;
	double y, m;
	if (utc.Month > 2)
	{
		y = utc.Year;
		m = utc.Month;
	}
	else
	{
		y = utc.Year - 1;
		m = utc.Month + 12;
	}
	double UT = (double)utc.Hour + (double)utc.Min / 60 + (double)utc.Sec / 3600;
	double JD = (int)(365.25 * y) + (int)(30.6001 * (m + 1)) + utc.Day + UT / 24 + 1720981.5;
	mjd.FracDay = JD - 2400000.5;
	mjd.Days = (int)mjd.FracDay;
	mjd.FracDay = mjd.FracDay - mjd.Days;
	return mjd;
}

UTC MJD2UTC(MJD mjd)//简化儒略日转通用时
{
	UTC utc;
	double JD = (double)mjd.Days + mjd.FracDay + 2400000.5;
	int a = (int)(JD + 0.5);
	int b = a + 1537;
	int c = (int)((b - 122.1) / 365.25);
	int d = (int)(365.25 * c);
	int e = (int)((b - d) / 30.6001);
	double D = b - d - (int)(30.6001 * e) + JD + 0.5 - a;
	utc.Day = (int)D;
	double h = 24 * (D - utc.Day);
	utc.Hour = (int)(h);
	double m = 60 * (h - utc.Hour);
	utc.Min = (int)(m);
	double s = 60 * (m - utc.Min);
	utc.Sec = (int)(s);
	utc.Month = e - 1 - 12 * (int)(e / 14);
	utc.Year = c - 4715 - (int)((7 + utc.Month) / 10);
	return utc;
}

GPSTIME MJD2GPSTIME(MJD mjd)
{
	GPSTIME gpst;
	int gpsw = int((mjd.Days + mjd.FracDay - 44244) / 7);
	double gpss = double((mjd.Days + mjd.FracDay - 44244 - gpsw * 7)) * 86400.0;
	gpst.Week = gpsw;
	gpst.SecOfWeek = gpss;
	return gpst;
}

MJD GPSTIME2MJD(GPSTIME gpst)
{
	MJD mjd;
	double Mjd = 44244 + gpst.Week * 7 + gpst.SecOfWeek / 86400.0;
	mjd.Days = int(Mjd);
	mjd.FracDay = double(Mjd - mjd.Days);
	return mjd;
}

GPSTIME GPSTIME2BDSTIME(GPSTIME gpst)
{
	GPSTIME bdst;
	bdst.Week = gpst.Week - 1356;
	bdst.SecOfWeek = gpst.SecOfWeek - 14;
	return bdst;
}


GPSTIME BDSTIME2GPSTIME(GPSTIME bdst)
{
	GPSTIME gpst;
	gpst.Week = bdst.Week + 1356;
	gpst.SecOfWeek = bdst.SecOfWeek + 14;
	return gpst;
}
/*****************************************/
// 此处考虑到，GPS时转化到MJD，再由MJD转化到UTC容易因为浮点数而丢失秒的数据，故对算法进行优化，在秒数上只进行整数运算
UTC GPSTIME2UTC(GPSTIME gpst)
{
	MJD mjd = GPSTIME2MJD(gpst);
	UTC u = MJD2UTC(mjd);
	UTC utc;
	utc.Year = u.Year;
	utc.Month = u.Month;
	int day = gpst.SecOfWeek / 86400;
	utc.Day = u.Day + day;
	int hour = (gpst.SecOfWeek - 86400 * day) / 3600;
	utc.Hour = hour;
	int min = (gpst.SecOfWeek - 86400 * day - 3600 * hour) / 60;
	utc.Min = min;
	utc.Sec = gpst.SecOfWeek - 86400 * day - 3600 * hour - 60 * min;
	return utc;
}

XYZ BLH2XYZ(BLH blh, double e2, double a)
{
	XYZ xyz;
	double B = blh.Lat * Pi / 180;
	double L = blh.Lon * Pi / 180;
	double N = a / sqrt(1 - e2 * sin(B) * sin(B));
	xyz.X = (blh.Height + N) * cos(B) * cos(L);
	xyz.Y = (blh.Height + N) * cos(B) * sin(L);
	xyz.Z = (N * (1 - e2) + blh.Height) * sin(B);
	return xyz;
}

BLH XYZ2BLH(XYZ xyz, double e2, double a)
{
	BLH blh;
	double L = 180 * atan2(xyz.Y, xyz.X) / Pi;
	while (L < 0 || L > 360)
	{
		if (L < 0)
		{
			L += 360;
		}
		if (L > 360)
		{
			L -= 360;
		}
	}
	blh.Lon = L;

	double B = atan(xyz.Z / sqrt(xyz.X * xyz.X + xyz.Y * xyz.Y));
	double B0 = 0;
	double n;
	int count = 0;
	do
	{
		B0 = B;
		n = a / sqrt(1 - e2 * sin(B0) * sin(B0));
		B = atan((xyz.Z + n * e2 * sin(B0)) / sqrt(xyz.X * xyz.X + xyz.Y * xyz.Y));
		count++;
	} while (abs(B - B0) > 1e-10 && count < 100);
	n = a / sqrt(1 - e2 * sin(B) * sin(B));

	blh.Lat = 180 * B / Pi;

	blh.Height = xyz.Z / sin(B) - n * (1 - e2);
	return blh;
}

//转化为当地坐标系
XYZ XYZ2ENU(XYZ xyz1, XYZ xyz2, int sys)
{
	XYZ enu;
	BLH blh;
	switch (sys)
	{
	case SYS_GPS:
		blh=XYZ2BLH(xyz1, WGS84_e2, WGS84_a);
		break;
	case SYS_BDS:
		blh=XYZ2BLH(xyz1, CGCS2000_e2, CGCS2000_a);
		break;
	default:
		break;
	}
	double B = degree2rad(blh.Lat);
	double L = degree2rad(blh.Lon);
	MatrixXd R = get_Rot(B, L);

	MatrixXd x0(3, 1);
	MatrixXd x_local(3, 1);
	x0(0, 0) = xyz2.X - xyz1.X;
	x0(1, 0) = xyz2.Y - xyz1.Y;
	x0(2, 0) = xyz2.Z - xyz1.Z;
	x_local = R * x0;
	enu.X = x_local(0, 0);
	enu.Y = x_local(1, 0);
	enu.Z = x_local(2, 0);
	return enu;
}

MatrixXd get_Rot(double B, double L)
{
	MatrixXd R(3, 3);
	R(0, 0) = -sin(L);
	R(0, 1) = cos(L);
	R(0, 2) = 0;
	R(1, 0) = -sin(B) * cos(L);
	R(1, 1) = -sin(B) * sin(L);
	R(1, 2) = cos(B);
	R(2, 0) = cos(B) * cos(L);
	R(2, 1) = cos(B) * sin(L);
	R(2, 2) = sin(B);
	return R;
}

MatrixXd get_XYZ(XYZ xyz)
{
	MatrixXd x(3, 1);
	x << xyz.X, xyz.Y, xyz.Z;
	return x;
}

MatrixXd get_BLH(BLH blh)
{
	MatrixXd b(3, 1);
	b << blh.Lat, blh.Lon, blh.Height;
	return b;
}

XYZ get_XYZ(MatrixXd xyz)
{
	return XYZ(xyz(0, 0), xyz(1, 0), xyz(2, 0));
}

BLH get_BLH(MatrixXd blh)
{
	return BLH(blh(0, 0), blh(1, 0), blh(2, 0));
}