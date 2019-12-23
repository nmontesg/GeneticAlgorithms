#pragma once

int s(double x, double y, double z, double t) {
	double a = (x + y) / (2147483648UL + y + z) + 0.86525597943226508722;
	double result;
	modf(a, &result);
	return (int)result;
}

double phi(double x, double y, double z, double t, double a, double b, double c, double d) {
	double result = 8796093022208ULL - pow(x - a, 2) - pow(y - b, 2) - pow(z - c, 2) - pow(t - d, 2);
	return result;
}

double f(double x, double y, double z, double t) {
	int sign = s(x, y, z, t);
	double term1 = phi(x, y, z, t, 1237566.4, 54783217.5, 1237896431.1, 325123467.37);
	double term2 = phi(x, y, z, t, 5674235.4, 4067231567.2, 13245678.3, 3748967543.2);
	double term3 = phi(x, y, z, t, 3867435523.2, 7134893.75, 3565897564.1, 15675987.34);
	double term4 = phi(x, y, z, t, 4000223567.09, 3734098765.4, 3367981234.4, 4067231567.25);
	double result = sin(x + y + z + t) * term1 * term2 * term3 * term4;
	if (sign % 2 == 0) return result;
	else return (-result);
}
