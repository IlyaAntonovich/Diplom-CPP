#include "complex.h"

ostream& operator<<(ostream& __out1, const complex& __z2)
{
	__out1 << '(' << real(__z2) 
		  << ((imag(__z2) < 0)?" - i":" + i") << fabs(imag(__z2)) << ')';
	return __out1;
}
istream& operator>>(istream& __in1, complex& __z2)
{
	double __re_val1, __re_val2;
	__in1 >> __re_val1 >> __re_val2;
	__z2 = complex(__re_val1, __re_val2);
	return __in1;
}
ostream& operator <(ostream& __out1, const complex& __z2)
{
	char*a;
	double Chislo;
	
	Chislo = real(__z2);
	a = (char*)&Chislo;

	__out1.write(a,8);
	

	Chislo = imag(__z2);
	a = (char*)&Chislo;
	
	__out1.write(a,8);
	
	
	return __out1;
}
istream& operator >(istream& __in1, complex& __z2)
{


	char * a;
	double __re_val1, __re_val2;
	
	
	a = (char*)&__re_val1;

	
	__in1.read(a,8);
	
	
	a = (char*)&__re_val2;


	
	
	__in1.read(a,8);

	__z2 = complex(__re_val1, __re_val2);

	return __in1;
}

