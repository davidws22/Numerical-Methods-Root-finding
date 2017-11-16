// ============================================================================
// file: main.cpp
// ============================================================================
// Programmer: David Shin	
// Date: 7/13/2017	
// Class: CS301 ("Numerical Methods")
// Time: TTH 5:00-8:50PM
// Instructor: Dr. Raheja	
// Project: 1
//
// Description: Implementation of 5 root finding algorithms:
//		Bisection, NewtonRaphson, Secant, Modified Secant, False Position.
//		Calculates roots of two functions defined in the function 'CalcFunction'
//      
//
// ============================================================================

#include <cmath>
#include <iostream>
#define NMAX 100
#define TOL 0.01
#define DELTA 0.01
using namespace std;


double CalcFunction(double x)
{
	//first function (a) to be calculated
	//double result = 2*x*x*x - 11.7*x*x + 17.7*x - 5;
	double result = x + 10 - x*cosh(50/x);
	//cout << result;
	return result;
}

double CalcDerivative(double x)
{
	//double result = 6 * x*x - 23.4*x + 17.7;
	double result = 1 + ((50*sinh(50/x)/x) - cosh(50/x));
	return result;

}
void Bisection(double a, double b)
{
	//variable declarations
	double c;
	double fa, fb, fc, error = 100.0;
	double cold;
	fa = CalcFunction(a);
	fb = CalcFunction(b);
	
	//check statement
	if ((fa > 0 && fb > 0) ||
		(fa < 0 && fb < 0))
	{
		cout << "Function has same signs at a and b. Terminating program." << endl;
		return;
	}

	for (int i = 0; i < NMAX; i++)
	{
		c = (a + b) / 2.0;
		fc = CalcFunction(c);
		if (i > 0)
		{
			error = fabs(c - cold) / fabs(c);
			cout << "Iteration " << i << " error: " << error << endl;
		}
		if (fabs(error) < TOL)
		{
			cout << "Convergence at: " << c << endl;
			return;
		}
		if (fa*fc < 0)
		{
			b = c;
			fb = fc;
			cold = c;
		}
		else
		{
			a = c;
			fa = fc;
			cold = c;
		}
	}
}

void FalsePosition(double a, double b)
{
	//variable declarations
	double c = a;
	double fa, fb, fc, error = 100.0;
	double cold;
	fa = CalcFunction(a);
	fb = CalcFunction(b);

	//check statement
	if ((fa > 0 && fb > 0) ||
		(fa < 0 && fb < 0))
	{
		cout << "Function has same signs at a and b. Terminating program." << endl;
		return;
	}
	
	for (int i = 0; i < NMAX; i++)
	{
		c = (a*fb - b*fa) / (fb - fa);
		fc = CalcFunction(c);
		if (i > 0)
		{
			error = fabs(c - cold) / fabs(c);
			cout << "Iteration " << i << " error: " << error << endl;
		}
		if (fabs(error) < TOL)
		{
			cout << "Convergence at: " << c << endl;
			return;
		}
		if (fc == 0)
		{
			cout << "Convergence at: " << c << endl;
			return;
		}
		else if (fc*fa < 0)
		{
			b = c;
			cold = c;
		}
		else
		{
			a = c;
			cold = c;
		}	
	}
}

void NewtonRaphson(double x)
{
	//variable declarations
	double fx, fp;
	double d;
	fx = CalcFunction(x);
	for (int i = 0; i < NMAX; i++)
	{
		fp = CalcDerivative(x);
		if (fabs(fp) < DELTA)
		{
			cout << "Small derivative" << endl;
			return;
		}
		d = fx / fp;
		x = x - d;
		fx = CalcFunction(x);
		if (i > 0)
		{
			cout << "Iteration " << i << " error: " << fabs(d) << endl;
		}
		if (fabs(d) < TOL)
		{
			cout << "Convergence to " << x << endl;
			return;
		}
	}
}

void Secant(double a, double b)
{
	//variable declarations
	double fa, fb, d;
	double temp;
	fa = CalcFunction(a);
	fb = CalcFunction(b);
	double old;

	if (fabs(fa) > fabs(fb))
	{
		temp = a;
		a = b;
		b = temp;
		temp = fa;
		fa = fb;
		fb = temp;
	}

	cout << "0 " << a << " " << fa << endl;
	
	cout << "1 " << b << " " << fb << endl;
	cout << "Iteration 1 " << "error: " << fabs(b - a) / fabs(b) << endl;
	old = b;
	
	for (int i = 1; i < NMAX; i++)
	{
		if (fabs(fa) > fabs(fb))
		{
			temp = a;
			a = b;
			b = temp;
			temp = fa;
			fa = fb;
			fb = temp;
		}
		d = (b - a) / (fb - fa);
		b = a;
		fb = fa;
		d = d*fa;
		
		if (fabs(d) < TOL)
		{
			cout << "Convergence to " << fabs(a) << endl;
			cout << "Iteration " << i + 1 << " error: " << (fabs((a-d) - a)) / fabs(a-d) << endl;
			
			return;
		}
		a = a - d;
		fa = CalcFunction(a);
		cout << i + 1 << " " << a << " " <<  fa << endl;
		cout << "Iteration " << i + 1 << " error: " << (fabs(a - old)) / fabs(a) << endl;
		old = a;
	}
}

void ModifiedSecant(double x)
{
	//variable declarations
	double fx, fb, d, temp;
	double old;
	fx = CalcFunction(x);
	fb = CalcFunction((x + DELTA*x));
	old = x - ((DELTA*x*fx) / (fb - fx));

	for (int i = 1; i < NMAX; i++)
	{
		x = old;
		fx = CalcFunction(x);
		fb = CalcFunction(x + DELTA*x);
		d = x - ((DELTA*x*fx) / (fb - fx));

		cout << "Iteration " << i << " error: " << fabs(d - old)/fabs(d)  << endl;
		
		if ((fabs(d-old)/fabs(d)) < TOL)
		{
			cout << "Convergence to " << fabs(d) << endl;
			return;
		}
		old = d;
	}
}

int main()
{
	double a, b;

	cout.setf(ios::fixed);
	cout.setf(ios::showpoint);
	cout.precision(8);
	/*
	Bisection(0, 1);
	cout << "Second bisection" << endl;
	Bisection(1, 2);
	cout << "Third Bisection" << endl;
	Bisection(3, 4);
	//Bisection(100, 200);
	
	FalsePosition(0, 1);
	FalsePosition(1, 2);
	FalsePosition(3, 4);
	

	NewtonRaphson(1);
	NewtonRaphson(2);
	NewtonRaphson(4);
	

	Secant(-1, 1);
	Secant(1, 2);
	Secant(3, 4);

	ModifiedSecant(1);
	ModifiedSecant(2);
	ModifiedSecant(4);


	*/
	//Bisection(120, 130);
	//FalsePosition(120, 130);
	//NewtonRaphson(120);
	//Secant(120, 130);
	//ModifiedSecant(120);
	return 0;
}
