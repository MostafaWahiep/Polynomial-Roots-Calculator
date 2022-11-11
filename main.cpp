#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

using namespace std;

const double epsilon = 1e-15; //stopping criteria when relative less than epsilon

//functions' prototype

complex<double> SolveFor(complex<double> x, double coefficients[], int n);
/*
to solve in the input equation with x
returns F(x)
*/

complex<double> *FillIntialGuesses(double coefficients[], int n);
/*
It builds a dynamic array filled with the starting points to iterate over
The basic idea is to evenly distribute the points on a circle
*/

void iteration(double coefficients[], complex<double> *initial, int n);
/*
It applies Durand-kerner method, iterating over the initials till it stop have the solutions
*/

bool CheckValidEquation(string s);
/*
checks the validity of the input string whether it follow the criteria or not to operate with the equation correctly
*/

vector<pair<double, int>> Operatecoefficients(string s, int *o);
/*
returns a vector of pairs; first is the coefficient of order of the second which is the power
example (2.3,7) is 2.3x^7
*/

void PrintSolutions(complex<double> *initial, int n);

int main()
{
	int n = 0; //assume that the degree is zero
	string s;
	getline(cin, s, '\n');

	if ((!CheckValidEquation(s))) //check the input equation
	{
		cout << "error";
		return 0;
	}

	vector<pair<double, int>> poly = Operatecoefficients(s, &n);
	double *coefficients = new double[n + 1]{ 0 };

	//filling the coefficients
	for (int i = 0; i < poly.size(); i++)
		coefficients[n - poly[i].second] += poly[i].first;

    //printing the coefficients
	for (int i = 0; i < n + 1; i++)
		cout << coefficients[i] << " ";
	cout << endl;

	//constant equations have no solution
	if (n == 0)
		cout << "no solution";
	else
	{
	    //normalizing the coefficients with respect to the highest coefficient
	    //essential step to avoid divergence in solutions
		if (coefficients[0] != 1)
		{
			double y = coefficients[0];
			for (int i = 0; i < n + 1; i++)
			{
				coefficients[i] /= y;
			}
		}

        //reducing the equation: x^3 + 2x^2 = x^2(x+2) and solve for x+2
		while (coefficients[n] == 0)
		{
			cout << 0 << endl;
			n--;
		}

		if (n != 0)
		{
			complex<double> *intial = FillIntialGuesses(coefficients, n);
			iteration(coefficients, intial, n);
			delete[] intial;
		}
	}

	delete[] coefficients;
}

complex<double> SolveFor(complex<double> x, double coefficients[], int n)
{
	complex<double> temp = coefficients[0];
	for (int i = 1; i <= n; i++)
	{
		temp = x *temp + coefficients[i];
	}

	return temp;
}

complex<double> *FillIntialGuesses(double coefficients[], int n)
{
	complex<double> *intial = new complex < double>[n];
	double a0 = coefficients[n], an = coefficients[0];
	double r = pow(abs(a0 / an), 1.0 / n);
	double theta = (2 *M_PI) / n;
	for (int i = 0; i < n; i++)
	{
		complex<double> temp(r* cos(i *theta + theta / 4), r* sin(i *theta + theta / 4)); //
		intial[i] = temp;
	}

	return intial;
}

void iteration(double coefficients[], complex<double> *initial, int n)
{
	int itr = 0;
	for (; itr < 300; itr++)
	{
		for (int i = 0; i < n; i++)
		{
			complex<double> temp = 1;
			for (int j = 0; j < n; j++)
			{
				if (j == i)
					continue;
				temp *= initial[i] - initial[j];
			}

			initial[i] = initial[i] - SolveFor(initial[i], coefficients, n) / temp;
		}
	}

	cout << "It takes "<< itr << " Iterations\n";

	//some solutions would have small imaginary due to floating point imprecision so make them real
	for (int i = 0; i < n; i++)
	{
		if (initial[i].imag() < 1e-5 && initial[i].imag() > (-1 * 1e-5))
		{
			complex<double> temp(initial[i].real(), 0);
			initial[i] = temp;
		}
	}

    //sorting with respect to the absolute of the imaginary part so the real solution are in the begining
	sort(initial, initial + n, [](complex < double>
		const &a, complex < double>const &b)->bool{

		return abs(a.imag()) < abs(b.imag());

	});

	PrintSolutions(initial,n);
}

void PrintSolutions(complex<double> *initial, int n)
{
    for (int i = 0; i < n; i++)
	{
		cout << "X" << setw(2) << left << i + 1 << " = " << setw(10) << initial[i].real();
		if (initial[i].imag() > 0)
			cout << "+ " << initial[i].imag() << 'i';
		else if (initial[i].imag() < 0)
		{
			cout << "- " << initial[i].imag() *-1 << 'i';
		}

		cout << endl;
	}
}

bool CheckValidEquation(string s)
{
	int i = 0;
	if (s[0] == '-')
		i++;

	for (; i < s.size(); i++)
	{
	    //contains only allowed characters
		if (!(s[i] == 'x' || s[i] == '+' || s[i] == '-' || s[i] == '^' || (s[i] >= '0' && s[i] <= '9') || s[i] == ' ' || s[i] == '.'))
			return false;
        //after the power is only a number
		if (s[i] == '^' && !(s[i + 1] >= '0' && s[i + 1] <= '9'))
			return false;
        //after x is space or power or an end line
		if (s[i] == 'x' && !(s[i + 1] == ' ' || s[i + 1] == '^' || s[i + 1] == '\0'))
			return false;
        //after the sign is space and after the space is x or a number
		if ((s[i] == '+' || s[i] == '-') && ((s[i + 1] != ' ' || s[i - 1] != ' ') || !(s[i + 2] == 'x' || (s[i + 2] >= '0' && s[i + 2] <= '9'))))
			return false;
	}

	return true;
}

vector<pair<double, int>> Operatecoefficients(string s, int *o)
{
	vector<pair<double, int>> poly;
	int terms = count(s.begin(), s.end(), '+') + count(s.begin(), s.end(), '-') + 1;
	int index = 0;
	bool sign = false;	// true for -ve false for +ve
	if (s[0] == '-')
	{
		index += 1;
		sign = true;
		terms--;
	}

	while (terms-- && index < s.size())
	{
		double coefficient = 1; //a default for coefficient is 1: x^3 is 1
		int power = 0;//a default for pwoer is 0: 2 is 2x^0

		//check for the sign and jump
		if (s[index] == '-')
		{
			sign = true;
			index += 2;
		}
		else if (s[index] == '+')
		{
			sign = false;
			index += 2;
		}

        //extract the coefficients
		if ((s[index] >= '0' && s[index] <= '9') && index < s.size())
		{
			int y = index;
			while ((s[index + 1] >= '0' && s[index + 1] <= '9') || s[index + 1] == '.') //find the last number
			{
				index++;
			}

			coefficient = stod(s.substr(y, index - y + 1));
			index++;
		}

        //the same for extracting the power if exists
		if (s[index] == 'x' && index < s.size())
		{
			if (s[index + 1] == '^')
			{
				index += 2;
				int y = index;
				while ((s[index + 1] >= '0' && s[index + 1] <= '9'))
				{
					index++;
				}

				power = stoi(s.substr(y, index - y + 1));
			}
			else
				power = 1;
		}

        //keeping track of the greatest power to know the degree
		if (power > *o)
			*
			o = power;

        //jump to the next term in the equation
		index += 2;

		if (sign)
			coefficient *= -1;

		poly.push_back(make_pair(coefficient, power));
		sign = false;	// true for -ve false for +ve
	}

	return poly;
}
