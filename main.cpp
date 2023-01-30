#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

using namespace std;

const double epsilon = 1e-15; //stopping criteria when relative less than epsilon

int counter; // To maintain a count of the number of solutions that have been printed

//functions' prototype

complex<double> SolveFor(complex<double> x, double coefficients[], int degree);
/*
to solve in the input equation with x
returns F(x)
*/

complex<double> *FillIntialGuesses(double coefficients[], int degree);
/*
It builds a dynamic array filled with the starting points to iterate over
The basic idea is to evenly distribute the points on a circle
*/

void iteration(double coefficients[], complex<double> *initial, int degree);
/*
It applies Durand-kerner method, iterating over the initials till it stop have the solutions
*/

bool CheckValidEquation(string &input_polynomial);
/*
checks the validity of the input string whether it follow the criteria or not to operate with the equation correctly
*/

vector<pair<double, int>> Operatecoefficients(string &input_polynomial, int *degree);
/*
returns a vector of pairs; first is the coefficient of order of the second which is the power
example (2.3,7) is 2.3x^7
*/

void PrintSolutions(complex<double> *initial, int degree);

int main()
{
	int degree = 0; //assume that the degree is zero
	string input_polynomial;
	getline(cin, input_polynomial, '\n');

	if ((!CheckValidEquation(input_polynomial))) //check the input equation
	{
		cout << "error";
		return 0;
	}

	vector<pair<double, int>> poly = Operatecoefficients(input_polynomial, &degree);
	double *coefficients = new double[degree + 1]{ 0 };

	//filling the coefficients
	for (int i = 0; i < poly.size(); i++)
		coefficients[degree - poly[i].second] += poly[i].first;


    //printing the coefficients
    cout << "Coefficients: ";
	for (int i = 0; i < degree + 1; i++)
		cout << coefficients[i] << " ";
	cout << endl;

	//constant equations have no solution
	if (degree == 0)
		cout << "no solution";
	else
	{
	    //normalizing the coefficients with respect to the highest coefficient
	    //essential step to avoid divergence in solutions
		if (coefficients[0] != 1)
		{
			double y = coefficients[0];
			for (int i = 0; i < degree + 1; i++)
			{
				coefficients[i] /= y;
			}
		}

        //reducing the equation: x^3 + 2x^2 = x^2(x+2) and solve for x+2
		while (coefficients[degree] == 0)
		{
			cout << "X" << setw(2) << left << ++counter << " = 0" << endl;
			degree--;
		}

		if (degree != 0)
		{
			complex<double> *intial = FillIntialGuesses(coefficients, degree);
			iteration(coefficients, intial, degree);
			delete[] intial;
		}
	}

	delete[] coefficients;
}

complex<double> SolveFor(complex<double> x, double coefficients[], int degree)
{
	complex<double> temp = coefficients[0];
	for (int i = 1; i <= degree; i++)
	{
		temp = x *temp + coefficients[i];
	}

	return temp;
}

complex<double> *FillIntialGuesses(double coefficients[], int degree)
{
	complex<double> *intial = new complex < double>[degree];
	double a0 = coefficients[degree], an = coefficients[0];
	double r = pow(abs(a0 / an), 1.0 / degree);
	double theta = (2 *M_PI) / degree;
	for (int i = 0; i < degree; i++)
	{
		complex<double> temp(r* cos(i *theta + theta / 4), r* sin(i *theta + theta / 4)); //
		intial[i] = temp;
	}

	return intial;
}

void iteration(double coefficients[], complex<double> *initial, int degree)
{
	int itr = 0;
	for (; itr < 300; itr++)
	{
		for (int i = 0; i < degree; i++)
		{
			complex<double> temp = 1;
			for (int j = 0; j < degree; j++)
			{
				if (j == i)
					continue;
				temp *= initial[i] - initial[j];
			}

			initial[i] = initial[i] - SolveFor(initial[i], coefficients, degree) / temp;
		}
	}

	//cout << "It takes "<< itr << " Iterations\n";

	//some solutions would have small imaginary due to floating point imprecision so make them real
	for (int i = 0; i < degree; i++)
	{
		if (initial[i].imag() < 1e-5 && initial[i].imag() > (-1 * 1e-5))
		{
			complex<double> temp(initial[i].real(), 0);
			initial[i] = temp;
		}
	}

    //sorting with respect to the absolute of the imaginary part so the real solution are in the begining
	sort(initial, initial + degree, [](complex < double>
		const &a, complex < double>const &b)->bool{

		return abs(a.imag()) < abs(b.imag());

	});

	PrintSolutions(initial, degree);
}

void PrintSolutions(complex<double> *initial, int degree)
{
    for (int i = 0; i < degree; i++)
	{
		cout << "X" << setw(2) << left << ++counter << " = " << setw(10) << fixed <<setprecision(5) << initial[i].real();
		if (initial[i].imag() > 0)
			cout << "+ " << initial[i].imag() << 'i';
		else if (initial[i].imag() < 0)
		{
			cout << "- " << initial[i].imag() *-1 << 'i';
		}

		cout << endl;
	}
}

bool CheckValidEquation(string &input_polynomial)
{
	int i = 0;
	if (input_polynomial[0] == '-')
		i++;

	for (; i < input_polynomial.size(); i++)
	{
	    //contains only allowed characters
		if (!(input_polynomial[i] == 'x' || input_polynomial[i] == '+' ||
                input_polynomial[i] == '-' || input_polynomial[i] == '^' ||
                (input_polynomial[i] >= '0' && input_polynomial[i] <= '9') || input_polynomial[i] == '.'))
			return false;
        //after the power is only a number
		if (input_polynomial[i] == '^' && !(input_polynomial[i + 1] >= '0' && input_polynomial[i + 1] <= '9'))
			return false;
        //after x is space or power or an end line
		if (input_polynomial[i] == 'x' && !(input_polynomial[i] == '+' || input_polynomial[i] == '-'
                || input_polynomial[i + 1] == '^' || input_polynomial[i + 1] == '\0'))
			return false;
        //after the sign is space and after the space is x or a number
		if ((input_polynomial[i] == '+' || input_polynomial[i] == '-')
            && !(input_polynomial[i + 1] == 'x' || (input_polynomial[i + 1] >= '0' && input_polynomial[i + 1] <= '9'))
            && !(input_polynomial[i - 1] == 'x' || (input_polynomial[i - 1] >= '0' && input_polynomial[i - 1] <= '9')))
			return false;
	}

	return true;
}

vector<pair<double, int>> Operatecoefficients(string &input_polynomial, int *degree)
{
	vector<pair<double, int>> poly;
	int terms = count(input_polynomial.begin(), input_polynomial.end(), '+') + count(input_polynomial.begin(), input_polynomial.end(), '-') + 1;
	int index = 0;
	bool sign = false;	// true for -ve false for +ve
	if (input_polynomial[0] == '-')
	{
		index += 1;
		sign = true;
		terms--;
	}

	while (terms-- && index < input_polynomial.size())
	{
		double coefficient = 1; //a default for coefficient is 1: x^3 is 1
		int power = 0;//a default for pwoer is 0: 2 is 2x^0

		//check for the sign and jump
		if (input_polynomial[index] == '-')
		{
			sign = true;
			index += 2;
		}
		else if (input_polynomial[index] == '+')
		{
			sign = false;
			index += 2;
		}

        //extract the coefficients
		if ((input_polynomial[index] >= '0' && input_polynomial[index] <= '9') && index < input_polynomial.size())
		{
			int y = index;
			while ((input_polynomial[index + 1] >= '0' && input_polynomial[index + 1] <= '9') || input_polynomial[index + 1] == '.') //find the last number
			{
				index++;
			}

			coefficient = stod(input_polynomial.substr(y, index - y + 1));
			index++;
		}

        //the same for extracting the power if exists
		if (input_polynomial[index] == 'x' && index < input_polynomial.size())
		{
			if (input_polynomial[index + 1] == '^')
			{
				index += 2;
				int y = index;
				while ((input_polynomial[index + 1] >= '0' && input_polynomial[index + 1] <= '9'))
				{
					index++;
				}

				power = stoi(input_polynomial.substr(y, index - y + 1));
			}
			else
				power = 1;
		}

        //keeping track of the greatest power to know the degree
		if (power > *degree)
			*degree = power;

        //jump to the next term in the equation
		index += 2;

		if (sign)
			coefficient *= -1;

		poly.push_back(make_pair(coefficient, power));
		sign = false;	// true for -ve false for +ve
	}

	return poly;
}
