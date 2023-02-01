#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

using namespace std;

const double epsilon = 1e-15;  // stopping criteria when relative less than epsilon

int counter;  // To maintain a count of the number of solutions that have been printed


/*******************************************************************************
 * @brief Solves a polynomial equation for a given complex value of x.
 *
 * Given an array of coefficients and a degree, this method calculates the result
 * of the polynomial equation using Horner's method, which evaluates a polynomial
 * in a numerically stable manner and in linear time.
 *
 * @param[in] x The complex value for which the equation is solved.
 * @param[in] coefficients The array of coefficients for the polynomial equation,
 *                         ordered by increasing degree.
 * @param[in] degree The degree of the polynomial equation.
 * @return The result of the polynomial equation for the given value of x.
 ******************************************************************************/
complex<double> EvaluatePolynomial(complex<double> x, double coefficients[], int degree);


/*******************************************************************************
 * @brief Calculates the initial guesses for iterating on a polynomial equation.
 *
 * Given an array of coefficients and a degree, this function calculates a set
 * of initial guesses that can be used to iterate on the solution of
 * a polynomial equation.
 * The initial guesses are calculated using a formula based on the coefficients
 * and degree of the equation.
 *
 * @param[in] coefficients The array of coefficients for the polynomial equation,
 *                         ordered by decreasing degree.
 * @param[in] degree The degree of the polynomial equation.
 * @return A pointer to an array of complex numbers representing the initial guesses.
 ******************************************************************************/
complex<double>* GetInitialGuessesForIteration(double coefficients[], int degree);


/*******************************************************************************
 * @brief Perform iterations to find the roots of a polynomial equation using
 *        Durand-Kerner method.
 *
 * Given an array of coefficients, a set of initial guesses, and degree of the polynomial,
 * this function performs a series of iterations to find the roots of the polynomial equation.
 * Currently, The number of iterations is fixed to 300, earlier update would have a stopping criteria.
 * The function has a time complexity of O(n) divisions and O(n^2) multiplications.
 * The final solutions are sorted based on the magnitude of their imaginary component.
 *
 * @param[in] coefficients The array of coefficients for the polynomial equation,
 *                         ordered by decreasing degree.
 * @param[in,out] initial A pointer to an array of complex numbers representing the initial guesses.
 *                        The function updates these values with the roots of the polynomial equation.
 * @param[in] degree The degree of the polynomial equation.
 ******************************************************************************/
void FindPolynomialRoots(double coefficients[], complex<double>* initial, int degree);


/*******************************************************************************
 * @brief Validates the input polynomial expression
 *
 * The function checks if the input polynomial expression contains only allowed characters such as
 * 'x', '+', '-', '^', numbers, and '.' (decimal point). It also checks if the power symbol '^'
 * is followed by a number, if the character 'x' is followed by a sign, power, or end line, and if
 * the sign is followed by 'x' or a number, or preceded by 'x' or a number. White spaces are removed
 * from the input expression.
 *
 * @param[in] input_polynomial The input polynomial expression to be validated
 * @return true if the input polynomial expression is valid
 * @return false otherwise
 ******************************************************************************/
bool IsValidPolynomialExpression(string& input_polynomial);


/*******************************************************************************
 * @brief Extracts the terms of a polynomial given as a string
 *
 * Given a polynomial in string format, this function extracts the terms and stores
 * them as pairs of coefficients and powers in a vector. The degree of the polynomial
 * is also determined and stored in a reference variable.
 *
 * @param[in] input_polynomial A reference to the string representing the polynomial
 * @param[in] degree A reference to the variable storing the degree of the polynomial
 * @return A vector of pairs, where each pair stores a coefficient and a power
 *
 * @note
 * The polynomial string should have the following format:
 *     [sign][coefficient][x][^power]
 * where:
 *     [sign] can be '+' or '-'
 *     [coefficient] is a double
 *     [x] is a variable
 *     [power] is an integer, only required if '^' is present
 *
 * Examples of valid polynomial strings:
 *     "3x^2-2x+1"
 *     "2.5x^3+3x^2-2x-1"
 *     "3"
 *
 * Examples of invalid polynomial strings:
 *     "3x^"
 *     "3x^2.5"
 *
 ******************************************************************************/
vector<pair<double, int>> ExtractPolynomialTerms(string &input_polynomial, int &degree);


/*******************************************************************************
 * @brief Prints the solutions to the equation in a formatted manner
 *
 * Given an array of complex numbers and the degree of the equation, this function
 * prints the solutions to the equation in a formatted manner. The output will include
 * the solution index, real part, and imaginary part.
 *
 * @param[in] initial Array of complex numbers representing the solutions to the equation
 * @param[in] degree Degree of the equation
 ******************************************************************************/
void PrintSolutions(complex<double>* initial, int degree);


int main()
{
	int degree = 0;	// Initialize degree to zero
	string input_polynomial;
	getline(cin, input_polynomial, '\n');

	if ((!IsValidPolynomialExpression(input_polynomial)))
	{
		cout << "error";
		return 0;
	}

	vector<pair<double, int>> terms = ExtractPolynomialTerms(input_polynomial, degree);

	// Allocate memory for the coefficients of the polynomial and initialize them to zero
	double *coefficients = new double[degree + 1]
	{ 0 };

	// Fill the coefficients array with values
	for (int i = 0; i < terms.size(); i++)
		coefficients[degree - terms[i].second] += terms[i].first;

	// Print the coefficients
	cout << "Coefficients: ";
	for (int i = 0; i < degree + 1; i++) cout << coefficients[i] << " ";
	cout << endl;

	// If the degree is zero, there is no solution for the polynomial
	if (degree == 0)
		cout << "no solution";
	else
	{
		// Normalize the coefficients with respect to the highest coefficient
		// This is essential to avoid divergence while iterating on the solutions
		if (coefficients[0] != 1)
		{
			double y = coefficients[0];
			for (int i = 0; i < degree + 1; i++)
			{
				coefficients[i] /= y;
			}
		}

		// Reduce the equation: x^3 + 2x^2 = x^2(x+2) and solve for x+2
		while (coefficients[degree] == 0)
		{
			cout << "X" << setw(2) << left << ++counter << " = 0" << endl;
			degree--;
		}

		// If the degree is not zero, find the roots of the polynomial
		if (degree != 0)
		{
			// Get initial guesses for the roots of the polynomial
			complex<double> *intial = GetInitialGuessesForIteration(coefficients, degree);

			// Find the roots of the polynomial
			FindPolynomialRoots(coefficients, intial, degree);

			// Deallocate memory
			delete[] intial;
		}
	}

	// Deallocate memory
	delete[] coefficients;
}

complex<double> EvaluatePolynomial(complex<double> x, double coefficients[], int degree)
{
	complex<double> result = coefficients[0];
	for (int i = 1; i <= degree; i++)
	{
		result = x *result + coefficients[i];
	}

	return result;
}

complex<double> *GetInitialGuessesForIteration(double coefficients[], int degree)
{
	complex<double> *intial = new complex < double>[degree];

	double a0 = coefficients[degree];
	double an = coefficients[0];
	double r = pow(abs(a0 / an), 1.0 / degree);
	double theta = (2 *M_PI) / degree;

	for (int i = 0; i < degree; i++)
	{
		double x = r* cos(i *theta + theta / 4);
		double y = r* sin(i *theta + theta / 4);
		complex<double> temp(x, y);
		intial[i] = temp;
	}

	return intial;
}

void FindPolynomialRoots(double coefficients[], complex<double> *initial, int degree)
{
	int iterations = 0;
	for (; iterations < 300; iterations++)
	{
		for (int i = 0; i < degree; i++)
		{
			complex<double> temp = 1;
			for (int j = 0; j < degree; j++)
			{
				if (j == i) continue;
				temp *= initial[i] - initial[j];
			}

			initial[i] = initial[i] - EvaluatePolynomial(initial[i], coefficients, degree) / temp;
		}
	}

	// Make solutions that are close to real into real solutions
	// To ensure that the imaginary component of each root is within a certain tolerance
	for (int i = 0; i < degree; i++)
	{
		if (abs(initial[i].imag()) < 1e-5)
		{
			initial[i] = complex<double> (initial[i].real(), 0);
		}
	}

	// Sort solutions based on magnitude of imaginary component
	sort(initial, initial + degree,
        [](const complex<double> &a, const complex<double> &b)
		{
			return abs(a.imag()) < abs(b.imag());
	});

	PrintSolutions(initial, degree);
}

void PrintSolutions(complex<double> *initial, int degree)
{
	for (int i = 0; i < degree; i++)
	{
		cout << "X" << setw(2) << left << ++counter << " = " << setw(10) << fixed <<
			setprecision(5) << initial[i].real();
		if (initial[i].imag() > 0)
			cout << "+ " << initial[i].imag() << 'i';
		else if (initial[i].imag() < 0)
		{
			cout << "- " << initial[i].imag() *-1 << 'i';
		}

		cout << endl;
	}
}

bool IsValidPolynomialExpression(string & input_polynomial)
{
	int start = 0;
	if (input_polynomial[0] == '-') start++;

	// Remove white spaces
	input_polynomial.erase(remove_if(input_polynomial.begin(), input_polynomial.end(),
            [](unsigned char x)
			{
				return isspace(x);
            }),input_polynomial.end());

	for (int i = start; i < input_polynomial.size(); i++)
	{
		char c = input_polynomial[i];
		// Check if character is allowed
		if (!(c == 'x' || c == '+' || c == '-' || c == '^' ||
				(c >= '0' && c <= '9') || c == '.'))
			return false;

		// Check if power is followed by a number
		if (c == '^' && !(input_polynomial[i + 1] >= '0' && input_polynomial[i + 1] <= '9'))
			return false;

		// Check if 'x' is followed by a sign, power or end line
		if (c == 'x' && !(input_polynomial[i + 1] == '+' || input_polynomial[i + 1] == '-' ||
				input_polynomial[i + 1] == '^' || input_polynomial[i + 1] == '\0'))
			return false;

		// Check if sign is followed by x or a number, or preceded by x or a number
		if ((c == '+' || c == '-') &&
			!(input_polynomial[i + 1] == 'x' ||
				(input_polynomial[i + 1] >= '0' && input_polynomial[i + 1] <= '9')) &&
			!(input_polynomial[i - 1] == 'x' ||
				(input_polynomial[i - 1] >= '0' && input_polynomial[i - 1] <= '9')))
			return false;
	}

	return true;
}

vector<pair<double, int>> ExtractPolynomialTerms(string &input_polynomial, int &degree)
{
	vector<pair<double, int>> terms;

	// count the number of terms in the polynomial
	int term_count = count(input_polynomial.begin(), input_polynomial.end(), '+') +
		count(input_polynomial.begin(), input_polynomial.end(), '-') + 1;

	int index = 0;
	bool is_negative = false;	// stores the sign of the term

	// Check if the first term is negative
	if (input_polynomial[0] == '-')
	{
		index++;
		is_negative = true;
		term_count--;
	}

	// loop through the terms and extract the coefficients and powers
	while (term_count-- && index < input_polynomial.size())
	{
		double coefficient = 1;	// default coefficient is 1
		int power = 0;	// default power is 0

		// check for the sign of the term
		if (input_polynomial[index] == '-')
		{
			is_negative = true;
			index++;
		}
		else if (input_polynomial[index] == '+')
		{
			is_negative = false;
			index++;
		}

		// extract the coefficient
		if ((input_polynomial[index] >= '0' && input_polynomial[index] <= '9') &&
			index < input_polynomial.size())
		{
			int start_index = index;

			// find the last digit of the coefficient
			while ((input_polynomial[index + 1] >= '0' &&
					input_polynomial[index + 1] <= '9') ||
				input_polynomial[index + 1] == '.')
			{
				index++;
			}

			coefficient = stod(input_polynomial.substr(start_index, index - start_index + 1));
			index++;
		}

		// extract the power
		if (input_polynomial[index] == 'x' && index < input_polynomial.size())
		{
			if (input_polynomial[index + 1] == '^')
			{
				index += 2;
				int start_index = index;

				// find the last digit of the power
				while ((input_polynomial[index + 1] >= '0' &&
						input_polynomial[index + 1] <= '9'))
				{
					index++;
				}

				power = stoi(input_polynomial.substr(start_index, index - start_index + 1));
			}
			else
			{
				power = 1;
			}
		}

		// store the greatest power to determine the degree of the polynomial
		if (power > degree)
		{
			degree = power;
		}

		// move to the next term in the polynomial
		index++;

		// store the term in the polynomial as a pair of coefficient and power
		if (is_negative)
		{
			coefficient *= -1;
		}

		terms.push_back(make_pair(coefficient, power));

		is_negative = false;	// reset the sign for the next term
	}

	return terms;
}
