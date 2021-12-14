#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>

const double PI = 3.141592653589793238;

using namespace std;

double function(double x) {
	return x * sin(x) * cos(2.0 * x);
}
double simpleIterFunction(double x) { // how
	return sin(x) * cos(2.0 * x)*2;
}
//Метод дихотомии
double BisectionMethod(double a, double b, double eps, int& iterations)
{
	double x;
	iterations = 0;
	while (b - a > eps) {
		x = (a + b) / 2;
		if ((function(x) + eps * eps) * (function(a)+eps*eps) > 0.0) {
			a = x;
		}
		else {
			b = x;
		}
		iterations++;
	}
	return (a + b) / 2;
}
//Метод хорд
double ChordMethod(double a, double b, double epsilon, int& iterations) {
	a += 0.1;
	double x0 = b;
	double x1 = x0 - function(x0) * (x0 - a) / (function(x0) - function(a));
	while (abs(x1 - x0) > epsilon)
	{
		x0 = x1;
		x1 = x0 - function(x0) * (x0 - b) / (function(x0) - function(b));
		iterations++;
	}
	return x0;
}
//Метод Простых итераций
double SimpleIterMethod(double a, double epsilon, int& iterations) {
	double b;
	iterations = 0;
	for (;;) {
		b = simpleIterFunction(a);
		if (abs(a - b) < epsilon || iterations > 1000) { break; }
		a = b;
		iterations++;
	}
	return b;
}

double derivativeFunction(double x) {
	return (sin(x) + x * cos(x)) * cos(2 * x) - 2 * x * sin(x) * sin(2 * x);
}
double NewtonMethod(double a, double epsilon, int& iterations) {
	double b;
	iterations = 0;
	for (;;) {
		b = a - function(a) / derivativeFunction(a);
		if (abs(a - b) < epsilon) break;
		a = b;
		iterations++;
	}
	return -b;
}
//Метод секущих
double SecantMethod(double a, double c, double epsilon, int& iterations) {
	a = 0.4;
	c = 2;
	double b;
	iterations = 0;
	for (;;) {
		b = a - function(a) / (function(a) - function(c)) * (a - c);
		if (abs(a - b) < epsilon) break;
		c = a;
		a = b;
		iterations++;
	}
	return -b;
}

void task_1() {
	double error = 1e-7;
	double exactAnswer = PI/4;
	std::vector<int> iterations(5);
	std::vector<double> results(5);
	std::cout << "Method name \t\tResult \t\t\tError margin \t\t Number of iterarions" << std::endl;
	results[0] = BisectionMethod(0, 2, error, iterations[0]);
	results[1] = ChordMethod(0, 2, error, iterations[1]);
	results[2] = SimpleIterMethod(1, error, iterations[2]);
	results[3] = NewtonMethod(0.5, error, iterations[3]);
	results[4] = SecantMethod(0.01, 2, error, iterations[4]);
	std::cout << std::fixed << std::setprecision(20) << "Bisection         \t" << results[0] << "\t" << abs(exactAnswer - results[0]) << "\t " << iterations[0] << std::endl;
	std::cout << std::fixed << std::setprecision(20) << "Chord             \t" << results[1] << "\t" << abs(exactAnswer - results[1]) << "\t " << iterations[1] << std::endl;
	std::cout << std::fixed << std::setprecision(20) << "Simple Iteration  \t" << results[2] << "\t" << abs(exactAnswer - results[2]) << "\t " << iterations[2] << std::endl;
	std::cout << std::fixed << std::setprecision(20) << "Newton            \t" << results[3] << "\t" << abs(exactAnswer - results[3]) << "\t " << iterations[3] << std::endl;
	std::cout << std::fixed << std::setprecision(20) << "Secant(sekushaya) \t" << results[4] << "\t" << abs(exactAnswer - results[4]) << "\t " << iterations[4] << std::endl;

}

#pragma endregion
#pragma region TASK2
double function2x(double x, double y) {
	return (y * y * exp(1 - x) - 1) / (2 * exp(1 - y));

}
double function2y(double x, double y) {
	return (x * x * exp(1 - y)) / (2 * exp(1 - x));
}


std::vector<double> SimpleIterMethod2(double a, double b, double epsilon, int& iterations) {
	double x, y;
	x = a;
	y = b;
	iterations = 0;
	std::vector<double> answer(2);
	do
	{
		iterations += 1;
		a = x;
		b = y;
		x = function2x(a, b);
		y = function2y(a, b);

	} while (abs(x - a) + abs(y - b) > 0.1);
	answer[0] = x;
	answer[1] = y;
	return answer;
}
void task_2() {

	double error = 1e-6;
	double a = 0.1; double b = 0.1;
	int iterations;
	std::vector<double> answer = SimpleIterMethod2(a, b, error, iterations);
	std::cout << std::fixed << std::setprecision(1) << "Initial approximation : x = " << a << ", y = " << b << std::endl;
	std::cout << std::fixed << std::setprecision(10) << "Simple iteration method \t" << "x = " << answer[0] << "\t" << "y = " << answer[1] << std::endl;
	std::cout << "Needed iterations: " << iterations << std::endl;

}

#pragma endregion

/*#pragma region TASK3


void function3(std::vector<double>& fx, std::vector<double>& x)
{
	fx[0] = f1(x[0], x[1]);
	fx[1] = f2(x[0], x[1]);
}

void createJacobi(std::vector<double>& jacobiM, std::vector<double>& x)
{
	jacobiM[0] = f1x(x[0], x[1]);
	jacobiM[1] = f1y(x[0], x[1]);
	jacobiM[2] = f2x(x[0], x[1]);
	jacobiM[3] = f2y(x[0], x[1]);
}

void inverseJacobi(std::vector<double>& jacobiM)
{
	double det = jacobiM[0] * jacobiM[3] - jacobiM[2] * jacobiM[1];
	double temp = jacobiM[0];
	jacobiM[0] = jacobiM[3];
	jacobiM[3] = temp;
	jacobiM[1] = -jacobiM[1];
	jacobiM[2] = -jacobiM[2];
	for (int i = 0; i < 4; i++) {
		jacobiM[i] = jacobiM[i] / det;
	}
}

std::vector<double>  productM(std::vector<double>& jacobiM, std::vector<double>& x) {
	std::vector <double> result(x.size());
	result[0] = jacobiM[0] * x[0] + jacobiM[1] * x[1];
	result[1] = jacobiM[2] * x[0] + jacobiM[3] * x[1];
	return result;
}

std::vector<double> minusV(std::vector<double>& x, std::vector<double>& p) {
	std::vector<double> result(x.size());
	for (int i = 0; i < x.size(); i++) {
		result[i] = x[i] - p[i];
	}

	return result;
}
std::vector<double> newtIter(std::vector<double> x) {
	std::vector<double> result(2);
	std::vector<double> product(2);
	std::vector<double> jacobiM(4);
	std::vector<double> fx(2);
	function3(fx, x);
	createJacobi(jacobiM, x);
	inverseJacobi(jacobiM);
	product = productM(jacobiM, fx);
	result = minusV(x, product);
	return result;
}
void task_3() {
	double a = 0.1, b = 0.1;
	std::vector<double> x = { a,b };
	double eps = 1e-6;
	double error = 1;
	double preverror = 0;
	std::vector<double>  currX;
	int iterations = 0;
	while (1) {
		currX = newtIter(x);
		if (abs(x[0] - currX[0]) + abs(x[1] - currX[1]) < eps)
			break;
		else
			x = currX;
		iterations++;
	}

	std::cout << "Initial approximation: x = " << std::fixed << std::setprecision(1) << a << ", y = " << b << std::endl;
	std::cout << "Newton method: " << "x = " << std::fixed << std::setprecision(10) << x[0] << "\t" << "y = " << x[1] << std::endl;
	std::cout << "Needed iterations: " << iterations << std::endl;
}

#pragma endregion */


double f(double x, double y)
{
	return pow(x, 2) * exp(1 - y) + pow(y, 2) * exp(1 - x) + x;
}

double f1(double x, double y)
{
	return 2 * x * exp(1 - y) - pow(y, 2) * exp(1 - x) + 1;
}

double f2(double x, double y)
{
	return -pow(x, 2) * exp(1 - y) + 2 * y * exp(1 - x);
}

double df1dx(double x, double y)
{
	return 2 * exp(1 - y) + pow(y, 2) * exp(1 - x);
}

double df1dy(double x, double y)
{
	return -2 * x * exp(1 - y) - 2 * y * exp(1 - x);
}

double df2dx(double x, double y)
{
	return -2 * x * exp(1 - y) - 2 * y * exp(1 - x);
}

double df2dy(double x, double y)
{
	return pow(x, 2) * exp(1 - y) + 2 * exp(1 - x);
}

void System_Method_Newton(double a, double b)
{
	int count_itr = 1;
	double x0 = a, x1, y0 = b, y1, delta = 0.00001, dmax;
	x1 = x0 - (df2dy(x0, y0) * f1(x0, y0) - df2dx(x0, y0) * f2(x0, y0)) / (df1dx(x0, y0) * df2dy(x0, y0) - df1dy(x0, y0) * df2dx(x0, y0));
	y1 = y0 - (-df1dy(x0, y0) * f1(x0, y0) + df1dx(x0, y0) * f2(x0, y0)) / (df1dx(x0, y0) * df2dy(x0, y0) - df1dy(x0, y0) * df2dx(x0, y0));
	if (abs(x1 - x0) > abs(y1 - y0))
	{
		dmax = abs(x1 - x0);
	}
	else
	{
		dmax = abs(y1 - y0);
	}
	while (dmax > delta)
	{
		x0 = x1;
		y0 = y1;
		x1 = x0 - (df2dy(x0, y0) * f1(x0, y0) - df2dx(x0, y0) * f2(x0, y0)) / (df1dx(x0, y0) * df2dy(x0, y0) - df1dy(x0, y0) * df2dx(x0, y0));
		y1 = y0 - (-df1dy(x0, y0) * f1(x0, y0) + df1dx(x0, y0) * f2(x0, y0)) / (df1dx(x0, y0) * df2dy(x0, y0) - df1dy(x0, y0) * df2dx(x0, y0));
		if (abs(x1 - x0) > abs(y1 - y0))
		{
			dmax = abs(x1 - x0);
		}
		else
		{
			dmax = abs(y1 - y0);
		}
		count_itr++;
	}
	cout << "Newton Method: " << endl;
	cout << "[x, y]: " << fixed << setprecision(10) << x1 << " " << y1 << endl;
	cout << "Answer: " << fixed << setprecision(10) << f(x1, y1) << endl;
	cout << "Number of iteration: " << count_itr << endl;
	cout << "-----------------------------------------------------------------" << endl;
}

double fi(double x, double y)
{
	return (pow(y, 2) * exp(1 - x) - 1) / (2 * exp(1 - y));
}

double phi(double x, double y)
{
	return (pow(x, 2) * exp(1 - y)) / (2 * exp(1 - x));
}

void System_Simple_Iteration_Method(double x0, double y0)
{
	int itr_count = 1;
	double delta = 0.000001, dmax;
	double x1 = fi(x0, y0);
	double y1 = phi(x0, y0);
	if (abs(x1 - x0) > abs(y1 - y0))
	{
		dmax = abs(x1 - x0);
	}
	else
	{
		dmax = abs(y1 - y0);
	}
	while (dmax > delta)
	{
		itr_count++;
		x0 = x1;
		y0 = y1;
		x1 = fi(x0, y0);
		y1 = phi(x0, y0);
		if (abs(x1 - x0) > abs(y1 - y0))
		{
			dmax = abs(x1 - x0);
		}
		else
		{
			dmax = abs(y1 - y0);
		}
	}
	cout << "Method Simple Iterations: " << endl;
	cout << "[x, y]: " << fixed << setprecision(10) << x1 << " " << y1 << endl;
	cout << "Answer: " << fixed << setprecision(10) << f(x1, y1) << endl;
	cout << "Number of iteration: " << itr_count << endl;
	cout << "-----------------------------------------------------------------" << endl;
	cout << endl;
}





int main() {
	int c = -1;
	while (c != 0) {
		cout << "enter number for task, 0 for exit" << endl;
		cin >> c;
		switch (c) {
		case 1:
			task_1();
			break;
		case 2:
			//task_2();
			System_Simple_Iteration_Method(-0.3, 0.3);
			break;
		case 3:
			//task_3();
			System_Method_Newton(-0.3, 0.3);

			break;
		default:
			break;
		}
	}
}