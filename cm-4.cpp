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
	//return sin(x) * cos(2.0 * x)*2;
	return x + function(x);
}

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
int main() {
	int c = 1;
	//while (c != 0) {
		//cout << "enter number for task, 0 for exit" << endl;
		//cin >> c;
		switch (c) {
		case 1:
			task_1();
			break;
		default:
			break;
		}
	//}
		cout << endl;
		cin >> c;
}