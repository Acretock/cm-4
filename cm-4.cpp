#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>

const double PI = 3.141592653589793238;


double function(double x) {
	return x * sin(x) * cos(2.0 * x);
}
double simpleIterFunction(double x) {
	return (exp(x) - x * x * x) / 2;
}
//Метод дихотомии
double BisectionMethod(double a, double b, double eps, int& iterations)
{
	double x;
	iterations = 0;
	while (b - a > eps) {
		x = (a + b) / 2;
		if (function(x) * function(a) > 0) {
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
	iterations = 0;
	while (fabs(b - a) > epsilon) {
		a = b - (b - a) * function(b) / (function(b) - function(a));
		b = a - (a - b) * function(a) / (function(a) - function(b));
		iterations++;
	}
	// a, b — (i - 1)-й и i-й члены

	return b;
}
//Метод Простых итераций
double SimpleIterMethod(double a, double epsilon, int& iterations) {
	double b;
	iterations = 0;
	for (;;) {
		b = simpleIterFunction(a);
		if (abs(a - b) < epsilon) break;
		a = b;
		iterations++;
	}
	return b;
}
//Метод Ньютона(касательных)
double derivativeFunction(double x) {
	return exp(x) - 3 * x * x - 2;
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
	return b;
}
//Метод секущих
double SecantMethod(double a, double c, double epsilon, int& iterations) {
	double b;
	iterations = 0;
	for (;;) {
		b = a - function(a) / (function(a) - function(c)) * (a - c);
		if (abs(a - b) < epsilon) break;
		c = a;
		a = b;
		iterations++;
	}
	return b;
}

void task_1() {
	double error = 1e-7;
	double exactAnswer = PI/4;
	std::vector<int> iterations(5);
	std::vector<double> results(5);
	std::cout << "Имя метода \tРезультат \t\tОшибка \t\t\t Количество итераций" << std::endl;
	results[0] = BisectionMethod(0, 2, error, iterations[0]);
	results[1] = ChordMethod(0, 2, error, iterations[1]);
	results[2] = SimpleIterMethod(1, error, iterations[2]);
	results[3] = NewtonMethod(1, error, iterations[3]);
	results[4] = SecantMethod(0, 2, error, iterations[4]);
	std::cout << std::fixed << std::setprecision(20) << "Метод бисекции \t" << results[0] << "\t" << abs(exactAnswer - results[0]) << "\t " << iterations[0] << std::endl;
	std::cout << std::fixed << std::setprecision(20) << "Метод хорд \t" << results[1] << "\t" << abs(exactAnswer - results[1]) << "\t " << iterations[1] << std::endl;
	std::cout << std::fixed << std::setprecision(20) << "Метод ПИ \t" << results[2] << "\t" << abs(exactAnswer - results[2]) << "\t " << iterations[2] << std::endl;
	std::cout << std::fixed << std::setprecision(20) << "Метод Ньютона \t" << results[3] << "\t" << abs(exactAnswer - results[3]) << "\t " << iterations[3] << std::endl;
	std::cout << std::fixed << std::setprecision(20) << "Метод секущих \t" << results[4] << "\t" << abs(exactAnswer - results[4]) << "\t " << iterations[4] << std::endl;

}