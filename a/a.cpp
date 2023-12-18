//#pragma once
#include <iostream>
#include <array>
#include <iomanip>

double alpha = 1.0;
double betta = 2.0;
double gamma = 2.0;
const int n = 21; /////////////
double h = 1.0 / double(n-1); ///////
double EPS = h * h * h;

double p(double x)
{
	return 1 + std::pow(x, gamma);
}
double p_(double x)
{
	return gamma * std::pow(x, gamma - 1);
}
double q(double x)
{
	return x + 1;
}
double u(double x)
{
	return std::pow(x, alpha) * std::pow(1 - x, betta);
}
double u_(double x)
{
	return alpha * betta * std::pow(x, alpha - 1) * std::pow(1 - x, betta - 1)
		- betta * std::pow(x, alpha) * std::pow(1 - x, betta - 1);
}
double u__(double x)
{
	return alpha * (alpha - 1) * std::pow(x, alpha - 2) * std::pow(1 - x, betta)
		- 2 * alpha * betta * std::pow(x, alpha - 1) * std::pow(1 - x, betta - 1)
		+ betta * (betta - 1) * std::pow(x, alpha) * std::pow(1 - x, betta - 2);
}
double f(double x)
{
	return -(p_(x) * u_(x) + p(x) * u__(x)) + q(x) * u(x);
}
double f(int i)
{
	return f(double(i) * h);
}
double a(int i)
{
	return p(double(i) * h);
}
double g(int i)
{
	return q(double(i) * h);
}

class Vector
{
	std::array<double, n> v;
public:
	Vector() { for (int i = 0; i < n; i++) v[i] = 0.0; }
	double& operator[](int i) { return v[i]; }
	const double& operator[](int i) const { return v[i]; }
	void printVert() const {
		for (int i = 0; i < n; i++)
			std::cout << v[i] << '\n';
	}
	friend std::ostream& operator<<(std::ostream& os, const Vector& v)
	{
		for (int i = 0; i < n; i++)
			os << v[i] << ' ';
		//os << '\n';
		return os;
	}
};
class Matrix
{
	std::array<Vector, n> mtx;
public:
	Vector& operator[](int i) { return mtx[i]; }
	const Vector& operator[](int i) const { return mtx[i]; }
	friend std::ostream& operator<<(std::ostream& os, const Matrix& mtx)
	{
		os << "Matrix<" << n << ">:\n";
		for (int i = 0; i < n - 1; i++)
			os << mtx[i] << '\n';
		os << mtx[n - 1];
		return os;
	}
};
double operator*(const Vector& v1, const Vector& v2)
{
	double res = 0.0;
	for (int i = 0; i < n; i++)
		res += v1[i] * v2[i];
	return res;
}
Vector operator*(double val, Vector v)
{
	for (int i = 0; i < n; i++)
		v[i] = v[i] * val;
	return v;
}
Vector operator*(const Matrix& m, const Vector& v)
{
	Vector res;
	for (int i = 0; i < n; i++)
	{
		res[i] = m[i] * v;
	}
	return res;
}
Vector operator-(Vector v1, const Vector& v2)
{
	for (int i = 0; i < n; i++)
		v1[i] = v1[i] - v2[i];
	return v1;
}
Vector abs(Vector v)
{
	for (int i = 0; i < n; i++)
		v[i] = std::abs(v[i]);
	return v;
}
bool operator>(const Vector& v, double val)
{
	for (int i = 1; i < n - 1; i++)
		if (v[i] < val) return false;
	return true;
}

class Tridiagonal
{
protected:
	Matrix mtx;
	Vector y;

	Vector alphas;
	Vector bettas;

public:

	void InitializeMatrix()
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				mtx[i][j] = 0.0;
		double ai = 0.0;
		double ci = a(0) + a(1) + h * h * g(0);
		double bi = a(1);
		if (n > 2)
		{
			alphas[0] = bi / ci;
			bettas[0] = 0.0;
			mtx[0][0] = ai;
			mtx[0][1] = ci;
			mtx[0][2] = bi;
		}

		for (int i = 1; i < n - 1; i++)
		{
			ai = alphas[i - 1] * a(i);
			ci = a(i) + a(i + 1) + h * h * g(i);
			bi = a(i + 1);

			mtx[i][i - 1] = ai;
			mtx[i][i] = ci;
			mtx[i][i + 1] = bi;

			alphas[i] = bi / (ci - ai);
			bettas[i] = (f(i) * h * h - bettas[i - 1] * a(i)) / (ci - ai);
		}

		ai = alphas[n - 2] * a(n - 1);
		ci = a(n - 1) + a(n) + h * h * g(n - 1);
		bi = 0.0;
		alphas[n - 1] = bi / (ci - ai);
		bettas[n - 1] = (f(n - 1) * h * h - bettas[n - 2] * a(n - 1)) / (ci - ai);
		mtx[n - 1][n - 3] = ai;
		mtx[n - 1][n - 2] = ci;
		mtx[n - 1][n - 1] = bi;
	}

	virtual const Matrix& getMatrix() const = 0;

	virtual const Vector& getAlphas() const = 0;
	virtual const Vector& getBettas() const = 0;

	virtual const Vector& getSolutions() = 0;
};
class Thomas : public Tridiagonal
{
public:
	Thomas()
	{
		InitializeMatrix();
	}

	const Matrix& getMatrix() const override
	{
		return mtx;
	};

	const Vector& getAlphas() const override
	{
		return alphas;
	}
	const Vector& getBettas() const override
	{
		return bettas;
	}

	const Vector& getSolutions() override
	{
		y[n - 1] = 0.0;
		y[0] = 0.0;
		for (int i = n - 2; i > 0; i--) {
			y[i] = alphas[i] * y[i + 1] + bettas[i];
		}
		return y;
	}
};
class Seidel : public Tridiagonal
{
public:
	Seidel()
	{
		InitializeMatrix();
	}

	const Matrix& getMatrix() const override
	{
		return mtx;
	};

	const Vector& getAlphas() const override
	{
		return alphas;
	}
	const Vector& getBettas() const override
	{
		return bettas;
	}

	const Vector& getSolutions() override
	{
		Vector p;
		int k = 0;
		do {
			for (int i = 0; i < n; i++)
				p[i] = y[i];
			for (int i = 1; i < n - 1; i++)
			{
				double var = 0.0;
				for (int j = 1; j < n - 1; j++)
					if (j != i) var += mtx[i][j] * y[j];
				y[i] = (f(i) * h * h) / mtx[i][i];
			}
			k++;
		} while (!converge(p));
		std::cout << "iterations: " << k << '\n';

		return y;
	}
	bool converge(const Vector& xkp)
	{
		double norm = 0.0;
		for (int i = 1; i < n - 1; i++)
			norm += (y[i] - xkp[i]) * (y[i] - xkp[i]);
		return (std::sqrt(norm) < EPS);
	}
	double okr(double x)
	{
		int i = 0;
		double neweps = EPS;
		while (neweps < 1)
		{
			i++;
			neweps *= n;
		}
		int okr = std::pow(double(n), i);
		x = int(x * okr + 0.5) / double(okr);

		return x;
	}
};
class Jacobi : public Tridiagonal
{
public:
	Jacobi()
	{
		InitializeMatrix();
	}

	const Matrix& getMatrix() const override
	{
		return mtx;
	};

	const Vector& getAlphas() const override
	{
		return alphas;
	}
	const Vector& getBettas() const override
	{
		return bettas;
	}

	const Vector& getSolutions() override
	{
		Vector temp;
		double norm = 0.0;
		int k = 0;
		do {
			for (int i = 1; i < n - 1; i++) {
				temp[i] = f(i) * h * h;
				for (int j = 1; j < n - 1; j++) {
					if (j != i) temp[i] -= mtx[i][j] * y[j];
				}
				temp[i] /= mtx[i][i];
			}
			norm = fabs(y[0] - temp[0]);
			for (int k = 1; k < n - 1; k++) {
				if (fabs(y[k] - temp[k]) > norm)
					norm = fabs(y[k] - temp[k]);
				y[k] = temp[k];
			}
			k++;
		} while (norm > EPS);
		std::cout << "iterations: " << k << '\n';

		return y;
	}
};
class DownRelaxation : public Tridiagonal
{
public:
	DownRelaxation()
	{
		InitializeMatrix();
	}
	const Matrix& getMatrix() const override
	{
		return mtx;
	};

	const Vector& getAlphas() const override
	{
		return alphas;
	}
	const Vector& getBettas() const override
	{
		return bettas;
	}

	const Vector& getSolutions() override
	{
		double w = 0.0;
		double norm = 0.0;
		Vector temp;
		int k = 0;
		do
		{
			for (int i = 1; i < n - 1; i++)
			{
				y[i] = f(i) * h * h;
				for (int j = 1; j < n - 1; j++)
				{
					if (i != j)
						y[i] = y[i] - mtx[i][j] * y[j];
				}
				y[i] /= mtx[i][i];
				y[i] = w * y[i] + (1 - w) * temp[i];
				w += 1.0 / n;
				if (fabs(y[i] - temp[i]) > norm)
					norm = fabs(y[i] - temp[i]);
				temp[i] = y[i];
			}
			k++;
		} while (norm < EPS);
		std::cout << "iterations: " << k << '\n';

		return y;
	}
};
class UpRelaxation :public Tridiagonal
{
public:
	UpRelaxation()
	{
		InitializeMatrix();
	}
	const Matrix& getMatrix() const override
	{
		return mtx;
	};

	const Vector& getAlphas() const override
	{
		return alphas;
	}
	const Vector& getBettas() const override
	{
		return bettas;
	}

	const Vector& getSolutions() override
	{
		double norm = 0.0;
		double w = 1.0;
		Vector temp;
		int k = 0;
		do
		{
			for (int i = 1; i < n - 1; i++)
			{
				y[i] = f(i) * h * h;
				for (int j = 1; j < n - 1; j++)
				{
					if (i != j)
						y[i] = y[i] - mtx[i][j] * y[j];
				}
				y[i] /= mtx[i][i];
				y[i] = w * y[i] + (1 - w) * temp[i];
				w += 1.0 / n;
				if (fabs(y[i] - temp[i]) > norm)
					norm = fabs(y[i] - temp[i]);
				temp[i] = y[i];
			}
			k++;
		} while (norm < EPS);
		std::cout << "iterations: " << k << " w:" << w <<'\n';

		return y;
	}
};
class GradientDescent : public Tridiagonal
{
public:
	GradientDescent()
	{
		InitializeMatrix();
	}
	const Matrix& getMatrix() const override
	{
		return mtx;
	};

	const Vector& getAlphas() const override
	{
		return alphas;
	}
	const Vector& getBettas() const override
	{
		return bettas;
	}

	const Vector& getSolutions() override
	{
		Vector r;
		double t = 0.0;
		int k = 0;
		Vector temp = y;
		do
		{
			for (int i = 1; i < n - 1; i++)
			{
				double ai = mtx[i][i - 1];
				double ci = mtx[i][i];
				double bi = mtx[i][i + 1];

				r[i] = -ai * y[i - 1] + ci * y[i] - bi * y[i + 1] - f(i) * h * h;
			}
			t = r * r / ((mtx * r) * r);
			y = temp - t * r;
			k++;
			temp = y;
		} while (abs(r) > EPS);
		std::cout << "iterations: " << k << '\n';
		return y;
	}
};

int main()
{
	Thomas method1;
	Matrix mtx = method1.getMatrix();
	std::cout << mtx << '\n';
	std::cout << "******Thomas test start**********\n";
	Vector y1 = method1.getSolutions();

	std::array<double, n> uy;

	for (int i = 0; i < n; i++) {
		uy[i] = u(y1[i]);
	}


	std::cout << "ih \t\t yi \t\t uy(i) \t\t |yi-uy[i]|\n";
	for (int i = 0; i < n-1; i++) ////////////////////
		std::cout << std::fixed << std::setprecision(15) << i * h << "\t\t" << y1[i] << "\t\t" << uy[i] << "\t\t" << fabs(y1[i] - uy[i]) << '\n';

	std::cout << "******Thomas test end************\n\n";


	std::cout << "******Jacobi test start**********\n";
	Jacobi method6;
	Vector y6 = method6.getSolutions();
	std::cout << "ih \t\t yi \t\t uy[i] \t\t |yi-uy[i]|\n";
	for (int i = 0; i < n-1; i++) ////////////////////
		std::cout << i * h << "\t\t" << y6[i] << "\t\t" << uy[i] << "\t\t" << fabs(y6[i] - uy[i]) << '\n';

	std::cout << "******Jacobi test end************\n\n";




	std::cout << "******Seidel test start**********\n";
	Seidel method2;
	Vector y2 = method2.getSolutions();
	std::cout << "ih \t\t yi \t\t uy[i] \t\t |yi-uy[i]|\n";
	for (int i = 0; i < n-1; i++) ////////////////////
		std::cout << i * h << "\t\t" << y2[i] << "\t\t" << uy[i] << "\t\t" << fabs(y2[i] - uy[i]) << '\n';

	std::cout << "******Seidel test end************\n\n";


	std::cout << "******UpRelaxation test start******\n";
	UpRelaxation method4;
	Vector y4 = method4.getSolutions();
	std::cout << "ih \t\t yi \t\t uy[i] \t\t |yi-uy[i]|\n";
	for (int i = 0; i < n-1; i++) ////////////////////
		std::cout << i * h << "\t\t" << y4[i] << "\t\t" << uy[i] << "\t\t" << fabs(y4[i] - uy[i]) << '\n';

	std::cout << "******UpRelaxation test end********\n\n";


	return 0;
}