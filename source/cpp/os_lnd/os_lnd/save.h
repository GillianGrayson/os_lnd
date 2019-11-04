#pragma once
#include <vector>
#include <iomanip>
#include <fstream>

template <class T>
void save_vector(const std::vector<T>& v, std::string fn, int precision = 16)
{
	std::ofstream f_out(fn);
	f_out << std::setprecision(precision) << std::scientific;

	for (auto const& x : v)
	{
		f_out << x << std::endl;
	}

	f_out.close();
}

template <class T>
void save_sp_mtx(Eigen::SparseMatrix<T> m, std::string fn, int precision = 16)
{
	std::ofstream f_out(fn);
	f_out << std::setprecision(precision) << std::scientific;

	for (int k = 0; k < m.outerSize(); ++k)
	{
		for (typename Eigen::SparseMatrix<T>::InnerIterator it(m, k); it; ++it)
		{
			f_out << it.row() << "\t"; // row index
			f_out << it.col() << "\t"; // col index (here it is equal to k)
			f_out << it.value() << std::endl;
		}
	}

	f_out.close();
}
