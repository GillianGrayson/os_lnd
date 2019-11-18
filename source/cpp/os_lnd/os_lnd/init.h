#pragma once
#include "third_party/eigen-git-mirror/Eigen/SparseCore"
#include "third_party/eigen-git-mirror/Eigen/Core"

typedef Eigen::SparseMatrix<std::complex<double>> sp_mtx;
typedef Eigen::Triplet<std::complex<double>> triplet;