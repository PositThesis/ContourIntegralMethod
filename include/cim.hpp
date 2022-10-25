#include <EigenIntegration/std_integration.hpp>
#include <Eigen/Dense>
#include <complex>
#include <functional>
#include <omp.h>
#include <random>
#include <tuple>
#include <vector>

template <typename Scalar> class ComplexMatrixFunction {
public:
  virtual Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> get(Scalar z);
  virtual int size();
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> operator()(Scalar z) {
    return this->get(z);
  }
};

template <typename Scalar>
std::vector<Scalar> gen_circular_contour(Scalar center, typename Scalar::value_type radius, size_t resolution) {
  using Real = typename Scalar::value_type;
  std::vector<Scalar> points;

  for (size_t idx = 0; idx < resolution; idx++) {
    Real x = ((Real)idx) / ((Real)resolution) * 2 * M_PI;

    points.push_back(
        center +
        radius *
            std::complex<typename Scalar::value_type>(
                std::cos(x),
                std::sin(x)));
  }

  return points;
}

template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
gen_random_full_rank(size_t &rows, size_t &cols, size_t seed) {
  using RealScalar = typename Scalar::value_type;
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix(rows, cols);

  std::default_random_engine e(seed);
  std::normal_distribution<double> dist(0, 1);

  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      double random_num_r = dist(e);
      double random_num_i = dist(e);
      matrix(i, j) =
          std::complex(RealScalar(random_num_r), RealScalar(random_num_i));
    }
  }

  return matrix;
}

template <typename MatrixType, typename Fn>
std::pair<MatrixType, MatrixType>
integrate(Fn &fn, std::vector<typename MatrixType::Scalar> points, bool invert,
          bool closed) {
  using Scalar = typename MatrixType::Scalar;
  assert(fn(points[0]).rows() == fn(points[0]).cols());
  assert(points.size() > 0);

  Eigen::Index size = fn(points[0]).rows();
  MatrixType acc0 = MatrixType::Zero(size, size);
  MatrixType acc1 = MatrixType::Zero(size, size);

  Eigen::Index max_index = closed ? points.size() : points.size() - 1;

#pragma omp parallel for
  for (Eigen::Index idx = 0; idx < max_index; idx++) {
    size_t upper_idx = (idx + 1) % points.size();

    Scalar average = (points[idx] + points[upper_idx]) /
                     Scalar(2.0); // if scalar is real, this works. if it is
                                  // complex, it still works
    // approximate the differential by difference
    Scalar delta = points[upper_idx] - points[idx];

    MatrixType fn_at_point = fn(average);
    if (invert)
      fn_at_point = fn_at_point.inverse();

#pragma omp critical
    acc0 += fn_at_point * delta;

#pragma omp critical
    acc1 += average * fn_at_point * delta;
  }

  return std::pair(acc0, acc1);
}

template <typename MatrixType>
std::vector<typename MatrixType::Scalar>
get_singular_values(Eigen::JacobiSVD<MatrixType> &svd,
                    typename MatrixType::Scalar::value_type gamma) {
  MatrixType singularValues = svd.singularValues();
  std::vector<typename MatrixType::Scalar> values_above_threshold;
  for (Eigen::Index idx = 0; idx < singularValues.rows(); idx++) {
    typename MatrixType::Scalar value = singularValues(idx, 0);

    if (value.real() * value.real() + value.imag() * value.imag() >
        gamma * gamma) {
      values_above_threshold.push_back(value);
    } else {
      // Felix noted that the algorithm depends on the singular values being
      // ordered. Therefore, as soon as any value is below the threshold, all
      // following values will be below it aswell.
      break;
    }
  }
  return values_above_threshold;
}

template <typename MatrixType>
MatrixType
build_inverse_diagonal(std::vector<typename MatrixType::Scalar> &values) {
  using RealScalar = typename MatrixType::Scalar::value_type;
  size_t k = values.size();
  MatrixType mat = MatrixType::Zero(k, k);
  for (size_t idx = 0; idx < k; idx++) {
    mat(idx, idx) = std::complex<RealScalar>(1.0, 0.0) / values[idx];
  }
  return mat;
}

template <typename MatrixType, typename T>
requires HasComplexScalar<MatrixType>
    std::pair<Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, 1>,
              MatrixType>
    cim_solve(T &X, std::vector<typename MatrixType::Scalar> points,
              typename MatrixType::Scalar::value_type gamma, size_t l) {
  using ComplexScalar = typename MatrixType::Scalar;
  using RealScalar = typename MatrixType::Scalar::value_type;

  assert(l > 0);
  assert(points.size() > 0);

  bool kfound = false;
  size_t m = X(points[0]).rows();
  size_t N = points.size();

  // The sum parts of A0 and A1 can be precomputed
  std::pair<MatrixType, MatrixType> accumulators =
      integrate<MatrixType>(X, points, true, true);
  MatrixType X_inv_sum = accumulators.first;
  MatrixType point_X_inv_sum = accumulators.second;

  while (!kfound) {
    std::cerr << "starting run with l = " << l << std::endl;
    MatrixType V_hat = gen_random_full_rank<ComplexScalar>(m, l, m + l);

    // A0 = 1 / iN * sum(X_inv) * V_hat
    MatrixType A0 =
        X_inv_sum * V_hat / std::complex<RealScalar>(0, (RealScalar)N);

    // Compute U, Sigma and V
    Eigen::JacobiSVD<MatrixType> svd(A0,
                                  Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Figure out the amount of nonzero singular values
    std::vector<ComplexScalar> values =
        get_singular_values<MatrixType>(svd, gamma);

    // If there are less than l nonzero singular values, i.e. if there are some
    // singular values that are 0
    size_t k = values.size();
    if (k < l || l == m) {
      if (k == 0) {
        return std::pair(MatrixType::Zero(0, 0), MatrixType::Zero(0, 0));
      }
      kfound = true;

      // Get The V0h and W0 matrices as subsets of the U and V matrices from the
      // SVD. caution: the V Matrix returned from the SVD is used for W0.
      MatrixType V0 = svd.matrixU().block(0, 0, m, k);
      MatrixType V0h = V0.conjugate().transpose();
      MatrixType W0 = svd.matrixV().block(0, 0, l, k);

      // build inverse sigma from the nonzero values
      MatrixType sigma_inv = build_inverse_diagonal<MatrixType>(values);

      // A1 = 1/iN * sum(point * X_inv) * V_hat
      MatrixType A1 = point_X_inv_sum * V_hat / std::complex<RealScalar>(0, N);

      // V0h:       kxm
      // A1:        mxl
      // W0:        lxk
      // sigma_inv: kxk

      MatrixType B = V0h * A1 * W0 * sigma_inv;
      Eigen::ComplexEigenSolver<MatrixType> ces(B);
      Eigen::Matrix<ComplexScalar, Eigen::Dynamic, 1> eigs = ces.eigenvalues();
      MatrixType vectors = V0 * ces.eigenvectors();

      return std::pair(eigs, vectors);
    } else {
      l += 1;
    }
  }

  throw std::exception();
}