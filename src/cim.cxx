#if USE_FDP == 1
#pragma message ("FDP is ON")
#include <EigenIntegration/Overrides.hpp>
#else
// make sure that it is defined at all
#pragma message ("FDP is OFF")
static_assert(USE_FDP == 0);
#endif

#include "cim.hpp"
#include <EigenIntegration/MtxIO.hpp>
#include <Eigen/Dense>
#include <chrono>
#include <complex>
#include <iostream>
#include <string.h>
#include <string>
#include <universal/number/posit/posit.hpp>
#include <chrono>

template <typename Scalar> using Vec = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
template <typename Scalar>
using Mat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
template <size_t nbits, size_t es>
using complex_posit = std::complex<sw::universal::posit<nbits, es>>;

template <typename Scalar>
Mat<Scalar> eval_poly(std::vector<Mat<Scalar>> matrices, Scalar z) {
  return matrices[0] * z * z + matrices[1] * z + matrices[2];
}

template <typename Scalar>
class ComplexPoly : public ComplexMatrixFunction<Scalar> {
private:
  int size_ = -1;
  std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> coeffs;

public:
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> get(Scalar z) {
    assert(this->size_ > 0);
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> result =
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(
            this->size_, this->size_);

    for (size_t idx = 0; idx < this->coeffs.size(); idx++) {
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> inner_result(
          this->size_, this->size_);

      inner_result = this->coeffs[idx];
      for (size_t idx_inner = 0; idx_inner < idx; idx_inner++) {
        inner_result *= z;
      }

      result += inner_result;
    }
    return result;
  }

  void push(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> M) {
    if (this->size_ == -1) {
      this->size_ = M.rows();
    }

    assert(M.rows() == M.cols());
    assert(M.rows() == this->size_);
    this->coeffs.push_back(M);
  }

  int size() { return this->size_; }
};
/**
 * Arguments:
 *  - Int: Number of points
 *  - Floating: Radius
 *  - Floating: Gamma
 *  - Int: L-start
 *  - [String]: Input Matrix name in polynomial order: [0] + [1]*z + [2]*z²
 */
int main(int argc, char **argv) {
#ifdef USE_Float
  using Scalar = std::complex<float>;
#endif
#ifdef USE_Double
  using Scalar = std::complex<double>;
#endif
#ifdef USE_LongDouble
  using Scalar = std::complex<long double>;
#endif
#ifdef USE_Posit16
  using Scalar = std::complex<sw::universal::posit<16, 2>>;
#endif
#ifdef USE_Posit32
  using Scalar = std::complex<sw::universal::posit<32, 2>>;
#endif
#ifdef USE_Posit64
  using Scalar = std::complex<sw::universal::posit<64, 2>>;
#endif
#ifdef USE_Posit643
  using Scalar = std::complex<sw::universal::posit<64, 3>>;
#endif
#ifdef USE_Posit644
  using Scalar = std::complex<sw::universal::posit<64, 4>>;
#endif
#ifdef USE_Posit645
  using Scalar = std::complex<sw::universal::posit<64, 5>>;
#endif
#ifdef USE_Posit323
  using Scalar = std::complex<sw::universal::posit<32, 3>>;
#endif
#ifdef USE_Posit324
  using Scalar = std::complex<sw::universal::posit<32, 4>>;
#endif
#ifdef USE_Posit325
  using Scalar = std::complex<sw::universal::posit<32, 5>>;
#endif
#ifdef USE_Posit163
  using Scalar = std::complex<sw::universal::posit<16, 3>>;
#endif
#ifdef USE_Posit164
  using Scalar = std::complex<sw::universal::posit<16, 4>>;
#endif
#ifdef USE_Posit165
  using Scalar = std::complex<sw::universal::posit<16, 5>>;
#endif
  using Real = Scalar::value_type;

  enum class Analysis {
    Unknown,
    Poly,

  };

  // configuration options
  Analysis analysis = Analysis::Unknown;
  std::vector<std::string> input_files;
  std::string output_file("-");
  Scalar center = Scalar(0);
  Real gamma = 1e-4;
  Real radius = 1.0;
  size_t l_start = 1;
  size_t num_points = 0;

  for (int i = 1; i < argc; i++) {
    if (i == argc - 1) {
      std::cerr << "Trailing argument" << std::endl;
      return -1;
    }
    if (std::strcmp(argv[i], "-cr") == 0) {
      i += 1;
      center = Scalar((Real)std::stod(argv[i]), center.imag());
      continue;
    }
    if (std::strcmp(argv[i], "-ci") == 0) {
      i += 1;
      center = Scalar(center.real(), (Real)std::stod(argv[i]));
      continue;
    }
    if (std::strcmp(argv[i], "-i") == 0) {
      i += 1;
      input_files.push_back(std::string(argv[i]));
      continue;
    }
    if (std::strcmp(argv[i], "-g") == 0) {
      i += 1;
      gamma = (Real)(std::stod(argv[i]));
      continue;
    }
    if (std::strcmp(argv[i], "-l") == 0) {
      i += 1;
      l_start = std::stoi(argv[i]);
      continue;
    }
    if (std::strcmp(argv[i], "-p") == 0) {
      i += 1;
      num_points = std::stoi(argv[i]);
      continue;
    }
    if (std::strcmp(argv[i], "-o") == 0) {
      i += 1;
      output_file = std::string(argv[i]);
      continue;
    }
    if (std::strcmp(argv[i], "-r") == 0) {
      i += 1;
      radius = (Real)std::stod(argv[i]);
      continue;
    }
    if (std::strcmp(argv[i], "-t") == 0) {
      i += 1;
      if (std::strcmp(argv[i], "poly") == 0)
        analysis = Analysis::Poly;
      continue;
    }
    std::cerr << "unrecognized option: " << argv[i] << std::endl;
    return -1;
  }

  assert(analysis != Analysis::Unknown);
  // check parameter validity
  assert(l_start > 0);
  if (analysis == Analysis::Poly)
    assert(input_files.size() > 0);
  assert(gamma > 0);
  assert(num_points > 0);
  assert(radius > 0);

  std::vector<Scalar> points = gen_circular_contour(center, radius, num_points);

  for (long unsigned int i = 0; i < points.size(); i++) {
    assert(std::abs(std::abs(points[i]) - radius) < 0.01);
  }

  ComplexPoly<Scalar> poly;
  for (std::string &f : input_files) {
    poly.push(get_matrix_from_mtx_file<Scalar>(f));
  }

  std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
  std::pair<Mat<Scalar>, Mat<Scalar>> result =
      cim_solve<Mat<Scalar>>(poly, points, gamma, l_start);
  std::chrono::time_point<std::chrono::high_resolution_clock> later = std::chrono::high_resolution_clock::now();

  if (output_file != "-")
    write_matrix(result.first, output_file + ".mtx");
  else
    std::cout << "Result: " << result.first << std::endl;

  Mat<Real> residuals = Mat<Real>::Zero(result.first.rows(), 1);
  for (Eigen::Index idx = 0; idx < result.first.rows(); idx++) {
    Scalar eigenvalue = result.first(idx, 0);
    Mat<Scalar> func_at_eigen = poly(eigenvalue);
    Mat<Scalar> residual_vec = func_at_eigen * result.second.col(idx);
    residuals(idx) = residual_vec.norm();
  }

  if (output_file != "-")
    write_matrix(residuals, output_file + "_residuals.mtx");
  else
    std::cout << "Residuals: " << residuals << std::endl;

  if (output_file != "-") {
    std::ofstream file(output_file + "_residuals.json");
    file << "{\n\t\"max\": " << residuals.maxCoeff()
         << ",\n\t\"min\": " << residuals.minCoeff()
         << ",\n\t\"microseconds\": " << std::chrono::duration_cast<std::chrono::microseconds>(later - start).count()
         << "\n}" << std::endl;
    file.close();
  } else {
    std::cout << "Max Residual: " << residuals.maxCoeff() << std::endl;
    std::cout << "Min Residual: " << residuals.minCoeff() << std::endl;
  }
}