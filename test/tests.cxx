#include <Eigen/Dense>
#include <EigenIntegration/Overrides.hpp>
#include "cim.hpp"
#include <complex>
#include <vector>
#include <iostream>
#include <string>
#include <tuple>

template <typename Scalar>
bool test_circular_contour_creation(Scalar center, typename Scalar::value_type radius, size_t resolution) {
    std::vector<Scalar> points = gen_circular_contour(center, radius, resolution);

    Scalar sum = Scalar(0, 0);
    typename Scalar::value_type angle_fraction = 2 * M_PI / (typename Scalar::value_type)resolution;;

    for (size_t idx = 0; idx < resolution; idx++) {
        Scalar point = points[idx];
        Scalar next_point = points[(idx+1) % resolution];
        sum += point;
        if (std::abs(std::abs(point - center) - std::abs(radius)) > 1e-4) {
            std::cerr << "Point: " << point << " is not on a circle around " << center << " with radius " << radius << std::endl;
            return false;
        }

        typename Scalar::value_type actual_angle = std::arg(next_point - center) - std::arg(point - center);
        typename Scalar::value_type angle = std::abs(actual_angle);
        angle = 2 * M_PI - angle > angle ? angle : 2 * M_PI - angle;
        if (std::abs(angle - angle_fraction) > 1e-4) {
            std::cerr << "Point: " << point << " and it's follow up " << next_point << " don't lie at an angle of " << angle_fraction << "; they are " << angle << " rad apart" << std::endl;
            return false;
        }

    }
    if (std::abs(sum/(typename Scalar::value_type)resolution - center) > 1e-4) {
        std::cerr << "The center of the contour " << sum/(typename Scalar::value_type)resolution << " does not coincide with the desired center " << center << std::endl;
        return false;
    }

    return true;
}

template <typename Scalar, typename Func>
bool test_integration(std::vector<Scalar> contour, Func f, Scalar expected) {
    std::pair<Eigen::Matrix<Scalar, -1, -1>, Eigen::Matrix<Scalar, -1, -1>> result = integrate<Eigen::Matrix<Scalar, -1, -1>, Func>(f, contour, false, true);

    if (std::abs(result.first(0, 0) - expected) > 1e-4) {
        std::cerr << "Integration yielded " << result.first(0, 0) << " but expected " << expected << std::endl;
        return false;
    }
    return true;
}

template <typename Func>
bool run_test(Func f, std::string name) {
    bool ok = f();

    if (!ok) {
        std::cerr << name << "FAIL" << std::endl;
    } else {
        std::cout << name << "OK" << std::endl;
    }

    return ok;
}

int main() {

    bool ok = true;
    ok &= run_test([]() -> bool {
        return test_circular_contour_creation(std::complex<double>(0, 0), 1, 50);
    }, "origin, unit radius, 50 points:                    ");

    ok &= run_test([]() -> bool {
        return test_circular_contour_creation(std::complex<double>(5, -3), 0.25, 29);
    }, "(5 - 3i), quarter radius, 29 points:               ");

    ok &= run_test([]() -> bool {
        return test_circular_contour_creation(std::complex<double>(5, -3), -1.25, 29);
    }, "(5 - 3i), -1.25 radius, 29 points:                 ");

    ok &= run_test([]() -> bool {
        std::vector<std::complex<double>> contour = gen_circular_contour(std::complex<double>(0, 0), 1.2, 50);
        Eigen::Matrix<std::complex<double>, -1, -1> mat = Eigen::Matrix<std::complex<double>, -1, -1>::Zero(1, 1);
        mat(0, 0) = std::complex<double>(1, 0);
        return test_integration(contour, [&](std::complex<double> z) {
            return z * z * mat;
        }, std::complex<double>(0, 0));
    }, "complex integration of holomorphic Scalar function: ");

    if (!ok) return -1;
    return 0;

}
