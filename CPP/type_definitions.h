
#ifndef TYPE_DEFINITIONS_H
#define TYPE_DEFINITIONS_H

#include <complex>

typedef unsigned int uint;

typedef Eigen::Matrix<unsigned long int, Eigen::Dynamic, Eigen::Dynamic> MatrixXuli;
typedef Eigen::Matrix<short int, Eigen::Dynamic, Eigen::Dynamic>         MatrixXsi;
typedef double variable_type;
typedef Eigen::Matrix<variable_type, Eigen::Dynamic, Eigen::Dynamic> MatrixX;

typedef double floatType;
typedef std::complex<floatType> cfloatType;
typedef std::complex<double> cdouble;

#endif // TYPE_DEFINITIONS_H

