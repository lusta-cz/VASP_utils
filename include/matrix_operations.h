#ifndef MATRIX_OPERATIONS_H_INCLUDED
#define MATRIX_OPERATIONS_H_INCLUDED

#include <array>

struct Matrix3x3 {
    double m[3][3];

    void multiplyByVector(const std::array<double,3>& in, std::array<double,3>& out) const;
    double determinant() const;
    Matrix3x3 inverse() const;
};

#endif // MATRIX_OPERATIONS_H_INCLUDED
