#include "matrix_operations.h"

#include <stdexcept>
#include <cmath>
#include <array>


void Matrix3x3::multiplyByVector(const std::array<double,3>& in, std::array<double,3>& out) const
{
    for (int i = 0; i < 3; ++i) {
        out[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
            out[i] += m[i][j] * in[j];
        }
    }
}


double Matrix3x3::determinant() const
{
    return
        m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) -
        m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0]) +
        m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);
}

Matrix3x3 Matrix3x3::inverse() const
{
    Matrix3x3 inv;
    double det = determinant();

    if (std::abs(det) < 1e-12)
        throw std::runtime_error("Matrix3x3::inverse(): singular matrix");

    inv.m[0][0] =  (m[1][1]*m[2][2] - m[1][2]*m[2][1]) / det;
    inv.m[0][1] = -(m[0][1]*m[2][2] - m[0][2]*m[2][1]) / det;
    inv.m[0][2] =  (m[0][1]*m[1][2] - m[0][2]*m[1][1]) / det;

    inv.m[1][0] = -(m[1][0]*m[2][2] - m[1][2]*m[2][0]) / det;
    inv.m[1][1] =  (m[0][0]*m[2][2] - m[0][2]*m[2][0]) / det;
    inv.m[1][2] = -(m[0][0]*m[1][2] - m[0][2]*m[1][0]) / det;

    inv.m[2][0] =  (m[1][0]*m[2][1] - m[1][1]*m[2][0]) / det;
    inv.m[2][1] = -(m[0][0]*m[2][1] - m[0][1]*m[2][0]) / det;
    inv.m[2][2] =  (m[0][0]*m[1][1] - m[0][1]*m[1][0]) / det;

    return inv;
}
