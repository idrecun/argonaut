#ifndef MATRIX_H
#define MATRIX_H

#include "array.h"

namespace morphi {

template<typename T>
class Matrix {
public:
    Matrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols), m_data(rows * cols) { }

    const T& at(size_t r, size_t c) const {
        return m_data[r * m_cols + c];
    }

    const T& set(size_t r, size_t c, const T& val) {
        return m_data[r * m_cols + c] = val;
    }


    Array<T> m_data;
    size_t m_rows = 0;
    size_t m_cols = 0;
};

} // namespace

#endif // MATRIX_H
