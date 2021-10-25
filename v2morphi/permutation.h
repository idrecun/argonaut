#ifndef PERMUTATION_H
#define PERMUTATION_H

#include "array.h"

namespace morphi {

template<typename T>
class Permutation {
public:

    Permutation() {}
    Permutation(size_t size) : m_forward(size), m_inverse(size) { }

    size_t size() const {
        return m_forward.m_size;
    }

    void swap(T a, T b) {
        if(a != b) {
            std::swap(m_forward[a], m_forward[b]);
            std::swap(m_inverse[m_forward[a]], m_inverse[m_forward[b]]);
        }
    }

    void set(T idx, T val) {
        m_forward[idx] = val;
        m_inverse[val] = idx;
    }

    const T& operator[](T idx) const {
        return m_forward[idx];
    }

    template<typename V>
    void copyFwd(const Permutation<V>& oth) {
        std::copy(oth.m_forward.m_data, oth.m_forward.m_end, m_forward.m_data);
        std::copy(oth.m_inverse.m_data, oth.m_inverse.m_end, m_inverse.m_data);
    }

    template<typename V>
    void copyInv(const Permutation<V>& oth) {
        std::copy(oth.m_forward.m_data, oth.m_forward.m_end, m_inverse.m_data);
        std::copy(oth.m_inverse.m_data, oth.m_inverse.m_end, m_forward.m_data);
    }

    Permutation<T> inverse() const {
        Permutation<T> ret(m_forward.m_size);
        ret.copyInv(*this);
        return ret;
    }

    Array<T> m_forward;
    Array<T> m_inverse;
};

} // namespace

#endif // PERMUTATION_H
