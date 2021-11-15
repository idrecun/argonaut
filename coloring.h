#ifndef COLORING_H
#define COLORING_H

#include "permutation.h"
#include "hash.h"

namespace morphi {

template<typename T>
class Coloring {
public:

    // colors: [vertex, color, vertex, color, ...]
    Coloring(size_t size, uint32_t* colors) : m_permutation(size), m_cell_end(size, 0), m_cell_level(size, 0) {
        std::sort((uint64_t*) colors, (uint64_t*) colors + size); // BYTE ORDER?
        size_t cell_prev = 0;
        m_cell_end[0] = size;
        for(size_t idx = 0; idx < size; idx++) {
            m_permutation.m_forward[idx] = colors[2 * idx];
            m_permutation.m_inverse[colors[2 * idx]] = idx;
            if(idx > 0 && colors[2 * idx + 1] != colors[2 * idx - 1]) {
                m_cell_end[cell_prev] = idx;
                m_cell_end[idx] = size;
                cell_prev = idx;
            }
        }
    }

    size_t size() const {
        return m_permutation.size();
    }

    const T& operator[](size_t idx) const {
        return m_permutation[idx];
    }

    T indexOf(T vertex) const {
        return m_permutation.m_inverse[vertex];
    }

    size_t cellSize(size_t cell_idx) const {
        return m_cell_end[cell_idx] - cell_idx;
    }

    Permutation<T> m_permutation;
    Array<T> m_cell_end;
    Array<T> m_cell_level;
};

template<typename T>
std::ostream& operator<<(std::ostream& output, const Coloring<T>& c) {
    for(size_t cell = 0; cell < c.size(); cell = c.m_cell_end[cell]) {
        if(cell > 0)
            output << ']';
        output << /*(size_t)c.m_cell_level[cell] <<*/ "[ ";
        for(size_t idx = cell; idx < c.m_cell_end[cell]; idx++)
            output << (size_t) c.m_permutation[idx] << ' ';
    }
    output << ']';
    return output;
}

} // namespace

#endif // COLORING_H
