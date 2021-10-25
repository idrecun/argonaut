#ifndef GRAPH_H
#define GRAPH_H

#include <cstddef>

#include "array.h"
#include "bitarray.h"
#include "permutation.h"

namespace morphi {

template<typename T>
class Graph {
public:
    Graph(size_t vertices, size_t edges, uint32_t* edge_list) : m_vertices(vertices), m_adj_ptr(vertices + 1, 0), m_adj_list(2 * edges), m_adj_matrix(vertices * (vertices - 1) / 2) {
        uint32_t* end = edge_list + 2 * edges;
        for(auto ptr = edge_list; ptr != end; ptr += 2) {
            uint32_t a = *ptr, b = *(ptr + 1);
            m_adj_matrix.set(matrixAt(a, b));
            m_adj_ptr[a]++;
            m_adj_ptr[b]++;
        }

        for(auto ptr = m_adj_ptr.m_data + 1; ptr != m_adj_ptr.m_end; ptr++)
            *ptr += *(ptr - 1);

        for(auto ptr = edge_list; ptr != end; ptr += 2) {
            uint32_t a = (uint32_t)*ptr, b = (uint32_t)*(ptr + 1);
            m_adj_list[--m_adj_ptr[a]] = (T) b;
            m_adj_list[--m_adj_ptr[b]] = (T) a;
        }
    }

    bool adjacent(size_t a, size_t b) const {
        if(a == b)
            return false;
        return m_adj_matrix[matrixAt(a, b)];
    }

    T* begin(T vertex) {
        return m_adj_list.m_data + m_adj_ptr[vertex];
    }

    T* end(T vertex) {
        return begin(vertex + 1);
    }

    // Compares G^inv(a) and G^inv(b)
    // a, b are expected to be inverse permutations
    bool less(const Permutation<T>& a, const Permutation<T>& b) const {
        uint8_t aij, bij;
        for(size_t i = 0; i < m_vertices; i++)
            for(size_t j = 0; j < m_vertices; j++) {
                aij = adjacent(a.m_inverse[i], a.m_inverse[j]);
                bij = adjacent(b.m_inverse[i], b.m_inverse[j]);
                if(aij < bij) return true;
                if(aij > bij) return false;
            }
        return false;
    }

    /*bool aut(const Permutation<T>& p) {

    }*/

    inline size_t matrixAt(size_t a, size_t b) const {
        assert(a != b);
        if(a > b)
            std::swap(a, b);
        return a * (2 * m_vertices - a - 1) / 2 + b - a - 1;
    }

    size_t m_vertices;
    Array<size_t> m_adj_ptr;
    Array<T> m_adj_list;
    BitArray m_adj_matrix;
};

} // namespace

#endif // GRAPH_H
