#ifndef GROUP_H
#define GROUP_H

#include "vector.h"
#include "bitarray.h"
#include "partition.h"
#include "permutation.h"
#include "assertions.h"

namespace morphi {

template<typename T>
class Group {
public:

    Group(size_t points, size_t max_elements)
        : m_elements(0),
          m_max_elements(max_elements),
          m_elem_cycles_ptr(max_elements + 1, 0),
          m_elem_cycles(2 * points * max_elements, 0),
          m_elem_fixed_points(points * max_elements),
          m_points(points),
          m_orbit_partition(points)
    {}

    void push(const Array<T>& permutation) {
        if(m_elements == m_max_elements)
            return;

        assert(m_points == permutation.m_size);
#ifdef DEBUG_SLOW_ASSERTS
        assertArrayElementsDistinct(permutation);
#endif
        size_t cycle_idx = m_elem_cycles_ptr[m_elements + 1] = m_elem_cycles_ptr[m_elements];
        BitArray visited(m_points);
        for(size_t point = 0; point < m_points; point++) {
            if(permutation[point] == point) {
                m_elem_fixed_points.set(fixedPointIndex(m_elements, point));
            }
            else if(!visited[point]) {
                while(!visited[point]) {
                    m_elem_cycles[cycle_idx++] = point;
                    m_orbit_partition.merge(point, permutation[point]);
                    visited.set(point);
                    point = permutation[point];
                }
                m_elem_cycles[cycle_idx++] = m_points;
            }
        }
        m_elem_cycles_ptr[++m_elements] = cycle_idx;
        for(size_t idx = 0; idx < m_elements; idx++)
            assert(m_elem_cycles_ptr[idx] <= m_elem_cycles_ptr[idx + 1]);
    }

    void updatePartitionCaller(const Vector<T>& stabilized, Partition<T>& partition, size_t& elem_idx, std::function<void(size_t, T, T)> onMerge) {
        for(; elem_idx < m_elements; elem_idx++) {
            bool stabilizes = true;
            for(size_t idx = 0; idx < stabilized.m_size; idx++)
                if(!m_elem_fixed_points[fixedPointIndex(elem_idx, stabilized[idx])]) {
                    stabilizes = false;
                    break;
                }
            if(stabilizes) {
                auto ptr = elemCyclesBegin(elem_idx);
                auto end = elemCyclesEnd(elem_idx);
                assert(ptr <= end);
                T cycle_rep = m_points;
                while(ptr != end) {
                    if(cycle_rep == m_points)
                        cycle_rep = *ptr;
                    else if(*ptr == m_points)
                        cycle_rep = m_points;
                    else {
                        onMerge(elem_idx, *(ptr - 1), *ptr);
                        partition.merge(*(ptr - 1), *ptr);
                    }
                    ptr++;
                }
            }
        }
    }

    void updatePartition(const Vector<T>& stabilized, Partition<T>& partition, size_t& elem_idx) {
        updatePartitionCaller(stabilized, partition, elem_idx, [](size_t, T, T) {});
    }

    size_t fixedPointIndex(size_t element, size_t point) {
        return element * m_points + point;
    }

    //   [begin(elem), end(elem)) = cycle1, separator, cycle2, separator, ..., cycleN, separator gde je sep = m_points (jer uvek staje u T, a ne koristi se)
    //   npr. (0 3)(1 4)(5 6) ce biti 0 3 7 1 4 7 5 6 7
    T* elemCyclesBegin(size_t element) {
        return m_elem_cycles.m_data + m_elem_cycles_ptr[element];
    }

    T* elemCyclesEnd(size_t element) {
        return elemCyclesBegin(element + 1);
    }

    // Group elements
    size_t m_elements;
    size_t m_max_elements;
    Array<size_t> m_elem_cycles_ptr;
    Array<T> m_elem_cycles;
    BitArray m_elem_fixed_points;

    // Group orbits
    size_t m_points;
    Partition<T> m_orbit_partition;
};

} // namespace

#endif // GROUP_H
