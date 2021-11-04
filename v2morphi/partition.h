#ifndef PARTITION_H
#define PARTITION_H

#include <algorithm>

#include "array.h"
#include "vector.h"

namespace morphi {

template<typename T>
class Partition {
public:

    struct PartitionClass {
        T size;
        T parent;
    };

    Partition(size_t size) : m_partition(size) {
        for(size_t idx = 0; idx < size; idx++)
            m_partition[idx] = {.size = 1, .parent = (T) idx};
    }

    T representative(T elem) {
        Vector<T> path(m_partition.m_size);
        while(elem != m_partition[elem].parent) {
            path.push(elem);
            elem = m_partition[elem].parent;
        }

        // Union-Find path compression
        for(size_t idx = 0; idx < path.m_size; idx++)
            m_partition[path[idx]].parent = elem;

        return elem;
    }

    void merge(T a, T b) {
        T pa = representative(a);
        T pb = representative(b);
        if(pa == pb)
            return;
        if(m_partition[pa].size < m_partition[pb].size) {
            m_partition[pa].parent = pb;
            m_partition[pb].size += m_partition[pa].size;
        }
        else {
            m_partition[pb].parent = pa;
            m_partition[pa].size += m_partition[pa].size;
        }
    }

    Array<PartitionClass> m_partition;
};

} // namespace

#endif // PARTITION_H
