#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <cstdint>
#include <limits>
#include "vector.h"
#include "bitarray.h"
#include "coloring.h"
#include "graph.h"
#include "hash.h"

namespace morphi {

template<typename T, typename HashType>
class AlgorithmDFS {
public:

    enum PathType {
        MaxPath = 1,
        AutPath = 2,
    };

    struct NodePath {
        NodePath(size_t size) : permutation(size), invariants(size) {}

        bool is_leaf = false;
        T lca_level = 0;
        Permutation<T> permutation;
        Vector<HashType> invariants;
    };

    AlgorithmDFS(size_t vertices, size_t edges, Array<uint32_t>& edge_list, Array<uint32_t>& colors)
        : graph(vertices, edges, edge_list.m_data),
          coloring(vertices, colors.m_data),
          invariants(vertices),
          stabilized(vertices),
          max_node(vertices),
          fst_node(vertices)
    {}

    const Permutation<T>& solve() {
        solve(PathType::MaxPath | PathType::AutPath);
        return max_node.permutation;
    }

    T solve(uint8_t flags) {
        assert(invariants.m_size <= stabilized.m_size);
        T level = stabilized.m_size;
#ifdef DEBUG_OUT
        std::cout << std::string(level, ' ') << "I " << coloring << std::endl;
#endif
        refine();
#ifdef DEBUG_OUT
        std::cout << std::string(level, ' ') << "R " << coloring << std::endl;
#endif
        /*
        if(!refine()) // zameniti sa refine() == BAD_PATH
            return level;
        */

        // update max_path, aut_path nodes
        bool max_path = (flags & PathType::MaxPath) &&
                       (max_node.invariants.m_size <= level ||
                        invariants.back() >= max_node.invariants[level]);
        if(!max_path) {
            unrefine();
            return level;
        }
        else {
            max_node.lca_level = level;

            if(max_node.invariants.m_size > level && invariants.back() > max_node.invariants[level])
                max_node.invariants.pop(level);

            if(max_node.invariants.m_size <= level) {
                max_node.invariants.push(invariants.back());
                max_node.is_leaf = false;
            }
        }

        size_t cell_idx;
        Array<T> cell_content = targetCell(cell_idx);
        for(auto ptr = cell_content.m_data; ptr != cell_content.m_end; ptr++) {
            individualize(cell_idx, *ptr);
            T backjump = solve(flags);
            unindividualize(cell_idx);

            if(level > backjump) {
                unrefine();
                return level;
            }
        }

        if(cell_content.m_size == 0) {
            if(max_path && (!max_node.is_leaf || (max_node.invariants.m_size == (size_t) level + 1 && graph.less(max_node.permutation, coloring.m_permutation.inverse())))) {
                max_node.is_leaf = true;
                max_node.permutation.copyInv(coloring.m_permutation);
            }
        }

        unrefine();
        return level;
    }

    void individualize(size_t cell_idx, T vertex) {
#ifdef DEBUG_OUT
        std::cout << "individualize" << std::endl;
#endif
        assert(coloring.m_cell_end[cell_idx] - cell_idx > 1);

        stabilized.push(vertex);

        T vertex_idx = coloring.m_permutation.m_inverse[vertex];
        coloring.m_permutation.swap(cell_idx, vertex_idx);

        coloring.m_cell_end[cell_idx + 1] = coloring.m_cell_end[cell_idx];
        coloring.m_cell_end[cell_idx] = cell_idx + 1;

        coloring.m_cell_level[cell_idx + 1] = stabilized.m_size - 1;
    }

    void unindividualize(size_t cell_idx) {
        stabilized.pop();

        coloring.m_cell_end[cell_idx] = coloring.m_cell_end[cell_idx + 1];
    }

    template<bool RootLevel>
    void refine1(size_t work_cell, Vector<T>& active_cells, BitArray& is_active, HashType& invariant) {
#ifdef DEBUG_OUT
        std::cout << std::string(stabilized.m_size, ' ') << "refine1" << std::endl;
#endif
        if constexpr (!RootLevel)
            hash::sequential32u(invariant, work_cell);

        size_t cell_beg = 0, cell_end = coloring.m_cell_end[0];
        size_t split = 0; // [cell_beg, split), [split, idx)
        for(size_t idx = 0; idx < coloring.size(); idx++) {
            if(!graph.adjacent(coloring[idx], coloring[work_cell]))
                coloring.m_permutation.swap(idx, split++);

            if(idx == cell_end - 1) {
                if(split != cell_beg && split != cell_end) {
                    coloring.m_cell_end[cell_beg] = split;
                    coloring.m_cell_end[split] = cell_end;
                    coloring.m_cell_level[split] = stabilized.m_size;

                    if(is_active[cell_beg] || split - cell_beg >= cell_end - split) {
                        active_cells.push(split);
                        is_active.set(split);
                    }
                    else {
                        active_cells.push(cell_beg);
                        is_active.set(cell_beg);
                    }

                    if constexpr (!RootLevel) {
                        hash::sequential32u(invariant, cell_beg);
                        hash::sequential32u(invariant, 0);
                        hash::sequential32u(invariant, split);
                        hash::sequential32u(invariant, 1);
                    }
                }

                if(cell_end < coloring.size()) {
                    cell_beg = split = cell_end;
                    cell_end = coloring.m_cell_end[idx + 1];
                }
            }
        }

        for(size_t idx = 0; idx < coloring.size(); idx = coloring.m_cell_end[idx]) {
            assert(coloring.m_cell_end[idx] > idx);
            assert(coloring.m_cell_end[idx] <= coloring.size());
        }
    }

    template<bool RootLevel>
    void refine2(size_t work_cell, Vector<T>& active_cells, BitArray& is_active, HashType& invariant) {
#ifdef DEBUG_OUT
        std::cout << std::string(stabilized.m_size, ' ') << "refine2" << std::endl;
#endif
        if constexpr (!RootLevel)
            hash::sequential32u(invariant, work_cell);

        size_t cell_beg = 0, cell_end = coloring.m_cell_end[0];
        size_t split1 = 0, split2 = cell_end; // [cell_beg, split1), [split1, idx), [split2, cell_end)
        size_t idx = 0;
        while(idx < coloring.size()) {
            uint8_t adj_count = (uint8_t) graph.adjacent(coloring[idx], coloring[work_cell]) +
                                (uint8_t) graph.adjacent(coloring[idx], coloring[work_cell + 1]);
            if(adj_count == 0)
                coloring.m_permutation.swap(idx++, split1++);
            else if(adj_count == 1)
                idx++;
            else
                coloring.m_permutation.swap(idx, --split2);

            if(idx == split2) {
                if(cell_beg != split1 && split1 != split2 && split2 != cell_end) {
                    coloring.m_cell_end[cell_beg] = split1;
                    coloring.m_cell_end[split1] = split2;
                    coloring.m_cell_level[split1] = stabilized.m_size;
                    coloring.m_cell_end[split2] = cell_end;
                    coloring.m_cell_level[split2] = stabilized.m_size;

                    if(is_active[cell_beg] || (split1 - cell_beg >= split2 - split1 && split1 - cell_beg >= cell_end - split2)) {
                        active_cells.push(split1);
                        is_active.set(split1);

                        active_cells.push(split2);
                        is_active.set(split2);
                    }
                    else if(split2 - split1 >= cell_end - split2) {
                        active_cells.push(cell_beg);
                        is_active.set(cell_beg);

                        active_cells.push(split2);
                        is_active.set(split2);
                    }
                    else {
                        active_cells.push(cell_beg);
                        is_active.set(cell_beg);

                        active_cells.push(split1);
                        is_active.set(split1);
                    }

                    if constexpr (!RootLevel) {
                        hash::sequential32u(invariant, cell_beg);
                        hash::sequential32u(invariant, 0);
                        hash::sequential32u(invariant, split1);
                        hash::sequential32u(invariant, 1);
                        hash::sequential32u(invariant, split2);
                        hash::sequential32u(invariant, 2);
                    }
                }
                else if((cell_beg != split1 && split1 != split2) ||
                        (split1 != split2 && split2 != cell_end) ||
                        (cell_beg != split1 && split2 != cell_end)) {
                    size_t split = split1 != cell_beg ? split1 : split2;

                    coloring.m_cell_end[cell_beg] = split;
                    coloring.m_cell_end[split] = cell_end;
                    coloring.m_cell_level[split] = stabilized.m_size;

                    if(is_active[cell_beg] || split - cell_beg >= cell_end - split) {
                        active_cells.push(split);
                        is_active.set(split);
                    }
                    else {
                        active_cells.push(cell_beg);
                        is_active.set(cell_beg);
                    }

                    if constexpr (!RootLevel) {
                        if(cell_beg != split1) {
                            hash::sequential32u(invariant, cell_beg);
                            hash::sequential32u(invariant, 0);
                        }
                        if(split1 != split2) {
                            hash::sequential32u(invariant, split1);
                            hash::sequential32u(invariant, 1);
                        }
                        if(split2 != cell_end) {
                            hash::sequential32u(invariant, split2);
                            hash::sequential32u(invariant, 2);
                        }
                    }
                }

                if(cell_end < coloring.size()) {
                    idx = cell_beg = split1 = cell_end;
                    cell_end = split2 = coloring.m_cell_end[idx];
                }
            }
        }

        for(size_t idx = 0; idx < coloring.size(); idx = coloring.m_cell_end[idx]) {
            assert(coloring.m_cell_end[idx] > idx);
            assert(coloring.m_cell_end[idx] <= coloring.size());
        }
    }

    template<bool RootLevel>
    void refineN(size_t work_cell, Vector<T>& active_cells, BitArray& is_active, HashType& invariant) {
#ifdef DEBUG_OUT
        std::cout << std::string(stabilized.m_size, ' ') << "refineN" << std::endl;
#endif

        Array<T> adj_count(coloring.size(), 0);

        size_t work_end = coloring.m_cell_end[work_cell];
        for(size_t idx = work_cell; idx != work_end; idx++) {
            T vertex = coloring[idx];
            for(auto ptr = graph.begin(vertex); ptr != graph.end(vertex); ptr++)
                adj_count[*ptr]++;
        }

        Array<T> buckets(coloring.size(), 0);
        for(auto ptr = adj_count.m_data; ptr != adj_count.m_end; ptr++)
            buckets[*ptr]++;

        for(auto ptr = buckets.m_data + 1; ptr != buckets.m_end; ptr++)
            *ptr += *(ptr - 1);

        struct vcPair {
            T vertex;
            T cell;
        };
        Array<vcPair> sorted(coloring.size());
        Array<T> cell_buckets(coloring.size());
        size_t cell_beg = 0, cell_end = coloring.m_cell_end[0];
        for(size_t idx = 0; idx < coloring.size(); idx++) {
            T vertex = coloring[idx];
            sorted[--buckets[adj_count[vertex]]] = {.vertex = vertex, .cell = (T) cell_beg};
            if(idx == cell_end - 1) {
                cell_buckets[cell_beg] = cell_end;
                cell_beg = cell_end;
                if(cell_beg < coloring.size()) cell_end = coloring.m_cell_end[cell_beg];
            }
        }

        for(size_t idx = coloring.size(); idx != 0; idx--) {
            vcPair& vertex = sorted[idx - 1];
            coloring.m_permutation.set(--cell_buckets[vertex.cell], vertex.vertex);
        }

        if constexpr (!RootLevel)
            hash::sequential32u(invariant, work_cell);

        size_t cell_prev, max_cell;
        for(cell_beg = 0; cell_beg != coloring.size(); cell_beg = cell_end) {
            cell_end = coloring.m_cell_end[cell_beg];
            assert(cell_end > cell_beg);
            cell_prev = cell_beg;

            if(is_active[cell_beg]) {
                if constexpr(!RootLevel) {
                    hash::sequential32u(invariant, cell_beg);
                    hash::sequential32u(invariant, adj_count[coloring[cell_beg]]);
                }
                for(size_t idx = cell_beg + 1; idx != cell_end; idx++)
                    if(adj_count[coloring[idx]] != adj_count[coloring[idx - 1]]) {
                        coloring.m_cell_end[cell_prev] = idx;
                        coloring.m_cell_end[idx] = cell_end;
                        coloring.m_cell_level[idx] = stabilized.m_size;
                        if constexpr(!RootLevel) {
                            hash::sequential32u(invariant, idx);
                            hash::sequential32u(invariant, adj_count[coloring[idx]]);
                        }
                        cell_prev = idx;
            assert(cell_prev < coloring.size());

                        active_cells.push(idx);
                        is_active.set(idx);
                    }
            }
            else {
                if constexpr(!RootLevel) {
                    hash::sequential32u(invariant, cell_beg);
                    hash::sequential32u(invariant, adj_count[coloring[cell_beg]]);
                }
                max_cell = cell_beg;
                for(size_t idx = cell_beg + 1; idx != cell_end; idx++) {
                    if(adj_count[coloring[idx]] != adj_count[coloring[idx - 1]]) {
                        coloring.m_cell_end[cell_prev] = idx;
                        coloring.m_cell_end[idx] = cell_end;
                        coloring.m_cell_level[idx] = stabilized.m_size;
                        if constexpr(!RootLevel) {
                            hash::sequential32u(invariant, idx);
                            hash::sequential32u(invariant, adj_count[coloring[idx]]);
                        }

                        if(coloring.cellSize(cell_prev) > coloring.cellSize(max_cell))
                            max_cell = cell_prev;
                        cell_prev = idx;
            assert(cell_prev < coloring.size());
                    }
                }
                if(coloring.cellSize(cell_prev) > coloring.cellSize(max_cell))
                    max_cell = cell_prev;

                for(size_t cell_idx = cell_beg; cell_idx != cell_end; cell_idx = coloring.m_cell_end[cell_idx])
                    if(cell_idx != max_cell) {
                        active_cells.push(cell_idx);
                        is_active.set(cell_idx);
                    }
            }
        }

        for(size_t idx = 0; idx < coloring.size(); idx = coloring.m_cell_end[idx]) {
            assert(coloring.m_cell_end[idx] > idx);
            assert(coloring.m_cell_end[idx] <= coloring.size());
        }
    }

    void refineCells(size_t work_cell, size_t work_size, Vector<T>& active_cells, BitArray& is_active, HashType& invariant) {
        if(stabilized.m_size == 0) {
            if(work_size == 1)
                refine1<true>(work_cell, active_cells, is_active, invariant);
            else if(work_size == 2)
                refine2<true>(work_cell, active_cells, is_active, invariant);
            else
                refineN<true>(work_cell, active_cells, is_active, invariant);
        }
        else {
            if(work_size == 1)
                refine1<false>(work_cell, active_cells, is_active, invariant);
            else if(work_size == 2)
                refine2<false>(work_cell, active_cells, is_active, invariant);
            else
                refineN<false>(work_cell, active_cells, is_active, invariant);
        }
    }

    uint8_t refine() {
        Vector<T> active_cells(coloring.size());
        BitArray is_active(coloring.size());

        if(stabilized.m_size > 0) {
            T cell_idx = coloring.m_permutation.m_inverse[stabilized.back()];
            active_cells.push(cell_idx);
            is_active.set(cell_idx);
        }
        else {
            for(size_t cell_idx = 0; cell_idx != coloring.size(); cell_idx = coloring.m_cell_end[cell_idx]) {
                active_cells.push(cell_idx);
                is_active.set(cell_idx);
            }
        }

        HashType invariant = 0;
        while(active_cells.m_size > 0) {
            size_t work_cell = active_cells.back();
            active_cells.pop();
            is_active.clear(work_cell);

            size_t work_size = coloring.m_cell_end[work_cell] - work_cell;
            refineCells(work_cell, work_size, active_cells, is_active, invariant);
#ifdef DEBUG_OUT
            std::cout << std::string(stabilized.m_size, ' ') << coloring << " : " << work_cell << std::endl;
#endif
        }
        invariants.push(invariant);

        return 0;
    }

    void unrefine() {
        size_t cell_beg = 0, cell_end = coloring.m_cell_end[0];
        size_t level = stabilized.m_size;
        while(cell_end != coloring.size()) {
            if(coloring.m_cell_level[cell_end] >= level) {
                cell_end = coloring.m_cell_end[cell_end];
                coloring.m_cell_end[cell_beg] = cell_end;
            }
            else {
                cell_beg = cell_end;
                cell_end = coloring.m_cell_end[cell_beg];
            }
        }
        if(invariants.m_size > 0)
            invariants.pop();
    }

    Array<T> targetCell(size_t& cell_idx) {
        cell_idx = 0;
        while(cell_idx < coloring.size() && coloring.m_cell_end[cell_idx] - cell_idx == 1)
            cell_idx++;
        if(cell_idx == coloring.size())
            return Array<T>(0);

        auto beg_ptr = coloring.m_permutation.m_forward.m_data + cell_idx;
        auto end_ptr = coloring.m_permutation.m_forward.m_data + coloring.m_cell_end[cell_idx];
        Array<T> cell_content(end_ptr - beg_ptr);

        std::copy(beg_ptr, end_ptr, cell_content.m_data);

        return cell_content;
    }

#ifdef QT_QML_DEBUG
    void test() {
#ifdef DEBUG_OUT
        /*for(size_t cell = 0; cell < graph.m_vertices; cell = coloring.m_cell_end[cell]) {
            std::cout << (size_t)coloring.m_cell_level[cell] << "| ";
            for(size_t idx = cell; idx < coloring.m_cell_end[cell]; idx++)
                std::cout << (size_t) coloring.m_permutation[idx] << ' ';
        }
        std::cout << std::endl;

        for(size_t vert = 0; vert < graph.m_vertices; vert++) {
            std::cout << "Neighbors of " << vert << ": ";
            for(auto ptr = graph.begin(vert); ptr != graph.end(vert); ptr++)
                std::cout << (size_t)*ptr << ' ';
            std::cout << std::endl;
        }

        for(size_t u = 0; u < graph.m_vertices; u++) {
            for(size_t v = 0; v < graph.m_vertices; v++)
                std::cout << graph.adjacent(u, v) << ' ';
            std::cout << std::endl;
        }

        refine();
        std::cout << coloring << std::endl;

        individualize(5, 1);
        std::cout << coloring << std::endl;

        refine();
        std::cout << coloring << std::endl;

        unrefine();
        std::cout << coloring << std::endl;

        unindividualize(5);
        std::cout << coloring << std::endl;*/

        solve();

        for(size_t u = 0; u < graph.m_vertices; u++) {
            for(size_t v = 0; v < graph.m_vertices; v++)
                std::cout << graph.adjacent(max_node.permutation.m_inverse[u], max_node.permutation.m_inverse[v]) << ' ';
            std::cout << std::endl;
        }
#endif

    }
#endif

    // Algorithm input
    Graph<T> graph;
    Coloring<T> coloring;

    // Algorithm state
    Vector<T> invariants;
    Vector<T> stabilized;

    NodePath max_node;
    NodePath fst_node;
};

} // namespace

#endif // ALGORITHMS_H
