#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <cstdint>
#include <limits>
#include "vector.h"
#include "bitarray.h"
#include "coloring.h"
#include "graph.h"
#include "hash.h"
#include "partition.h"
#include "group.h"
#include "assertions.h"

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
          fst_node(vertices),
          automorphisms(vertices, 500)
    {}

    const Permutation<T>& solve() {
        solve(PathType::MaxPath | PathType::AutPath);
        std::cerr << "Tree size: " << statistics.tree_size << std::endl;
        std::cerr << "Bad nodes: " << statistics.bad_nodes << std::endl;
        std::cerr << "Aut size: " << automorphisms.m_elements << std::endl;
        std::cerr << "Orbit prunes: " << statistics.orbit_prunes << std::endl;
        std::cerr << "Max path nodes: " << statistics.max_nodes << std::endl;
        std::cerr << "Aut path nodes: " << statistics.aut_nodes << std::endl;
        std::cerr << "Max path length: " << max_node.invariants.m_size << std::endl;
        std::cerr << "Fst path length: " << fst_node.invariants.m_size << std::endl;
        return max_node.permutation;
    }

    T solve(uint8_t flags) {
        statistics.tree_size++;

        /*struct OnReturn {
            AlgorithmDFS* m_context;
            OnReturn(AlgorithmDFS* context) : m_context(context) {}
            ~OnReturn() { m_context->unrefine(); }
        } on_return(this);*/

        assert(invariants.m_size <= stabilized.m_size);
        T level = stabilized.m_size;
#ifdef DEBUG_OUT
        std::cerr << std::string(2 * level, ' ') << "I " << coloring << std::endl;
#endif
        refine();
#ifdef DEBUG_OUT
        std::cerr << std::string(2 * level, ' ') << "R " << coloring << std::endl;
        std::cerr << std::string(2 * level, ' ') << (size_t) invariants.back() << std::endl;
#endif
        /*
        if(!refine()) // zameniti sa refine() == BAD_PATH
            return level;
        */

        // update max_path, aut_path nodes
        if(!fst_node.is_leaf)
            fst_node.invariants.push(invariants.back());

        bool max_path = (bool)(flags & PathType::MaxPath) &&
                         (max_node.invariants.m_size <= level ||
                          invariants.back() >= max_node.invariants[level]);
        bool aut_path = (bool)(flags & PathType::AutPath) &&
                        fst_node.invariants.m_size > level &&
                        invariants.back() == fst_node.invariants[level];

        if(max_path)
            statistics.max_nodes++;
        if(aut_path)
            statistics.aut_nodes++;
#ifdef DEBUG_OUT
        std::cerr << std::string(2 * level, ' ');
        if(max_path) std::cerr << "MAX ";
        if(aut_path) std::cerr << "AUT ";
        if(!max_path && !aut_path) std::cerr << "BAD ";
        std::cerr << std::endl;
#endif

        if(!max_path && !aut_path) {
            statistics.bad_nodes++;
            unrefine();
            return level;
        }
        if(max_path) {
            max_node.lca_level = level;

            if(max_node.invariants.m_size > level && invariants.back() > max_node.invariants[level])
                max_node.invariants.pop(level);

            if(max_node.invariants.m_size <= level) {
                max_node.invariants.push(invariants.back());
                max_node.is_leaf = false;
            }
        }

        flags = 0;
        flags |= max_path ? PathType::MaxPath : 0;
        flags |= aut_path ? PathType::AutPath : 0;

        size_t cell_idx;
        Array<T> cell_content = targetCell(cell_idx);
#ifdef DEBUG_OUT
        std::cerr << std::string(2 * level, ' ') << "Target cell: ";
        for(auto ptr = cell_content.m_data; ptr != cell_content.m_end; ptr++)
            std::cerr << (size_t) *ptr << ' ';
        std::cerr << std::endl;
#endif
        Partition<T> orbit_partition(coloring.size());
        size_t aut_counter = 0;
        for(auto ptr = cell_content.m_data; ptr != cell_content.m_end; ptr++) {
            if(ptr > cell_content.m_data) {
                if(level != fst_node.lca_level) {
                    automorphisms.updatePartition(stabilized, orbit_partition, aut_counter);
                    if(orbit_partition.mcr(*ptr) != *ptr) {
                        statistics.orbit_prunes++;
                        continue;
                    }
                }
                else if(automorphisms.m_orbit_partition.mcr(*ptr) != *ptr) {
                    statistics.orbit_prunes++;
                    continue;
                }
            }

            individualize(cell_idx, *ptr);
            T backjump = solve(flags);
            unindividualize(cell_idx);

            if(level > backjump) {
                unrefine();
                return backjump;
            }

            fst_node.lca_level = std::min(fst_node.lca_level, level);
            max_node.lca_level = std::min(max_node.lca_level, level);
        }

        if(cell_content.m_size == 0) {
            if(max_path && (!max_node.is_leaf || (max_node.invariants.m_size == (size_t) level + 1 && graph.less(max_node.permutation.m_inverse, coloring.m_permutation.m_forward)))) {
                max_node.is_leaf = true;
                max_node.lca_level = level;
                max_node.permutation.copyInv(coloring.m_permutation);
            }

            if(!fst_node.is_leaf) {
                fst_node.is_leaf = true;
                fst_node.lca_level = level;
                fst_node.permutation.copyInv(coloring.m_permutation);
            }

#ifdef DEBUG_NO_AUT
            unrefine();
            return level;
#endif

            //if(fst_node.lca_level != level) {
                Array<T> leaf_quotient(coloring.size());
                for(size_t idx = 0; idx < leaf_quotient.m_size; idx++)
                    //leaf_quotient[idx] = fst_node.permutation.m_inverse[coloring.m_permutation.m_inverse[idx]];
                    leaf_quotient[idx] = coloring.m_permutation.m_forward[fst_node.permutation.m_forward[idx]];

                if(graph.isAutomorphism(leaf_quotient)) {
                    automorphisms.push(leaf_quotient);
                    unrefine();
                    return fst_node.lca_level;
                }
            //}

            //if(max_node.lca_level != level) {
                //Array<T> leaf_quotient(coloring.size());
                for(size_t idx = 0; idx < leaf_quotient.m_size; idx++)
                    //leaf_quotient[idx] = max_node.permutation.m_inverse[coloring.m_permutation.m_inverse[idx]];
                    leaf_quotient[idx] = coloring.m_permutation.m_forward[max_node.permutation.m_forward[idx]];

                if(graph.isAutomorphism(leaf_quotient)) {
                    automorphisms.push(leaf_quotient);
                    unrefine();
                    return automorphisms.m_orbit_partition.mcr(stabilized[fst_node.lca_level]) != stabilized[fst_node.lca_level] ? fst_node.lca_level : max_node.lca_level;
                }
            //}
        }

        unrefine();
        return level;
    }

    void individualize(size_t cell_idx, T vertex) {
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
        std::cerr << std::string(2 * stabilized.m_size, ' ') << "refine1" << std::endl;
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
                else if constexpr (!RootLevel) {
                    hash::sequential32u(invariant, cell_beg);
                    if(split == cell_end)
                        hash::sequential32u(invariant, 0);
                    if(split == cell_beg)
                        hash::sequential32u(invariant, 1);
                }

                if(cell_end < coloring.size()) {
                    cell_beg = split = cell_end;
                    cell_end = coloring.m_cell_end[idx + 1];
                }
            }
        }
#ifdef DEBUG_SLOW_ASSERTS
        assertValidColoring(coloring);
#endif
    }

    template<bool RootLevel>
    void refine2(size_t work_cell, Vector<T>& active_cells, BitArray& is_active, HashType& invariant) {
#ifdef DEBUG_OUT
        std::cerr << std::string(2 * stabilized.m_size, ' ') << "refine2" << std::endl;
#endif
        if constexpr (!RootLevel)
            hash::sequential32u(invariant, work_cell);

        size_t cell_beg = 0, cell_end = coloring.m_cell_end[0];
        size_t split1 = 0, split2 = cell_end; // [cell_beg, split1), [split1, idx), [split2, cell_end)
        size_t idx = 0;
        while(idx < coloring.size()) {
            assert(cell_beg <= split1);
            assert(split1 <= idx);
            assert(idx <= split2);
            assert(split2 <= cell_end);
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
                else if constexpr (!RootLevel) {
                    hash::sequential32u(invariant, cell_beg);
                    if(split1 == cell_end)
                        hash::sequential32u(invariant, 0);
                    if(split1 == cell_beg && split2 == cell_end)
                        hash::sequential32u(invariant, 1);
                    if(split2 == cell_beg)
                        hash::sequential32u(invariant, 2);
                }

                if(cell_end < coloring.size()) {
                    idx = cell_beg = split1 = cell_end;
                    cell_end = split2 = coloring.m_cell_end[idx];
                }
                else break;
            }
        }
#ifdef DEBUG_SLOW_ASSERTS
        assertValidColoring(coloring);
#endif
    }

    template<bool RootLevel>
    void refineN(size_t work_cell, Vector<T>& active_cells, BitArray& is_active, HashType& invariant) {
#ifdef DEBUG_OUT
        std::cerr << std::string(2 * stabilized.m_size, ' ') << "refineN" << std::endl;
#endif

#ifdef QT_QML_DEBUG
        Array<T> tmp(coloring.size());
        for(size_t cell = 0; cell < coloring.size(); cell = coloring.m_cell_end[cell])
            for(size_t idx = cell; idx < coloring.m_cell_end[cell]; idx++)
                tmp[coloring[idx]] = cell;
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

        /*Array<T> sorted(coloring.size());
        std::copy(coloring.m_permutation.m_forward.m_data, coloring.m_permutation.m_forward.m_end, sorted.m_data);
        for(size_t cell_idx = 0; cell_idx < coloring.size(); cell_idx = coloring.m_cell_end[cell_idx])
            std::sort(sorted.m_data + cell_idx, sorted.m_data + coloring.m_cell_end[cell_idx], [&adj_count](T a, T b) { return adj_count[a] == adj_count[b] ? a < b : adj_count[a] < adj_count[b]; });
        for(size_t idx = 0; idx < coloring.size(); idx++)
            coloring.m_permutation.set(idx, sorted[idx]);
        size_t cell_beg, cell_end;*/

#ifdef QT_QML_DEBUG
        for(size_t cell = 0; cell < coloring.size(); cell = coloring.m_cell_end[cell]) {
            for(size_t idx = cell; idx < coloring.m_cell_end[cell]; idx++)
                assert(tmp[coloring[idx]] == cell);
            for(size_t idx = cell + 1; idx < coloring.m_cell_end[cell]; idx++)
                assert(adj_count[coloring[idx - 1]] <= adj_count[coloring[idx]]);
        }
#endif

        if constexpr (!RootLevel)
            hash::sequential32u(invariant, work_cell);

        size_t cell_prev, max_cell;
        for(cell_beg = 0; cell_beg != coloring.size(); cell_beg = cell_end) {
            cell_end = coloring.m_cell_end[cell_beg];
            assert(cell_end > cell_beg);
            cell_prev = cell_beg;

            if constexpr(!RootLevel) {
                hash::sequential32u(invariant, cell_beg);
                hash::sequential32u(invariant, adj_count[coloring[cell_beg]]);
            }

            if(is_active[cell_beg]) {
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
#ifdef DEBUG_SLOW_ASSERTS
            assertCellSplittingValid(coloring, cell_beg, cell_end, adj_count, stabilized.m_size);
#endif
        }

#ifdef DEBUG_SLOW_ASSERTS
        assertValidColoring(coloring);
#endif
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

    void refine() {
        Vector<T> active_cells(coloring.size());
        BitArray is_active(coloring.size());

        if(stabilized.m_size > 0) {
            T cell_idx = coloring.m_permutation.m_inverse[stabilized.back()];
            active_cells.push(cell_idx);
            is_active.set(cell_idx);
            assert(coloring.m_cell_end[cell_idx] != 0);
        }
        else {
            for(size_t cell_idx = 0; cell_idx != coloring.size(); cell_idx = coloring.m_cell_end[cell_idx]) {
                active_cells.push(cell_idx);
                is_active.set(cell_idx);
                assert(coloring.m_cell_end[cell_idx] != 0);
            }
        }

        HashType invariant = 0;
        while(active_cells.m_size > 0) {
            size_t work_cell = active_cells.back();
            active_cells.pop();
            is_active.clear(work_cell);

            assert(coloring.m_cell_end[work_cell] != 0);

            size_t work_size = coloring.m_cell_end[work_cell] - work_cell;
            refineCells(work_cell, work_size, active_cells, is_active, invariant);

#ifdef DEBUG_SLOW_ASSERTS
            assertColoringSplittingValid(coloring, graph, work_cell, work_size);
#endif

            size_t cell_idx = 0;
            while(coloring.cellSize(cell_idx) == 1)
                cell_idx++;
            if(cell_idx == coloring.size())
                break;

#ifdef DEBUG_OUT
            std::cerr << std::string(2 * stabilized.m_size, ' ') << coloring << " : " << work_cell << std::endl;
#endif
        }

#ifdef DEBUG_SLOW_ASSERTS
        assertEquitableColoring(coloring, graph);
#endif
        //invariant = calculatePartialQuotientInvariant();
        invariants.push(invariant);
    }

    HashType calculateQuotientInvariant() {
        HashType invariant = 0;
        if(stabilized.m_size == 0)
            return invariant;
        for(size_t cell = 0; cell < coloring.size(); cell = coloring.m_cell_end[cell]) {
            hash::sequential32u(invariant, cell);
            hash::sequential32u(invariant, coloring.cellSize(cell));
            for(size_t oth_cell = 0; oth_cell < coloring.size(); oth_cell = coloring.m_cell_end[oth_cell]) {
                T adj_count = 0;
                for(size_t idx = cell; idx < coloring.m_cell_end[cell]; idx++)
                    adj_count += graph.adjacent(coloring.m_permutation[oth_cell], coloring.m_permutation[idx]);
                hash::sequential32u(invariant, adj_count);
            }
        }
        return invariant;
    }

    HashType calculateMorphiInvariant() {
        HashType invariant = 0;
        if(stabilized.m_size == 0)
            return invariant;
        struct pair { size_t cell; size_t color; };
        Vector<pair> new_cells(coloring.size());
        for(size_t cell_idx = 0, color = 0; cell_idx < coloring.size(); cell_idx = coloring.m_cell_end[cell_idx], color++) {
            if(coloring[cell_idx] == stabilized.back() || (cell_idx > 0 && coloring[cell_idx - 1] == stabilized.back())) {
                new_cells.push({.cell = cell_idx, .color = color });
            }
            else if(coloring.m_cell_level[cell_idx] == stabilized.m_size ||
               (coloring.m_cell_end[cell_idx] < coloring.size() && coloring.m_cell_level[coloring.m_cell_end[cell_idx]] == stabilized.m_size))
                new_cells.push({.cell = cell_idx, .color = color });
        }

        for(size_t idx = 0; idx < new_cells.m_size; idx++) {
            hash::sequential32u(invariant, new_cells[idx].color);
            hash::sequential32u(invariant, coloring.cellSize(new_cells[idx].cell));
            for(size_t jdx = 0; jdx < new_cells.m_size; jdx++) {
                T adj_count = 0;
                for(size_t kdx = new_cells[idx].cell; kdx < coloring.m_cell_end[new_cells[idx].cell]; kdx++)
                    adj_count += (uint8_t) graph.adjacent(coloring[new_cells[jdx].cell], coloring[kdx]);
                hash::sequential32u(invariant, adj_count);
            }
        }

        return invariant;
    }

    HashType calculateCellSequenceInvariant() {
        HashType invariant = 0;
        if(stabilized.m_size == 0)
            return invariant;
        for(size_t cell_idx = 0; cell_idx < coloring.size(); cell_idx++)
            hash::sequential32u(invariant, cell_idx);
        return invariant;
    }

    HashType calculateMultisetQuotientInvariant() {
        auto hashTriple = [](HashType x, HashType y, HashType z) {
            HashType hash = 0;
            hash::sequential32u(hash, x);
            hash::sequential32u(hash, y);
            hash::sequential32u(hash, z);
            return hash;
        };
        HashType invariant = 0;
        for(size_t cell_idx = 0; cell_idx < coloring.size(); cell_idx = coloring.m_cell_end[cell_idx]) {
            hash::multiset32add(invariant, cell_idx);
            for(size_t cell_jdx = 0; cell_jdx < cell_idx; cell_jdx = coloring.m_cell_end[cell_idx]) {
                HashType adj_count = 0;
                for(size_t idx = cell_jdx; idx < coloring.m_cell_end[cell_jdx]; idx++)
                    adj_count += graph.adjacent(coloring[cell_idx], coloring[idx]);
                hash::multiset32add(invariant, hashTriple(cell_idx, cell_jdx, adj_count));
            }
        }
        return invariant;
    }

    HashType calculatePartialQuotientInvariant() {
        HashType invariant = invariants.back();
        for(size_t cell_idx = 0; cell_idx < coloring.size(); cell_idx = coloring.m_cell_end[cell_idx]) {
            hash::sequential32u(invariant, cell_idx);
            for(size_t cell_jdx = 0; cell_jdx < cell_idx; cell_jdx = coloring.m_cell_end[cell_jdx]) {
                HashType adj_count = 0;
                for(size_t idx = cell_jdx; idx < coloring.m_cell_end[cell_jdx]; idx++)
                    adj_count += graph.adjacent(coloring[cell_idx], coloring[idx]);
                hash::sequential32u(invariant, adj_count);
                //hash::multiset32add(invariant, hashTriple(cell_idx, cell_jdx, adj_count));
            }
        }
        return invariant;
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

        auto beg_ptr = coloring.m_permutation.m_forward.m_data + (size_t) cell_idx;
        auto end_ptr = coloring.m_permutation.m_forward.m_data + (size_t) coloring.m_cell_end[cell_idx];
        Array<T> cell_content(end_ptr - beg_ptr);
        std::copy(beg_ptr, end_ptr, cell_content.m_data);
        std::sort(cell_content.m_data, cell_content.m_end);

        return cell_content;
    }

    // Algorithm input
    Graph<T> graph;
    Coloring<T> coloring;

    // Algorithm state
    Vector<HashType> invariants;
    Vector<T> stabilized;

    NodePath max_node;
    NodePath fst_node;

    Group<T> automorphisms;

    // Statistics
    struct Statistics {
        size_t tree_size = 0;
        size_t bad_nodes = 0;
        size_t orbit_prunes = 0;
        size_t max_nodes = 0;
        size_t aut_nodes = 0;
    } statistics;
};

} // namespace

#endif // ALGORITHMS_H
