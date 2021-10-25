#ifndef ALGORITHM_SELECTOR_H
#define ALGORITHM_SELECTOR_H

#include <cstdint>
#include <limits>
#include <string>
#include <iostream>

#include "array.h"
#include "vector.h"
#include "permutation.h"
#include "algorithms.h"

namespace morphi {

class AlgorithmSelector {
public:

    struct Options;

    // AlgorithmSelector(command line options, input stream, other stuff)
    AlgorithmSelector(std::istream& input, const Options& opt) : m_options(opt) {
        fromStream(input);
        m_canon = Permutation<uint32_t>(m_vertices);
    }

    void fromStream(std::istream& input) {
        std::string s;
        input >> s;
        while(s == "c") {
            input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            input >> s;
        }

        // trebalo bi da je s == "p"
        // procitaj i skipuj "edge"
        input >> s;

        size_t vertices, edges;
        input >> vertices >> edges;

        Array<uint32_t> colors(2 * vertices, 0);
        Vector<uint32_t> edge_list(2 * edges);

        for(size_t idx = 0; idx < vertices; idx++)
            colors.m_data[2 * idx] = idx;

        uint32_t a, b;
        while(input >> s) {
            input >> a >> b;
            if(s == "n")
                colors[2 * a - 1] = b;
            else if(s == "e") {
                edge_list.push(a - 1);
                edge_list.push(b - 1);
            }
        }

        m_vertices = vertices;
        m_edges = edges;
        m_colors = std::move(colors);
        m_edge_list = std::move(edge_list.m_array);
    }

    template<typename AlgorithmType>
    void runWith() {
        AlgorithmType solver(m_vertices, m_edges, m_edge_list, m_colors);
#ifdef DEBUG_OUT
        solver.test();
#else
        m_canon.copyFwd(solver.solve());
#endif
    }

    void run() {
        // gomilu if-else zavisno od opcija i velicine grafa
        if(m_vertices <= 0xff)
            runWith< AlgorithmDFS<uint8_t, uint32_t> >();
        else if(m_vertices <= 0xffff)
            runWith< AlgorithmDFS<uint16_t, uint32_t> >();
        else
            runWith< AlgorithmDFS<uint32_t, uint32_t> >();
    }

    // Command line options
    struct Options {
        bool RandomRelabel;
    } m_options;

    // Raw graph input
    size_t m_vertices;
    size_t m_edges;
    Array<uint32_t> m_colors;
    Array<uint32_t> m_edge_list;

    // Algorithm output
    Permutation<uint32_t> m_canon;
};

} // namespace

#endif // ALGORITHM_SELECTOR_H
