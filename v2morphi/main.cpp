#include <iostream>
#include <fstream>
#include <sstream>

#include "array.h"
#include "permutation.h"
#include "algorithms.h"
#include "algorithm_selector.h"

void test_permutation() {
    morphi::Permutation<uint32_t> p(4);
    p.set(0, 3);
    p.set(1, 2);
    p.set(2, 0);
    p.set(3, 1);

    std::cout << "Forward permutation (operator[]): ";
    for(size_t i = 0; i < 4; i++)
        std::cout << p[i] << ' ';
    std::cout << std::endl;

    std::cout << "Forward permutation (m_forward): ";
    for(size_t i = 0; i < 4; i++)
        std::cout << p.m_forward[i] << ' ';
    std::cout << std::endl;

    std::cout << "Inverse permutation (m_inverse): ";
    for(size_t i = 0; i < 4; i++)
        std::cout << p.m_inverse[i] << ' ';
    std::cout << std::endl;

    morphi::Permutation<uint32_t> q = p.inverse();

    std::cout << "Inverse permutation (p.inverse(), operator[]): ";
    for(size_t i = 0; i < 4; i++)
        std::cout << q[i] << ' ';
    std::cout << std::endl;

    p.swap(0, 2);
    std::cout << "Forward permutation (swapped idx 0 and 2): ";
    for(size_t i = 0; i < 4; i++)
        std::cout << p[i] << ' ';
    std::cout << std::endl;
    std::cout << "Inverse permutation (swapped idx 0 and 2): ";
    for(size_t i = 0; i < 4; i++)
        std::cout << p.m_inverse[i] << ' ';
    std::cout << std::endl;
}

morphi::Array<int> make_array() {
    morphi::Array<int> tmp(10, 3);
    tmp[2] = 4;
    return tmp;
}

void test_bitarray_onebyte() {
    morphi::BitArray arr(3);
    arr.set(0);
    std::cout << (size_t)arr.m_data[0] << std::endl;
    arr.clear(0);
    std::cout << (size_t)arr.m_data[0] << std::endl;
    arr.set(1);
    std::cout << (size_t)arr.m_data[0] << std::endl;
    arr.set(2);
    std::cout << (size_t)arr.m_data[0] << std::endl;
    arr.clear(1);
    std::cout << (size_t)arr.m_data[0] << std::endl;
}

void test_graph_single(std::string filename, unsigned num_passes) {
    std::ifstream input;
    input.open(filename);
    morphi::AlgorithmSelector selector(input, {.relabel = false});
    input.close();


    std::ofstream output;
    output.open(filename + ".log");
    auto coutBuffer = std::cout.rdbuf(output.rdbuf());

    selector.run();
    std::cerr << "Finished initial pass for " << filename << std::endl;
    while(--num_passes) {
        selector.relabel();
        clock_t start_time = clock();
        selector.run();
        std::cerr << "Time elapsed: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << std::endl;
        std::cerr << "Finished pass for " << filename << ". " << (num_passes - 1) << " passes left." << std::endl;
    }

    std::cout.rdbuf(coutBuffer);
    output.close();
}

void test_graphs() {
    std::string test_files[] = {
        /*"/home/idrecun/repos/morphi/tests/undirected_dim/mz/mz-4",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-2",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-3",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-4",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-5",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-6",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-7",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-8",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-9",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-10",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-5",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-10",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-15",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-20",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-25",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-30",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-35",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-40",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-45",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-50",
        "/home/idrecun/repos/morphi/tests/undirected_dim/lattice/lattice-4",
        "/home/idrecun/repos/morphi/tests/undirected_dim/lattice/lattice-5",
        "/home/idrecun/repos/morphi/tests/undirected_dim/lattice/lattice-10",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-7",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-9",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-13",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-15",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-19",*/
        "/home/idrecun/repos/morphi/tests/undirected_dim/mz/mz-8",
    };
    unsigned num_passes = 5;
    for(auto test_file : test_files) {
        test_graph_single(test_file, num_passes);

        std::stringstream sstream;
        sstream << "diff <(head -n " << (num_passes - 1) << " " << test_file << ".log) <(tail -n " << (num_passes - 1) << " " << test_file << ".log)";
        std::cerr << sstream.str() << std::endl;
        std::system(sstream.str().c_str());

        std::cerr << "rm " + test_file + ".log" << std::endl;
        std::system(("rm " + test_file + ".log").c_str());
    }
}

int main()
{
    srand(time(0));
    morphi::global_alloc.reserve(1ull << 28);
    // morphi::Array<int> arr = make_array();

    // test_permutation();
    // test_bitarray_onebyte();

    /*
    std::string test_input_nonorbit =
            "c lorem ipsum dolor sit amet\n"
            "p edge 8 9\n"
            "e 1 2\n"
            "e 2 3\n"
            "e 3 1\n"
            "e 3 4\n"
            "e 4 5\n"
            "e 5 6\n"
            "e 6 7\n"
            "e 6 8\n"
            "e 8 7\n";
    std::string test_input_moresteps =
            "p edge 7 8\n"
            "e 1 2\n"
            "e 3 2\n"
            "e 3 4\n"
            "e 5 4\n"
            "e 5 6\n"
            "e 7 6\n"
            "e 7 1\n"
            "e 2 7\n";
    std::stringstream sstream(test_input_moresteps);

    std::fstream infile;
    infile.open("/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-15");
    morphi::AlgorithmSelector selector(infile, {.relabel = true});
    infile.close();

    selector.run();

#ifndef DEBUG_OUT
    for(auto ptr = selector.m_canon.m_forward.m_data; ptr != selector.m_canon.m_forward.m_end; ptr++)
        std::cout << (size_t) *ptr << ' ';
    std::cout << std::endl;
#endif
    */

    test_graphs();

    return 0;
}
