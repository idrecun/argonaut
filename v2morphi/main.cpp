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

int main()
{
    srand(time(0));
    morphi::global_alloc.reserve(1ull << 20);
    // morphi::Array<int> arr = make_array();

    // test_permutation();
    // test_bitarray_onebyte();

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
    infile.open("/home/idrecun/repos/morphi/tests/undirected_dim/mz/mz-4");
    morphi::AlgorithmSelector selector(infile, {.relabel = true});
    infile.close();

    selector.run();

#ifndef DEBUG_OUT
    for(auto ptr = selector.m_canon.m_forward.m_data; ptr != selector.m_canon.m_forward.m_end; ptr++)
        std::cout << (size_t) *ptr << ' ';
    std::cout << std::endl;
#endif

    return 0;
}
