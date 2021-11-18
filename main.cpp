#include <iostream>
#include <fstream>
#include <sstream>

#include "array.h"
#include "permutation.h"
#include "algorithms.h"
#include "algorithm_selector.h"
#include "group.h"
#include "hash.h"

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

    morphi::Permutation<uint32_t> q(p.m_forward.m_size);
    q.copyInv(p);

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

void test_partition() {
    uint n = 10;
    morphi::Partition<uint> p(n);
    for(size_t i = 0; i < n; i++)
        std::cout << p.representative(i) << ' ';
    std::cout << std::endl;
    p.merge(0, 4);
    for(size_t i = 0; i < n; i++)
        std::cout << p.representative(i) << ' ';
    std::cout << std::endl;
    p.merge(1, 5);
    for(size_t i = 0; i < n; i++)
        std::cout << p.representative(i) << ' ';
    std::cout << std::endl;
    p.merge(2, 3);
    for(size_t i = 0; i < n; i++)
        std::cout << p.representative(i) << ' ';
    std::cout << std::endl;
    p.merge(3, 6);
    for(size_t i = 0; i < n; i++)
        std::cout << p.representative(i) << ' ';
    std::cout << std::endl;
    p.merge(4, 3);
    for(size_t i = 0; i < n; i++)
        std::cout << p.representative(i) << ' ';
    std::cout << std::endl;
    p.merge(9, 8);
    for(size_t i = 0; i < n; i++)
        std::cout << p.representative(i) << ' ';
    std::cout << std::endl;
    p.merge(9, 5);
    for(size_t i = 0; i < n; i++)
        std::cout << p.representative(i) << ' ';
    std::cout << std::endl;
    p.merge(7, 3);
    for(size_t i = 0; i < n; i++)
        std::cout << p.representative(i) << ' ';
    std::cout << std::endl;
    p.merge(6, 1);
    for(size_t i = 0; i < n; i++)
        std::cout << p.representative(i) << ' ';
    std::cout << std::endl;
}

void test_multiset_hash() {
    uint32_t h = 0;
    morphi::hash::multiset32add(h, 123);
    std::cerr << h << std::endl;
    morphi::hash::multiset32add(h, 15);
    std::cerr << h << std::endl;
    morphi::hash::multiset32sub(h, 123);
    std::cerr << h << std::endl;
    morphi::hash::multiset32sub(h, 15);
    std::cerr << h << std::endl;
    morphi::hash::multiset32add(h, 15);
    std::cerr << h << std::endl;
}

void test_group() {
    uint n = 6;
    morphi::Group<uint> g(n, 20);
    for(size_t i = 0; i < n; i++)
        std::cout << g.m_orbit_partition.representative(i) << ' ';
    std::cout << std::endl;

    morphi::Array<uint> p1(n, 0);
    std::iota(p1.m_data, p1.m_end, 0);
    std::swap(p1.m_data[0], p1.m_data[2]);
    g.push(p1);
    for(size_t i = 0; i < n; i++)
        std::cout << g.m_orbit_partition.representative(i) << ' ';
    std::cout << std::endl;
    for(size_t i = 0; i < g.m_elements; i++) {
        std::cout << "p" << i + 1 << std::endl;
        for(auto cycle_ptr = g.elemCyclesBegin(i); cycle_ptr != g.elemCyclesEnd(i); cycle_ptr++)
            std::cout << *cycle_ptr << ' ';
        std::cout << std::endl;
        for(size_t j = 0; j < n; j++)
            std::cout << g.m_elem_fixed_points[g.fixedPointIndex(i, j)];
        std::cout << std::endl;
    }

    morphi::Array<uint> p2(n, 0);
    std::iota(p2.m_data, p2.m_end, 0);
    std::swap(p2.m_data[2], p2.m_data[4]);
    g.push(p2);
    for(size_t i = 0; i < n; i++)
        std::cout << g.m_orbit_partition.representative(i) << ' ';
    std::cout << std::endl;
    for(size_t i = 0; i < g.m_elements; i++) {
        std::cout << "p" << i + 1 << std::endl;
        for(auto cycle_ptr = g.elemCyclesBegin(i); cycle_ptr != g.elemCyclesEnd(i); cycle_ptr++)
            std::cout << *cycle_ptr << ' ';
        std::cout << std::endl;
        for(size_t j = 0; j < n; j++)
            std::cout << g.m_elem_fixed_points[g.fixedPointIndex(i, j)];
        std::cout << std::endl;
    }

    morphi::Vector<uint> stabilized(n);
    stabilized.push(4);

    size_t count = 0;
    morphi::Partition<uint> stabilizer(n);
    g.updatePartition(stabilized, stabilizer, count);
    for(size_t i = 0; i < n; i++)
        std::cout << stabilizer.representative(i) << ' ';
    std::cout << std::endl;

    morphi::Array<uint> p3(n, 0);
    std::iota(p3.m_data, p3.m_end, 0);
    std::swap(p3.m_data[0], p3.m_data[2]);
    std::swap(p3.m_data[1], p3.m_data[3]);
    std::swap(p3.m_data[3], p3.m_data[5]);
    g.push(p3);
    for(size_t i = 0; i < n; i++)
        std::cout << g.m_orbit_partition.representative(i) << ' ';
    std::cout << std::endl;
    for(size_t i = 0; i < g.m_elements; i++) {
        std::cout << "p" << i + 1 << std::endl;
        for(auto cycle_ptr = g.elemCyclesBegin(i); cycle_ptr != g.elemCyclesEnd(i); cycle_ptr++)
            std::cout << *cycle_ptr << ' ';
        std::cout << std::endl;
        for(size_t j = 0; j < n; j++)
            std::cout << g.m_elem_fixed_points[g.fixedPointIndex(i, j)];
        std::cout << std::endl;
    }

    g.updatePartition(stabilized, stabilizer, count);
    for(size_t i = 0; i < n; i++)
        std::cout << stabilizer.representative(i) << ' ';
    std::cout << std::endl;
}

void test_graph_single(std::string filename, unsigned num_passes) {
    std::ifstream input;
    input.open(filename);
    morphi::AlgorithmSelector selector(input, {.relabel = false});
    input.close();


    std::ofstream output;
    output.open(filename + ".log");
    auto coutBuffer = std::cout.rdbuf(output.rdbuf());

    clock_t start_time = clock();
    selector.run();
    std::cerr << "Time elapsed: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << std::endl;
    std::cerr << "Finished initial pass for " << filename << std::endl;
    while(--num_passes) {
        selector.relabel();
        start_time = clock();
        selector.run();
        std::cerr << "Time elapsed: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << std::endl;
        std::cerr << "Finished pass for " << filename << ". " << (num_passes - 1) << " passes left." << std::endl;
    }

    std::cout.rdbuf(coutBuffer);
    output.close();
}

void test_graphs() {
    std::string test_files[] = {
        "/home/idrecun/repos/argonaut/graphs/milan.bliss",
        "/home/idrecun/repos/argonaut/graphs/regular.bliss",
        "/home/idrecun/repos/argonaut/graphs/k33.bliss",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-2",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-3",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-4",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-5",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-6",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-7",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-8",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-9",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-10",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-11",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-12",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-13",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-14",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-15",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-16",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-17",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-18",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-19",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-20",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-21",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-22",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-23",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-24",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-25",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-26",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-27",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-28",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-29",
        "/home/idrecun/repos/morphi/tests/undirected_dim/latin/latin-30",
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
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-55",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-60",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-65",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-70",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-75",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-80",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-85",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-90",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-95",
        "/home/idrecun/repos/morphi/tests/undirected_dim/grid/grid-2-100",
        "/home/idrecun/repos/morphi/tests/undirected_dim/lattice/lattice-4",
        "/home/idrecun/repos/morphi/tests/undirected_dim/lattice/lattice-5",
        "/home/idrecun/repos/morphi/tests/undirected_dim/lattice/lattice-30",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-7",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-9",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-13",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-15",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-19",
        "/home/idrecun/repos/morphi/tests/undirected_dim/mz/mz-4",
        "/home/idrecun/repos/morphi/tests/undirected_dim/mz/mz-6",
        "/home/idrecun/repos/morphi/tests/undirected_dim/mz/mz-8",
        "/home/idrecun/repos/morphi/tests/undirected_dim/mz/mz-10",
        "/home/idrecun/repos/morphi/tests/undirected_dim/mz/mz-20",
        "/home/idrecun/repos/morphi/tests/undirected_dim/mz/mz-30",
        "/home/idrecun/repos/morphi/tests/undirected_dim/mz/mz-40",
        "/home/idrecun/repos/morphi/tests/undirected_dim/mz/mz-50",
        "/home/idrecun/repos/morphi/tests/undirected_dim/ag/ag2-8",
        "/home/idrecun/repos/morphi/tests/undirected_dim/cfi/cfi-20",
        "/home/idrecun/repos/morphi/tests/undirected_dim/cfi/cfi-40",
        "/home/idrecun/repos/morphi/tests/undirected_dim/cfi/cfi-60",
        "/home/idrecun/repos/morphi/tests/undirected_dim/cfi/cfi-80",
        "/home/idrecun/repos/morphi/tests/undirected_dim/cfi/cfi-100",
        "/home/idrecun/repos/morphi/tests/undirected_dim/cfi/cfi-120",
        "/home/idrecun/repos/morphi/tests/undirected_dim/cfi/cfi-140",
        "/home/idrecun/repos/morphi/tests/undirected_dim/cfi/cfi-160",
        "/home/idrecun/repos/morphi/tests/undirected_dim/cfi/cfi-180",
        "/home/idrecun/repos/morphi/tests/undirected_dim/cfi/cfi-200",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-37",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-43",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-49",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-55",
        /*"/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-61",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-67",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-73",
        "/home/idrecun/repos/morphi/tests/undirected_dim/sts/sts-79",*/
    };
    unsigned num_passes = 3;
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
    morphi::global_alloc.reserve(1ull << 30);

    // morphi::Array<int> arr = make_array();
    // test_permutation();
    // test_bitarray_onebyte();
    test_graphs();
    // test_partition();
    // test_group();
    // test_multiset_hash();

    return 0;
}