TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

debug {
    #DEFINES += DEBUG_OUT
}

SOURCES += \
        main.cpp

HEADERS += \
    algorithm_selector.h \
    algorithms.h \
    array.h \
    bitarray.h \
    coloring.h \
    graph.h \
    hash.h \
    matrix.h \
    memory.h \
    permutation.h \
    vector.h
