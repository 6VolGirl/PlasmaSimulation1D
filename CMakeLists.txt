cmake_minimum_required(VERSION 3.30)
project(PlasmaSimulator1D)

set(CMAKE_CXX_STANDARD 20)

add_executable(PlasmaSimulator1D
        main.cpp
        FDTD1D.cpp
        DrudeADE.cpp
        PMLCoefficients.cpp
        SpectrumAnalyzer.cpp

)
