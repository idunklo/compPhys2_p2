cmake_minimum_required(VERSION 2.8)

set(GCC_COMPILE_FLAGS "--std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_FLAGS}")

project("compPhys2_proj2")

set(SOURCES 
  InitialState/InitialState.cpp
  InitialState/RandomUniform.cpp
  Sampler/ProgressBar.cpp
  Sampler/Timer.cpp
  Sampler/Sampler.cpp
  System/Particle.cpp
  System/System.cpp
  WaveFunction/WaveFunction.cpp
  WaveFunction/TrialWaveFunction.cpp
  Hamiltonian/HarmonicOscillator.cpp
  Hamiltonian/Hamiltonian.cpp
  )

set(HEADERS 
  InitialState/RandomUniform.h
  InitialState/InitialState.h
  Sampler/ProgressBar.h
  Sampler/Sampler.h
  Sampler/Timer.h
  System/System.h
  System/Particle.h
  WaveFunction/TrialWaveFunction.h
  WaveFunction/WaveFunction.h
  Hamiltonian/Hamiltonian.h
  Hamiltonian/HarmonicOscillator.h
  )

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR}/lib)
add_executable(mainMPI main/main.cpp ${SOURCES} ${HEADERS})
target_link_libraries(mainMPI ${MPI_LIBRARIES})

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
#add_executable(main main/main.cpp ${SOURCES} ${HEADERS})


