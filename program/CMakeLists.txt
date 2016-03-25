cmake_minimum_required(VERSION 2.8)

project("compPhys2_proj1" CXX)



set(GCC_COMPILE_FLAGS "-O3 --std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_FLAGS}")

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

add_executable(main main/main.cpp ${SOURCES} ${HEADERS})


