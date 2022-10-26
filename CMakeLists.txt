cmake_minimum_required(VERSION 3.10)

project(KrylovSolvers VERSION 0.1)
include(CTest)

find_package(universal CONFIG REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(EigenUniversalIntegration 0.1.0 CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(DEBUG   "Enable Debugging"               OFF)
option(OPENMP  "Enable OpenMP"                  ON)
option(TESTS   "Enable tests"                   OFF)

if(DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
else(DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ftree-vectorize") # add -DNDEBUG later when everything works
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "aarch64")
        set(CMAKE_BUILD_PARALLEL_LEVEL 1)
    else(${CMAKE_SYSTEM_PROCESSOR} MATCHES "aarch64")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse -mavx2 -msse3 -msse2")
    endif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "aarch64")

    if(OPENMP)
        message("using openmp")
        find_package(OpenMP REQUIRED)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fopenmp")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
    endif(OPENMP)
endif(DEBUG)

message("Flags: ${CMAKE_CXX_FLAGS}")

list(APPEND USE_Type_Flags USE_Float USE_Double USE_LongDouble USE_Posit16 USE_Posit32 USE_Posit64 USE_Posit643 USE_Posit644 USE_Posit645 USE_Posit323 USE_Posit324 USE_Posit325 USE_Posit163 USE_Posit164 USE_Posit165)

list(APPEND USE_Type_Flags_no_fdp USE_Posit16 USE_Posit32 USE_Posit64 USE_Posit643 USE_Posit644 USE_Posit645 USE_Posit323 USE_Posit324 USE_Posit325 USE_Posit163 USE_Posit164 USE_Posit165)

foreach(type IN LISTS USE_Type_Flags)
    add_executable(CIM${type} src/cim.cxx)
    target_link_libraries(CIM${type} Eigen3::Eigen OpenMP::OpenMP_CXX)
    target_include_directories(CIM${type} PUBLIC "${PROJECT_BINARY_DIR}" "include" "3rd_party" universal)
    target_compile_definitions(CIM${type} PUBLIC ${type} USE_FDP=1)
    install(TARGETS CIM${type} DESTINATION bin)
endforeach()

foreach(type IN LISTS USE_Type_Flags_no_fdp)
    add_executable(CIM${type}_no_fdp src/cim.cxx)
    target_link_libraries(CIM${type}_no_fdp Eigen3::Eigen OpenMP::OpenMP_CXX)
    target_include_directories(CIM${type}_no_fdp PUBLIC "${PROJECT_BINARY_DIR}" "include" "3rd_party" universal)
    target_compile_definitions(CIM${type}_no_fdp PUBLIC ${type} USE_FDP=0)
    install(TARGETS CIM${type}_no_fdp DESTINATION bin)
endforeach()

if(TESTS)

    add_executable(TestIEEE test/tests.cxx)
    target_link_libraries(TestIEEE Eigen3::Eigen OpenMP::OpenMP_CXX EigenUniversalIntegration::EigenUniversalIntegration)
    target_include_directories(TestIEEE PUBLIC "${PROJECT_BINARY_DIR}" "include" "include/3rd_party" universal)
    target_compile_definitions(TestIEEE PUBLIC ${type} USEIEEE)
    target_compile_features(TestIEEE INTERFACE cxx_std_20)
    add_test(NAME TestIEEE COMMAND TestIEEE)
    install(TARGETS TestIEEE DESTINATION bin)

    add_executable(TestPosit test/tests.cxx)
    target_link_libraries(TestPosit Eigen3::Eigen OpenMP::OpenMP_CXX EigenUniversalIntegration::EigenUniversalIntegration)
    target_include_directories(TestPosit PUBLIC "${PROJECT_BINARY_DIR}" "include" "include/3rd_party" universal)
    target_compile_definitions(TestPosit PUBLIC ${type} USEPosit)
    target_compile_features(TestPosit INTERFACE cxx_std_20)
    add_test(NAME TestPosit COMMAND TestPosit)
    install(TARGETS TestPosit DESTINATION bin)

endif(TESTS)