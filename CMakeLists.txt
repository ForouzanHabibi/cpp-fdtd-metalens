cmake_minimum_required(VERSION 3.16)
project(cpp_fdtd_metalens LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP)

add_executable(unit_cell_library_verification
    src/unit_cell_library_verification.cpp)

add_executable(full_device_subpixel_33nm
    src/full_device_subpixel_33nm.cpp)

add_executable(full_device_subpixel_25nm
    src/full_device_subpixel_25nm.cpp)

add_executable(spectral_band_localization
    src/spectral_band_localization.cpp)

if(OpenMP_CXX_FOUND)
    target_link_libraries(unit_cell_library_verification PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(full_device_subpixel_33nm PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(full_device_subpixel_25nm PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(spectral_band_localization PUBLIC OpenMP::OpenMP_CXX)
endif()
