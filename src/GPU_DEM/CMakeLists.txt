cmake_minimum_required(VERSION 3.18)

project(LIGGGHTS_GPU_DEM LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

option(LIGGGHTS_GPU_DEM_BUILD_TESTS "Build GPU_DEM standalone tests" ON)

add_library(liggghts_gpu_dem
  gpu_contact_pipeline.cu
  gpu_dem_context.cu
  gpu_integrator.cu
  gpu_neighbor_builder.cu
  gpu_particle_data.cu
  gpu_timestep.cpp
  gpu_wall_contact.cu
)

target_include_directories(liggghts_gpu_dem
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(liggghts_gpu_dem
  PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
)

set_target_properties(liggghts_gpu_dem PROPERTIES
  CUDA_SEPARABLE_COMPILATION ON
  POSITION_INDEPENDENT_CODE ON
)

if(LIGGGHTS_GPU_DEM_BUILD_TESTS)
  add_executable(gpu_particle_data_test
    tests/gpu_particle_data_test.cu
  )
  target_link_libraries(gpu_particle_data_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_integrator_test
    tests/gpu_integrator_test.cu
  )
  target_link_libraries(gpu_integrator_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_neighbor_builder_test
    tests/gpu_neighbor_builder_test.cpp
  )
  target_link_libraries(gpu_neighbor_builder_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_contact_pipeline_test
    tests/gpu_contact_pipeline_test.cpp
  )
  target_link_libraries(gpu_contact_pipeline_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_hooke_collision_test
    tests/gpu_hooke_collision_test.cpp
  )
  target_link_libraries(gpu_hooke_collision_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_hooke_damped_collision_test
    tests/gpu_hooke_damped_collision_test.cpp
  )
  target_link_libraries(gpu_hooke_damped_collision_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_wall_contact_test
    tests/gpu_wall_contact_test.cpp
  )
  target_link_libraries(gpu_wall_contact_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_timestep_driver_test
    tests/gpu_timestep_driver_test.cpp
  )
  target_link_libraries(gpu_timestep_driver_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_timestep_profiling_test
    tests/gpu_timestep_profiling_test.cpp
  )
  target_link_libraries(gpu_timestep_profiling_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_synthetic_benchmark
    tests/gpu_synthetic_benchmark.cpp
  )
  target_link_libraries(gpu_synthetic_benchmark PRIVATE liggghts_gpu_dem)
endif()
