cmake_minimum_required(VERSION 3.18)

project(LIGGGHTS_GPU_DEM LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

option(LIGGGHTS_GPU_DEM_BUILD_TESTS "Build GPU_DEM standalone tests" ON)

add_library(liggghts_gpu_dem
  gpu_contact_pipeline.cu
  gpu_contact_history.cu
  gpu_dem_context.cu
  gpu_dem_driver.cpp
  gpu_execution_policy.cpp
  gpu_external_forces.cu
  gpu_integrator.cu
  gpu_liggghts_bridge.cpp
  gpu_liggghts_feature_probe.cpp
  gpu_mesh_wall.cu
  gpu_multisphere.cpp
  gpu_neighbor_builder.cu
  gpu_particle_data.cu
  gpu_particle_insertion.cu
  gpu_runtime_facade.cpp
  gpu_timestep.cpp
  gpu_wall_contact.cu
)

target_include_directories(liggghts_gpu_dem
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
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

  add_executable(gpu_liggghts_bridge_test
    tests/gpu_liggghts_bridge_test.cpp
  )
  target_link_libraries(gpu_liggghts_bridge_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_liggghts_bridge_timestep_test
    tests/gpu_liggghts_bridge_timestep_test.cpp
  )
  target_link_libraries(gpu_liggghts_bridge_timestep_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_dem_driver_test
    tests/gpu_dem_driver_test.cpp
  )
  target_link_libraries(gpu_dem_driver_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_execution_policy_test
    tests/gpu_execution_policy_test.cpp
  )
  target_link_libraries(gpu_execution_policy_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_liggghts_feature_probe_test
    tests/gpu_liggghts_feature_probe_test.cpp
  )
  target_link_libraries(gpu_liggghts_feature_probe_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_runtime_facade_policy_test
    tests/gpu_runtime_facade_policy_test.cpp
  )
  target_link_libraries(gpu_runtime_facade_policy_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_runtime_facade_cuda_test
    tests/gpu_runtime_facade_cuda_test.cpp
  )
  target_link_libraries(gpu_runtime_facade_cuda_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_integrator_test
    tests/gpu_integrator_test.cu
  )
  target_link_libraries(gpu_integrator_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_gravity_freefall_test
    tests/gpu_gravity_freefall_test.cpp
  )
  target_link_libraries(gpu_gravity_freefall_test PRIVATE liggghts_gpu_dem)

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

  add_executable(gpu_hertz_contact_test
    tests/gpu_hertz_contact_test.cpp
  )
  target_link_libraries(gpu_hertz_contact_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_tangential_history_test
    tests/gpu_tangential_history_test.cpp
  )
  target_link_libraries(gpu_tangential_history_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_hertz_collision_test
    tests/gpu_hertz_collision_test.cpp
  )
  target_link_libraries(gpu_hertz_collision_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_hooke_damped_collision_test
    tests/gpu_hooke_damped_collision_test.cpp
  )
  target_link_libraries(gpu_hooke_damped_collision_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_wall_contact_test
    tests/gpu_wall_contact_test.cpp
  )
  target_link_libraries(gpu_wall_contact_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_mesh_wall_contact_test
    tests/gpu_mesh_wall_contact_test.cpp
  )
  target_link_libraries(gpu_mesh_wall_contact_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_stl_loader_test
    tests/gpu_stl_loader_test.cpp
  )
  target_link_libraries(gpu_stl_loader_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_multisphere_template_test
    tests/gpu_multisphere_template_test.cpp
  )
  target_link_libraries(gpu_multisphere_template_test PRIVATE liggghts_gpu_dem)

  add_executable(gpu_hertz_mesh_timestep_test
    tests/gpu_hertz_mesh_timestep_test.cpp
  )
  target_link_libraries(gpu_hertz_mesh_timestep_test PRIVATE liggghts_gpu_dem)

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
