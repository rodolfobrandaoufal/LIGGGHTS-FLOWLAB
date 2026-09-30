# ----------------------------------------------------------------------
#    This is the
#
#    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
#    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
#    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
#    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
#    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
#    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®
#
#    DEM simulation engine, released by
#    DCS Computing Gmbh, Linz, Austria
#    http://www.dcs-computing.com, office@dcs-computing.com
#
#    LIGGGHTS® is part of CFDEM®project:
#    http://www.liggghts.com | http://www.cfdem.com
#
#    Core developer and main author:
#    Christoph Kloss, christoph.kloss@dcs-computing.com
#
#    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
#    License, version 2 or later. It is distributed in the hope that it will
#    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
#    received a copy of the GNU General Public License along with LIGGGHTS®.
#    If not, see http://www.gnu.org/licenses . See also top-level README
#    and LICENSE files.
#
#    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
#    the producer of the LIGGGHTS® software and the CFDEM®coupling software
#    See http://www.cfdem.com/terms-trademark-policy for details.
#
#-------------------------------------------------------------------------
#    Contributing author:
#    Tóth János (MATE, Gödöllő)
#-------------------------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.5.0)
PROJECT(LIGGGHTS C CXX)

# Default to an optimised build; the optimisation level comes only from
# CMAKE_BUILD_TYPE (Release: -O3 -DNDEBUG), never from hard-coded flags (P0-08).
GET_PROPERTY(LIGGGHTS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
IF(NOT CMAKE_BUILD_TYPE AND NOT LIGGGHTS_MULTI_CONFIG)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Release, RelWithDebInfo, Debug, MinSizeRel)" FORCE)
  MESSAGE(STATUS "CMAKE_BUILD_TYPE not set, using Release")
ENDIF()

SET(LIGGGHTS_MAJOR_VERSION 3)
SET(LIGGGHTS_MINOR_VERSION 8)
SET(LIGGGHTS_PATCH_VERSION 0)
SET(LIGGGHTS_VERSION ${LIGGGHTS_MAJOR_VERSION}.${LIGGGHTS_MINOR_VERSION}.${LIGGGHTS_PATCH_VERSION})

# NOTE: Implement these?
# https://www.cfdem.com/media/DEM/docu/Section_start.html#making-liggghts-r-public
# -DLAMMPS_MEMALIGN
# -DLAMMPS_XDR
# -DPACK_ARRAY
# -DPACK_POINTER
# -DPACK_MEMCPY

#=======================================
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cMake")

# Generated headers (style_*.h, style_contact_model.h, version_liggghts.h) are
# written into the build tree, never into src/ (P0-06). Several build trees
# with different options can therefore share one source tree, and a CMake
# configure never touches the headers of the Make route.
SET(LIGGGHTS_GENERATED_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
FILE(MAKE_DIRECTORY "${LIGGGHTS_GENERATED_DIR}")
INCLUDE_DIRECTORIES(BEFORE "${LIGGGHTS_GENERATED_DIR}")

# '#include "style_fix.h"' is looked up in the directory of the including file
# (src/) BEFORE any -I path, so generated headers left in src/ by the Make
# route ("make <machine>" runs Make.sh style/models) or by an old CMake would
# silently replace the ones generated here (different whitelist, style list
# and version banner).
OPTION(LIGGGHTS_ALLOW_SOURCE_TREE_HEADERS
       "Only warn (instead of failing) when src/ contains generated style_*.h / version_liggghts.h that shadow the CMake-generated ones" OFF)
FILE(GLOB LIGGGHTS_STALE_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
     "${CMAKE_CURRENT_SOURCE_DIR}/style_*.h" "${CMAKE_CURRENT_SOURCE_DIR}/version_liggghts.h")
IF(LIGGGHTS_STALE_HEADERS)
  STRING(REPLACE ";" " " stale_list "${LIGGGHTS_STALE_HEADERS}")
  SET(stale_msg
    "The source directory ${CMAKE_CURRENT_SOURCE_DIR} contains generated headers "
    "(from the Make route or an older CMake build): ${stale_list}\n"
    "They take precedence over the headers CMake generates in ${LIGGGHTS_GENERATED_DIR}, "
    "so this build would silently use their style lists, contact-model whitelist and version banner.\n"
    "Remove them (they are untracked; 'make <machine>' re-creates them for the Make route):\n"
    "  rm ${CMAKE_CURRENT_SOURCE_DIR}/style_*.h ${CMAKE_CURRENT_SOURCE_DIR}/version_liggghts.h\n"
    "or configure with -DLIGGGHTS_ALLOW_SOURCE_TREE_HEADERS=ON to build with them anyway.")
  IF(LIGGGHTS_ALLOW_SOURCE_TREE_HEADERS)
    MESSAGE(WARNING ${stale_msg})
  ELSE()
    MESSAGE(FATAL_ERROR ${stale_msg})
  ENDIF()
ENDIF()

INCLUDE(Version)
INCLUDE(Style)
INCLUDE(Model)

#=======================================
SET(ENABLED_OPTIONS "")
SET(DISABLED_OPTIONS "")

SET(ENABLED_MODELS "")

OPTION(DISABLE_ALL "Disable all features/models (OFF by default)" OFF)
OPTION(ENABLE_ALL  "Enable all features/models (OFF by default)"  OFF)

SET(DEFAULT_ON ON CACHE INTERNAL "Default ON value")
SET(DEFAULT_OFF OFF CACHE INTERNAL "Default OFF value")

IF(DISABLE_ALL)
  SET(DEFAULT_ON OFF)
ENDIF()

IF(ENABLE_ALL)
  SET(DEFAULT_OFF ON)
ENDIF()

OPTION(ENABLE_MPI    "Use MPI"      ${DEFAULT_ON})
# ENABLE_VTK: ON = require VTK, OFF = no VTK dumps, AUTO (default) = use VTK
# when it is found, so that a plain configure works on machines without VTK.
IF(DEFAULT_ON)
  SET(ENABLE_VTK AUTO CACHE STRING "Use dump_vtk: ON (required), OFF, AUTO (if found)")
ELSE()
  SET(ENABLE_VTK OFF CACHE STRING "Use dump_vtk: ON (required), OFF, AUTO (if found)")
ENDIF()
SET_PROPERTY(CACHE ENABLE_VTK PROPERTY STRINGS ON OFF AUTO)
OPTION(ENABLE_JPEG   "Use libjpeg"  ${DEFAULT_OFF})
OPTION(ENABLE_PNG    "Use libpng"   ${DEFAULT_OFF})
OPTION(ENABLE_FFMPEG "Use ffmpeg"   ${DEFAULT_OFF})
OPTION(ENABLE_GZIP   "Use gzip"     ${DEFAULT_OFF})

OPTION(ENABLE_SQ  "Use Superquadrics" ${DEFAULT_OFF})

OPTION(ENABLE_SMALLBIG   "Use of 4-byte (small) vs 8-byte (big) integers."        ON)
OPTION(ENABLE_BIGBIG     "Huge molecular systems with more than 2 billion atoms." OFF)
OPTION(ENABLE_SMALLSMALL "Machine does not support 64-bit integers."              OFF)

OPTION(ENABLE_LONGLONG_TO_LONG "System does not recognize “long long” data types." OFF)

# Modernization-branch options (documented in doc/Section_start.txt)
OPTION(LIGGGHTS_ENABLE_HDF5 "Build dump hdf5 and dump mesh/hdf5 (requires a parallel HDF5)" OFF)
OPTION(LIGGGHTS_ENABLE_GPU_DEM "Also configure the experimental src/GPU_DEM subproject" OFF)
OPTION(LIGGGHTS_FAST_MATH "Compile with -ffast-math (changes floating-point results; not bitwise reproducible with the Make route)" OFF)
OPTION(LIGGGHTS_NATIVE_ARCH "Compile with -march=native (binary only runs on CPUs like the build host)" OFF)
OPTION(LIGGGHTS_MATH_ERRNO "Keep errno semantics for libm calls. OFF adds -fno-math-errno: no measured speed-up, and NOT bitwise identical for some contact models (hooke/hysteresis); see audit/fixes/phaseC/kernel/REPORT.md" ON)
OPTION(LIGGGHTS_NO_CONTACT_MODEL_FALLBACK "Strict mode: contact-model combinations outside the whitelist are an error instead of running through the runtime fallback" OFF)
OPTION(LIGGGHTS_ENABLE_OPENMP "OpenMP threading of pair gran and fix wall/gran (runtime: OMP_NUM_THREADS or 'package omp'; doc/package_omp.txt). OFF builds exactly the serial code" OFF)
OPTION(LIGGGHTS_OPENMP_SERIAL_BITWISE "With LIGGGHTS_ENABLE_OPENMP: compile the threaded pair kernel with -fno-inline so that its results are bitwise identical to the serial kernel (costs ~10% Pair time; OFF: identical for all thread counts >= 2, but some models then differ in the last bits from 1 thread)" ON)
SET(LIGGGHTS_SANITIZE "" CACHE STRING "Comma-separated -fsanitize= list, e.g. address,undefined (empty = none)")
OPTION(LIGGGHTS_ENABLE_TESTING "Register the tests in ../tests with CTest" ON)

# 9 normal models
OPTION(ENABLE_MODEL_HERTZ "Hertz Model" ${DEFAULT_ON})
OPTION(ENABLE_MODEL_HOOKE "Hooke Model" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_HOOKE_HYSTERESIS "Hooke Hysteresis Model" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_HOOKE_STIFFNESS "Hooke Stiffness Model" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_EDINBURGH "Edinburgh Model" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_EDINBURGH_STIFFNESS "Edinburgh Stiffness Model" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_HERTZ_STIFFNESS "Hertz Stiffness Model" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_LUDING "Luding Model" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_THORNTON_NING "Thornton Ning Model" ${DEFAULT_OFF})

# 4 tangential models
OPTION(ENABLE_MODEL_TANGENTIAL_OFF "Tangential Model Off" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_TANGENTIAL_HISTORY "Tangential Model History" ${DEFAULT_ON})
OPTION(ENABLE_MODEL_TANGENTIAL_NO_HISTORY "Tangential Model No History" ${DEFAULT_ON})
OPTION(ENABLE_MODEL_TANGENTIAL_LUDING "Tangential Model Luding" ${DEFAULT_OFF})

# 5 cohesion models
OPTION(ENABLE_MODEL_COHESION_OFF "Cohesion Model Off" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_COHESION_SJKR "Cohesion Model SJKR" ${DEFAULT_ON})
OPTION(ENABLE_MODEL_COHESION_SJKR2 "Cohesion Model SJKR2" ${DEFAULT_ON})
OPTION(ENABLE_MODEL_COHESION_WASHINO_CAPILLARY_VISCOUS "Cohesion Model Washino/capillary/viscous" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_COHESION_EASO_CAPILLARY_VISCOUS "Cohesion Model Easo/capillary/viscous" ${DEFAULT_OFF})

# 6 rolling models
OPTION(ENABLE_MODEL_ROLLING_OFF "Rolling Model Off" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_ROLLING_CDT "Rolling Model CDT" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_ROLLING_EPSD "Rolling Model EPSD" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_ROLLING_EPSD2 "Rolling Model EPSD2" ${DEFAULT_ON})
OPTION(ENABLE_MODEL_ROLLING_EPSD3 "Rolling Model ESPD3" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_ROLLING_LUDING "Rolling Model Luding" ${DEFAULT_OFF})

# 3 surface models
OPTION(ENABLE_MODEL_SURFACE_DEFAULT "Surface Model Default" ${DEFAULT_ON})
OPTION(ENABLE_MODEL_SURFACE_SUPERQUADRIC "Surface Model Superquadric" ${DEFAULT_OFF})
OPTION(ENABLE_MODEL_SURFACE_MULTICONTACT "Surface Model Multicontact" ${DEFAULT_OFF})

#=======================================
IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0.0)
    SET(CMAKE_CXX_STANDARD 14)
  ELSE()
    SET(CMAKE_CXX_STANDARD 17)
  ENDIF()
  SET(CMAKE_CXX_STANDARD_REQUIRED ON)  

  SET(OPT_CXX_FLAGS "")
  MESSAGE(STATUS "Compiler: clang")
ELSEIF(("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
  IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0.0)
    SET(CMAKE_CXX_STANDARD 14)
  ELSE()
    SET(CMAKE_CXX_STANDARD 17)
  ENDIF()
  SET(CMAKE_CXX_STANDARD_REQUIRED ON)

  # no -O level and no -ffast-math here (P0-08): the optimisation level comes
  # from CMAKE_BUILD_TYPE, fast math only from LIGGGHTS_FAST_MATH
  SET(OPT_CXX_FLAGS "-funroll-loops -fstrict-aliasing -Wall -Wextra -Wno-uninitialized -Wno-literal-suffix -fPIC")
  MESSAGE(STATUS "Compiler: gcc")
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	SET(OPT_CXX_FLAGS "-O3 -funroll-loops -fstrict-aliasing -W -Wno-uninitialized -march=native -xHost")
  MESSAGE(STATUS "Compiler: intel")
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  SET(OPT_CXX_FLAGS "/O2 /bigobj")

  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES -DNOMINMAX)
  ADD_DEFINITIONS("/wd4101 /wd4244 /wd4267")

  MESSAGE(STATUS "Compiler: msvc")
ENDIF()

IF(("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang"))
  IF(LIGGGHTS_FAST_MATH)
    SET(OPT_CXX_FLAGS "${OPT_CXX_FLAGS} -ffast-math")
    SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} FAST_MATH")
  ENDIF()
  IF(NOT LIGGGHTS_MATH_ERRNO)
    # roadmap B5 / PF-05: opt-in. sqrt is correctly rounded either way, but with
    # this flag 3 hooke/hysteresis combinations of tests/kernel/model_matrix
    # changed in the last bits (cause not isolated; likely different FMA
    # contraction once libm calls are treated as pure)
    SET(OPT_CXX_FLAGS "${OPT_CXX_FLAGS} -fno-math-errno")
    SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} NO_MATH_ERRNO")
  ENDIF()
  IF(LIGGGHTS_NATIVE_ARCH)
    SET(OPT_CXX_FLAGS "${OPT_CXX_FLAGS} -march=native")
    SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} NATIVE_ARCH")
  ENDIF()
  IF(LIGGGHTS_SANITIZE)
    SET(OPT_CXX_FLAGS "${OPT_CXX_FLAGS} -fsanitize=${LIGGGHTS_SANITIZE} -fno-omit-frame-pointer")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=${LIGGGHTS_SANITIZE}")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=${LIGGGHTS_SANITIZE}")
    SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} SANITIZE=${LIGGGHTS_SANITIZE}")
  ENDIF()
ELSEIF(LIGGGHTS_FAST_MATH OR LIGGGHTS_NATIVE_ARCH OR LIGGGHTS_SANITIZE)
  MESSAGE(WARNING "LIGGGHTS_FAST_MATH / LIGGGHTS_NATIVE_ARCH / LIGGGHTS_SANITIZE are only implemented for GCC and Clang")
ENDIF()

IF(LIGGGHTS_NO_CONTACT_MODEL_FALLBACK)
  ADD_DEFINITIONS(-DLIGGGHTS_NO_CONTACT_MODEL_FALLBACK)
  SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} NO_CONTACT_MODEL_FALLBACK")
ENDIF()

# roadmap C2/B8: OpenMP threading of the granular pair and wall kernels.
# With the option OFF no OpenMP flag and no LIGGGHTS_OMP define is added, so
# the preprocessed sources equal the serial build.
IF(LIGGGHTS_ENABLE_OPENMP)
  FIND_PACKAGE(OpenMP REQUIRED COMPONENTS CXX)
  SET(OPT_CXX_FLAGS "${OPT_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ADD_DEFINITIONS(-DLIGGGHTS_OMP)
  # threaded pair kernel: no inlining, so every contact sub-model call goes to
  # the same out-of-line code as the serial kernel in lammps.cpp (first in
  # link order); inlining would change FMA contraction and break the bitwise
  # equality of the deterministic mode with the serial kernel
  IF(LIGGGHTS_OPENMP_SERIAL_BITWISE)
    SET_SOURCE_FILES_PROPERTIES(pair_gran_omp.cpp PROPERTIES COMPILE_OPTIONS "-fno-inline")
  ELSE()
    SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} OPENMP_NOT_SERIAL_BITWISE")
  ENDIF()
  SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} OPENMP")
ELSE()
  SET(DISABLED_OPTIONS "${DISABLED_OPTIONS} OPENMP")
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPT_CXX_FLAGS}")

IF(WIN32)
  INCLUDE_DIRECTORIES(WINDOWS/extra)

  FIND_PACKAGE(Cygwin REQUIRED)

  FOREACH(flag_var
    CMAKE_C_FLAGS_RELEASE CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    SET(${flag_var} "${${flag_var}} /MP4")
  ENDFOREACH()
ENDIF()

#=======================================
# based on: https://github.com/schrummy14/LIGGGHTS_Flexible_Fibers/blob/master/src/WINDOWS/CMake_patch.zip

SCAN_STYLES()
GENERATE_STYLES()

#=======================================
# based on: https://github.com/schrummy14/LIGGGHTS_Flexible_Fibers/blob/master/src/WINDOWS/CMake_patch.zip

IF(ENABLE_MODEL_HERTZ)
  ADD_CONTACT_MODEL("HERTZ" "" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} HERTZ")
ENDIF()
IF(ENABLE_MODEL_HOOKE)
  ADD_CONTACT_MODEL("HOOKE" "" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} HOOKE")
ENDIF()
IF(ENABLE_MODEL_HOOKE_HYSTERESIS)
  ADD_CONTACT_MODEL("HOOKE_HYSTERESIS" "" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} HOOKE_HYSTERESIS")
ENDIF()
IF(ENABLE_MODEL_HOOKE_STIFFNESS)
  ADD_CONTACT_MODEL("HOOKE_STIFFNESS" "" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} HOOKE_STIFFNESS")
ENDIF()
IF(ENABLE_MODEL_EDINBURGH)
  ADD_CONTACT_MODEL("EDINBURGH" "" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} EDINBURGH")
ENDIF()
IF(ENABLE_MODEL_EDINBURGH_STIFFNESS)
  ADD_CONTACT_MODEL("EDINBURGH_STIFFNESS" "" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} EDINBURGH_STIFFNESS")
ENDIF()
IF(ENABLE_MODEL_HERTZ_STIFFNESS)
  ADD_CONTACT_MODEL("HERTZ_STIFFNESS" "" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} HERTZ_STIFFNESS")
ENDIF()
IF(ENABLE_MODEL_LUDING)
  ADD_CONTACT_MODEL("LUDING" "" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} LUDING")
ENDIF()
IF(ENABLE_MODEL_THORNTON_NING)
  ADD_CONTACT_MODEL("THORNTON_NING" "" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} THORNTON_NING")
ENDIF()

# 4 tangential models
IF(ENABLE_MODEL_TANGENTIAL_OFF)
  ADD_CONTACT_MODEL("" "TANGENTIAL_OFF" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} TANGENTIAL_OFF")
ENDIF()
IF(ENABLE_MODEL_TANGENTIAL_HISTORY)
  ADD_CONTACT_MODEL("" "TANGENTIAL_HISTORY" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} TANGENTIAL_HISTORY")
ENDIF()
IF(ENABLE_MODEL_TANGENTIAL_NO_HISTORY)
  ADD_CONTACT_MODEL("" "TANGENTIAL_NO_HISTORY" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} TANGENTIAL_NO_HISTORY")
ENDIF()
IF(ENABLE_MODEL_TANGENTIAL_LUDING)
  ADD_CONTACT_MODEL("" "TANGENTIAL_LUDING" "" "" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} TANGENTIAL_LUDING")
ENDIF()

# 5 cohesion models
IF(ENABLE_MODEL_COHESION_OFF)
  ADD_CONTACT_MODEL("" "" "COHESION_OFF" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} COHESION_OFF")
ENDIF()
IF(ENABLE_MODEL_COHESION_SJKR)
  ADD_CONTACT_MODEL("" "" "COHESION_SJKR" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} COHESION_SJKR")
ENDIF()
IF(ENABLE_MODEL_COHESION_SJKR2)
  ADD_CONTACT_MODEL("" "" "COHESION_SJKR2" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} COHESION_SJKR2")
ENDIF()
IF(ENABLE_MODEL_COHESION_WASHINO_CAPILLARY_VISCOUS)
  ADD_CONTACT_MODEL("" "" "COHESION_WASHINO_CAPILLARY_VISCOUS" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} COHESION_WASHINO_CAPILLARY_VISCOUS")
ENDIF()
IF(ENABLE_MODEL_COHESION_EASO_CAPILLARY_VISCOUS)
  ADD_CONTACT_MODEL("" "" "COHESION_EASO_CAPILLARY_VISCOUS" "" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} COHESION_EASO_CAPILLARY_VISCOUS")
ENDIF()

# 6 rolling models
IF(ENABLE_MODEL_ROLLING_OFF)
  ADD_CONTACT_MODEL("" "" "" "ROLLING_OFF" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} ROLLING_OFF")
ENDIF()
IF(ENABLE_MODEL_ROLLING_CDT)
  ADD_CONTACT_MODEL("" "" "" "ROLLING_CDT" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} ROLLING_CDT")
ENDIF()
IF(ENABLE_MODEL_ROLLING_EPSD)
  ADD_CONTACT_MODEL("" "" "" "ROLLING_EPSD" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} ROLLING_EPSD")
ENDIF()
IF(ENABLE_MODEL_ROLLING_EPSD2)
  ADD_CONTACT_MODEL("" "" "" "ROLLING_EPSD2" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} ROLLING_EPSD2")
ENDIF()
IF(ENABLE_MODEL_ROLLING_EPSD3)
  ADD_CONTACT_MODEL("" "" "" "ROLLING_EPSD3" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} ROLLING_EPSD3")
ENDIF()
IF(ENABLE_MODEL_ROLLING_LUDING)
  ADD_CONTACT_MODEL("" "" "" "ROLLING_LUDING" "")
  SET(ENABLED_MODELS "${ENABLED_MODELS} ROLLING_LUDING")
ENDIF()

# 3 surface models
IF(ENABLE_MODEL_SURFACE_DEFAULT)
  ADD_CONTACT_MODEL("" "" "" "" "SURFACE_DEFAULT")
  SET(ENABLED_MODELS "${ENABLED_MODELS} SURFACE_DEFAULT")
ENDIF()
IF(ENABLE_MODEL_SURFACE_SUPERQUADRIC)
  ADD_CONTACT_MODEL("" "" "" "" "SURFACE_SUPERQUADRIC")
  SET(ENABLED_MODELS "${ENABLED_MODELS} SURFACE_SUPERQUADRIC")
ENDIF()
IF(ENABLE_MODEL_SURFACE_MULTICONTACT)
  ADD_CONTACT_MODEL("" "" "" "" "SURFACE_MULTICONTACT")
  SET(ENABLED_MODELS "${ENABLED_MODELS} SURFACE_MULTICONTACT")
ENDIF()

WRITE_WHITELIST()

#========================================
# Sub-projects are opt-in; subdirectories are no longer added automatically
# (P0-07), so a CPU-only configure never touches src/GPU_DEM.
IF(LIGGGHTS_ENABLE_GPU_DEM)
  IF(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/GPU_DEM/CMakeLists.txt")
    MESSAGE(FATAL_ERROR "LIGGGHTS_ENABLE_GPU_DEM=ON but ${CMAKE_CURRENT_SOURCE_DIR}/GPU_DEM/CMakeLists.txt does not exist")
  ENDIF()
  ADD_SUBDIRECTORY(GPU_DEM)
  SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} GPU_DEM")
ENDIF()

#=======================================
FILE(GLOB SOURCES *.cpp)
FILE(GLOB HEADERS *.h)

ADD_LIBRARY(liggghts_obj OBJECT ${SOURCES})

ADD_LIBRARY(liggghts_static STATIC $<TARGET_OBJECTS:liggghts_obj>)
SET_TARGET_PROPERTIES(liggghts_static PROPERTIES OUTPUT_NAME liggghts)
SET_PROPERTY(TARGET liggghts_static PROPERTY PUBLIC_HEADER ${HEADERS})

ADD_LIBRARY(liggghts_shared SHARED $<TARGET_OBJECTS:liggghts_obj>)
SET_TARGET_PROPERTIES(liggghts_shared PROPERTIES OUTPUT_NAME liggghts)
SET_PROPERTY(TARGET liggghts_shared PROPERTY PUBLIC_HEADER ${HEADERS})

ADD_EXECUTABLE(liggghts_bin $<TARGET_OBJECTS:liggghts_obj>)
SET_TARGET_PROPERTIES(liggghts_bin PROPERTIES OUTPUT_NAME liggghts)

IF(LIGGGHTS_ENABLE_OPENMP)
  TARGET_LINK_LIBRARIES(liggghts_static PUBLIC OpenMP::OpenMP_CXX)
  TARGET_LINK_LIBRARIES(liggghts_shared PUBLIC OpenMP::OpenMP_CXX)
  TARGET_LINK_LIBRARIES(liggghts_bin PUBLIC OpenMP::OpenMP_CXX)
ENDIF()

#=======================================
IF(ENABLE_MPI)
  FIND_PACKAGE(MPI)

  IF(MPI_FOUND)
    INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(liggghts_static PUBLIC ${MPI_LIBRARIES})
    TARGET_LINK_LIBRARIES(liggghts_shared PUBLIC ${MPI_LIBRARIES})
    TARGET_LINK_LIBRARIES(liggghts_bin PUBLIC ${MPI_LIBRARIES})

    IF(MPI_COMPILE_FLAGS)
      SET_TARGET_PROPERTIES(liggghts_static PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
      SET_TARGET_PROPERTIES(liggghts_shared PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
      SET_TARGET_PROPERTIES(liggghts_bin PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    ENDIF()

    IF(MPI_LINK_FLAGS)
      SET_TARGET_PROPERTIES(liggghts_static PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
      SET_TARGET_PROPERTIES(liggghts_shared PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
      SET_TARGET_PROPERTIES(liggghts_bin PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    ENDIF()

    SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} MPI")
  ELSE()
    MESSAGE(FATAL_ERROR "MPI NOT found!")
  ENDIF()
ELSE()
  SET(DISABLED_OPTIONS "${DISABLED_OPTIONS} MPI")
  MESSAGE(STATUS "Using MPI stubs")
ENDIF()

#=======================================
IF("${ENABLE_VTK}" STREQUAL "AUTO")
  FIND_PACKAGE(VTK COMPONENTS NO_MODULE QUIET)
  IF(NOT VTK_FOUND)
    MESSAGE(STATUS "VTK not found: dump */vtk styles disabled (ENABLE_VTK=AUTO; set ENABLE_VTK=ON to require VTK)")
  ENDIF()
  SET(LIGGGHTS_USE_VTK ${VTK_FOUND})
ELSEIF(ENABLE_VTK)
  FIND_PACKAGE(VTK COMPONENTS NO_MODULE)
  SET(LIGGGHTS_USE_VTK ON)
ELSE()
  SET(LIGGGHTS_USE_VTK OFF)
ENDIF()

IF(LIGGGHTS_USE_VTK)

  IF(VTK_FOUND)
    IF(VTK_VERSION VERSION_LESS "8.90.0")
      INCLUDE(${VTK_USE_FILE})
    ELSE()
      INCLUDE_DIRECTORIES("${VTK_PREFIX_PATH}/include/vtk-${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}")
    ENDIF()
    
    TARGET_LINK_LIBRARIES(liggghts_static PUBLIC ${VTK_LIBRARIES})
    TARGET_LINK_LIBRARIES(liggghts_shared PUBLIC ${VTK_LIBRARIES})
    TARGET_LINK_LIBRARIES(liggghts_bin PUBLIC ${VTK_LIBRARIES})

    ADD_DEFINITIONS(-DLAMMPS_VTK)
    SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} VTK")

    MESSAGE(STATUS "Found VTK: version " ${VTK_VERSION})
  ELSE()
    MESSAGE(FATAL_ERROR "VTK NOT found!")
  ENDIF()
ELSE()
  SET(DISABLED_OPTIONS "${DISABLED_OPTIONS} VTK")
ENDIF()

#=======================================
IF(ENABLE_JPEG)
  FIND_PACKAGE(JPEG)

  IF(JPEG_FOUND)
    INCLUDE_DIRECTORIES(${JPEG_INCLUDE_DIR})
    ADD_DEFINITIONS(-DLAMMPS_JPEG)
    TARGET_LINK_LIBRARIES(liggghts_static PUBLIC jpeg)
    TARGET_LINK_LIBRARIES(liggghts_shared PUBLIC jpeg)
    TARGET_LINK_LIBRARIES(liggghts_bin PUBLIC jpeg)
    SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} JPEG")
  ELSE()
    MESSAGE(FATAL_ERROR "JPEG NOT found!")
  ENDIF()
ELSE()
  SET(DISABLED_OPTIONS "${DISABLED_OPTIONS} JPEG")
ENDIF()

#=======================================
IF(ENABLE_PNG)
  FIND_PACKAGE(PNG)

  IF(PNG_FOUND)
    INCLUDE_DIRECTORIES(${PNG_INCLUDE_DIR})
    ADD_DEFINITIONS(-DLAMMPS_PNG)
    TARGET_LINK_LIBRARIES(liggghts_static PUBLIC ${PNG_LIBRARY})
    TARGET_LINK_LIBRARIES(liggghts_shared PUBLIC ${PNG_LIBRARY})
    TARGET_LINK_LIBRARIES(liggghts_bin PUBLIC ${PNG_LIBRARY})
    SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} PNG")
  ELSE()
    MESSAGE(FATAL_ERROR "PNG NOT found!")
  ENDIF()
ELSE()
  SET(DISABLED_OPTIONS "${DISABLED_OPTIONS} PNG")
ENDIF()

#=======================================
IF(ENABLE_FFMPEG)
  FIND_PATH(AVCODEC_INCLUDE_DIR libavcodec/avcodec.h)
  FIND_LIBRARY(AVCODEC_LIBRARY avcodec)

  FIND_PATH(AVFORMAT_INCLUDE_DIR libavformat/avformat.h)
  FIND_LIBRARY(AVFORMAT_LIBRARY avformat)

  FIND_PATH(AVUTIL_INCLUDE_DIR libavutil/avutil.h)
  FIND_LIBRARY(AVUTIL_LIBRARY avutil)

  IF((NOT AVCODEC_INCLUDE_DIR-NOTFOUND)
    AND (NOT AVFORMAT_INCLUDE_DIR-NOTFOUND)
    AND (NOT AVUTIL_INCLUDE_DIR-NOTFOUND)
    AND (NOT AVCODEC_LIBRARY-NOTFOUND)
    AND (NOT AVFORMAT_LIBRARY-NOTFOUND)
    AND (NOT AVUTIL_LIBRARY-NOTFOUND)
  )
    INCLUDE_DIRECTORIES(${AVCODEC_INCLUDE_DIR} ${AVFORMAT_INCLUDE_DIR} ${AVUTIL_INCLUDE_DIR})
    ADD_DEFINITIONS(-DLAMMPS_FFMPEG)
    TARGET_LINK_LIBRARIES(liggghts_static PUBLIC ${AVCODEC_LIBRARY} ${AVFORMAT_LIBRARY} ${AVUTIL_LIBRARY})
    TARGET_LINK_LIBRARIES(liggghts_shared PUBLIC ${AVCODEC_LIBRARY} ${AVFORMAT_LIBRARY} ${AVUTIL_LIBRARY})
    TARGET_LINK_LIBRARIES(liggghts_bin PUBLIC ${AVCODEC_LIBRARY} ${AVFORMAT_LIBRARY} ${AVUTIL_LIBRARY})
    
    SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} FFMPEG")
    
    MESSAGE(STATUS "FFMPEG found!")
  ELSE()
    MESSAGE(FATAL_ERROR "FFMPEG NOT found!")
  ENDIF()
ELSE()
  SET(DISABLED_OPTIONS "${DISABLED_OPTIONS} FFMPEG")
ENDIF()

#=======================================
IF(ENABLE_GZIP)
  FIND_PROGRAM(GZIP "gzip")

  IF(GZIP)
    ADD_DEFINITIONS(-DLAMMPS_GZIP)

    SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} GZIP")

    MESSAGE(STATUS "GZIP found!")
  ELSE()
    MESSAGE(FATAL_ERROR "GZIP NOT found!")
  ENDIF()
ELSE()
  SET(DISABLED_OPTIONS "${DISABLED_OPTIONS} GZIP")
ENDIF()

#=======================================
IF(ENABLE_SQ)
  FIND_PACKAGE(Boost)

  IF(Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS(-DSUPERQUADRIC_ACTIVE_FLAG -DNONSPHERICAL_ACTIVE_FLAG)
    TARGET_LINK_LIBRARIES(liggghts_static PUBLIC ${Boost_LIBRARY})
    TARGET_LINK_LIBRARIES(liggghts_shared PUBLIC ${Boost_LIBRARY})
    TARGET_LINK_LIBRARIES(liggghts_bin PUBLIC ${Boost_LIBRARY})

    SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} SQ")

    MESSAGE(STATUS "Using Superquadrics")
  ELSE()
    MESSAGE(FATAL_ERROR "Boost NOT found!")
  ENDIF()
ELSE()
  SET(DISABLED_OPTIONS "${DISABLED_OPTIONS} SQ")
ENDIF()

#=======================================
# Parallel HDF5 for dump hdf5 / dump mesh/hdf5 (P0-09, F-26). Same define and
# library as src/MAKE/Makefile.hdf5mpi, but located with FindHDF5 instead of
# hard-coded OpenMPI paths. The dumps use MPI-IO collective writes, so a serial
# HDF5 is rejected here instead of failing obscurely at compile or link time.
IF(LIGGGHTS_ENABLE_HDF5)
  IF(NOT ENABLE_MPI)
    MESSAGE(FATAL_ERROR "LIGGGHTS_ENABLE_HDF5 requires ENABLE_MPI=ON (parallel HDF5)")
  ENDIF()
  IF(NOT DEFINED HDF5_PREFER_PARALLEL)
    SET(HDF5_PREFER_PARALLEL TRUE)   # prefer h5pcc when serial and parallel HDF5 are both installed
  ENDIF()
  FIND_PACKAGE(HDF5 COMPONENTS C)
  IF(NOT HDF5_FOUND)
    MESSAGE(FATAL_ERROR "LIGGGHTS_ENABLE_HDF5=ON but no HDF5 was found. Install a parallel HDF5 "
                        "(Debian/Ubuntu: libhdf5-openmpi-dev) or point CMake to it with "
                        "-DHDF5_C_COMPILER_EXECUTABLE=<path to h5pcc> or -DHDF5_ROOT=<prefix>.")
  ENDIF()
  IF(NOT HDF5_IS_PARALLEL)
    MESSAGE(FATAL_ERROR "LIGGGHTS_ENABLE_HDF5=ON but the HDF5 found (version ${HDF5_VERSION}, "
                        "include ${HDF5_C_INCLUDE_DIRS}) is a SERIAL build. dump hdf5 and dump mesh/hdf5 "
                        "need a parallel (MPI-IO) HDF5 built with the same MPI as LIGGGHTS. "
                        "Select it with -DHDF5_C_COMPILER_EXECUTABLE=<path to h5pcc> or -DHDF5_ROOT=<prefix> "
                        "(Debian/Ubuntu: libhdf5-openmpi-dev, h5pcc), or set LIGGGHTS_ENABLE_HDF5=OFF.")
  ENDIF()
  IF(HDF5_C_INCLUDE_DIRS)
    INCLUDE_DIRECTORIES(${HDF5_C_INCLUDE_DIRS})
  ELSE()
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
  ENDIF()
  ADD_DEFINITIONS(-DLIGGGHTS_HDF5 ${HDF5_C_DEFINITIONS})
  IF(HDF5_C_LIBRARIES)
    SET(LIGGGHTS_HDF5_LIBS ${HDF5_C_LIBRARIES})
  ELSE()
    SET(LIGGGHTS_HDF5_LIBS ${HDF5_LIBRARIES})
  ENDIF()
  TARGET_LINK_LIBRARIES(liggghts_static PUBLIC ${LIGGGHTS_HDF5_LIBS})
  TARGET_LINK_LIBRARIES(liggghts_shared PUBLIC ${LIGGGHTS_HDF5_LIBS})
  TARGET_LINK_LIBRARIES(liggghts_bin PUBLIC ${LIGGGHTS_HDF5_LIBS})
  SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} HDF5")
  MESSAGE(STATUS "Using parallel HDF5 ${HDF5_VERSION}: ${LIGGGHTS_HDF5_LIBS}")
ELSE()
  SET(DISABLED_OPTIONS "${DISABLED_OPTIONS} HDF5")
ENDIF()

#=======================================
IF(ENABLE_SMALLBIG)
  IF(ENABLE_BIGBIG OR ENABLE_SMALLSMALL)
    MESSAGE(FATAL_ERROR "Use only one of these options: SMALLBIG, BIGBIG, SMALLSMALL")
  ENDIF()

  ADD_DEFINITIONS(-DLAMMPS_SMALLBIG)
  SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} SMALLBIG")
ELSEIF(ENABLE_BIGBIG)
  IF(ENABLE_SMALLBIG OR ENABLE_SMALLSMALL)
    MESSAGE(FATAL_ERROR "Use only one of these options: SMALLBIG, BIGBIG, SMALLSMALL")
  ENDIF()

  ADD_DEFINITIONS(-DLAMMPS_BIGBIG)
  SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} BIGBIG")
ELSEIF(ENABLE_SMALLSMALL)
  IF(ENABLE_BIGBIG OR ENABLE_SMALLBIG)
    MESSAGE(FATAL_ERROR "Use only one of these options: SMALLBIG, BIGBIG, SMALLSMALL")
  ENDIF()
  
  ADD_DEFINITIONS(-DLAMMPS_SMALLSMALL)
  SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} SMALLSMALL")
ELSE()
  ADD_DEFINITIONS(-DLAMMPS_SMALLBIG)
  SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} SMALLBIG")
ENDIF()

IF(ENABLE_LONGLONG_TO_LONG)
  ADD_DEFINITIONS(-DLAMMPS_LONGLONG_TO_LONG)
  SET(ENABLED_OPTIONS "${ENABLED_OPTIONS} LONGLONG_TO_LONG")
ENDIF()

#=======================================
GENERATE_VERSION_H()

#=======================================
MESSAGE(STATUS "=======================================")
MESSAGE(STATUS "liggghts ${LIGGGHTS_VERSION}")
MESSAGE(STATUS "Enabled options:${ENABLED_OPTIONS}")
MESSAGE(STATUS "Disabled options:${DISABLED_OPTIONS}")
MESSAGE(STATUS "Enabled models:${ENABLED_MODELS}")
MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}; generated headers: ${LIGGGHTS_GENERATED_DIR}")
MESSAGE(STATUS "=======================================")

#=======================================
# CTest (roadmap A1): the regression suites in ../tests run against the
# liggghts_bin target of this build. See tests/CMakeLists.txt.
IF(LIGGGHTS_ENABLE_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../tests/CMakeLists.txt")
  ENABLE_TESTING()
  SET(LIGGGHTS_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  ADD_SUBDIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/../tests" "${CMAKE_CURRENT_BINARY_DIR}/tests")
ENDIF()

#=======================================
INSTALL(TARGETS liggghts_bin liggghts_static liggghts_shared
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include)
IF(WIN32)
  INSTALL(DIRECTORY WINDOWS DESTINATION include)
ENDIF()
