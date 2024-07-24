cmake_minimum_required (VERSION 3.5)
project (cint C)
set(cint_VERSION_MAJOR "6")
set(cint_VERSION_MINOR "1")
set(cint_VERSION_PATCH "1")
set(cint_VERSION_TWEAK "0")
set(cint_VERSION "${cint_VERSION_MAJOR}.${cint_VERSION_MINOR}.${cint_VERSION_PATCH}")
set(cint_SOVERSION "${cint_VERSION_MAJOR}")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if ("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flto") # --save-temps -fuse-ld=lld")
    message("Clang found")
endif ()


if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
  message("Build type RELWITHDEBINFO")
else()
  message("Build type ${CMAKE_BUILD_TYPE}")
endif()
set(CMAKE_VERBOSE_MAKEFILE OFF)

option(BUILD_MARCH_NATIVE "gcc flag -march=native" off)
if (BUILD_MARCH_NATIVE)
  include(CheckCCompilerFlag)
  CHECK_C_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
      message("Add CFLAGS -march=native -unroll-aggressive")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -unroll-aggressive -ipo")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-math-errno")
    else()
      message("Add CFLAGS -march=native -ftree-vectorize")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -funroll-loops -ftree-vectorize")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-math-errno")
    endif()
  endif()
endif()

cmake_policy(SET CMP0042 NEW)
include(GNUInstallDirs)


## ENZYME

find_package(Enzyme REQUIRED CONFIG NO_DEFAULT_PATH)
message("found dir ${Enzyme_DIR}")
get_property(importTargetsAfter DIRECTORY "${CMAKE_SOURCE_DIR}" PROPERTY IMPORTED_TARGETS)
message("imported targets ${importTargetsAfter}")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(QUADMATH)
if(QUADMATH_FOUND)
  set(HAVE_QUADMATH_H 1 CACHE INTERNAL "Have QUADMATH")
endif()

if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.12)
  cmake_policy(SET CMP0075 NEW)
endif()

include(CheckFunctionExists)
set(CMAKE_REQUIRED_INCLUDES math.h)
set(CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(expl HAVE_EXPL)
check_function_exists(sqrtl HAVE_SQRTL)

include_directories(
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_BINARY_DIR}/src)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_BINARY_DIR}/src/config.h")

configure_file(
  "${PROJECT_SOURCE_DIR}/include/cint.h.in"
  "${PROJECT_BINARY_DIR}/include/cint.h")

set(cintSrc
  src/c2f.c src/cart2sph.c src/cint1e.c src/cint2e.c src/cint_bas.c
  src/fblas.c src/g1e.c src/g2e.c src/misc.c src/optimizer.c
  src/fmt.c src/rys_wheeler.c src/eigh.c src/rys_roots.c src/find_roots.c
  src/cint2c2e.c src/g2c2e.c src/cint3c2e.c src/g3c2e.c
  src/cint3c1e.c src/g3c1e.c src/breit.c
  src/cint1e_a.c src/cint3c1e_a.c
  src/cint1e_grids.c src/g1e_grids.c
  src/autocode/breit1.c src/autocode/dkb.c src/autocode/gaunt1.c
  src/autocode/grad1.c src/autocode/grad2.c src/autocode/hess.c
  src/autocode/int3c1e.c src/autocode/int3c2e.c src/autocode/intor1.c
  src/autocode/intor2.c src/autocode/intor3.c src/autocode/intor4.c
  src/autocode/deriv3.c src/autocode/int1e_grids1.c
  src/autocode/deriv4.c src/autocode/lresc.c)

option(BUILD_SHARED_LIBS "build shared libraries" 1)
option(ENABLE_EXAMPLE "build examples" 0)
option(ENABLE_TEST "build tests" 0)
option(ENABLE_STATIC "Enforce static library build" 1)

if(ENABLE_STATIC)
  set(BUILD_SHARED_LIBS 0)
endif()

add_library(cint ${cintSrc})
target_link_libraries(cint -lm -flto) # LLDEnzymeFlags LLDEnzymeLooseTypeFlags)

option(SO "build so" 0)
option(GRAD "build grad" 0)
option(SCF "build scf" 0)

add_subdirectory(diagonal)

if (SO)
  set(gradSrc
    scf/scf.c)
  add_library(grad SHARED ${gradSrc} ${cintSrc})
  target_link_libraries(grad -lm -flto diag LLDEnzymeFlags LLDEnzymeLooseTypeFlags)
else()
  if (GRAD)
    add_subdirectory(grad)
  endif()

  if (SCF)
    add_subdirectory(scf)
  endif()
endif()

set(CintHeaders
  ${PROJECT_SOURCE_DIR}/include/cint_funcs.h
  ${PROJECT_BINARY_DIR}/include/cint.h)

install(TARGETS cint DESTINATION "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" COMPONENT "lib")
install(FILES ${CintHeaders} DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR} COMPONENT "dev")

include(CPackComponent)

set(CPACK_PACKAGE_CHECKSUM "SHA512")
set(CPACK_GENERATOR "DEB")
set(CPACK_DEBIAN_COMPRESSION_TYPE "xz")
set(CPACK_PACKAGE_VENDOR "Qiming @sunqm Sun")
set(CPACK_PACKAGE_CONTACT "https://github.com/sunqm/libcint")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/sunqm/libcint")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.rst")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/doc/libcint.CPack.txt")

set(CPACK_PACKAGE_VERSION_MAJOR "${cint_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${cint_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${cint_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION_TWEAK "${cint_VERSION_TWEAK}")

set(CPACK_COMPONENTS_GROUPING "IGNORE")

set(CPACK_MONOLITHIC_INSTALL OFF)
set(CPACK_DEB_COMPONENT_INSTALL ON)
set(CPACK_RPM_COMPONENT_INSTALL ON)
set(CPACK_NSIS_COMPONENT_INSTALL ON)

set("CPACK_DEBIAN_LIB_PACKAGE_NAME" "libcint")
set("CPACK_DEBIAN_DEV_PACKAGE_NAME" "${CPACK_DEBIAN_LIB_PACKAGE_NAME}-dev")

set(CPACK_DEBIAN_LIB_PACKAGE_DEPENDS "libopenblas0-pthread")
set(CPACK_DEBIAN_DEV_PACKAGE_DEPENDS "${CPACK_DEBIAN_LIB_PACKAGE_NAME},libopenblas-pthread-dev")
if(HAVE_QUADMATH)
  set(CPACK_DEBIAN_LIB_PACKAGE_DEPENDS "${CPACK_DEBIAN_LIB_PACKAGE_DEPENDS},libquadmath0")
endif()

set(CPACK_COMPONENTS_GROUPING "IGNORE")