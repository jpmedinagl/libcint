cmake_minimum_required (VERSION 3.5)
project (library C)


set(src
  src/file.c
  src/file.h)

set(scfSrc
  grad/scf.c)

include_directories(src)

add_library(scf SHARED ${src} ${scfSrc})

#add_library(scf SHARED ${src})
#add_executable(sub ${subSrc})

#target_link_libraries(sub scf)