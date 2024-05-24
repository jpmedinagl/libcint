cmake_minimum_required(VERSION 3.13.0)

add_library(diag "")

target_sources(diag PUBLIC blaswrap.h f2c.h dlamch.c  dlansy.c  dlascl.c  dorgtr.c  dscal.c  dsteqr.c  dsterf.c  dsyev.c  dsytrd.c  ilaenv.c  lsame.c  xerbla.c)

target_include_directories(diag PUBLIC ${CMAKE_CURRENT_LIST_DIR})
