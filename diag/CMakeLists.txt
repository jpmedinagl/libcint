add_library(diag "")

target_sources(diag PUBLIC blaswrap.h dlamch.c  dlansy.c  dlascl.c  dorgtr.c  dscal.c  dsteqr.c  dsterf.c  dsyev.c  dsytrd.c  ilaenv.c  lsame.c  xerbla.c)

target_include_directories(diag PUBLIC ${CMAKE_CURRENT_LIST_DIR})

