cmake_minimum_required(VERSION 3.13.0)

add_library(diag "")

target_sources(diag PUBLIC blaswrap.h  ddot.c   dger.c    dlaisnan.c  dlansy.c  dlarf.c   dlartg.c  dlatrd.c  dorg2r.c  dorgtr.c  dsteqr.c  dsyev.c  dsyr2k.c  dtrmm.c  fio.h     iladlc.c  i_nint.c  pow_di.c  sysdep1.h
daxpy.c     dgemm.c  disnan.c  dlamch.c    dlapy2.c  dlarfg.c  dlascl.c  dnrm2.c   dorgql.c  dscal.c   dsterf.c  dsymv.c  dsytd2.c  dtrmv.c  fmt.h     iladlr.c  iparmq.c  s_cmp.c   wsfe.c
dcopy.c     dgemv.c  dlaev2.c  dlanst.c    dlarfb.c  dlarft.c  dlasr.c   dorg2l.c  dorgqr.c  d_sign.c  dswap.c   dsyr2.c  dsytrd.c  f2c.h    ieeeck.c  ilaenv.c  lsame.c   s_copy.c  xerbla.c
)

target_include_directories(diag PUBLIC ${CMAKE_CURRENT_LIST_DIR})
