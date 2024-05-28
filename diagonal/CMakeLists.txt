cmake_minimum_required(VERSION 3.13.0)

add_library(diag "")

target_sources(diag PUBLIC 
backspac.c  dcopy.c  dgemv.c   dlaev2.c    dlansy.c  dlarfg.c  dlaset.c  dnrm2.c   dorgqr.c  dsteqr.c  dsymv.c   dsytrd.c   err.c  fmt.h      ieeeck.c  ilaenv.c   lio.h     pow_di.c  sfe.c      wref.c
blaswrap.h  ddot.c   dger.c    dlaisnan.c  dlapy2.c  dlarft.c  dlasr.c   dorg2l.c  dorgtr.c  dsterf.c  dsyr2.c   dtrmm.c    f2c.h  fmtlib.c   iio.c     i_nint.c   lsame.c   rewind.c  sig_die.c  wrtfmt.c
close.c     dfe.c    disnan.c  dlamch.c    dlarfb.c  dlartg.c  dlasrt.c  dorg2r.c  dscal.c   dswap.c   dsyr2k.c  dtrmv.c    fio.h  fp.h       iladlc.c  inquire.c  lwrite.c  s_cmp.c   sysdep1.h  wsfe.c
daxpy.c     dgemm.c  dlae2.c   dlanst.c    dlarf.c   dlascl.c  dlatrd.c  dorgql.c  d_sign.c  dsyev.c   dsytd2.c  endfile.c  fmt.c  i77vers.c  iladlr.c  iparmq.c   open.c    s_copy.c  util.c     xerbla.c
)

target_include_directories(diag PUBLIC ${CMAKE_CURRENT_LIST_DIR})
