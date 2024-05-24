cmake_minimum_required(VERSION 3.13.0)

add_library(diag "")

target_sources(diag PUBLIC 
backspac.c  dgemv.c     dlapy2.c  dlatrd.c  d_sign.c  dsytd2.c   fmt.c      iladlr.c   open.c     util.c
blaswrap.h  dger.c      dlarfb.c  dnrm2.c   dsteqr.c  dsytrd.c   fmt.h      ilaenv.c   pow_di.c   wref.c
close.c     disnan.c    dlarf.c   dorg2l.c  dsterf.c  dtrmm.c    fmtlib.c   i_nint.c   rewind.c   wrtfmt.c
daxpy.c     dlaev2.c    dlarfg.c  dorg2r.c  dswap.c   dtrmv.c    fp.h       inquire.c  s_cmp.c    wsfe.c
dcopy.c     dlaisnan.c  dlarft.c  dorgql.c  dsyev.c   endfile.c  i77vers.c  iparmq.c   s_copy.c   xerbla.c
ddot.c      dlamch.c    dlartg.c  dorgqr.c  dsymv.c   err.c      ieeeck.c   lio.h      sfe.c
dfe.c       dlanst.c    dlascl.c  dorgtr.c  dsyr2.c   f2c.h      iio.c      lsame.c    sig_die.c
dgemm.c     dlansy.c    dlasr.c   dscal.c   dsyr2k.c  fio.h      iladlc.c   lwrite.c   sysdep1.h
)

target_include_directories(diag PUBLIC ${CMAKE_CURRENT_LIST_DIR})
