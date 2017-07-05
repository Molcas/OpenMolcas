!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2015, Ignacio Fdez. Galvan                             *
!***********************************************************************
module link_blas
  use iso_c_binding
  use iso_c_utilities
  use dlfcn
  use blas_mod, &
      int_dasum=>dasum, &
      int_daxpy=>daxpy, &
      int_dcabs1=>dcabs1, &
      int_dcopy=>dcopy, &
      int_ddot=>ddot, &
      int_dgemm=>dgemm, &
      int_dgemv=>dgemv, &
      int_dger=>dger, &
      int_dnrm2=>dnrm2, &
      int_drot=>drot, &
      int_dscal=>dscal, &
      int_dspmv=>dspmv, &
      int_dspr2=>dspr2, &
      int_dspr=>dspr, &
      int_dswap=>dswap, &
      int_dsymm=>dsymm, &
      int_dsymv=>dsymv, &
      int_dsyr2=>dsyr2, &
      int_dsyr2k=>dsyr2k, &
      int_dsyrk=>dsyrk, &
      int_dtpmv=>dtpmv, &
      int_dtpsv=>dtpsv, &
      int_dtrmm=>dtrmm, &
      int_dtrmv=>dtrmv, &
      int_dtrsm=>dtrsm, &
      int_dtrsv=>dtrsv, &
      int_dznrm2=>dznrm2, &
      int_idamax=>idamax, &
      int_lsame=>lsame, &
      int_scopy=>scopy, &
      int_xerbla=>xerbla, &
      int_zaxpy=>zaxpy, &
      int_zcopy=>zcopy, &
      int_zdotc=>zdotc, &
      int_zdscal=>zdscal, &
      int_zgemm=>zgemm, &
      int_zgemv=>zgemv, &
      int_zgerc=>zgerc, &
      int_zhemv=>zhemv, &
      int_zher2=>zher2, &
      int_zher2k=>zher2k, &
      int_zhpmv=>zhpmv, &
      int_zhpr2=>zhpr2, &
      int_zscal=>zscal, &
      int_zswap=>zswap, &
      int_ztrmm=>ztrmm, &
      int_ztrmv=>ztrmv
  use lapack_mod, &
      int_dbdsqr=>dbdsqr, &
      int_dgebd2=>dgebd2, &
      int_dgebrd=>dgebrd, &
      int_dgelq2=>dgelq2, &
      int_dgelqf=>dgelqf, &
      int_dgels=>dgels, &
      int_dgeqr2=>dgeqr2, &
      int_dgeqrf=>dgeqrf, &
      int_dgesvd=>dgesvd, &
      int_dgesv=>dgesv, &
      int_dgetrf=>dgetrf, &
      int_dgetrf2=>dgetrf2, &
      int_dgetri=>dgetri, &
      int_dgetrs=>dgetrs, &
      int_disnan=>disnan, &
      int_dlabad=>dlabad, &
      int_dlabrd=>dlabrd, &
      int_dlacpy=>dlacpy, &
      int_dladiv=>dladiv, &
      int_dladiv1=>dladiv1, &
      int_dladiv2=>dladiv2, &
      int_dlae2=>dlae2, &
      int_dlaebz=>dlaebz, &
      int_dlaev2=>dlaev2, &
      int_dlagtf=>dlagtf, &
      int_dlagts=>dlagts, &
      int_dlaisnan=>dlaisnan, &
      int_dlamch=>dlamch, &
      int_dlaneg=>dlaneg, &
      int_dlange=>dlange, &
      int_dlansp=>dlansp, &
      int_dlanst=>dlanst, &
      int_dlansy=>dlansy, &
      int_dlapy2=>dlapy2, &
      int_dlapy3=>dlapy3, &
      int_dlar1v=>dlar1v, &
      int_dlarfb=>dlarfb, &
      int_dlarf=>dlarf, &
      int_dlarfg=>dlarfg, &
      int_dlarft=>dlarft, &
      int_dlarnv=>dlarnv, &
      int_dlarra=>dlarra, &
      int_dlarrb=>dlarrb, &
      int_dlarrc=>dlarrc, &
      int_dlarrd=>dlarrd, &
      int_dlarre=>dlarre, &
      int_dlarrf=>dlarrf, &
      int_dlarrj=>dlarrj, &
      int_dlarrk=>dlarrk, &
      int_dlarrr=>dlarrr, &
      int_dlarrv=>dlarrv, &
      int_dlartg=>dlartg, &
      int_dlaruv=>dlaruv, &
      int_dlas2=>dlas2, &
      int_dlascl=>dlascl, &
      int_dlaset=>dlaset, &
      int_dlasq1=>dlasq1, &
      int_dlasq2=>dlasq2, &
      int_dlasq3=>dlasq3, &
      int_dlasq4=>dlasq4, &
      int_dlasq5=>dlasq5, &
      int_dlasq6=>dlasq6, &
      int_dlasr=>dlasr, &
      int_dlasrt=>dlasrt, &
      int_dlassq=>dlassq, &
      int_dlasv2=>dlasv2, &
      int_dlaswp=>dlaswp, &
      int_dlatrd=>dlatrd, &
      int_dopgtr=>dopgtr, &
      int_dopmtr=>dopmtr, &
      int_dorg2l=>dorg2l, &
      int_dorg2r=>dorg2r, &
      int_dorgbr=>dorgbr, &
      int_dorgl2=>dorgl2, &
      int_dorglq=>dorglq, &
      int_dorgql=>dorgql, &
      int_dorgqr=>dorgqr, &
      int_dorgtr=>dorgtr, &
      int_dorm2l=>dorm2l, &
      int_dorm2r=>dorm2r, &
      int_dormbr=>dormbr, &
      int_dorml2=>dorml2, &
      int_dormlq=>dormlq, &
      int_dormql=>dormql, &
      int_dormqr=>dormqr, &
      int_dormtr=>dormtr, &
      int_dposv=>dposv, &
      int_dpotrf=>dpotrf, &
      int_dpotrf2=>dpotrf2, &
      int_dpotrs=>dpotrs, &
      int_dpptrf=>dpptrf, &
      int_dspev=>dspev, &
      int_dspgst=>dspgst, &
      int_dspgv=>dspgv, &
      int_dsptrd=>dsptrd, &
      int_dstebz=>dstebz, &
      int_dstein=>dstein, &
      int_dstemr=>dstemr, &
      int_dsteqr=>dsteqr, &
      int_dsterf=>dsterf, &
      int_dstevr=>dstevr, &
      int_dsyev=>dsyev, &
      int_dsyevr=>dsyevr, &
      int_dsygs2=>dsygs2, &
      int_dsygst=>dsygst, &
      int_dsygv=>dsygv, &
      int_dsytd2=>dsytd2, &
      int_dsytrd=>dsytrd, &
      int_dtrti2=>dtrti2, &
      int_dtrtri=>dtrtri, &
      int_dtrtrs=>dtrtrs, &
      int_ieeeck=>ieeeck, &
      int_iladlc=>iladlc, &
      int_iladlr=>iladlr, &
      int_ilaenv=>ilaenv, &
      int_ilazlc=>ilazlc, &
      int_ilazlr=>ilazlr, &
      int_iparmq=>iparmq, &
      int_iparam2stage=>iparam2stage, &
      int_zheev=>zheev, &
      int_zhetd2=>zhetd2, &
      int_zhetrd=>zhetrd, &
      int_zhpev=>zhpev, &
      int_zhptrd=>zhptrd, &
      int_zlacgv=>zlacgv, &
      int_zladiv=>zladiv, &
      int_zlanhe=>zlanhe, &
      int_zlanhp=>zlanhp, &
      int_zlarf=>zlarf, &
      int_zlarfb=>zlarfb, &
      int_zlarfg=>zlarfg, &
      int_zlarft=>zlarft, &
      int_zlascl=>zlascl, &
      int_zlaset=>zlaset, &
      int_zlasr=>zlasr, &
      int_zlassq=>zlassq, &
      int_zlatrd=>zlatrd, &
      int_zsteqr=>zsteqr, &
      int_zung2l=>zung2l, &
      int_zung2r=>zung2r, &
      int_zungql=>zungql, &
      int_zungqr=>zungqr, &
      int_zungtr=>zungtr, &
      int_zupgtr=>zupgtr
  use legacy_mod, &
      int_dgetf2=>dgetf2, &
      int_dpotf2=>dpotf2

  implicit none

  type(c_ptr), dimension(:), allocatable, private :: handles
!
! Initializing procedure pointers is a F2008 feature, not supported by all compilers.
! When it is implemented, thhe pointers below should look like:
!
!   procedure(int_dasum), pointer :: lb_dasum=>int_dasum
!
! and then the exact placement of the initialization call (in start.f) is not critical
!
! BLAS procedures
  procedure(int_dasum), pointer :: lb_dasum
  procedure(int_daxpy), pointer :: lb_daxpy
  procedure(int_dcabs1), pointer :: lb_dcabs1
  procedure(int_dcopy), pointer :: lb_dcopy
  procedure(int_ddot), pointer :: lb_ddot
  procedure(int_dgemm), pointer :: lb_dgemm
  procedure(int_dgemv), pointer :: lb_dgemv
  procedure(int_dger), pointer :: lb_dger
  procedure(int_dnrm2), pointer :: lb_dnrm2
  procedure(int_drot), pointer :: lb_drot
  procedure(int_dscal), pointer :: lb_dscal
  procedure(int_dspmv), pointer :: lb_dspmv
  procedure(int_dspr2), pointer :: lb_dspr2
  procedure(int_dspr), pointer :: lb_dspr
  procedure(int_dswap), pointer :: lb_dswap
  procedure(int_dsymm), pointer :: lb_dsymm
  procedure(int_dsymv), pointer :: lb_dsymv
  procedure(int_dsyr2), pointer :: lb_dsyr2
  procedure(int_dsyr2k), pointer :: lb_dsyr2k
  procedure(int_dsyrk), pointer :: lb_dsyrk
  procedure(int_dtpmv), pointer :: lb_dtpmv
  procedure(int_dtpsv), pointer :: lb_dtpsv
  procedure(int_dtrmm), pointer :: lb_dtrmm
  procedure(int_dtrmv), pointer :: lb_dtrmv
  procedure(int_dtrsm), pointer :: lb_dtrsm
  procedure(int_dtrsv), pointer :: lb_dtrsv
  procedure(int_dznrm2), pointer :: lb_dznrm2
  procedure(int_idamax), pointer :: lb_idamax
  procedure(int_lsame), pointer :: lb_lsame
  procedure(int_scopy), pointer :: lb_scopy
  procedure(int_xerbla), pointer :: lb_xerbla
  procedure(int_zaxpy), pointer :: lb_zaxpy
  procedure(int_zcopy), pointer :: lb_zcopy
  procedure(int_zdotc), pointer :: lb_zdotc
  procedure(int_zdscal), pointer :: lb_zdscal
  procedure(int_zgemm), pointer :: lb_zgemm
  procedure(int_zgemv), pointer :: lb_zgemv
  procedure(int_zgerc), pointer :: lb_zgerc
  procedure(int_zhemv), pointer :: lb_zhemv
  procedure(int_zher2), pointer :: lb_zher2
  procedure(int_zher2k), pointer :: lb_zher2k
  procedure(int_zhpmv), pointer :: lb_zhpmv
  procedure(int_zhpr2), pointer :: lb_zhpr2
  procedure(int_zscal), pointer :: lb_zscal
  procedure(int_zswap), pointer :: lb_zswap
  procedure(int_ztrmm), pointer :: lb_ztrmm
  procedure(int_ztrmv), pointer :: lb_ztrmv
! LAPACK procedures
  procedure(int_dbdsqr), pointer :: lb_dbdsqr
  procedure(int_dgebd2), pointer :: lb_dgebd2
  procedure(int_dgebrd), pointer :: lb_dgebrd
  procedure(int_dgelq2), pointer :: lb_dgelq2
  procedure(int_dgelqf), pointer :: lb_dgelqf
  procedure(int_dgels), pointer :: lb_dgels
  procedure(int_dgeqr2), pointer :: lb_dgeqr2
  procedure(int_dgeqrf), pointer :: lb_dgeqrf
  procedure(int_dgesvd), pointer :: lb_dgesvd
  procedure(int_dgesv), pointer :: lb_dgesv
  procedure(int_dgetrf), pointer :: lb_dgetrf
  procedure(int_dgetrf2), pointer :: lb_dgetrf2
  procedure(int_dgetri), pointer :: lb_dgetri
  procedure(int_dgetrs), pointer :: lb_dgetrs
  procedure(int_disnan), pointer :: lb_disnan
  procedure(int_dlabad), pointer :: lb_dlabad
  procedure(int_dlabrd), pointer :: lb_dlabrd
  procedure(int_dlacpy), pointer :: lb_dlacpy
  procedure(int_dladiv), pointer :: lb_dladiv
  procedure(int_dladiv1), pointer :: lb_dladiv1
  procedure(int_dladiv2), pointer :: lb_dladiv2
  procedure(int_dlae2), pointer :: lb_dlae2
  procedure(int_dlaebz), pointer :: lb_dlaebz
  procedure(int_dlaev2), pointer :: lb_dlaev2
  procedure(int_dlagtf), pointer :: lb_dlagtf
  procedure(int_dlagts), pointer :: lb_dlagts
  procedure(int_dlaisnan), pointer :: lb_dlaisnan
  procedure(int_dlamch), pointer :: lb_dlamch
  procedure(int_dlaneg), pointer :: lb_dlaneg
  procedure(int_dlange), pointer :: lb_dlange
  procedure(int_dlansp), pointer :: lb_dlansp
  procedure(int_dlanst), pointer :: lb_dlanst
  procedure(int_dlansy), pointer :: lb_dlansy
  procedure(int_dlapy2), pointer :: lb_dlapy2
  procedure(int_dlapy3), pointer :: lb_dlapy3
  procedure(int_dlar1v), pointer :: lb_dlar1v
  procedure(int_dlarfb), pointer :: lb_dlarfb
  procedure(int_dlarf), pointer :: lb_dlarf
  procedure(int_dlarfg), pointer :: lb_dlarfg
  procedure(int_dlarft), pointer :: lb_dlarft
  procedure(int_dlarnv), pointer :: lb_dlarnv
  procedure(int_dlarra), pointer :: lb_dlarra
  procedure(int_dlarrb), pointer :: lb_dlarrb
  procedure(int_dlarrc), pointer :: lb_dlarrc
  procedure(int_dlarrd), pointer :: lb_dlarrd
  procedure(int_dlarre), pointer :: lb_dlarre
  procedure(int_dlarrf), pointer :: lb_dlarrf
  procedure(int_dlarrj), pointer :: lb_dlarrj
  procedure(int_dlarrk), pointer :: lb_dlarrk
  procedure(int_dlarrr), pointer :: lb_dlarrr
  procedure(int_dlarrv), pointer :: lb_dlarrv
  procedure(int_dlartg), pointer :: lb_dlartg
  procedure(int_dlaruv), pointer :: lb_dlaruv
  procedure(int_dlas2), pointer :: lb_dlas2
  procedure(int_dlascl), pointer :: lb_dlascl
  procedure(int_dlaset), pointer :: lb_dlaset
  procedure(int_dlasq1), pointer :: lb_dlasq1
  procedure(int_dlasq2), pointer :: lb_dlasq2
  procedure(int_dlasq3), pointer :: lb_dlasq3
  procedure(int_dlasq4), pointer :: lb_dlasq4
  procedure(int_dlasq5), pointer :: lb_dlasq5
  procedure(int_dlasq6), pointer :: lb_dlasq6
  procedure(int_dlasr), pointer :: lb_dlasr
  procedure(int_dlasrt), pointer :: lb_dlasrt
  procedure(int_dlassq), pointer :: lb_dlassq
  procedure(int_dlasv2), pointer :: lb_dlasv2
  procedure(int_dlaswp), pointer :: lb_dlaswp
  procedure(int_dlatrd), pointer :: lb_dlatrd
  procedure(int_dopgtr), pointer :: lb_dopgtr
  procedure(int_dopmtr), pointer :: lb_dopmtr
  procedure(int_dorg2l), pointer :: lb_dorg2l
  procedure(int_dorg2r), pointer :: lb_dorg2r
  procedure(int_dorgbr), pointer :: lb_dorgbr
  procedure(int_dorgl2), pointer :: lb_dorgl2
  procedure(int_dorglq), pointer :: lb_dorglq
  procedure(int_dorgql), pointer :: lb_dorgql
  procedure(int_dorgqr), pointer :: lb_dorgqr
  procedure(int_dorgtr), pointer :: lb_dorgtr
  procedure(int_dorm2l), pointer :: lb_dorm2l
  procedure(int_dorm2r), pointer :: lb_dorm2r
  procedure(int_dormbr), pointer :: lb_dormbr
  procedure(int_dorml2), pointer :: lb_dorml2
  procedure(int_dormlq), pointer :: lb_dormlq
  procedure(int_dormql), pointer :: lb_dormql
  procedure(int_dormqr), pointer :: lb_dormqr
  procedure(int_dormtr), pointer :: lb_dormtr
  procedure(int_dposv), pointer :: lb_dposv
  procedure(int_dpotrf), pointer :: lb_dpotrf
  procedure(int_dpotrf2), pointer :: lb_dpotrf2
  procedure(int_dpotrs), pointer :: lb_dpotrs
  procedure(int_dpptrf), pointer :: lb_dpptrf
  procedure(int_dspev), pointer :: lb_dspev
  procedure(int_dspgst), pointer :: lb_dspgst
  procedure(int_dspgv), pointer :: lb_dspgv
  procedure(int_dsptrd), pointer :: lb_dsptrd
  procedure(int_dstebz), pointer :: lb_dstebz
  procedure(int_dstein), pointer :: lb_dstein
  procedure(int_dstemr), pointer :: lb_dstemr
  procedure(int_dsteqr), pointer :: lb_dsteqr
  procedure(int_dsterf), pointer :: lb_dsterf
  procedure(int_dstevr), pointer :: lb_dstevr
  procedure(int_dsyev), pointer :: lb_dsyev
  procedure(int_dsyevr), pointer :: lb_dsyevr
  procedure(int_dsygs2), pointer :: lb_dsygs2
  procedure(int_dsygst), pointer :: lb_dsygst
  procedure(int_dsygv), pointer :: lb_dsygv
  procedure(int_dsytd2), pointer :: lb_dsytd2
  procedure(int_dsytrd), pointer :: lb_dsytrd
  procedure(int_dtrti2), pointer :: lb_dtrti2
  procedure(int_dtrtri), pointer :: lb_dtrtri
  procedure(int_dtrtrs), pointer :: lb_dtrtrs
  procedure(int_ieeeck), pointer :: lb_ieeeck
  procedure(int_iladlc), pointer :: lb_iladlc
  procedure(int_iladlr), pointer :: lb_iladlr
  procedure(int_ilaenv), pointer :: lb_ilaenv
  procedure(int_ilazlc), pointer :: lb_ilazlc
  procedure(int_ilazlr), pointer :: lb_ilazlr
  procedure(int_iparmq), pointer :: lb_iparmq
  procedure(int_iparam2stage), pointer :: lb_iparam2stage
  procedure(int_zheev), pointer :: lb_zheev
  procedure(int_zhetd2), pointer :: lb_zhetd2
  procedure(int_zhetrd), pointer :: lb_zhetrd
  procedure(int_zhpev), pointer :: lb_zhpev
  procedure(int_zhptrd), pointer :: lb_zhptrd
  procedure(int_zlacgv), pointer :: lb_zlacgv
  procedure(int_zladiv), pointer :: lb_zladiv
  procedure(int_zlanhe), pointer :: lb_zlanhe
  procedure(int_zlanhp), pointer :: lb_zlanhp
  procedure(int_zlarf), pointer :: lb_zlarf
  procedure(int_zlarfb), pointer :: lb_zlarfb
  procedure(int_zlarfg), pointer :: lb_zlarfg
  procedure(int_zlarft), pointer :: lb_zlarft
  procedure(int_zlascl), pointer :: lb_zlascl
  procedure(int_zlaset), pointer :: lb_zlaset
  procedure(int_zlasr), pointer :: lb_zlasr
  procedure(int_zlassq), pointer :: lb_zlassq
  procedure(int_zlatrd), pointer :: lb_zlatrd
  procedure(int_zsteqr), pointer :: lb_zsteqr
  procedure(int_zung2l), pointer :: lb_zung2l
  procedure(int_zung2r), pointer :: lb_zung2r
  procedure(int_zungql), pointer :: lb_zungql
  procedure(int_zungqr), pointer :: lb_zungqr
  procedure(int_zungtr), pointer :: lb_zungtr
  procedure(int_zupgtr), pointer :: lb_zupgtr
! Legacy procedures
  procedure(int_dgetf2), pointer :: lb_dgetf2
  procedure(int_dpotf2), pointer :: lb_dpotf2

contains

!define _DEBUG_

!===============================================================================

  subroutine lb_initialize(lib,prlev)
    character(len=*), intent(in) :: lib
    integer, intent(in) :: prlev
    character(len=1024), dimension(:), allocatable :: libs
    character(kind=c_char,len=1024) :: libname
    type(c_funptr) :: funptr=c_null_funptr
    type(DL_Info), target :: info
    integer :: i
    logical :: try_load,loaded

    try_load = (lib /= 'Internal')

    call lb_close()

!***************************************************
!   Try to load all the libraries specified by "lib"
!***************************************************
    call split_string(lib,libs)

    allocate(handles(size(libs)))
    handles(:)=c_null_ptr

    loaded=.false.
    if (try_load) then
      do i=1,size(libs)
        libname=libs(i)
        if (len_trim(libs(i)) > 0) then
          handles(i)=dlopen(trim(libname)//c_null_char, int(ior(rtld_global,rtld_lazy),kind=c_int))
          loaded=(loaded .or. c_associated(handles(i)))
          if (prlev > 0) then
            if (c_associated(handles(i))) then
              write (6,*) trim(libname),' loaded'
            else
              write(6,*) c_f_string(dlerror())
            end if
          end if
        end if
      end do
    end if

    deallocate(libs)

!********************************************************************
!   Associate all BLAS and LAPACK routines with the library functions
!********************************************************************
    if (loaded) then
!
!     BLAS procedures
!
      funptr=link_func('dasum')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dasum)
      end if
!
      funptr=link_func('daxpy')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_daxpy)
      end if
!
      funptr=link_func('dcabs1')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dcabs1)
      end if
!
      funptr=link_func('dcopy')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dcopy)
      end if
!
      funptr=link_func('ddot')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_ddot)
      end if
!
      funptr=link_func('dgemm')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgemm)
      end if
!
      funptr=link_func('dgemv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgemv)
      end if
!
      funptr=link_func('dger')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dger)
      end if
!
      funptr=link_func('dnrm2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dnrm2)
      end if
!
      funptr=link_func('drot')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_drot)
      end if
!
      funptr=link_func('dscal')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dscal)
      end if
!
      funptr=link_func('dspmv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dspmv)
      end if
!
      funptr=link_func('dspr2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dspr2)
      end if
!
      funptr=link_func('dspr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dspr)
      end if
!
      funptr=link_func('dswap')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dswap)
      end if
!
      funptr=link_func('dsymm')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsymm)
      end if
!
      funptr=link_func('dsymv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsymv)
      end if
!
      funptr=link_func('dsyr2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsyr2)
      end if
!
      funptr=link_func('dsyr2k')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsyr2k)
      end if
!
      funptr=link_func('dsyrk')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsyrk)
      end if
!
      funptr=link_func('dtpmv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dtpmv)
      end if
!
      funptr=link_func('dtpsv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dtpsv)
      end if
!
      funptr=link_func('dtrmm')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dtrmm)
      end if
!
      funptr=link_func('dtrmv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dtrmv)
      end if
!
      funptr=link_func('dtrsm')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dtrsm)
      end if
!
      funptr=link_func('dtrsv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dtrsv)
      end if
!
      funptr=link_func('dznrm2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dznrm2)
      end if
!
      funptr=link_func('idamax')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_idamax)
      end if
!
      funptr=link_func('lsame')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_lsame)
      end if
!
      funptr=link_func('scopy')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_scopy)
      end if
!
      funptr=link_func('xerbla')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_xerbla)
      end if
!
      funptr=link_func('zaxpy')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zaxpy)
      end if
!
      funptr=link_func('zcopy')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zcopy)
      end if
!
      funptr=link_func('zdotc')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zdotc)
      end if
!
      funptr=link_func('zdscal')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zdscal)
      end if
!
      funptr=link_func('zgemm')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zgemm)
      end if
!
      funptr=link_func('zgemv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zgemv)
      end if
!
      funptr=link_func('zgerc')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zgerc)
      end if
!
      funptr=link_func('zhemv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zhemv)
      end if
!
      funptr=link_func('zher2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zher2)
      end if
!
      funptr=link_func('zher2k')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zher2k)
      end if
!
      funptr=link_func('zhpmv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zhpmv)
      end if
!
      funptr=link_func('zhpr2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zhpr2)
      end if
!
      funptr=link_func('zscal')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zscal)
      end if
!
      funptr=link_func('zswap')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zswap)
      end if
!
      funptr=link_func('ztrmm')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_ztrmm)
      end if
!
      funptr=link_func('ztrmv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_ztrmv)
      end if
!
!     LAPACK procedures
!
      funptr=link_func('dbdsqr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dbdsqr)
      end if
!
      funptr=link_func('dgebd2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgebd2)
      end if
!
      funptr=link_func('dgebrd')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgebrd)
      end if
!
      funptr=link_func('dgelq2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgelq2)
      end if
!
      funptr=link_func('dgelqf')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgelqf)
      end if
!
      funptr=link_func('dgels')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgels)
      end if
!
      funptr=link_func('dgeqr2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgeqr2)
      end if
!
      funptr=link_func('dgeqrf')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgeqrf)
      end if
!
      funptr=link_func('dgesvd')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgesvd)
      end if
!
      funptr=link_func('dgesv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgesv)
      end if
!
      funptr=link_func('dgetrf')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgetrf)
      end if
!
      funptr=link_func('dgetrf2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgetrf2)
      end if
!
      funptr=link_func('dgetri')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgetri)
      end if
!
      funptr=link_func('dgetrs')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgetrs)
      end if
!
      funptr=link_func('disnan')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_disnan)
      end if
!
      funptr=link_func('dlabad')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlabad)
      end if
!
      funptr=link_func('dlabrd')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlabrd)
      end if
!
      funptr=link_func('dlacpy')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlacpy)
      end if
!
      funptr=link_func('dladiv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dladiv)
      end if
!
      funptr=link_func('dladiv1')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dladiv1)
      end if
!
      funptr=link_func('dladiv2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dladiv2)
      end if
!
      funptr=link_func('dlae2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlae2)
      end if
!
      funptr=link_func('dlaebz')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlaebz)
      end if
!
      funptr=link_func('dlaev2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlaev2)
      end if
!
      funptr=link_func('dlagtf')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlagtf)
      end if
!
      funptr=link_func('dlagts')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlagts)
      end if
!
      funptr=link_func('dlaisnan')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlaisnan)
      end if
!
      funptr=link_func('dlamch')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlamch)
      end if
!
      funptr=link_func('dlaneg')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlaneg)
      end if
!
      funptr=link_func('dlange')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlange)
      end if
!
      funptr=link_func('dlansp')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlansp)
      end if
!
      funptr=link_func('dlanst')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlanst)
      end if
!
      funptr=link_func('dlansy')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlansy)
      end if
!
      funptr=link_func('dlapy2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlapy2)
      end if
!
      funptr=link_func('dlapy3')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlapy3)
      end if
!
      funptr=link_func('dlar1v')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlar1v)
      end if
!
      funptr=link_func('dlarfb')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarfb)
      end if
!
      funptr=link_func('dlarf')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarf)
      end if
!
      funptr=link_func('dlarfg')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarfg)
      end if
!
      funptr=link_func('dlarft')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarft)
      end if
!
      funptr=link_func('dlarnv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarnv)
      end if
!
      funptr=link_func('dlarra')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarra)
      end if
!
      funptr=link_func('dlarrb')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarrb)
      end if
!
      funptr=link_func('dlarrc')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarrc)
      end if
!
      funptr=link_func('dlarrd')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarrd)
      end if
!
      funptr=link_func('dlarre')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarre)
      end if
!
      funptr=link_func('dlarrf')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarrf)
      end if
!
      funptr=link_func('dlarrj')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarrj)
      end if
!
      funptr=link_func('dlarrk')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarrk)
      end if
!
      funptr=link_func('dlarrr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarrr)
      end if
!
      funptr=link_func('dlarrv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlarrv)
      end if
!
      funptr=link_func('dlartg')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlartg)
      end if
!
      funptr=link_func('dlaruv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlaruv)
      end if
!
      funptr=link_func('dlas2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlas2)
      end if
!
      funptr=link_func('dlascl')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlascl)
      end if
!
      funptr=link_func('dlaset')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlaset)
      end if
!
      funptr=link_func('dlasq1')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlasq1)
      end if
!
      funptr=link_func('dlasq2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlasq2)
      end if
!
      funptr=link_func('dlasq3')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlasq3)
      end if
!
      funptr=link_func('dlasq4')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlasq4)
      end if
!
      funptr=link_func('dlasq5')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlasq5)
      end if
!
      funptr=link_func('dlasq6')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlasq6)
      end if
!
      funptr=link_func('dlasr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlasr)
      end if
!
      funptr=link_func('dlasrt')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlasrt)
      end if
!
      funptr=link_func('dlassq')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlassq)
      end if
!
      funptr=link_func('dlasv2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlasv2)
      end if
!
      funptr=link_func('dlaswp')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlaswp)
      end if
!
      funptr=link_func('dlatrd')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dlatrd)
      end if
!
      funptr=link_func('dopgtr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dopgtr)
      end if
!
      funptr=link_func('dopmtr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dopmtr)
      end if
!
      funptr=link_func('dorg2l')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dorg2l)
      end if
!
      funptr=link_func('dorg2r')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dorg2r)
      end if
!
      funptr=link_func('dorgbr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dorgbr)
      end if
!
      funptr=link_func('dorgl2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dorgl2)
      end if
!
      funptr=link_func('dorglq')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dorglq)
      end if
!
      funptr=link_func('dorgql')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dorgql)
      end if
!
      funptr=link_func('dorgqr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dorgqr)
      end if
!
      funptr=link_func('dorgtr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dorgtr)
      end if
!
      funptr=link_func('dorm2l')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dorm2l)
      end if
!
      funptr=link_func('dorm2r')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dorm2r)
      end if
!
      funptr=link_func('dormbr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dormbr)
      end if
!
      funptr=link_func('dorml2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dorml2)
      end if
!
      funptr=link_func('dormlq')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dormlq)
      end if
!
      funptr=link_func('dormql')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dormql)
      end if
!
      funptr=link_func('dormqr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dormqr)
      end if
!
      funptr=link_func('dormtr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dormtr)
      end if
!
      funptr=link_func('dposv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dposv)
      end if
!
      funptr=link_func('dpotrf')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dpotrf)
      end if
!
      funptr=link_func('dpotrf2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dpotrf2)
      end if
!
      funptr=link_func('dpotrs')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dpotrs)
      end if
!
      funptr=link_func('dpptrf')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dpptrf)
      end if
!
      funptr=link_func('dspev')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dspev)
      end if
!
      funptr=link_func('dspgst')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dspgst)
      end if
!
      funptr=link_func('dspgv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dspgv)
      end if
!
      funptr=link_func('dsptrd')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsptrd)
      end if
!
      funptr=link_func('dstebz')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dstebz)
      end if
!
      funptr=link_func('dstein')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dstein)
      end if
!
      funptr=link_func('dstemr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dstemr)
      end if
!
      funptr=link_func('dsteqr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsteqr)
      end if
!
      funptr=link_func('dsterf')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsterf)
      end if
!
      funptr=link_func('dstevr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dstevr)
      end if
!
      funptr=link_func('dsyev')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsyev)
      end if
!
      funptr=link_func('dsyevr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsyevr)
      end if
!
      funptr=link_func('dsygs2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsygs2)
      end if
!
      funptr=link_func('dsygst')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsygst)
      end if
!
      funptr=link_func('dsygv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsygv)
      end if
!
      funptr=link_func('dsytd2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsytd2)
      end if
!
      funptr=link_func('dsytrd')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dsytrd)
      end if
!
      funptr=link_func('dtrti2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dtrti2)
      end if
!
      funptr=link_func('dtrtri')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dtrtri)
      end if
!
      funptr=link_func('dtrtrs')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dtrtrs)
      end if
!
      funptr=link_func('ieeeck')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_ieeeck)
      end if
!
      funptr=link_func('iladlc')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_iladlc)
      end if
!
      funptr=link_func('iladlr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_iladlr)
      end if
!
      funptr=link_func('ilaenv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_ilaenv)
      end if
!
      funptr=link_func('ilazlc')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_ilazlc)
      end if
!
      funptr=link_func('ilazlr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_ilazlr)
      end if
!
      funptr=link_func('iparmq')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_iparmq)
      end if
!
      funptr=link_func('iparam2stage')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_iparam2stage)
      end if
!
      funptr=link_func('zheev')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zheev)
      end if
!
      funptr=link_func('zhetd2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zhetd2)
      end if
!
      funptr=link_func('zhetrd')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zhetrd)
      end if
!
      funptr=link_func('zhpev')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zhpev)
      end if
!
      funptr=link_func('zhptrd')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zhptrd)
      end if
!
      funptr=link_func('zlacgv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zlacgv)
      end if
!
      funptr=link_func('zladiv')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zladiv)
      end if
!
      funptr=link_func('zlanhe')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zlanhe)
      end if
!
      funptr=link_func('zlanhp')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zlanhp)
      end if
!
      funptr=link_func('zlarf')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zlarf)
      end if
!
      funptr=link_func('zlarfb')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zlarfb)
      end if
!
      funptr=link_func('zlarfg')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zlarfg)
      end if
!
      funptr=link_func('zlarft')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zlarft)
      end if
!
      funptr=link_func('zlascl')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zlascl)
      end if
!
      funptr=link_func('zlaset')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zlaset)
      end if
!
      funptr=link_func('zlasr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zlasr)
      end if
!
      funptr=link_func('zlassq')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zlassq)
      end if
!
      funptr=link_func('zlatrd')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zlatrd)
      end if
!
      funptr=link_func('zsteqr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zsteqr)
      end if
!
      funptr=link_func('zung2l')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zung2l)
      end if
!
      funptr=link_func('zung2r')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zung2r)
      end if
!
      funptr=link_func('zungql')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zungql)
      end if
!
      funptr=link_func('zungqr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zungqr)
      end if
!
      funptr=link_func('zungtr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zungtr)
      end if
!
      funptr=link_func('zupgtr')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_zupgtr)
      end if
!
!     Legacy procedures
!
      funptr=link_func('dgetf2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dgetf2)
      end if
!
      funptr=link_func('dpotf2')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_dpotf2)
      end if
!
!****************************************
!   Or use the fallback internal routines
!****************************************
    else
      if (prlev > 0) then
        write(6,*) 'Using internal BLAS+LAPACK'
      end if
      !
      ! BLAS
      lb_dasum=>int_dasum
      lb_daxpy=>int_daxpy
      lb_dcabs1=>int_dcabs1
      lb_dcopy=>int_dcopy
      lb_ddot=>int_ddot
      lb_dgemm=>int_dgemm
      lb_dgemv=>int_dgemv
      lb_dger=>int_dger
      lb_dnrm2=>int_dnrm2
      lb_drot=>int_drot
      lb_dscal=>int_dscal
      lb_dspmv=>int_dspmv
      lb_dspr2=>int_dspr2
      lb_dspr=>int_dspr
      lb_dswap=>int_dswap
      lb_dsymm=>int_dsymm
      lb_dsymv=>int_dsymv
      lb_dsyr2=>int_dsyr2
      lb_dsyr2k=>int_dsyr2k
      lb_dsyrk=>int_dsyrk
      lb_dtpmv=>int_dtpmv
      lb_dtpsv=>int_dtpsv
      lb_dtrmm=>int_dtrmm
      lb_dtrmv=>int_dtrmv
      lb_dtrsm=>int_dtrsm
      lb_dtrsv=>int_dtrsv
      lb_dznrm2=>int_dznrm2
      lb_idamax=>int_idamax
      lb_lsame=>int_lsame
      lb_xerbla=>int_xerbla
      lb_scopy=>int_scopy
      lb_zaxpy=>int_zaxpy
      lb_zcopy=>int_zcopy
      lb_zdotc=>int_zdotc
      lb_zdscal=>int_zdscal
      lb_zgemm=>int_zgemm
      lb_zgemv=>int_zgemv
      lb_zgerc=>int_zgerc
      lb_zhemv=>int_zhemv
      lb_zher2=>int_zher2
      lb_zher2k=>int_zher2k
      lb_zhpmv=>int_zhpmv
      lb_zhpr2=>int_zhpr2
      lb_zscal=>int_zscal
      lb_zswap=>int_zswap
      lb_ztrmm=>int_ztrmm
      lb_ztrmv=>int_ztrmv
      !
      ! LAPACK
      lb_dbdsqr=>int_dbdsqr
      lb_dgebd2=>int_dgebd2
      lb_dgebrd=>int_dgebrd
      lb_dgelq2=>int_dgelq2
      lb_dgelqf=>int_dgelqf
      lb_dgels=>int_dgels
      lb_dgeqr2=>int_dgeqr2
      lb_dgeqrf=>int_dgeqrf
      lb_dgesvd=>int_dgesvd
      lb_dgesv=>int_dgesv
      lb_dgetrf=>int_dgetrf
      lb_dgetrf2=>int_dgetrf2
      lb_dgetri=>int_dgetri
      lb_dgetrs=>int_dgetrs
      lb_disnan=>int_disnan
      lb_dlabad=>int_dlabad
      lb_dlabrd=>int_dlabrd
      lb_dlacpy=>int_dlacpy
      lb_dladiv=>int_dladiv
      lb_dladiv1=>int_dladiv1
      lb_dladiv2=>int_dladiv2
      lb_dlae2=>int_dlae2
      lb_dlaebz=>int_dlaebz
      lb_dlaev2=>int_dlaev2
      lb_dlagtf=>int_dlagtf
      lb_dlagts=>int_dlagts
      lb_dlaisnan=>int_dlaisnan
      lb_dlamch=>int_dlamch
      lb_dlaneg=>int_dlaneg
      lb_dlange=>int_dlange
      lb_dlansp=>int_dlansp
      lb_dlanst=>int_dlanst
      lb_dlansy=>int_dlansy
      lb_dlapy2=>int_dlapy2
      lb_dlapy3=>int_dlapy3
      lb_dlar1v=>int_dlar1v
      lb_dlarfb=>int_dlarfb
      lb_dlarf=>int_dlarf
      lb_dlarfg=>int_dlarfg
      lb_dlarft=>int_dlarft
      lb_dlarnv=>int_dlarnv
      lb_dlarra=>int_dlarra
      lb_dlarrb=>int_dlarrb
      lb_dlarrc=>int_dlarrc
      lb_dlarrd=>int_dlarrd
      lb_dlarre=>int_dlarre
      lb_dlarrf=>int_dlarrf
      lb_dlarrj=>int_dlarrj
      lb_dlarrk=>int_dlarrk
      lb_dlarrr=>int_dlarrr
      lb_dlarrv=>int_dlarrv
      lb_dlartg=>int_dlartg
      lb_dlaruv=>int_dlaruv
      lb_dlas2=>int_dlas2
      lb_dlascl=>int_dlascl
      lb_dlaset=>int_dlaset
      lb_dlasq1=>int_dlasq1
      lb_dlasq2=>int_dlasq2
      lb_dlasq3=>int_dlasq3
      lb_dlasq4=>int_dlasq4
      lb_dlasq5=>int_dlasq5
      lb_dlasq6=>int_dlasq6
      lb_dlasr=>int_dlasr
      lb_dlasrt=>int_dlasrt
      lb_dlassq=>int_dlassq
      lb_dlasv2=>int_dlasv2
      lb_dlaswp=>int_dlaswp
      lb_dlatrd=>int_dlatrd
      lb_dopgtr=>int_dopgtr
      lb_dopmtr=>int_dopmtr
      lb_dorg2l=>int_dorg2l
      lb_dorg2r=>int_dorg2r
      lb_dorgbr=>int_dorgbr
      lb_dorgl2=>int_dorgl2
      lb_dorglq=>int_dorglq
      lb_dorgql=>int_dorgql
      lb_dorgqr=>int_dorgqr
      lb_dorgtr=>int_dorgtr
      lb_dorm2l=>int_dorm2l
      lb_dorm2r=>int_dorm2r
      lb_dormbr=>int_dormbr
      lb_dorml2=>int_dorml2
      lb_dormlq=>int_dormlq
      lb_dormql=>int_dormql
      lb_dormqr=>int_dormqr
      lb_dormtr=>int_dormtr
      lb_dposv=>int_dposv
      lb_dpotrf=>int_dpotrf
      lb_dpotrf2=>int_dpotrf2
      lb_dpotrs=>int_dpotrs
      lb_dpptrf=>int_dpptrf
      lb_dspev=>int_dspev
      lb_dspgst=>int_dspgst
      lb_dspgv=>int_dspgv
      lb_dsptrd=>int_dsptrd
      lb_dstebz=>int_dstebz
      lb_dstein=>int_dstein
      lb_dstemr=>int_dstemr
      lb_dsteqr=>int_dsteqr
      lb_dsterf=>int_dsterf
      lb_dstevr=>int_dstevr
      lb_dsyev=>int_dsyev
      lb_dsyevr=>int_dsyevr
      lb_dsygs2=>int_dsygs2
      lb_dsygst=>int_dsygst
      lb_dsygv=>int_dsygv
      lb_dsytd2=>int_dsytd2
      lb_dsytrd=>int_dsytrd
      lb_dtrti2=>int_dtrti2
      lb_dtrtri=>int_dtrtri
      lb_dtrtrs=>int_dtrtrs
      lb_ieeeck=>int_ieeeck
      lb_iladlc=>int_iladlc
      lb_iladlr=>int_iladlr
      lb_ilaenv=>int_ilaenv
      lb_ilazlc=>int_ilazlc
      lb_ilazlr=>int_ilazlr
      lb_iparmq=>int_iparmq
      lb_iparam2stage=>int_iparam2stage
      lb_zheev=>int_zheev
      lb_zhetd2=>int_zhetd2
      lb_zhetrd=>int_zhetrd
      lb_zhpev=>int_zhpev
      lb_zhptrd=>int_zhptrd
      lb_zlacgv=>int_zlacgv
      lb_zladiv=>int_zladiv
      lb_zlanhe=>int_zlanhe
      lb_zlanhp=>int_zlanhp
      lb_zlarf=>int_zlarf
      lb_zlarfb=>int_zlarfb
      lb_zlarfg=>int_zlarfg
      lb_zlarft=>int_zlarft
      lb_zlascl=>int_zlascl
      lb_zlaset=>int_zlaset
      lb_zlasr=>int_zlasr
      lb_zlassq=>int_zlassq
      lb_zlatrd=>int_zlatrd
      lb_zsteqr=>int_zsteqr
      lb_zung2l=>int_zung2l
      lb_zung2r=>int_zung2r
      lb_zungql=>int_zungql
      lb_zungqr=>int_zungqr
      lb_zungtr=>int_zungtr
      lb_zupgtr=>int_zupgtr
      !
      ! Legacy
      lb_dgetf2=>int_dgetf2
      lb_dpotf2=>int_dpotf2
    end if

    if (prlev > 0) then
      ! BLAS
      !
      if (DLAddr(c_funloc(lb_dasum),c_loc(info)) /= 0) then
        write(6,*) 'dasum from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dasum found!'
      end if
      if (DLAddr(c_funloc(lb_daxpy),c_loc(info)) /= 0) then
        write(6,*) 'daxpy from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no daxpy found!'
      end if
      if (DLAddr(c_funloc(lb_dcabs1),c_loc(info)) /= 0) then
        write(6,*) 'dcabs1 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dcabs1 found!'
      end if
      if (DLAddr(c_funloc(lb_dcopy),c_loc(info)) /= 0) then
        write(6,*) 'dcopy from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dcopy found!'
      end if
      if (DLAddr(c_funloc(lb_ddot),c_loc(info)) /= 0) then
        write(6,*) 'ddot from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no ddot found!'
      end if
      if (DLAddr(c_funloc(lb_dgemm),c_loc(info)) /= 0) then
        write(6,*) 'dgemm from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgemm found!'
      end if
      if (DLAddr(c_funloc(lb_dgemv),c_loc(info)) /= 0) then
        write(6,*) 'dgemv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgemv found!'
      end if
      if (DLAddr(c_funloc(lb_dger),c_loc(info)) /= 0) then
        write(6,*) 'dger from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dger found!'
      end if
      if (DLAddr(c_funloc(lb_dnrm2),c_loc(info)) /= 0) then
        write(6,*) 'dnrm2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dnrm2 found!'
      end if
      if (DLAddr(c_funloc(lb_drot),c_loc(info)) /= 0) then
        write(6,*) 'drot from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no drot found!'
      end if
      if (DLAddr(c_funloc(lb_dscal),c_loc(info)) /= 0) then
        write(6,*) 'dscal from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dscal found!'
      end if
      if (DLAddr(c_funloc(lb_dspmv),c_loc(info)) /= 0) then
        write(6,*) 'dspmv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dspmv found!'
      end if
      if (DLAddr(c_funloc(lb_dspr2),c_loc(info)) /= 0) then
        write(6,*) 'dspr2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dspr2 found!'
      end if
      if (DLAddr(c_funloc(lb_dspr),c_loc(info)) /= 0) then
        write(6,*) 'dspr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dspr found!'
      end if
      if (DLAddr(c_funloc(lb_dswap),c_loc(info)) /= 0) then
        write(6,*) 'dswap from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dswap found!'
      end if
      if (DLAddr(c_funloc(lb_dsymm),c_loc(info)) /= 0) then
        write(6,*) 'dsymm from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsymm found!'
      end if
      if (DLAddr(c_funloc(lb_dsymv),c_loc(info)) /= 0) then
        write(6,*) 'dsymv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsymv found!'
      end if
      if (DLAddr(c_funloc(lb_dsyr2),c_loc(info)) /= 0) then
        write(6,*) 'dsyr2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsyr2 found!'
      end if
      if (DLAddr(c_funloc(lb_dsyr2k),c_loc(info)) /= 0) then
        write(6,*) 'dsyr2k from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsyr2k found!'
      end if
      if (DLAddr(c_funloc(lb_dsyrk),c_loc(info)) /= 0) then
        write(6,*) 'dsyrk from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsyrk found!'
      end if
      if (DLAddr(c_funloc(lb_dtpmv),c_loc(info)) /= 0) then
        write(6,*) 'dtpmv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dtpmv found!'
      end if
      if (DLAddr(c_funloc(lb_dtpsv),c_loc(info)) /= 0) then
        write(6,*) 'dtpsv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dtpsv found!'
      end if
      if (DLAddr(c_funloc(lb_dtrmm),c_loc(info)) /= 0) then
        write(6,*) 'dtrmm from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dtrmm found!'
      end if
      if (DLAddr(c_funloc(lb_dtrmv),c_loc(info)) /= 0) then
        write(6,*) 'dtrmv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dtrmv found!'
      end if
      if (DLAddr(c_funloc(lb_dtrsm),c_loc(info)) /= 0) then
        write(6,*) 'dtrsm from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dtrsm found!'
      end if
      if (DLAddr(c_funloc(lb_dtrsv),c_loc(info)) /= 0) then
        write(6,*) 'dtrsv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dtrsv found!'
      end if
      if (DLAddr(c_funloc(lb_dznrm2),c_loc(info)) /= 0) then
        write(6,*) 'dznrm2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dznrm2 found!'
      end if
      if (DLAddr(c_funloc(lb_idamax),c_loc(info)) /= 0) then
        write(6,*) 'idamax from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no idamax found!'
      end if
      if (DLAddr(c_funloc(lb_lsame),c_loc(info)) /= 0) then
        write(6,*) 'lsame from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no lsame found!'
      end if
      if (DLAddr(c_funloc(lb_scopy),c_loc(info)) /= 0) then
        write(6,*) 'scopy from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no scopy found!'
      end if
      if (DLAddr(c_funloc(lb_xerbla),c_loc(info)) /= 0) then
        write(6,*) 'xerbla from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no xerbla found!'
      end if
      if (DLAddr(c_funloc(lb_zaxpy),c_loc(info)) /= 0) then
        write(6,*) 'zaxpy from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zaxpy found!'
      end if
      if (DLAddr(c_funloc(lb_zcopy),c_loc(info)) /= 0) then
        write(6,*) 'zcopy from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zcopy found!'
      end if
      if (DLAddr(c_funloc(lb_zdotc),c_loc(info)) /= 0) then
        write(6,*) 'zdotc from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zdotc found!'
      end if
      if (DLAddr(c_funloc(lb_zdscal),c_loc(info)) /= 0) then
        write(6,*) 'zdscal from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zdscal found!'
      end if
      if (DLAddr(c_funloc(lb_zgemm),c_loc(info)) /= 0) then
        write(6,*) 'zgemm from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zgemm found!'
      end if
      if (DLAddr(c_funloc(lb_zgemv),c_loc(info)) /= 0) then
        write(6,*) 'zgemv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zgemv found!'
      end if
      if (DLAddr(c_funloc(lb_zgerc),c_loc(info)) /= 0) then
        write(6,*) 'zgerc from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zgerc found!'
      end if
      if (DLAddr(c_funloc(lb_zhemv),c_loc(info)) /= 0) then
        write(6,*) 'zhemv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zhemv found!'
      end if
      if (DLAddr(c_funloc(lb_zher2),c_loc(info)) /= 0) then
        write(6,*) 'zher2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zher2 found!'
      end if
      if (DLAddr(c_funloc(lb_zher2k),c_loc(info)) /= 0) then
        write(6,*) 'zher2k from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zher2k found!'
      end if
      if (DLAddr(c_funloc(lb_zhpmv),c_loc(info)) /= 0) then
        write(6,*) 'zhpmv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zhpmv found!'
      end if
      if (DLAddr(c_funloc(lb_zhpr2),c_loc(info)) /= 0) then
        write(6,*) 'zhpr2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zhpr2 found!'
      end if
      if (DLAddr(c_funloc(lb_zscal),c_loc(info)) /= 0) then
        write(6,*) 'zscal from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zscal found!'
      end if
      if (DLAddr(c_funloc(lb_zswap),c_loc(info)) /= 0) then
        write(6,*) 'zswap from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zswap found!'
      end if
      if (DLAddr(c_funloc(lb_ztrmm),c_loc(info)) /= 0) then
        write(6,*) 'ztrmm from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no ztrmm found!'
      end if
      if (DLAddr(c_funloc(lb_ztrmv),c_loc(info)) /= 0) then
        write(6,*) 'ztrmv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no ztrmv found!'
      end if

      ! LAPACK
      !
      if (DLAddr(c_funloc(lb_dbdsqr),c_loc(info)) /= 0) then
        write(6,*) 'dbdsqr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dbdsqr found!'
      end if
      if (DLAddr(c_funloc(lb_dgebd2),c_loc(info)) /= 0) then
        write(6,*) 'dgebd2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgebd2 found!'
      end if
      if (DLAddr(c_funloc(lb_dgebrd),c_loc(info)) /= 0) then
        write(6,*) 'dgebrd from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgebrd found!'
      end if
      if (DLAddr(c_funloc(lb_dgelq2),c_loc(info)) /= 0) then
        write(6,*) 'dgelq2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgelq2 found!'
      end if
      if (DLAddr(c_funloc(lb_dgelqf),c_loc(info)) /= 0) then
        write(6,*) 'dgelqf from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgelqf found!'
      end if
      if (DLAddr(c_funloc(lb_dgels),c_loc(info)) /= 0) then
        write(6,*) 'dgels from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgels found!'
      end if
      if (DLAddr(c_funloc(lb_dgeqr2),c_loc(info)) /= 0) then
        write(6,*) 'dgeqr2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgeqr2 found!'
      end if
      if (DLAddr(c_funloc(lb_dgeqrf),c_loc(info)) /= 0) then
        write(6,*) 'dgeqrf from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgeqrf found!'
      end if
      if (DLAddr(c_funloc(lb_dgesvd),c_loc(info)) /= 0) then
        write(6,*) 'dgesvd from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgesvd found!'
      end if
      if (DLAddr(c_funloc(lb_dgesv),c_loc(info)) /= 0) then
        write(6,*) 'dgesv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgesv found!'
      end if
      if (DLAddr(c_funloc(lb_dgetrf),c_loc(info)) /= 0) then
        write(6,*) 'dgetrf from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgetrf found!'
      end if
      if (DLAddr(c_funloc(lb_dgetrf2),c_loc(info)) /= 0) then
        write(6,*) 'dgetrf2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgetrf2 found!'
      end if
      if (DLAddr(c_funloc(lb_dgetri),c_loc(info)) /= 0) then
        write(6,*) 'dgetri from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgetri found!'
      end if
      if (DLAddr(c_funloc(lb_dgetrs),c_loc(info)) /= 0) then
        write(6,*) 'dgetrs from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgetrs found!'
      end if
      if (DLAddr(c_funloc(lb_disnan),c_loc(info)) /= 0) then
        write(6,*) 'disnan from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no disnan found!'
      end if
      if (DLAddr(c_funloc(lb_dlabad),c_loc(info)) /= 0) then
        write(6,*) 'dlabad from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlabad found!'
      end if
      if (DLAddr(c_funloc(lb_dlabrd),c_loc(info)) /= 0) then
        write(6,*) 'dlabrd from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlabrd found!'
      end if
      if (DLAddr(c_funloc(lb_dlacpy),c_loc(info)) /= 0) then
        write(6,*) 'dlacpy from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlacpy found!'
      end if
      if (DLAddr(c_funloc(lb_dladiv),c_loc(info)) /= 0) then
        write(6,*) 'dladiv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dladiv found!'
      end if
      if (DLAddr(c_funloc(lb_dladiv1),c_loc(info)) /= 0) then
        write(6,*) 'dladiv1 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dladiv1 found!'
      end if
      if (DLAddr(c_funloc(lb_dladiv2),c_loc(info)) /= 0) then
        write(6,*) 'dladiv2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dladiv2 found!'
      end if
      if (DLAddr(c_funloc(lb_dlae2),c_loc(info)) /= 0) then
        write(6,*) 'dlae2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlae2 found!'
      end if
      if (DLAddr(c_funloc(lb_dlaebz),c_loc(info)) /= 0) then
        write(6,*) 'dlaebz from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlaebz found!'
      end if
      if (DLAddr(c_funloc(lb_dlaev2),c_loc(info)) /= 0) then
        write(6,*) 'dlaev2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlaev2 found!'
      end if
      if (DLAddr(c_funloc(lb_dlagtf),c_loc(info)) /= 0) then
        write(6,*) 'dlagtf from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlagtf found!'
      end if
      if (DLAddr(c_funloc(lb_dlagts),c_loc(info)) /= 0) then
        write(6,*) 'dlagts from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlagts found!'
      end if
      if (DLAddr(c_funloc(lb_dlaisnan),c_loc(info)) /= 0) then
        write(6,*) 'dlaisnan from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlaisnan found!'
      end if
      if (DLAddr(c_funloc(lb_dlamch),c_loc(info)) /= 0) then
        write(6,*) 'dlamch from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlamch found!'
      end if
      if (DLAddr(c_funloc(lb_dlaneg),c_loc(info)) /= 0) then
        write(6,*) 'dlaneg from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlaneg found!'
      end if
      if (DLAddr(c_funloc(lb_dlange),c_loc(info)) /= 0) then
        write(6,*) 'dlange from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlange found!'
      end if
      if (DLAddr(c_funloc(lb_dlansp),c_loc(info)) /= 0) then
        write(6,*) 'dlansp from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlansp found!'
      end if
      if (DLAddr(c_funloc(lb_dlanst),c_loc(info)) /= 0) then
        write(6,*) 'dlanst from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlanst found!'
      end if
      if (DLAddr(c_funloc(lb_dlansy),c_loc(info)) /= 0) then
        write(6,*) 'dlansy from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlansy found!'
      end if
      if (DLAddr(c_funloc(lb_dlapy2),c_loc(info)) /= 0) then
        write(6,*) 'dlapy2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlapy2 found!'
      end if
      if (DLAddr(c_funloc(lb_dlapy3),c_loc(info)) /= 0) then
        write(6,*) 'dlapy3 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlapy3 found!'
      end if
      if (DLAddr(c_funloc(lb_dlar1v),c_loc(info)) /= 0) then
        write(6,*) 'dlar1v from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlar1v found!'
      end if
      if (DLAddr(c_funloc(lb_dlarfb),c_loc(info)) /= 0) then
        write(6,*) 'dlarfb from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarfb found!'
      end if
      if (DLAddr(c_funloc(lb_dlarf),c_loc(info)) /= 0) then
        write(6,*) 'dlarf from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarf found!'
      end if
      if (DLAddr(c_funloc(lb_dlarfg),c_loc(info)) /= 0) then
        write(6,*) 'dlarfg from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarfg found!'
      end if
      if (DLAddr(c_funloc(lb_dlarft),c_loc(info)) /= 0) then
        write(6,*) 'dlarft from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarft found!'
      end if
      if (DLAddr(c_funloc(lb_dlarnv),c_loc(info)) /= 0) then
        write(6,*) 'dlarnv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarnv found!'
      end if
      if (DLAddr(c_funloc(lb_dlarra),c_loc(info)) /= 0) then
        write(6,*) 'dlarra from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarra found!'
      end if
      if (DLAddr(c_funloc(lb_dlarrb),c_loc(info)) /= 0) then
        write(6,*) 'dlarrb from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarrb found!'
      end if
      if (DLAddr(c_funloc(lb_dlarrc),c_loc(info)) /= 0) then
        write(6,*) 'dlarrc from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarrc found!'
      end if
      if (DLAddr(c_funloc(lb_dlarrd),c_loc(info)) /= 0) then
        write(6,*) 'dlarrd from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarrd found!'
      end if
      if (DLAddr(c_funloc(lb_dlarre),c_loc(info)) /= 0) then
        write(6,*) 'dlarre from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarre found!'
      end if
      if (DLAddr(c_funloc(lb_dlarrf),c_loc(info)) /= 0) then
        write(6,*) 'dlarrf from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarrf found!'
      end if
      if (DLAddr(c_funloc(lb_dlarrj),c_loc(info)) /= 0) then
        write(6,*) 'dlarrj from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarrj found!'
      end if
      if (DLAddr(c_funloc(lb_dlarrk),c_loc(info)) /= 0) then
        write(6,*) 'dlarrk from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarrk found!'
      end if
      if (DLAddr(c_funloc(lb_dlarrr),c_loc(info)) /= 0) then
        write(6,*) 'dlarrr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarrr found!'
      end if
      if (DLAddr(c_funloc(lb_dlarrv),c_loc(info)) /= 0) then
        write(6,*) 'dlarrv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlarrv found!'
      end if
      if (DLAddr(c_funloc(lb_dlartg),c_loc(info)) /= 0) then
        write(6,*) 'dlartg from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlartg found!'
      end if
      if (DLAddr(c_funloc(lb_dlaruv),c_loc(info)) /= 0) then
        write(6,*) 'dlaruv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlaruv found!'
      end if
      if (DLAddr(c_funloc(lb_dlas2),c_loc(info)) /= 0) then
        write(6,*) 'dlas2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlas2 found!'
      end if
      if (DLAddr(c_funloc(lb_dlascl),c_loc(info)) /= 0) then
        write(6,*) 'dlascl from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlascl found!'
      end if
      if (DLAddr(c_funloc(lb_dlaset),c_loc(info)) /= 0) then
        write(6,*) 'dlaset from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlaset found!'
      end if
      if (DLAddr(c_funloc(lb_dlasq1),c_loc(info)) /= 0) then
        write(6,*) 'dlasq1 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlasq1 found!'
      end if
      if (DLAddr(c_funloc(lb_dlasq2),c_loc(info)) /= 0) then
        write(6,*) 'dlasq2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlasq2 found!'
      end if
      if (DLAddr(c_funloc(lb_dlasq3),c_loc(info)) /= 0) then
        write(6,*) 'dlasq3 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlasq3 found!'
      end if
      if (DLAddr(c_funloc(lb_dlasq4),c_loc(info)) /= 0) then
        write(6,*) 'dlasq4 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlasq4 found!'
      end if
      if (DLAddr(c_funloc(lb_dlasq5),c_loc(info)) /= 0) then
        write(6,*) 'dlasq5 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlasq5 found!'
      end if
      if (DLAddr(c_funloc(lb_dlasq6),c_loc(info)) /= 0) then
        write(6,*) 'dlasq6 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlasq6 found!'
      end if
      if (DLAddr(c_funloc(lb_dlasr),c_loc(info)) /= 0) then
        write(6,*) 'dlasr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlasr found!'
      end if
      if (DLAddr(c_funloc(lb_dlasrt),c_loc(info)) /= 0) then
        write(6,*) 'dlasrt from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlasrt found!'
      end if
      if (DLAddr(c_funloc(lb_dlassq),c_loc(info)) /= 0) then
        write(6,*) 'dlassq from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlassq found!'
      end if
      if (DLAddr(c_funloc(lb_dlasv2),c_loc(info)) /= 0) then
        write(6,*) 'dlasv2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlasv2 found!'
      end if
      if (DLAddr(c_funloc(lb_dlaswp),c_loc(info)) /= 0) then
        write(6,*) 'dlaswp from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlaswp found!'
      end if
      if (DLAddr(c_funloc(lb_dlatrd),c_loc(info)) /= 0) then
        write(6,*) 'dlatrd from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dlatrd found!'
      end if
      if (DLAddr(c_funloc(lb_dopgtr),c_loc(info)) /= 0) then
        write(6,*) 'dopgtr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dopgtr found!'
      end if
      if (DLAddr(c_funloc(lb_dopmtr),c_loc(info)) /= 0) then
        write(6,*) 'dopmtr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dopmtr found!'
      end if
      if (DLAddr(c_funloc(lb_dorg2l),c_loc(info)) /= 0) then
        write(6,*) 'dorg2l from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dorg2l found!'
      end if
      if (DLAddr(c_funloc(lb_dorg2r),c_loc(info)) /= 0) then
        write(6,*) 'dorg2r from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dorg2r found!'
      end if
      if (DLAddr(c_funloc(lb_dorgbr),c_loc(info)) /= 0) then
        write(6,*) 'dorgbr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dorgbr found!'
      end if
      if (DLAddr(c_funloc(lb_dorgl2),c_loc(info)) /= 0) then
        write(6,*) 'dorgl2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dorgl2 found!'
      end if
      if (DLAddr(c_funloc(lb_dorglq),c_loc(info)) /= 0) then
        write(6,*) 'dorglq from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dorglq found!'
      end if
      if (DLAddr(c_funloc(lb_dorgql),c_loc(info)) /= 0) then
        write(6,*) 'dorgql from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dorgql found!'
      end if
      if (DLAddr(c_funloc(lb_dorgqr),c_loc(info)) /= 0) then
        write(6,*) 'dorgqr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dorgqr found!'
      end if
      if (DLAddr(c_funloc(lb_dorgtr),c_loc(info)) /= 0) then
        write(6,*) 'dorgtr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dorgtr found!'
      end if
      if (DLAddr(c_funloc(lb_dorm2l),c_loc(info)) /= 0) then
        write(6,*) 'dorm2l from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dorm2l found!'
      end if
      if (DLAddr(c_funloc(lb_dorm2r),c_loc(info)) /= 0) then
        write(6,*) 'dorm2r from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dorm2r found!'
      end if
      if (DLAddr(c_funloc(lb_dormbr),c_loc(info)) /= 0) then
        write(6,*) 'dormbr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dormbr found!'
      end if
      if (DLAddr(c_funloc(lb_dorml2),c_loc(info)) /= 0) then
        write(6,*) 'dorml2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dorml2 found!'
      end if
      if (DLAddr(c_funloc(lb_dormlq),c_loc(info)) /= 0) then
        write(6,*) 'dormlq from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dormlq found!'
      end if
      if (DLAddr(c_funloc(lb_dormql),c_loc(info)) /= 0) then
        write(6,*) 'dormql from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dormql found!'
      end if
      if (DLAddr(c_funloc(lb_dormqr),c_loc(info)) /= 0) then
        write(6,*) 'dormqr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dormqr found!'
      end if
      if (DLAddr(c_funloc(lb_dormtr),c_loc(info)) /= 0) then
        write(6,*) 'dormtr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dormtr found!'
      end if
      if (DLAddr(c_funloc(lb_dposv),c_loc(info)) /= 0) then
        write(6,*) 'dposv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dposv found!'
      end if
      if (DLAddr(c_funloc(lb_dpotrf),c_loc(info)) /= 0) then
        write(6,*) 'dpotrf from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dpotrf found!'
      end if
      if (DLAddr(c_funloc(lb_dpotrf2),c_loc(info)) /= 0) then
        write(6,*) 'dpotrf2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dpotrf2 found!'
      end if
      if (DLAddr(c_funloc(lb_dpotrs),c_loc(info)) /= 0) then
        write(6,*) 'dpotrs from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dpotrs found!'
      end if
      if (DLAddr(c_funloc(lb_dpptrf),c_loc(info)) /= 0) then
        write(6,*) 'dpptrf from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dpptrf found!'
      end if
      if (DLAddr(c_funloc(lb_dspev),c_loc(info)) /= 0) then
        write(6,*) 'dspev from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dspev found!'
      end if
      if (DLAddr(c_funloc(lb_dspgst),c_loc(info)) /= 0) then
        write(6,*) 'dspgst from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dspgst found!'
      end if
      if (DLAddr(c_funloc(lb_dspgv),c_loc(info)) /= 0) then
        write(6,*) 'dspgv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dspgv found!'
      end if
      if (DLAddr(c_funloc(lb_dsptrd),c_loc(info)) /= 0) then
        write(6,*) 'dsptrd from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsptrd found!'
      end if
      if (DLAddr(c_funloc(lb_dstebz),c_loc(info)) /= 0) then
        write(6,*) 'dstebz from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dstebz found!'
      end if
      if (DLAddr(c_funloc(lb_dstein),c_loc(info)) /= 0) then
        write(6,*) 'dstein from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dstein found!'
      end if
      if (DLAddr(c_funloc(lb_dstemr),c_loc(info)) /= 0) then
        write(6,*) 'dstemr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dstemr found!'
      end if
      if (DLAddr(c_funloc(lb_dsteqr),c_loc(info)) /= 0) then
        write(6,*) 'dsteqr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsteqr found!'
      end if
      if (DLAddr(c_funloc(lb_dsterf),c_loc(info)) /= 0) then
        write(6,*) 'dsterf from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsterf found!'
      end if
      if (DLAddr(c_funloc(lb_dstevr),c_loc(info)) /= 0) then
        write(6,*) 'dstevr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dstevr found!'
      end if
      if (DLAddr(c_funloc(lb_dsyev),c_loc(info)) /= 0) then
        write(6,*) 'dsyev from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsyev found!'
      end if
      if (DLAddr(c_funloc(lb_dsyevr),c_loc(info)) /= 0) then
        write(6,*) 'dsyevr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsyevr found!'
      end if
      if (DLAddr(c_funloc(lb_dsygs2),c_loc(info)) /= 0) then
        write(6,*) 'dsygs2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsygs2 found!'
      end if
      if (DLAddr(c_funloc(lb_dsygst),c_loc(info)) /= 0) then
        write(6,*) 'dsygst from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsygst found!'
      end if
      if (DLAddr(c_funloc(lb_dsygv),c_loc(info)) /= 0) then
        write(6,*) 'dsygv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsygv found!'
      end if
      if (DLAddr(c_funloc(lb_dsytd2),c_loc(info)) /= 0) then
        write(6,*) 'dsytd2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsytd2 found!'
      end if
      if (DLAddr(c_funloc(lb_dsytrd),c_loc(info)) /= 0) then
        write(6,*) 'dsytrd from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dsytrd found!'
      end if
      if (DLAddr(c_funloc(lb_dtrti2),c_loc(info)) /= 0) then
        write(6,*) 'dtrti2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dtrti2 found!'
      end if
      if (DLAddr(c_funloc(lb_dtrtri),c_loc(info)) /= 0) then
        write(6,*) 'dtrtri from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dtrtri found!'
      end if
      if (DLAddr(c_funloc(lb_dtrtrs),c_loc(info)) /= 0) then
        write(6,*) 'dtrtrs from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dtrtrs found!'
      end if
      if (DLAddr(c_funloc(lb_ieeeck),c_loc(info)) /= 0) then
        write(6,*) 'ieeeck from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no ieeeck found!'
      end if
      if (DLAddr(c_funloc(lb_iladlc),c_loc(info)) /= 0) then
        write(6,*) 'iladlc from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no iladlc found!'
      end if
      if (DLAddr(c_funloc(lb_iladlr),c_loc(info)) /= 0) then
        write(6,*) 'iladlr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no iladlr found!'
      end if
      if (DLAddr(c_funloc(lb_ilaenv),c_loc(info)) /= 0) then
        write(6,*) 'ilaenv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no ilaenv found!'
      end if
      if (DLAddr(c_funloc(lb_ilazlc),c_loc(info)) /= 0) then
        write(6,*) 'ilazlc from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no ilazlc found!'
      end if
      if (DLAddr(c_funloc(lb_ilazlr),c_loc(info)) /= 0) then
        write(6,*) 'ilazlr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no ilazlr found!'
      end if
      if (DLAddr(c_funloc(lb_iparmq),c_loc(info)) /= 0) then
        write(6,*) 'iparmq from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no iparmq found!'
      end if
      if (DLAddr(c_funloc(lb_iparam2stage),c_loc(info)) /= 0) then
        write(6,*) 'iparam2stage from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no iparam2stage found!'
      end if
      if (DLAddr(c_funloc(lb_zheev),c_loc(info)) /= 0) then
        write(6,*) 'zheev from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zheev found!'
      end if
      if (DLAddr(c_funloc(lb_zhetd2),c_loc(info)) /= 0) then
        write(6,*) 'zhetd2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zhetd2 found!'
      end if
      if (DLAddr(c_funloc(lb_zhetrd),c_loc(info)) /= 0) then
        write(6,*) 'zhetrd from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zhetrd found!'
      end if
      if (DLAddr(c_funloc(lb_zhpev),c_loc(info)) /= 0) then
        write(6,*) 'zhpev from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zhpev found!'
      end if
      if (DLAddr(c_funloc(lb_zhptrd),c_loc(info)) /= 0) then
        write(6,*) 'zhptrd from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zhptrd found!'
      end if
      if (DLAddr(c_funloc(lb_zlacgv),c_loc(info)) /= 0) then
        write(6,*) 'zlacgv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zlacgv found!'
      end if
      if (DLAddr(c_funloc(lb_zladiv),c_loc(info)) /= 0) then
        write(6,*) 'zladiv from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zladiv found!'
      end if
      if (DLAddr(c_funloc(lb_zlanhe),c_loc(info)) /= 0) then
        write(6,*) 'zlanhe from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zlanhe found!'
      end if
      if (DLAddr(c_funloc(lb_zlanhp),c_loc(info)) /= 0) then
        write(6,*) 'zlanhp from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zlanhp found!'
      end if
      if (DLAddr(c_funloc(lb_zlarf),c_loc(info)) /= 0) then
        write(6,*) 'zlarf from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zlarf found!'
      end if
      if (DLAddr(c_funloc(lb_zlarfb),c_loc(info)) /= 0) then
        write(6,*) 'zlarfb from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zlarfb found!'
      end if
      if (DLAddr(c_funloc(lb_zlarfg),c_loc(info)) /= 0) then
        write(6,*) 'zlarfg from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zlarfg found!'
      end if
      if (DLAddr(c_funloc(lb_zlarft),c_loc(info)) /= 0) then
        write(6,*) 'zlarft from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zlarft found!'
      end if
      if (DLAddr(c_funloc(lb_zlascl),c_loc(info)) /= 0) then
        write(6,*) 'zlascl from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zlascl found!'
      end if
      if (DLAddr(c_funloc(lb_zlaset),c_loc(info)) /= 0) then
        write(6,*) 'zlaset from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zlaset found!'
      end if
      if (DLAddr(c_funloc(lb_zlasr),c_loc(info)) /= 0) then
        write(6,*) 'zlasr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zlasr found!'
      end if
      if (DLAddr(c_funloc(lb_zlassq),c_loc(info)) /= 0) then
        write(6,*) 'zlassq from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zlassq found!'
      end if
      if (DLAddr(c_funloc(lb_zlatrd),c_loc(info)) /= 0) then
        write(6,*) 'zlatrd from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zlatrd found!'
      end if
      if (DLAddr(c_funloc(lb_zsteqr),c_loc(info)) /= 0) then
        write(6,*) 'zsteqr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zsteqr found!'
      end if
      if (DLAddr(c_funloc(lb_zung2l),c_loc(info)) /= 0) then
        write(6,*) 'zung2l from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zung2l found!'
      end if
      if (DLAddr(c_funloc(lb_zung2r),c_loc(info)) /= 0) then
        write(6,*) 'zung2r from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zung2r found!'
      end if
      if (DLAddr(c_funloc(lb_zungql),c_loc(info)) /= 0) then
        write(6,*) 'zungql from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zungql found!'
      end if
      if (DLAddr(c_funloc(lb_zungqr),c_loc(info)) /= 0) then
        write(6,*) 'zungqr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zungqr found!'
      end if
      if (DLAddr(c_funloc(lb_zungtr),c_loc(info)) /= 0) then
        write(6,*) 'zungtr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zungtr found!'
      end if
      if (DLAddr(c_funloc(lb_zupgtr),c_loc(info)) /= 0) then
        write(6,*) 'zupgtr from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no zupgtr found!'
      end if

      ! Legacy
      !
      if (DLAddr(c_funloc(lb_dgetf2),c_loc(info)) /= 0) then
        write(6,*) 'dgetf2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dgetf2 found!'
      end if
      if (DLAddr(c_funloc(lb_dpotf2),c_loc(info)) /= 0) then
        write(6,*) 'dpotf2 from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no dpotf2 found!'
      end if
    end if

  end subroutine lb_initialize

!===============================================================================

  subroutine lb_close()
    integer :: i,rc

    if (allocated(handles)) then
      do i=1,size(handles)
        if (c_associated(handles(i))) then
          rc = dlclose(handles(i))
          if (rc /= 0) then
            write(6,*) c_f_string(dlerror())
          end if
        end if
      end do
      deallocate(handles)
    end if

  end subroutine lb_close

!===============================================================================

  subroutine split_string(string,array)
    character(len=*) :: string
    character(len=*), dimension(:), allocatable :: array
    integer :: i,n,offset,lim

    n=1
    do i=1,len(string)
      if (string(i:i) == ':') n=n+1
    end do
    allocate(array(n))
    offset=0
    do i=1,n-1
      lim=offset+index(string(offset+1:),':')
      array(i)=string(offset+1:lim-1)
      offset=lim
    end do
    array(n)=string(offset+1:)

  end subroutine split_string

!===============================================================================

  function link_func(funname)
    character(kind=c_char,len=*) :: funname
    type(c_funptr) :: link_func
    integer :: i
    logical :: success

!******************************************************
!   Try to associate a function on all loaded libraries
!******************************************************
    success=.false.
    link_func=c_null_funptr
    do i=1,size(handles)
      link_func=dlsym(handles(i),trim(funname)//'_'//c_null_char)
      if (c_associated(link_func)) then
        success=.true.
        exit
#ifdef _DEBUG_
      else
        write(6,*) c_f_string(dlerror())
#endif
      end if
    end do
#ifdef _DEBUG_
    if (.not. success) then
      write(6,*) 'no ',trim(funname),' found'
    end if
#endif

  end function link_func

end module link_blas
