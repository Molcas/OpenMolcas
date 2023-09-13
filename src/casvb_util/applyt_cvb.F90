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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

!************************************
!** Routines involving CI and ORBS **
!************************************
!***********************************************************************
!*                                                                     *
!*  APPLYT := transform CI vector IVEC according to orbital            *
!*         := transformation in IGJORB (from GAUSSJ)                   *
!*                                                                     *
!*  :> IVEC:      starting CI vector                                   *
!*  :> IGJORB:    simple orbital updates                               *
!*  :> N_APPLYT   stats                                                *
!*  :> I1ALF:     string handling                                      *
!*  :> I1BET:           do.                                            *
!*  :> IATO:            do.                                            *
!*  :> IBTO:            do.                                            *
!*  :> PHATO:           do.                                            *
!*  :> PHBTO:           do.                                            *
!*                                                                     *
!*  :< IVEC:      transformed CI vector                                *
!*  :< N_APPLYT   stats                                                *
!*                                                                     *
!***********************************************************************
subroutine applyt_cvb(cvec,gjorb)

implicit real*8(a-h,o-z)
dimension gjorb(*), cvec(*)

call applyt_cvb_internal(gjorb)

! This is to allow type punning without an explicit interface
contains

subroutine applyt_cvb_internal(gjorb)
  use iso_c_binding
  real*8, target :: gjorb(*)
  integer, pointer :: igjorb(:)
  call c_f_pointer(c_loc(gjorb(1)),igjorb,[1])
  call iapplyt_cvb(cvec,igjorb)
  nullify(igjorb)
end subroutine applyt_cvb_internal

end subroutine applyt_cvb
