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

use casvb_global, only: gjorb_type, i1alf, i1bet, iato, ibto, icnt_ci, iform_ci, n_applyt, ndet, phato, phbto
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: cvec(0:ndet)
type(gjorb_type), intent(in) :: gjorb
integer(kind=iwp) :: ivec

ivec = nint(cvec(0))
n_applyt = n_applyt+1
if (iform_ci(ivec) == 0) then
  call permci_cvb(cvec,gjorb%i1)
  call applyt2_cvb(cvec(1:),gjorb%r,gjorb%i2,i1alf,i1bet,iato,ibto,phato,phbto)
else
  write(u6,*) ' Unsupported format in APPLYT :',iform_ci(ivec)
  call abend_cvb()
end if
icnt_ci(ivec) = 0

return

end subroutine applyt_cvb
