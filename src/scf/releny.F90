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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine RelEny(ERelMV,ERelDC,D,VMs,Dar,lth)
!***********************************************************************
!                                                                      *
!     purpose: Compute relativistic corrections to the energy          *
!                                                                      *
!     input:                                                           *
!       D       : density matrix of length lth                         *
!       VMs     : mass velocity integrals of length lth                *
!       Dar     : Darwin integrals of length lth                       *
!                                                                      *
!     output:                                                          *
!       ERelMV  : mass velocity correction                             *
!       ERelDC  : Darwin correction                                    *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: ERelMV, ERelDC
integer(kind=iwp), intent(in) :: lth
real(kind=wp), intent(in) :: D(lth), VMs(lth), Dar(lth)
real(kind=wp), external :: DDot_

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

ERelMV = DDot_(lth,D,1,VMs,1)
ERelDC = DDot_(lth,D,1,Dar,1)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine RelEny
