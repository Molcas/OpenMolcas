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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

function Cho_LK_ScreeningThreshold(delta)
! Thomas Bondo Pedersen, May 2013.
!
! Return the basic LK screening threshold. Input is a damping
! parameter, f.ex. the max. Fock matrix element in the previous
! iteration. Use delta < 0 or delta > 1 to avoid damping.
!
! The Cholesky environment must have been set up prior to calling
! this function [by calling Cho_X_Init(..)]

use Constants, only: Zero, One
use Definitions, only: wp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp) :: Cho_LK_ScreeningThreshold
real(kind=wp), intent(in) :: delta
#include "cholesky.fh"
real(kind=wp) :: thr
real(kind=wp), parameter :: thr0 = 1.0e-6_wp

thr = min(ThrCom,thr0)
if ((delta >= Zero) .and. (delta <= One)) then
  thr = thr*delta
end if
Cho_LK_ScreeningThreshold = max(thr,1.0e-15_wp)
#ifdef _DEBUGPRINT_
write(u6,'(1P,4(A,ES15.6))') 'ThrCom=',ThrCom,' thr0=',thr0,' delta=',delta,' Cho_LK_ScreeningThreshold=',Cho_LK_ScreeningThreshold
#endif

end function Cho_LK_ScreeningThreshold
