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

real*8 function Cho_LK_ScreeningThreshold(delta)
! Thomas Bondo Pedersen, May 2013.
!
! Return the basic LK screening threshold. Input is a damping
! parameter, f.ex. the max. Fock matrix element in the previous
! iteration. Use delta<0.0d0 or delta>1.0d0 to avoid damping.
!
! The Cholesky environment must have been set up prior to calling
! this function [by calling Cho_X_Init(..)]

implicit none
real*8 delta
#include "cholesky.fh"
real*8 thr0
parameter(thr0=1.0d-6)
real*8 thr

thr = min(ThrCom,thr0)
if ((delta >= 0.0d0) .and. (delta <= 1.0d0)) then
  thr = thr*delta
end if
Cho_LK_ScreeningThreshold = max(thr,1.0d-15)
#ifdef _DEBUGPRINT_
write(6,'(1P,4(A,D15.6))') 'ThrCom=',ThrCom,' thr0=',thr0,' delta=',delta,' Cho_LK_ScreeningThreshold=',Cho_LK_ScreeningThreshold
#endif

end function Cho_LK_ScreeningThreshold
