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
!               2003-2005, Valera Veryazov                             *
!               2017, Roland Lindh                                     *
!***********************************************************************

subroutine Reduce_Thresholds(EThr_,SIntTh)

use Gateway_Info, only: ThrInt
use InfSCF, only: DltNTh, DltNTh_old, DThr, DThr_old, EThr, EThr_old, FThr, FThr_old, ThrInt_old
use Constants, only: Zero, One, Ten
use Definitions, only: wp, u6

implicit none
real(kind=wp), intent(in) :: EThr_
real(kind=wp), intent(inout) :: SIntTh
real(kind=wp) :: Relax

write(u6,*)
write(u6,*) 'Temporary increase of thresholds...'
write(u6,*)
EThr_old = EThr
DThr_old = DThr
DltNTh_old = DltNTh
FThr_old = FThr

! Get threshold used in connection of products of integrals and densities

ThrInt_old = ThrInt

EThr = EThr_
if (EThr_old == Zero) then
  Relax = One
else
  Relax = EThr/EThr_old
end if
SIntTh = SIntTh*Relax
DThr = DThr*Relax
DltNTh = Ten**2*EThr
FThr = FThr*Relax
ThrInt = ThrInt*Relax

return

end subroutine Reduce_Thresholds
