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

use InfSO, only: DltNTh
use InfSCF, only: EThr, DThr, FThr
use Save_Stuff, only: DltNTh_old, DThr_Old, EThr_old, FThr_Old, SIntTh_old, ThrInt_old
use Constants, only: Zero, One, Ten

implicit none
real*8 EThr_, SIntTh, Relax
real*8, external :: Get_ThrInt

write(6,*)
write(6,*) 'Temporary increase of thresholds...'
write(6,*)
SIntTh_old = SIntTh
EThr_old = EThr
DThr_old = DThr
DltNTh_old = DltNTh
FThr_old = FThr

! Get threshold used in connection of products of integrals and densities

ThrInt_Old = Get_ThrInt()

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
call xSet_ThrInt(ThrInt_Old*Relax)

return

end subroutine Reduce_Thresholds
