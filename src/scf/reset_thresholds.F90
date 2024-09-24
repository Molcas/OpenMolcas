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

subroutine Reset_Thresholds()

use InfSO, only: DltNTh
use InfSCF, only: EThr, DThr, FThr
use Save_Stuff, only: DltNTh_old, DThr_Old, EThr_old, FThr_Old, ThrInt_Old

implicit none

write(6,*)
write(6,*) 'Restore thresholds...'
write(6,*)
EThr = EThr_old
DThr = DThr_old
DltNTh = DltNTh_old
FThr = FThr_old
call xSet_ThrInt(ThrInt_Old)

return

end subroutine Reset_Thresholds
