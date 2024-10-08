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

use Gateway_Info, only: ThrInt
use InfSCF, only: DltNTh, DltNTh_old, DThr, DThr_Old, EThr, EThr_old, FThr, FThr_Old, ThrInt_Old
use Definitions, only: u6

implicit none

write(u6,*)
write(u6,*) 'Restore thresholds...'
write(u6,*)
EThr = EThr_old
DThr = DThr_old
DltNTh = DltNTh_old
FThr = FThr_old
ThrInt = ThrInt_Old

return

end subroutine Reset_Thresholds
