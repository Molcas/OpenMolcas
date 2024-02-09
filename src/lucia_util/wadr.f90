!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
Module wadr
Private
INTEGER, Public:: LDMAT,LPMAT,LPA,LDIA,LSXN,LBM,LF1,LF2,LG,LH,LHD, NLX,nPWXY, LFI, ipDens, lDSPN
Real*8, Allocatable, Public:: TUVX(:), FockOcc(:)
End Module wadr
