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
!#define _DEBUGPRINT_
      SUBROUTINE MKLTV(NVERT,NLEV,IDRT,LTV)
!     PURPOSE: FIND THE MIDLEVEL
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      IMPLICIT None
!
      Integer NVERT, NLEV
      Integer IDRT(NVERT,5)
      Integer LTV(-1:NLEV)

      Integer, Parameter:: LTAB=1
      Integer IV, LEV
!
!     SET UP A LEVEL-TO-VERTEX TABLE, LTV, AND IDENTIFY MIDVERTICES:
!
      LTV(:)=0
!
      DO IV=1,NVERT
        LEV=IDRT(IV,LTAB)
        LTV(LEV)=LTV(LEV)+1
      End Do
!
      DO LEV=NLEV,0,-1
        LTV(LEV-1)=LTV(LEV-1)+LTV(LEV)
      End Do
!
      DO LEV=-1,NLEV-1
        LTV(LEV)=1+LTV(LEV+1)
      End Do

      END SUBROUTINE MKLTV
