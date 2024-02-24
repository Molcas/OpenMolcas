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
SUBROUTINE MKMAW_RASSI(NLEV,NVERT,IDOWN,IDAW,IUP,IRAW,IMAW,LTV,MIDLEV)
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION IDOWN(NVERT,0:3),IDAW(NVERT,0:4)
DIMENSION IUP(NVERT,0:3),IRAW(NVERT,0:4)
DIMENSION IMAW(NVERT,0:3),LTV(-1:NLEV)
! Purpose: Construct the Modified Arc Weights

Call MKRAW(NVERT,IDOWN,IUP,IRAW)

! USE IDAW,IRAW TABLES TO DETERMINE MIDLEV.
! ASSUMPTION IS THAT A BALANCED NUMBER OF UPPER/LOWER WALKS
! IS THE BEST CHOICE.
MIDLEV=1
MINW=1000000
DO IL=1,NLEV-1
  NW=0
  DO IV=LTV(IL),LTV(IL-1)-1
    NW=NW+IRAW(IV,4)-IDAW(IV,4)
  END DO
  NW=ABS(NW)
  IF(NW.LT.MINW) THEN
    MIDLEV=IL
    MINW=NW
  END IF
END DO
MVSta=LTV(MIDLEV)
MVEnd=LTV(MIDLEV-1)-1

CALL MKMAW(IDOWN,IDAW,IUP,IRAW,IMAW,NVERT, MVSta, MVEnd)


END SUBROUTINE MKMAW_RASSI
