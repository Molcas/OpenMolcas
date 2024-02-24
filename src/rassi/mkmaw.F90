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

Call MKRAW(NVERT,IDOWN,IUP,IRAW)

Call MKMID(NVERT,NLEV,IDAW,IRAW,LTV,MIDLEV, NMIDV, MVSta, MVEnd, MXUP, MXDWN)

CALL MKMAW(IDOWN,IDAW,IUP,IRAW,IMAW,NVERT, MVSta, MVEnd)


END SUBROUTINE MKMAW_RASSI
