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
! Copyright (C) 1991, Jeppe Olsen                                      *
!***********************************************************************

subroutine SMOST(NSMST,NSMCI,MXPCSM,ISMOST)
! ISMOST(ISYM,ITOTSM) : Symmetry of an internal state if ITOTSM
!                       if symmetry of 1 string is ISYM, the
!                       symmetry of the other string is
!                       ISMOST(ISYM,ITOTSM)
!
! Jeppe Olsen, Spring of 1991

implicit real*8(A-H,O-Z)
dimension ISMOST(MXPCSM,MXPCSM)

do ITOTSM=1,NSMCI
  do ISTSM=1,NSMST
    !    SYMCOM(ITASK,IOBJ,I1,I2,I12)
    call SYMCOM(2,1,ISTSM,JSTSM,ITOTSM)
    ISMOST(ISTSM,ITOTSM) = JSTSM
  end do
end do

NTEST = 0
if (NTEST /= 0) then
  write(6,*) ' ==============='
  write(6,*) ' Info from SMOST'
  write(6,*) ' ==============='
  do ITOTSM=1,NSMCI
    write(6,*) ' ISMOST array for ITOTSM = ',ITOTSM
    call IWRTMA(ISMOST(1,ITOTSM),1,NSMST,1,NSMST)
  end do
end if

end subroutine SMOST
