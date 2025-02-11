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

use Symmetry_Info, only: Mul
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: NSMST, NSMCI, MXPCSM, ISMOST(MXPCSM,MXPCSM)
integer(kind=iwp) :: ISTSM, ITOTSM, JSTSM, NTEST

do ITOTSM=1,NSMCI
  do ISTSM=1,NSMST
    JSTSM = Mul(ISTSM,ITOTSM)
    ISMOST(ISTSM,ITOTSM) = JSTSM
  end do
end do

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) ' ==============='
  write(u6,*) ' Info from SMOST'
  write(u6,*) ' ==============='
  do ITOTSM=1,NSMCI
    write(u6,*) ' ISMOST array for ITOTSM = ',ITOTSM
    call IWRTMA(ISMOST(1,ITOTSM),1,NSMST,1,NSMST)
  end do
end if

end subroutine SMOST
