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

!#define _DEBUGPRINT_
subroutine SMOST(NSMST,NSMCI,MXPCSM,ISMOST)
! ISMOST(ISYM,ITOTSM) : Symmetry of an internal state if ITOTSM
!                       if symmetry of 1 string is ISYM, the
!                       symmetry of the other string is
!                       ISMOST(ISYM,ITOTSM)
!
! Jeppe Olsen, Spring of 1991

use Symmetry_Info, only: Mul
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NSMST, NSMCI, MXPCSM
integer(kind=iwp), intent(out) :: ISMOST(MXPCSM,MXPCSM)
integer(kind=iwp) :: ISTSM, ITOTSM, JSTSM

do ITOTSM=1,NSMCI
  do ISTSM=1,NSMST
    JSTSM = Mul(ISTSM,ITOTSM)
    ISMOST(ISTSM,ITOTSM) = JSTSM
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' ==============='
write(u6,*) ' Info from SMOST'
write(u6,*) ' ==============='
do ITOTSM=1,NSMCI
  write(u6,*) ' ISMOST array for ITOTSM = ',ITOTSM
  call IWRTMA(ISMOST(:,ITOTSM),1,NSMST,1,NSMST)
end do
#endif

end subroutine SMOST
