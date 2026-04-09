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
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine DETSTR_MCLR(IDET,IASTR,IBSTR,NEL,NAEL,NBEL,ISGN,IWORK)
! A DETERMINANT,IDET,IS GIVEN AS A SET OF OCCUPIED SPIN ORBITALS,
! POSITIVE NUMBER INDICATES ALPHA ORBITAL AND NEGATIVE NUMBER
! INDICATES BETA ORBITAL.
!
! FIND CORRESPONDING ALPHA STRING AND BETA STRING,
! AND DETERMINE SIGN NEEDED TO CHANGE DETERMINANT
! INTO PRODUCT OF ORDERED ALPHA STRING AND
! BETA STRING
!
! JEPPE OLSEN NOVEMBER 1988

use Index_Functions, only: nTri_Elem
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NEL, NAEL, NBEL
integer(kind=iwp), intent(inout) :: IDET(NEL)
integer(kind=iwp), intent(out) :: IASTR(NAEL), IBSTR(NBEL), ISGN, IWORK(NEL)
integer(kind=iwp) :: ITMP

! FIRST REORDER SPIN ORBITALS IN ASCENDING SEQUENCE
! THIS WILL AUTOMATICALLY SPLIT ALPHA AND BETASTRING

call ORDSTR_MCLR(IDET,IWORK,NEL,ISGN)

! ALPHA STRING IS LAST NAEL ORBITALS
IASTR(:) = IWORK(NBEL+1:NBEL+NAEL)

! BETA  STRING MUST BE COMPLETELY TURNED AROUND
IBSTR(:) = -IWORK(NBEL:1:-1)
! SIGN CHANGE FOR SWITCH OF BETA ORBITALS
iTmp = nTri_Elem(NBEL)
ISGN = ISGN*(-1)**iTmp

#ifdef _DEBUGPRINT_
write(u6,*) ' INPUT DETERMINANT'
call IWRTMA(IDET,1,NEL,1,NEL)
write(u6,*) ' CORRESPONDING ALPHA STRING'
call IWRTMA(IASTR,1,NAEL,1,NAEL)
write(u6,*) ' CORRESPONDING BETA STRING'
call IWRTMA(IBSTR,1,NBEL,1,NBEL)
write(u6,*) ' ISGN FOR SWITCH ',ISGN
#endif

!if (doDMRG .and. doMCLR) then ! yma
!  do I=1,NEL
!    write(117,'(1X,I5)',advance='no') IDET(I)
!  end do
!  write(117,'(A,1X,I2)',advance='no') ' SIGN',ISGN
!end if

end subroutine DETSTR_MCLR
