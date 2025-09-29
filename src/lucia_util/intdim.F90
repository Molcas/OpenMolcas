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

subroutine INTDIM()
! Number of integrals and storage mode

use lucia_data, only: I1234S, I12S, I34S, NBINT1, NBINT2, NSMOB

implicit none

! 1 : Number of one-electron integrals
!NINT1 = NSXFSM(NSMOB,NTOOBS,NTOOBS,1,1)
! 2 : Number of two-electron integrals
!if (PNTGRP == 1) then
!  ! Full eightfold symmetry can be used
I12S = 1
I34S = 1
I1234S = 1
!else
!  ! Only symmetry between 12 and 34
!  I12S = 0
!  I34S = 0
!  I1234S = 1
!end if
!NINT2 = NDXFSM(NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,1,I12S,I34S,I1234S)
! Number of integrals without complex conjugation symmetry
! (used for T1 transformed Hamiltonian)
!I12 = 0
!I34 = 0
!I1234 = 1
!NINT2_NO_CCSYM = NDXFSM(NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,1,I12,I34,I1234)
!if (ISIMTRH == 1) write(u6,*) ' Number of two-electron integrals in exp(-T1)Hexp(T1) ',NINT2_NO_CCSYM
! Number of symmetry blocks of one- and two-electron integrals
NBINT1 = NSMOB
NBINT2 = NSMOB**3

end subroutine INTDIM
