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
! Copyright (C) 2000, Jeppe Olsen                                      *
!***********************************************************************

subroutine INTPNT(IPNT1,ISL1,IPNT2,ISL2)
! Pointers to symmetry blocks of integrals
! IPNT1 : Pointer to given one-electron block, total symmetric
! ISL1  : Symmetry of last index for given first index, 1 e-
! IPNT2 : Pointer to given two-electron block
! ISL1  : Symmetry of last index for given first index, 1 e-
!
! In addition pointers to one-electron integrals with general
! symmetry is generated in PGINT1(ISM)%A
!
! Pointers for similarity transformed Hamiltonian may also be
! generated
!
! Jeppe Olsen, Last Update : August 2000

use lucia_data, only: I1234S, I12S, I34S, NSMOB, NTOOBS, PGINT1, PGINT1A
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: IPNT1(NSMOB), ISL1(NSMOB), IPNT2(NSMOB,NSMOB,NSMOB), ISL2(NSMOB,NSMOB,NSMOB)
integer(kind=iwp) :: ISM

! 0 : Pointers to one-integrals, all symmetries, Lower half matrices
do ISM=1,NSMOB
  call PNT2DM(1,NSMOB,NTOOBS,NTOOBS,ISM,ISL1,PGINT1(ISM)%A)
end do
! 0.5 : Pointers to one-electron integrals, all symmetries, complete form
do ISM=1,NSMOB
  call PNT2DM(0,NSMOB,NTOOBS,NTOOBS,ISM,ISL1,PGINT1A(ISM)%A)
end do
! 1 : Number of one-electron integrals
call PNT2DM(1,NSMOB,NTOOBS,NTOOBS,1,ISL1,IPNT1)
! 2 : two-electron integrals
call PNT4DM(NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,1,I12S,I34S,I1234S,IPNT2,ISL2)

end subroutine INTPNT
