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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine W1TW2(IVEC,JVEC,CI,SGM,nCI)
! Given contravariant indices of two wave operators W1 and W2,
! in the vectors numbered IVEC and JVEC on file (unit LUSOLV),
! compute the vector in CAS space
!   | SGM > := | SGM > + (W1 conj)*(W2)*| CI >

use Index_Functions, only: nTri_Elem, nTri3_Elem
use caspt2_module, only: nAshT, STSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC, nCI
real(kind=wp), intent(in) :: ci(nCI)
real(kind=wp), intent(inout) :: sgm(nCI)
integer(kind=iwp) :: nOp1, nOp2, nOp3
real(kind=wp) :: OP0
real(kind=wp), allocatable :: OP1(:), OP2(:), OP3(:)

! (1): Compute a representation of the operator PCAS*W1T*W2
NOP1 = NASHT**2
NOP2 = nTri_Elem(NOP1)
NOP3 = nTri3_Elem(NOP1)
call mma_allocate(OP1,NOP1,Label='OP1')
call mma_allocate(OP2,NOP2,Label='OP2')
call mma_allocate(OP3,NOP3,Label='OP3')

call MKWWOP(IVEC,JVEC,OP0,OP1,NOP2,OP2,NOP3,OP3)

! Modify the coefficients, see subroutine MODOP.

call MODOP(OP1,NOP2,OP2,NOP3,OP3)

! (2) Apply the operators:
call HAM3(OP0,OP1,NOP2,OP2,NOP3,OP3,STSYM,CI,SGM,nCI)

call mma_deallocate(OP1)
call mma_deallocate(OP2)
call mma_deallocate(OP3)

end subroutine W1TW2
