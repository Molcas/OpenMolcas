************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1998  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE W1TW2(IVEC,JVEC,CI,SGM)
      use definitions, only: iwp, wp
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: nAshT, STSym
      implicit none


      integer(kind=iwp), intent(in):: IVEC, JVEC
      Real(kind=wp), intent(in):: ci(*)
      Real(kind=wp), intent(inout)::  sgm(*)

      integer(kind=iwp) :: nOp1, nOp2, nOp3
      Real(kind=wp), Allocatable:: OP1(:), OP2(:), OP3(:)
      Real(kind=wp):: OP0

C Given contravariant indices of two wave operators W1 and W2,
C in the vectors numbered IVEC and JVEC on file (unit LUSOLV),
C compute the vector in CAS space
C   | SGM > := | SGM > + (W1 conj)*(W2)*| CI >


C (1): Compute a representation of the operator PCAS*W1T*W2
      NOP1=NASHT**2
      NOP2=(NOP1*(NOP1+1))/2
      NOP3=(NOP2*(NOP1+2))/3
      CALL mma_allocate(OP1,NOP1,Label='OP1')
      CALL mma_allocate(OP2,NOP2,Label='OP2')
      CALL mma_allocate(OP3,NOP3,Label='OP3')

      CALL MKWWOP(IVEC,JVEC,OP0,OP1,NOP2,OP2,NOP3,OP3)

C Modify the coefficients, see subroutine MODOP.

      CALL MODOP(OP1,NOP2,OP2,NOP3,OP3)

C (2) Apply the operators:
      CALL HAM3(OP0,OP1,NOP2,OP2,NOP3,OP3,STSYM,CI,SGM)

      CALL mma_deallocate(OP1)
      CALL mma_deallocate(OP2)
      CALL mma_deallocate(OP3)

      END SUBROUTINE W1TW2
