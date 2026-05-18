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
      SUBROUTINE MKWWOP(IVEC,JVEC,OP0,OP1,NOP2,OP2,NOP3,OP3)
      use caspt2_module, only: NASHT
      use constants, only: Zero
      use definitions, only: iwp, wp
      IMPLICIT None

! Presently symmetry blocking is disregarded for OP2, OP3, but
! index pair C permutation symmetry is used.
! NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
! NOP3=(NASHT**2+2 over 3)  (Binomial coefficient)
      integer(kind=iwp), intent(in):: IVEC, JVEC, NOP2, NOP3
      real(kind=wp), intent(out):: OP0, OP1(NASHT,NASHT), OP2(NOP2),    &
     &                             OP3(NOP3)

! Given the coefficients for two excitation operators in the
! vectors numbered IVEC and C JVEC on file, construct the
! zero-, one-, two-, and three-body
! expansions of the product (Op in IVEC conjugated)(Op in JVEC)
! as operating on the CASSCF space.

      OP0=Zero
      OP1(:,:)=Zero
      OP2(:)=Zero
      OP3(:)=Zero
      CALL MKWWOPA(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
      CALL MKWWOPB(IVEC,JVEC,OP0,OP1,NOP2,OP2)
      CALL MKWWOPC(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
      CALL MKWWOPD(IVEC,JVEC,OP1,NOP2,OP2)
      CALL MKWWOPE(IVEC,JVEC,OP0,OP1)
      CALL MKWWOPF(IVEC,JVEC,NOP2,OP2)
      CALL MKWWOPG(IVEC,JVEC,OP1)
      CALL MKWWOPH(IVEC,JVEC,OP0)

      END SUBROUTINE MKWWOP
