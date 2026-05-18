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
      SUBROUTINE MKWWOPH(IVEC,JVEC,OP0)
      use definitions, only: iwp, wp
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM, NASUP, NISUP, NINDEP
      IMPLICIT None

      integer(kind=iwp), intent(in):: IVEC, JVEC
      real(kind=wp), intent(inout):: OP0

      real(kind=wp), ALLOCATABLE:: W1(:), W2(:)
      integer(kind=iwp) ICASE, ISYM, NAS, NIS, MDVEC, ISCT, IISTA,      &
     &                  IIEND, NCOL, NSCT
      real(kind=wp), external:: DDot_

! Given the coefficients for two excitation operators, available in
! vectors numbered IVEC and JVEC on file, use the blocks for
! excitation case BJAI, i.e. cases 12 and 13, to express the
! product (Op in IVEC conjugated)(Op in JVEC) as a zero-body
! operator, i.e. a scalar factor, in the CASSCF space.
! Formula used:
!  (W1(ij,ab) conj)*(W2(kl,cd)) = dik*djl*dac*dbd
! the same for both cases 10 and 11.

! Loop over cases
      DO ICASE=12,13
! Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
! Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        CALL mma_allocate(W2,NAS*MDVEC,Label='W2')
! Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         NSCT=NAS*NCOL
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1,NAS*MDVEC)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2,NAS*MDVEC)
! Pick up a symmetry block of W1 and W2
         OP0=OP0+DDOT_(NSCT,W1,1,W2,1)
        END DO
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2)
      END DO
! End of loop over cases.
      END DO
      END SUBROUTINE MKWWOPH
