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
      SUBROUTINE MKWWOPG(IVEC,JVEC,OP1)
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM,NASUP,NISUP,NINDEP,NASHT,NAES
      IMPLICIT None

! Presently symmetry blocking is disregarded.
      integer(kind=iwp), intent(in):: IVEC, JVEC
      real(kind=wp), intent(inout):: OP1(NASHT,NASHT)

      real(kind=wp), ALLOCATABLE:: W1(:), W2(:), WPROD(:)
      integer(kind=iwp) ICASE, ISYM, NAS, NIS, MDVEC, NWPROD, ISCT,     &
     &                  IISTA, IIEND, NCOL, IW1, ITABS, IW2, IT, IX,    &
     &                  IWPROD, IXABS
      real(kind=wp) W_PROD

! Given the coefficients for two excitation operators, available in
! vectors numbered IVEC and JVEC on file, use the blocks for
! excitation case BJAT, i.e. cases 10 and 11, to express the
! product (Op in IVEC conjugated)(Op in JVEC) as a one-body
! operator on the CASSCF space.
! Formula used:
!  (W1(t,aij) conj)*(W2(x,bkl)) = dik*djl*dab* Etx
! the same for both cases 10 and 11.

! Loop over cases
      DO ICASE=10,11
! Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
! Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,LABEL='W1')
        CALL mma_allocate(W2,NAS*MDVEC,LABEL='W2')
        NWPROD=NAS**2
! Allocate space for the contraction:
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
        WPROD(:)=Zero
! Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1,NAS*MDVEC)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2,NAS*MDVEC)
! Multiply WProd = (W1)*(W2 transpose)
         CALL DGEMM_('N','T',                                           &
     &              NAS,NAS,NCOL,                                       &
     &              One,W1,NAS,                                         &
     &              W2,NAS,                                             &
     &              One,WPROD,NAS)
        END DO
! Deallocate space for this block of excitation amplitudes:
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2)

! Loop over (T)
          DO IT=1,NAS
            IW1=IT
            ITABS=IT+NAES(ISYM)
! Loop over (X)
          DO IX=1,NAS
            IW2=IX
            IXABS=IX+NAES(ISYM)
            IWPROD=IW1+NAS*(IW2-1)
            W_PROD=WPROD(IWPROD)
            OP1(ITABS,IXABS)=OP1(ITABS,IXABS)+W_PROD
          END DO
        END DO
! Deallocate matrix product
        CALL mma_deallocate(WPROD)
      END DO
! End of loop over cases.
      END DO
      END SUBROUTINE MKWWOPG
