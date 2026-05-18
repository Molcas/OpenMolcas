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
      SUBROUTINE MKWWOPF(IVEC,JVEC,NOP2,OP2)
      use definitions, only: iwp, wp
      use constants, only: Zero, One, Two
      USE SUPERINDEX, only: MTGEU, MTGTU
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT, NSYM, NASUP, NISUP, NINDEP,       &
     &                         NTGEUES, NTGTUES
      IMPLICIT None

! Presently symmetry blocking is disregarded, but index pair
! permutation symmetry is used.
! NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
      integer(kind=iwp), intent(in):: IVEC, JVEC, NOP2
      real(kind=wp), intent(inout):: OP2(NOP2)

      real(kind=wp), ALLOCATABLE, TARGET:: W1(:), W2_H(:)
      real(kind=wp), ALLOCATABLE:: WPROD(:)
      real(kind=wp), POINTER:: W2(:)
      integer(kind=iwp) ICASE, ISYM, NAS, NIS, MDVEC, NWPROD, ISCT,     &
     &                  IISTA, IIEND, NCOL, IW1, ITABS, IW2,            &
     &                  IWPROD, IXABS, ITU, ITUABS, ITX, ITY, IUABS,    &
     &                  IUX, IUY, IXY, IXYABS, IYABS, JTXUY, JTYUX
      real(kind=wp) W_PROD

! Given the coefficients for two excitation operators, available in
! vectors numbered IVEC and JVEC on file, use the blocks for
! excitation cases BVAT(+) and BVAT(-), i.e. cases 8 and 9, to
! construct the zero-, one-, and two-body
! expansions of the product (Op in IVEC conjugated)(Op in JVEC)
! as operating on the CASSCF space.
! Formulae used:
! For the F+ case (i.e. case 8)
! W1(tu,ab)(conj)*W2(xy,cd) = (dac*dbd)*(2 Etxuy + 2 Etyux)
! For the F- case (i.e. case 9)
! W1(tu,ab)(conj)*W2(xy,cd) = (dac*dbd)*(2 Etxuy - 2 Etyux)

! FIRST THE F+ i.e. BVAT+ i.e. CASE 8 -----------------------------
      ICASE=8
! Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
! Allocate space for one section of excitation amplitudes:
! Pick up a symmetry block of W1 and W2
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        IF(JVEC.EQ.IVEC) THEN
          W2=>W1
        ELSE
          CALL mma_allocate(W2_H,NAS*MDVEC,Label='W2_H')
          W2=>W2_H
        END IF
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
! Multiply WProd = (W1 sect )*(W2 sect transpose)
         CALL DGEMM_('N','T',                                           &
     &              NAS,NAS,NCOL,                                       &
     &              One,W1,NAS,                                         &
     &              W2,NAS,                                             &
     &              One,WPROD,NAS)
         END DO
! Deallocate W1 and W2
        CALL mma_deallocate(W1)
        IF(JVEC.NE.IVEC) CALL mma_deallocate(W2_H)
        W2=>Null()

! Loop over (TU)
        DO ITU=1,NAS
          IW1=ITU
          ITUABS=ITU+NTGEUES(ISYM)
          ITABS=MTGEU(1,ITUABS)
          IUABS=MTGEU(2,ITUABS)
! Loop over (XY)
          DO IXY=1,NAS
            IW2=IXY
            IXYABS=IXY+NTGEUES(ISYM)
            IXABS=MTGEU(1,IXYABS)
            IYABS=MTGEU(2,IXYABS)
            ITX=ITABS+NASHT*(IXABS-1)
            IUY=IUABS+NASHT*(IYABS-1)
            ITY=ITABS+NASHT*(IYABS-1)
            IUX=IUABS+NASHT*(IXABS-1)
            IWPROD=IW1+NAS*(IW2-1)
            W_PROD=WPROD(IWPROD)
! Remember: C For the F+ case (i.e. case 8)
! W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Etxuy + 2 Etyux)
! Contrib to 2-particle operator, from 2 Etxuy:
            IF(ITX.GE.IUY) THEN
              JTXUY=(ITX*(ITX-1))/2+IUY
            ELSE
              JTXUY=(IUY*(IUY-1))/2+ITX
            END IF
            OP2(JTXUY)=OP2(JTXUY)+Two*W_PROD
! Contrib to 2-particle operator, from 2 Etyux:
            IF(ITY.GE.IUX) THEN
              JTYUX=(ITY*(ITY-1))/2+IUX
            ELSE
              JTYUX=(IUX*(IUX-1))/2+ITY
            END IF
            OP2(JTYUX)=OP2(JTYUX)+Two*W_PROD
          END DO
        END DO
! Deallocate matrix product:
        CALL mma_deallocate(WPROD)
      END DO
! THEN THE F- i.e. BVAT- i.e. CASE 9 -----------------------------
      ICASE=9
! Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
! Allocate space for one section of excitation amplitudes:
! Pick up a symmetry block of W1 and W2
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        IF(JVEC.EQ.IVEC) THEN
          W2=>W1
        ELSE
          CALL mma_allocate(W2_H,NAS*MDVEC,Label='W2_H')
          W2=>W2_H
        END IF
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
! Multiply WProd = (W1 sect )*(W2 sect transpose)
         CALL DGEMM_('N','T',                                           &
     &              NAS,NAS,NCOL,                                       &
     &              One,W1,NAS,                                         &
     &              W2,NAS,                                             &
     &              One,WPROD,NAS)
        END DO
! Deallocate W1 and W2
        CALL mma_deallocate(W1)
        IF(JVEC.NE.IVEC) CALL mma_deallocate(W2_H)
        W2=>Null()

! Loop over (TU)
        DO ITU=1,NAS
          IW1=ITU
          ITUABS=ITU+NTGTUES(ISYM)
          ITABS=MTGTU(1,ITUABS)
          IUABS=MTGTU(2,ITUABS)
! Loop over (XY)
          DO IXY=1,NAS
            IW2=IXY
            IXYABS=IXY+NTGTUES(ISYM)
            IXABS=MTGTU(1,IXYABS)
            IYABS=MTGTU(2,IXYABS)
            ITX=ITABS+NASHT*(IXABS-1)
            IUY=IUABS+NASHT*(IYABS-1)
            ITY=ITABS+NASHT*(IYABS-1)
            IUX=IUABS+NASHT*(IXABS-1)
            IWPROD=IW1+NAS*(IW2-1)
            W_PROD=WPROD(IWPROD)
! Remember: C For the F- case (i.e. case 9)
! W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Etxuy - 2 Etyux)
! Contrib to 2-particle operator, from 2 Etxuy:
            IF(ITX.GE.IUY) THEN
              JTXUY=(ITX*(ITX-1))/2+IUY
            ELSE
              JTXUY=(IUY*(IUY-1))/2+ITX
            END IF
            OP2(JTXUY)=OP2(JTXUY)+Two*W_PROD
! Contrib to 2-particle operator, from -2 Etyux:
            IF(ITY.GE.IUX) THEN
              JTYUX=(ITY*(ITY-1))/2+IUX
            ELSE
              JTYUX=(IUX*(IUX-1))/2+ITY
            END IF
            OP2(JTYUX)=OP2(JTYUX)-Two*W_PROD
          END DO
        END DO
! Deallocate matrix product:
        CALL mma_deallocate(WPROD)
      END DO
      END SUBROUTINE MKWWOPF
