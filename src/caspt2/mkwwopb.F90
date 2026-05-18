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
      SUBROUTINE MKWWOPB(IVEC,JVEC,OP0,OP1,NOP2,OP2)
      use definitions, only: iwp, wp
      use Constants, only: Zero, One, Two, Four, Six, Twelve
      USE SUPERINDEX, only: MTGEU, MTGTU
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT, NSYM, NASUP, NISUP, NINDEP,       &
     &                         NTGEUES, NTGTUES
      IMPLICIT None

! Presently symmetry blocking is disregarded, but index pair
! permutation symmetry is used.
! NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
      integer(kind=iwp), Intent(in):: IVEC, JVEC, NOP2
      real(kind=wp), intent(inout):: OP0, OP1(NASHT,NASHT),OP2(NOP2)

      real(kind=wp), Allocatable, TARGET:: W1(:), W2_H(:), WPROD(:)
      real(kind=wp), POINTER::  W2(:)
      integer(kind=iwp) ICASE, IIEND, IISTA, ISCT, ISYM, ITABS, ITU,    &
     &                  ITUABS, IUABS, IW1, IW2, IWPROD, IXABS, IXT,    &
     &                  IXU, IXY, IXYABS, IYABS, IYT, IYU, JXTYU,       &
     &                  JYTXU, MDVEC, NAS, NCOL, NIS, NWPROD
      real(kind=wp) W_PROD
! Given the coefficients for two excitation operators, available in
! vectors numbered IVEC and JVEC on file, use the blocks for
! excitation cases VJTI(+) and VJTI(-), i.e. cases 2 and 3, to
! construct the zero-, one-, and two-body
! expansions of the product (Op in IVEC conjugated)(Op in JVEC)
! as operating on the CASSCF space.
! Formulae used:
! For the B+ case (i.e. case 2)
! W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Extyu + 2 Eytxu -2dxt Eyu
!           -2dyu Ext -2dyt Exu -2dxu Eyt + 4 dxt dyu + 4 dxu dyt)
! For the B- case (i.e. case 3)
! W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Extyu - 2 Eytxu -6dxt Eyu
!           -6dyu Ext +6dyt Exu +6dxu Eyt +12 dxt dyu -12 dxu dyt)

! FIRST THE B+ i.e. VJTI+ i.e. CASE 2 -----------------------------
      ICASE=2
! Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
! Allocate space for one section of excitation amplitudes:
! Pick up a symmetry block of W1 and W2
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        IF(IVEC.EQ.JVEC) THEN
          W2=>W1
        ELSE
          CALL mma_allocate(W2_H,NAS*MDVEC,Label='W2_H')
          W2=>W2_H
        END IF
        NWPROD=NAS**2
! Allocate space for the contraction:
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
        WPROD(:)=Zero
! Loop over sections:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA+MDVEC-1,MDVEC)
         NCOL=IIEND-IISTA+1
        CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1,NAS*MDVEC)
        IF (IVEC.NE.JVEC) CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2,NAS*MDVEC)
! Multiply WProd = (W1 sect )*(W2 sect transpose)
        CALL DGEMM_('N','T',                                            &
     &              NAS,NAS,NCOL,                                       &
     &              One,W1,NAS,                                         &
     &              W2,NAS,                                             &
     &              One,WPROD,NAS)
        END DO
! Deallocate W1 and W2
        CALL mma_deallocate(W1)
        IF(IVEC.NE.JVEC) CALL mma_deallocate(W2_H)
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
            IXT=IXABS+NASHT*(ITABS-1)
            IYU=IYABS+NASHT*(IUABS-1)
            IYT=IYABS+NASHT*(ITABS-1)
            IXU=IXABS+NASHT*(IUABS-1)
            IWPROD=IW1+NAS*(IW2-1)
            W_PROD=WPROD(IWPROD)
! Remember:
! W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Extyu + 2 Eytxu -2dxt Eyu
!           -2dyu Ext -2dyt Exu -2dxu Eyt + 4 dxt dyu + 4 dxu dyt)
! Contrib to 2-particle operator, from 2 Extyu:
            IF(IXT.GE.IYU) THEN
              JXTYU=(IXT*(IXT-1))/2+IYU
            ELSE
              JXTYU=(IYU*(IYU-1))/2+IXT
            END IF
            OP2(JXTYU)=OP2(JXTYU)+Two*W_PROD
! Contrib to 1-particle operator, from -2dxt Eyu
            IF(IXABS.EQ.ITABS) THEN
              OP1(IYABS,IUABS)=OP1(IYABS,IUABS)-Two*W_PROD
            END IF
! Contrib to 1-particle operator, from -2dyu Ext
            IF(IYABS.EQ.IUABS) THEN
              OP1(IXABS,ITABS)=OP1(IXABS,ITABS)-Two*W_PROD
! Contrib to 0-particle operator, from +4 dxt dyu
              IF(IXABS.EQ.ITABS) OP0=OP0 + Four*W_PROD
            END IF
! Contrib to 2-particle operator, from 2 Eytxu:
            IF(IYT.GT.IXU) THEN
              JYTXU=(IYT*(IYT-1))/2+IXU
            ELSE
              JYTXU=(IXU*(IXU-1))/2+IYT
            END IF
            OP2(JYTXU)=OP2(JYTXU)+Two*W_PROD
! Contrib to 1-particle operator, from -2dyt Exu
            IF(IYABS.EQ.ITABS) THEN
              OP1(IXABS,IUABS)=OP1(IXABS,IUABS)-Two*W_PROD
            END IF
! Contrib to 1-particle operator, from -2dxu Eyt
            IF(IXABS.EQ.IUABS) THEN
              OP1(IYABS,ITABS)=OP1(IYABS,ITABS)-Two*W_PROD
! Contrib to 0-particle operator, from +4 dyt dxu
              IF(IYABS.EQ.ITABS) OP0=OP0 + Four*W_PROD
            END IF
          END DO
        END DO
! Deallocate matrix product:
        CALL mma_deallocate(WPROD)
      END DO
! Then THE B- i.e. VJTI- i.e. CASE 3 -----------------------------
      ICASE=3
! Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
! Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        CALL mma_allocate(W2_H,NAS*MDVEC,Label='W2_H')
        W2=>W2_H
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
        CALL DGEMM_('N','T',                                            &
     &              NAS,NAS,NCOL,                                       &
     &              One,W1,NAS,                                         &
     &              W2,NAS,                                             &
     &              One,WPROD,NAS)
        END DO
! Deallocate W1, W2
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2_H)
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
            IXT=IXABS+NASHT*(ITABS-1)
            IYU=IYABS+NASHT*(IUABS-1)
            IYT=IYABS+NASHT*(ITABS-1)
            IXU=IXABS+NASHT*(IUABS-1)
            IWPROD=IW1+NAS*(IW2-1)
            W_PROD=WPROD(IWPROD)
! Remember:
! W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Extyu - 2 Eytxu -6dxt Eyu
!           -6dyu Ext +6dyt Exu +6dxu Eyt +12 dxt dyu -12 dxu dyt)
! Contrib to 2-particle operator, from 2 Extyu:
            IF(IXT.GE.IYU) THEN
              JXTYU=(IXT*(IXT-1))/2+IYU
            ELSE
              JXTYU=(IYU*(IYU-1))/2+IXT
            END IF
            OP2(JXTYU)=OP2(JXTYU)+Two*W_PROD
! Contrib to 1-particle operator, from -6dxt Eyu
            IF(IXABS.EQ.ITABS) THEN
              OP1(IYABS,IUABS)=OP1(IYABS,IUABS)-Six*W_PROD
            END IF
! Contrib to 1-particle operator, from -6dyu Ext
            IF(IYABS.EQ.IUABS) THEN
              OP1(IXABS,ITABS)=OP1(IXABS,ITABS)-Six*W_PROD
! Contrib to 0-particle operator, from +12 dxt dyu
              IF(IXABS.EQ.ITABS) OP0=OP0 + Twelve*W_PROD
            END IF
! Contrib to 2-particle operator, from -2 Eytxu:
            IF(IYT.GE.IXU) THEN
              JYTXU=(IYT*(IYT-1))/2+IXU
            ELSE
              JYTXU=(IXU*(IXU-1))/2+IYT
            END IF
            OP2(JYTXU)=OP2(JYTXU)-Two*W_PROD
! Contrib to 1-particle operator, from +6dyt Exu
            IF(IYABS.EQ.ITABS) THEN
              OP1(IXABS,IUABS)=OP1(IXABS,IUABS)+Six*W_PROD
            END IF
! Contrib to 1-particle operator, from +6dxu Eyt
            IF(IXABS.EQ.IUABS) THEN
              OP1(IYABS,ITABS)=OP1(IYABS,ITABS)+Six*W_PROD
! Contrib to 0-particle operator, from -12 dyt dxu
              IF(IYABS.EQ.ITABS) OP0=OP0 -Twelve*W_PROD
            END IF
          END DO
        END DO
! Deallocate matrix product
        CALL mma_deallocate(WPROD)
      END DO
      END SUBROUTINE MKWWOPB
