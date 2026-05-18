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
!*******************************************************************************
! Case B (ICASE=2,3)
!*******************************************************************************
      SUBROUTINE MKSB(DREF,NDREF,PREF,NPREF)
      use definitions, only: iwp, wp
      use constants, only: Two, Four, Eight
      USE SUPERINDEX, only: MTU,MTGEU,KTU,KTGTU
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDSMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT,NSYM,NINDEP,NTU,NTUES,NTGEU,NTGTU, &
     &                         NTGEUES,NTGTUES
      IMPLICIT None

      INTEGER(kind=iwp), intent(in):: NDREF,NPREF
      REAL(kind=wp), intent(in)::  DREF(NDREF),PREF(NPREF)

      REAL(kind=wp), ALLOCATABLE:: SB(:), SBP(:), SBM(:)

      integer(kind=iwp) ISYM,NINP,NAS,NSB,ITUABS,ITABS,IUABS,IXY,IXYABS,&
     &                  IXABS,IYABS,ISADR,IXT,IYU,IP1,IP2,IP,ID1,ID2,ID,&
     &                  IDISK,ISMADR,ISPADR,ITGEU,ITGEUABS,ITGTU,ITU,   &
     &                  IXGEY,IXGEYABS,IXGTY,IYX,NASM,NASP,NSBM,NSBP
      REAL(kind=wp) VALUE,STUXY,STUYX

! Set up the matrices SBP(tu,xy) and SBM(tu,xy)
! Formulae used:
!    SB(tu,xy)=
!    = 4 Pxtyu -4dxt Dyu -4dyu Dxt +2dyt Dxu + 8 dxt dyu
!      -4dxu dyt + 2dxu Dyt
!    SBP(tu,xy)=SB(tu,xy)+SB(tu,yx)
!    SBM(tu,xy)=SB(tu,xy)-SB(tu,yx)


! Loop over superindex symmetry.
      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,2)
        IF(NINP.EQ.0) CYCLE
        NAS=NTU(ISYM)
        NSB=(NAS*(NAS+1))/2
        IF(NSB.GT.0) CALL mma_allocate(SB,NSB,Label='SB')
        DO ITU=1,NAS
          ITUABS=ITU+NTUES(ISYM)
          ITABS=MTU(1,ITUABS)
          IUABS=MTU(2,ITUABS)
          DO IXY=1,ITU
            IXYABS=IXY+NTUES(ISYM)
            IXABS=MTU(1,IXYABS)
            IYABS=MTU(2,IXYABS)
            ISADR=(ITU*(ITU-1))/2+IXY
            IXT=IXABS+NASHT*(ITABS-1)
            IYU=IYABS+NASHT*(IUABS-1)
            IP1=MAX(IXT,IYU)
            IP2=MIN(IXT,IYU)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=Four*PREF(IP)
! Add  -4 dxt Dyu + 8dxt dyu
            IF(IXABS.EQ.ITABS) THEN
              ID1=MAX(IYABS,IUABS)
              ID2=MIN(IYABS,IUABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE-Four*DREF(ID)
              IF(IYABS.EQ.IUABS) VALUE=VALUE+Eight
            END IF
! Add  -4 dyu Dxt
            IF(IYABS.EQ.IUABS) THEN
              ID1=MAX(IXABS,ITABS)
              ID2=MIN(IXABS,ITABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE-Four*DREF(ID)
            END IF
! Add  +2 dyt Dxu
            IF(IYABS.EQ.ITABS) THEN
              ID1=MAX(IXABS,IUABS)
              ID2=MIN(IXABS,IUABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE+Two*DREF(ID)
            END IF
! Add  -4dxu dyt + 2dxu Dyt
            IF(IXABS.EQ.IUABS) THEN
              ID1=MAX(IYABS,ITABS)
              ID2=MIN(IYABS,ITABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE+Two*DREF(ID)
              IF(IYABS.EQ.ITABS) VALUE=VALUE-Four
            END IF
            ISADR=(ITU*(ITU-1))/2+IXY
            SB(ISADR)=VALUE
          END DO
        END DO
        NASP=NTGEU(ISYM)
        NSBP=(NASP*(NASP+1))/2
        IF(NSBP.GT.0) CALL mma_allocate(SBP,NSBP,Label='SBP')
        NASM=NTGTU(ISYM)
        NSBM=(NASM*(NASM+1))/2
        IF(NSBM.GT.0) CALL mma_allocate(SBM,NSBM,Label='SBM')
        DO ITGEU=1,NASP
          ITGEUABS=ITGEU+NTGEUES(ISYM)
          ITABS=MTGEU(1,ITGEUABS)
          IUABS=MTGEU(2,ITGEUABS)
          ITU=KTU(ITABS,IUABS)-NTUES(ISYM)
          DO IXGEY=1,ITGEU
            IXGEYABS=IXGEY+NTGEUES(ISYM)
            IXABS=MTGEU(1,IXGEYABS)
            IYABS=MTGEU(2,IXGEYABS)
            IXY=KTU(IXABS,IYABS)-NTUES(ISYM)
            IYX=KTU(IYABS,IXABS)-NTUES(ISYM)
            IF(ITU.GE.IXY) THEN
              ISADR=(ITU*(ITU-1))/2+IXY
            ELSE
              ISADR=(IXY*(IXY-1))/2+ITU
            END IF
            STUXY=SB(ISADR)
            IF(ITU.GE.IYX) THEN
              ISADR=(ITU*(ITU-1))/2+IYX
            ELSE
              ISADR=(IYX*(IYX-1))/2+ITU
            END IF
            STUYX=SB(ISADR)
            ISPADR=(ITGEU*(ITGEU-1))/2+IXGEY
            SBP(ISPADR)=STUXY+STUYX
            IF(ITABS.EQ.IUABS) CYCLE
            IF(IXABS.EQ.IYABS) CYCLE
            ITGTU=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
            IXGTY=KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
            ISMADR=(ITGTU*(ITGTU-1))/2+IXGTY
            SBM(ISMADR)=STUXY-STUYX
          END DO
        END DO
        IF(NSB.GT.0) CALL mma_deallocate(SB)

! Write to disk, and save size and address.
        IF(NSBP.GT.0) THEN
          IDISK=IDSMAT(ISYM,2)
          CALL DDAFILE(LUSBT,1,SBP,NSBP,IDISK)
          CALL mma_deallocate(SBP)
        END IF
        IF(NSBM.GT.0) THEN
          IF(NINDEP(ISYM,3).GT.0) THEN
            IDISK=IDSMAT(ISYM,3)
            CALL DDAFILE(LUSBT,1,SBM,NSBM,IDISK)
          END IF
          CALL mma_deallocate(SBM)
        END IF
      END DO

      END SUBROUTINE MKSB
