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
      SUBROUTINE MKSF(PREF,NPREF)
      use definitions, only: iwp, wp
      use constants, only: Four
      USE SUPERINDEX, only: MTU,MTGEU,KTU,MTGEU,KTGTU
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDSMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM,NINDEP,NTU,NTUES,NASHT,NTGEU,       &
     &                         NTGEUES,NTGTU,NTGTUES
      IMPLICIT NONE

      INTEGER(kind=iwp), intent(in)::  NPREF
      REAL(kind=wp), intent(in)::  PREF(NPREF)

      REAL(kind=wp), ALLOCATABLE:: SF(:), SFP(:), SFM(:)
      INTEGER(kind=iwp) ISYM,NINP,NAS,NSF,ITU,ITUABS,ITABS,IUABS,IXY,   &
     &                  IXYABS,IXABS,IYABS,ISADR,ITX,IUY,IP1,IP2,IP,    &
     &                  IDISK,ISMADR,ISPADR,ITGEU,ITGEUABS,ITGTU,IXGEY, &
     &                  IXGEYABS,IXGTY,IYX,NASM,NASP,NSFM,NSFP
      REAL(kind=wp) VALUE,STUXY,STUYX
! Set up the matrices SFP(tu,xy) and SFM(tu,xy)
! Formulae used:
!    SF(tu,xy)= 4 Ptxuy
!    SFP(tu,xy)=SF(tu,xy)+SF(tu,yx)
!    SFM(tu,xy)=SF(tu,xy)-SF(tu,yx)



! Loop over superindex symmetry.
      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,8)
        IF(NINP.EQ.0) CYCLE
        NAS=NTU(ISYM)
        NSF=(NAS*(NAS+1))/2
        IF(NSF.GT.0) CALL mma_allocate(SF,NSF,Label='SF')
        DO ITU=1,NAS
          ITUABS=ITU+NTUES(ISYM)
          ITABS=MTU(1,ITUABS)
          IUABS=MTU(2,ITUABS)
          DO IXY=1,ITU
            IXYABS=IXY+NTUES(ISYM)
            IXABS=MTU(1,IXYABS)
            IYABS=MTU(2,IXYABS)
            ISADR=(ITU*(ITU-1))/2+IXY
            ITX=ITABS+NASHT*(IXABS-1)
            IUY=IUABS+NASHT*(IYABS-1)
            IP1=MAX(ITX,IUY)
            IP2=MIN(ITX,IUY)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=Four*PREF(IP)
            SF(ISADR)=VALUE
          END DO
        END DO
        NASP=NTGEU(ISYM)
        NSFP=(NASP*(NASP+1))/2
        IF(NSFP.GT.0) CALL mma_allocate(SFP,NSFP,Label='SFP')
        NASM=NTGTU(ISYM)
        NSFM=(NASM*(NASM+1))/2
        IF(NSFM.GT.0) THEN
          CALL mma_allocate(SFM,NSFM,Label='SFM')
        END IF
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
            STUXY=SF(ISADR)
            IF(ITU.GE.IYX) THEN
              ISADR=(ITU*(ITU-1))/2+IYX
            ELSE
              ISADR=(IYX*(IYX-1))/2+ITU
            END IF
            STUYX=SF(ISADR)
            ISPADR=(ITGEU*(ITGEU-1))/2+IXGEY
            SFP(ISPADR)=STUXY+STUYX
            IF(ITABS.EQ.IUABS) CYCLE
            IF(IXABS.EQ.IYABS) CYCLE
            ITGTU=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
            IXGTY=KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
            ISMADR=(ITGTU*(ITGTU-1))/2+IXGTY
            SFM(ISMADR)=STUXY-STUYX
          END DO
        END DO
        IF(NSF.GT.0) CALL mma_deallocate(SF)

! Write to disk
        IF(NSFP.GT.0.and.NINDEP(ISYM,8).GT.0) THEN
          IDISK=IDSMAT(ISYM,8)
          CALL DDAFILE(LUSBT,1,SFP,NSFP,IDISK)
          CALL mma_deallocate(SFP)
        END IF
        IF(NSFM.GT.0) THEN
          IF(NINDEP(ISYM,9).GT.0) THEN
           IDISK=IDSMAT(ISYM,9)
           CALL DDAFILE(LUSBT,1,SFM,NSFM,IDISK)
          END IF
          CALL mma_deallocate(SFM)
        END IF
      END DO

      END SUBROUTINE MKSF
