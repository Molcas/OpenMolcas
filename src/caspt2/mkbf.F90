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
      SUBROUTINE MKBF(DREF,NDREF,PREF,NPREF,FP)
      use definitions, only: iwp, wp
      use constants, only: Half, Four
      USE SUPERINDEX, only: MTU, MTGEU,KTU,KTGTU
      use caspt2_global, only:ipea_shift
      use caspt2_global, only:LUSBT
      use EQSOLV, only: IDSMAT,IDBMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM,NINDEP,NTU,NTUES,NASHT,NTGEU,EASUM, &
     &                         NTGTU,NTGEUES,NTGTUES
      IMPLICIT NONE

      integer(kind=iwp), intent(in):: NDREF,NPREF
      REAL(KIND=WP) PREF(NPREF),FP(NPREF),DREF(NDREF)

      REAL(KIND=WP), ALLOCATABLE:: BF(:), BFP(:), SDP(:),               &
     &                             SP(:), BFM(:), SDM(:),               &
     &                             SM(:)
      INTEGER(KIND=IWP) ISYM,NINP,NAS,NBF,ITU,ITUABS,ITABS,IUABS,IXY,   &
     &                  IXYABS,IXABS,IYABS,IBADR,ITX,IUY,IP1,IP2,IP,    &
     &                  NASP,NBFP,NSP,IDSP,IDIAG,I,NASM,NBFM,NSM,IDSM,  &
     &                  IBMADR,IBPADR,IDISK,IDT,IDU,INSM,ITGEU,ITGEUABS,&
     &                  ITGTU,IXGEY,IXGEYABS,IXGTY,IYX
      REAL(KIND=WP) BTUXY,BTUYX

! Set up the matrices BFP(tu,xy) and BFM(tu,xy)
! Formulae used:
!    BF(tu,xy)= 2*(Ftxuy - EASUM*Gtxuy)
!    BFP(tu,xy)=BF(tu,xy)+BF(tu,yx)
!    BFM(tu,xy)=BF(tu,xy)-BF(tu,yx)


! Loop over superindex symmetry.
      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,8)
        IF(NINP.EQ.0) CYCLE
        NAS=NTU(ISYM)
        NBF=(NAS*(NAS+1))/2
        IF(NBF.GT.0) THEN
          CALL mma_allocate(BF,NBF,LABEL='BF')
        END IF
        DO ITU=1,NAS
          ITUABS=ITU+NTUES(ISYM)
          ITABS=MTU(1,ITUABS)
          IUABS=MTU(2,ITUABS)
          DO IXY=1,ITU
            IXYABS=IXY+NTUES(ISYM)
            IXABS=MTU(1,IXYABS)
            IYABS=MTU(2,IXYABS)
            IBADR=(ITU*(ITU-1))/2+IXY
            ITX=ITABS+NASHT*(IXABS-1)
            IUY=IUABS+NASHT*(IYABS-1)
            IP1=MAX(ITX,IUY)
            IP2=MIN(ITX,IUY)
            IP=(IP1*(IP1-1))/2+IP2
            BF(IBADR)=Four*(FP(IP)-EASUM*PREF(IP))
          END DO
        END DO
        NASP=NTGEU(ISYM)
        NBFP=(NASP*(NASP+1))/2
        IF(NBFP.GT.0) THEN
          CALL mma_allocate(BFP,NBFP,Label='BFP')
!GG.Nov03  Load in SDP the diagonal elements of SFP matrix:
          NSP=(NASP*(NASP+1))/2
          CALL mma_allocate(SP,NSP,Label='SP')
          CALL mma_allocate(SDP,NASP,Label='SDP')
          IDSP=IDSMAT(ISYM,8)
          CALL DDAFILE(LUSBT,2,SP,NSP,IDSP)
          IDIAG=0
          DO I=1,NASP
            IDIAG=IDIAG+I
            SDP(I)=SP(IDIAG)
          END DO
          CALL mma_deallocate(SP)
!GG End
        END IF
        NASM=NTGTU(ISYM)
        NBFM=(NASM*(NASM+1))/2
        IF(NBFM.GT.0) THEN
          CALL mma_allocate(BFM,NBFM,Label='BFM')
!GG.Nov03  Load in SDM the diagonal elements of SFM matrix:
          NSM=(NASM*(NASM+1))/2
          CALL mma_allocate(SM,NSM,Label='SM')
          CALL mma_allocate(SDM,NASM,Label='SDM')
          IDSM=IDSMAT(ISYM,9)
          CALL DDAFILE(LUSBT,2,SM,NSM,IDSM)
          IDIAG=0
          DO I=1,NASM
            IDIAG=IDIAG+I
            SDM(I)=SM(IDIAG)
          END DO
          CALL mma_deallocate(SM)
!GG End
        END IF
        INSM=1
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
              IBADR=(ITU*(ITU-1))/2+IXY
            ELSE
              IBADR=(IXY*(IXY-1))/2+ITU
            END IF
            BTUXY=BF(IBADR)
            IF(ITU.GE.IYX) THEN
              IBADR=(ITU*(ITU-1))/2+IYX
            ELSE
              IBADR=(IYX*(IYX-1))/2+ITU
            END IF
            BTUYX=BF(IBADR)
            IBPADR=(ITGEU*(ITGEU-1))/2+IXGEY
            BFP(IBPADR)=BTUXY+BTUYX
!GG.Nov03
            IF (ITGEU.eq.IXGEY) THEN
              IDT=(ITABS*(ITABS+1))/2
              IDU=(IUABS*(IUABS+1))/2
              BFP(IBPADR)=BFP(IBPADR)+ipea_shift*Half*                  &
     &                    (Four-DREF(IDT)-DREF(IDU))*SDP(ITGEU)
            ENDIF
!GG End
            IF(ITABS.EQ.IUABS) CYCLE
            IF(IXABS.EQ.IYABS) CYCLE
            ITGTU=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
            IXGTY=KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
            IBMADR=(ITGTU*(ITGTU-1))/2+IXGTY
            BFM(IBMADR)=BTUXY-BTUYX
!GG.Nov03
            IF (ITGEU.eq.IXGEY) THEN
              IDT=(ITABS*(ITABS+1))/2
              IDU=(IUABS*(IUABS+1))/2
              BFM(IBMADR)=BFM(IBMADR)+ipea_shift*Half*                  &
     &                    (Four-DREF(IDT)-DREF(IDU))*SDM(INSM)
              INSM=INSM+1
            ENDIF

          END DO
        END DO
        IF(NBF.GT.0) CALL mma_deallocate(BF)

! Write to disk
        IF(NBFP.GT.0.and.NINDEP(ISYM,8).GT.0) THEN
          IDISK=IDBMAT(ISYM,8)
          CALL DDAFILE(LUSBT,1,BFP,NBFP,IDISK)
          CALL mma_deallocate(BFP)
!GG.Nov03 DisAlloc SDP
          CALL mma_deallocate(SDP)
!GG End
        END IF
        IF(NBFM.GT.0) THEN
         IF(NINDEP(ISYM,9).GT.0) THEN
          IDISK=IDBMAT(ISYM,9)
          CALL DDAFILE(LUSBT,1,BFM,NBFM,IDISK)
         END IF
         CALL mma_deallocate(BFM)
!GG.Nov03 DisAlloc SDM
         CALL mma_deallocate(SDM)
!GG End
        END IF
      END DO

      END SUBROUTINE MKBF
