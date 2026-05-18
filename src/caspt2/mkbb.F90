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
      SUBROUTINE MKBB(DREF,NDREF,PREF,NPREF,FD,FP)
      use definitions, only: iwp, wp
      use constants, only: Half, Two, Four, Eight
      USE SUPERINDEX, only: MTU, MTGEU, KTU, KTGTU
      use caspt2_global, only:ipea_shift
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDSMAT,IDBMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM,NINDEP,NTU,NTUES,NASHT,EASUM,       &
     &                         EPSA,NTGEU,NTGEUES,NTGTU,NTGTUES
      IMPLICIT None

      INTEGER(KIND=IWP), INTENT(IN):: NDREF,NPREF
      REAL(KIND=WP), INTENT(IN):: DREF(NDREF),PREF(NPREF)
      REAL(KIND=WP), INTENT(IN):: FD(NDREF),FP(NPREF)

      REAL(KIND=WP), ALLOCATABLE:: BB(:), BBP(:), SP(:), SDP(:),        &
     &                             BBM(:), SM(:), SDM(:)

      INTEGER(KIND=IWP) ISYM,NINP,NAS,NBB,ITU,ITUABS,ITABS,IUABS,IXY,   &
     &                  IXYABS,IXABS,IYABS,IBADR,IXT,IYU,IP1,IP2,IP,I,  &
     &                  IBMADR,IBPADR,ID,ID1,ID2,IDIAG,IDISK,IDSM,IDSP, &
     &                  IDT,IDU,INSM,ITGEU,ITGEUABS,ITGTU,IXGEY,        &
     &                  IXGEYABS,IXGTY,IYX,NASM,NASP,NBBM,NBBP,NSM,NSP
      REAL(KIND=WP) ET,EU,EX,EY,ATUXY,ATUX,ATYU,ATUY,ATYX,BTUYX,        &
     &              BTUXY,VALUE
! Set up the matrices BBP(tu,xy) and BBM(tu,xy)
! Formulae used:
!    BB(tu,xy)= 2*( Fyuxt - (A-Et-Eu-Ex-Ey)*Gyuxt )
!      + 4*dxt ( (A-Et-Ey-Eu)*Dyu - Fyu)
!      + 4*dyu ( (A-Et-Ey-Ex)*Dxt - Fxt)
!      - 2*dyt ( (A-Et-Eu-Ex)*Dxu - Fxu)
!      - 2*dxu ( (A-Et-Eu-Ey)*Dyt - Fyt)
!      + 8*dxt*dyu (Et+Ey)
!      - 4*dxu*dyt (Et+Ex)
! where A= EASUM= sum over active w of (Ew*Dww).
!    BBP(tu,xy)=BB(tu,xy)+BB(tu,yx)
!    BBM(tu,xy)=BB(tu,xy)-BB(tu,yx)


! Loop over superindex symmetry.
      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,2)
        IF(NINP.EQ.0) CYCLE
        NAS=NTU(ISYM)
        NBB=(NAS*(NAS+1))/2
        IF(NBB.GT.0) THEN
          CALL mma_allocate(BB,NBB,Label='BB')
        END IF
        DO ITU=1,NAS
          ITUABS=ITU+NTUES(ISYM)
          ITABS=MTU(1,ITUABS)
          IUABS=MTU(2,ITUABS)
          ET=EPSA(ITABS)
          EU=EPSA(IUABS)
          DO IXY=1,ITU
            IXYABS=IXY+NTUES(ISYM)
            IXABS=MTU(1,IXYABS)
            IYABS=MTU(2,IXYABS)
            EX=EPSA(IXABS)
            EY=EPSA(IYABS)
            IBADR=(ITU*(ITU-1))/2+IXY
            IXT=IXABS+NASHT*(ITABS-1)
            IYU=IYABS+NASHT*(IUABS-1)
            IP1=MAX(IXT,IYU)
            IP2=MIN(IXT,IYU)
            IP=(IP1*(IP1-1))/2+IP2
            ATUXY=EASUM-ET-EU-EX-EY
            VALUE=Four*(FP(IP)-ATUXY*PREF(IP))
! Add  + 4*dxt ( (A-Et-Ey-Eu)*Dyu - Fyu)
            IF(IXABS.EQ.ITABS) THEN
              ID1=MAX(IYABS,IUABS)
              ID2=MIN(IYABS,IUABS)
              ID=(ID1*(ID1-1))/2+ID2
              ATYU=EASUM-ET-EY-EU
              VALUE=VALUE+Four*(ATYU*DREF(ID)-FD(ID))
! Add  + 8*dxt*dyu (Et+Ey)
              IF(IYABS.EQ.IUABS) THEN
                VALUE=VALUE+Eight*(ET+EY)
              END IF
            END IF
! Add  + 4*dyu ( (A-Et-Ey-Ex)*Dxt - Fxt)
            IF(IYABS.EQ.IUABS) THEN
              ID1=MAX(IXABS,ITABS)
              ID2=MIN(IXABS,ITABS)
              ID=(ID1*(ID1-1))/2+ID2
              ATYX=EASUM-ET-EY-EX
              VALUE=VALUE+Four*(ATYX*DREF(ID)-FD(ID))
            END IF
! Add  - 2*dyt ( (A-Et-Eu-Ex)*Dxu - Fxu)
            IF(IYABS.EQ.ITABS) THEN
              ID1=MAX(IXABS,IUABS)
              ID2=MIN(IXABS,IUABS)
              ID=(ID1*(ID1-1))/2+ID2
              ATUX=EASUM-ET-EU-EX
              VALUE=VALUE-Two*(ATUX*DREF(ID)-FD(ID))
! Add  - 4*dxu*dyt (Et+Ex)
              IF(IXABS.EQ.IUABS) THEN
                VALUE=VALUE-Four*(ET+EX)
              END IF
            END IF
! Add  - 2*dxu ( (A-Et-Eu-Ey)*Dyt - Fyt)
            IF(IXABS.EQ.IUABS) THEN
              ID1=MAX(IYABS,ITABS)
              ID2=MIN(IYABS,ITABS)
              ID=(ID1*(ID1-1))/2+ID2
              ATUY=EASUM-ET-EU-EY
              VALUE=VALUE-Two*(ATUY*DREF(ID)-FD(ID))
            END IF
            BB(IBADR)=VALUE
          END DO
        END DO
        NASP=NTGEU(ISYM)
        NBBP=(NASP*(NASP+1))/2
        IF(NBBP.GT.0) THEN
          CALL mma_allocate(BBP,NBBP,Label='BBP')
!GG.Nov03  Load in SDP the diagonal elements of SBP matrix:
          NSP=(NASP*(NASP+1))/2
          CALL mma_allocate(SP,NSP,Label='SP')
          CALL mma_allocate(SDP,NASP,Label='SDP')
          IDSP=IDSMAT(ISYM,2)
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
        NBBM=(NASM*(NASM+1))/2
        IF(NBBM.GT.0) THEN
          CALL mma_allocate(BBM,NBBM,Label='BBM')
!GG.Nov03  Load in SDM the diagonal elements of SBM matrix:
          NSM=(NASM*(NASM+1))/2
          CALL mma_allocate(SDM,NASM,Label='SDM')
          IF (NINDEP(ISYM,3)>0) THEN
             CALL mma_allocate(SM,NSM,Label='SM')
             IDSM=IDSMAT(ISYM,3)
             CALL DDAFILE(LUSBT,2,SM,NSM,IDSM)
             IDIAG=0
             DO I=1,NASM
               IDIAG=IDIAG+I
               SDM(I)=SM(IDIAG)
             END DO
             CALL mma_deallocate(SM)
          END IF
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
            BTUXY=BB(IBADR)
            IF(ITU.GE.IYX) THEN
              IBADR=(ITU*(ITU-1))/2+IYX
            ELSE
              IBADR=(IYX*(IYX-1))/2+ITU
            END IF
            BTUYX=BB(IBADR)
            IBPADR=(ITGEU*(ITGEU-1))/2+IXGEY
            BBP(IBPADR)=BTUXY+BTUYX
!GG.Nov03
            IF (ITGEU.eq.IXGEY) THEN
              IDT=(ITABS*(ITABS+1))/2
              IDU=(IUABS*(IUABS+1))/2
              BBP(IBPADR)=BBP(IBPADR)+ipea_shift*half*                  &
     &                          (DREF(IDT)+DREF(IDU))*SDP(ITGEU)
            ENDIF
!GG End
            IF (NINDEP(ISYM,3)<1) CYCLE
            IF(ITABS.EQ.IUABS) CYCLE
            IF(IXABS.EQ.IYABS) CYCLE
            ITGTU=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
            IXGTY=KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
            IBMADR=(ITGTU*(ITGTU-1))/2+IXGTY
            BBM(IBMADR)=BTUXY-BTUYX
!GG.Nov03
            IF (ITGEU.eq.IXGEY) THEN
              IDT=(ITABS*(ITABS+1))/2
              IDU=(IUABS*(IUABS+1))/2
              BBM(IBMADR)=BBM(IBMADR)+ipea_shift*half*                  &
     &                           (DREF(IDT)+DREF(IDU))*SDM(INSM)
              INSM=INSM+1
            ENDIF
!GG.End
          END DO
        END DO
        IF(NBB.GT.0) CALL mma_deallocate(BB)

! Write to disk, and save size and address.
        IF(NBBP.GT.0.and.NINDEP(ISYM,2).GT.0) THEN
          IDISK=IDBMAT(ISYM,2)
          CALL DDAFILE(LUSBT,1,BBP,NBBP,IDISK)
          CALL mma_deallocate(BBP)
!GG.Nov03 DisAlloc SDP
          CALL mma_deallocate(SDP)
!GG End
        END IF
        IF(NBBM.GT.0) THEN
          IF(NINDEP(ISYM,3).GT.0) THEN
            IDISK=IDBMAT(ISYM,3)
            CALL DDAFILE(LUSBT,1,BBM,NBBM,IDISK)
          END IF
          CALL mma_deallocate(BBM)
!GG.Nov03 DisAlloc SDM
          CALL mma_deallocate(SDM)
!GG End
        END IF
      END DO

      END SUBROUTINE MKBB
