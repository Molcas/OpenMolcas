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
      SUBROUTINE MKBD(DREF,NDREF,PREF,NPREF,FD,FP)
      use definitions, only: iwp, wp
      use constants, only: Half, Two, Four
      USE SUPERINDEX, only: MTU
      use caspt2_global, only:ipea_shift
      use caspt2_global, only:LUSBT
      use EQSOLV, only: IDSMAT, IDBMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM,NINDEP,NTU,EASUM,NASHT,NTUES,EPSA
      IMPLICIT NONE

      INTEGER(KIND=IWP), INTENT(IN):: NDREF,NPREF
      REAL(KIND=WP), INTENT(IN):: DREF(NDREF),PREF(NPREF)
      REAL(KIND=WP), INTENT(IN):: FD(NDREF),FP(NPREF)

      REAL(KIND=WP), ALLOCATABLE:: BD(:), S(:), SD(:)
      INTEGER(KIND=IWP) ISYM,NIN,NAS,NBD,NS2,NAS2,IDS,IDIAG,I,ITU,ITU2, &
     &                  ITUABS,ITABS,IUABS,IXY,IXY2,IXYABS,IXABS,IYABS, &
     &                  IB11,IB21,IB12,IB22,IUTP,IXYP,IP1,IP2,IP,ID,ID1,&
     &                  ID2,IDISK,IDT,IDU,IUYP,IXTP
      REAL(KIND=WP) ET,EX,B11,B22,ETX,FUY,DUY

! Set up the matrix BD(tuP,xyQ),P and Q are 1 or 2,
! Formulae used:
!    BD(tu1,xy1)=
!      2*Futxy + 2*(Ex+Et-A)*Gutxy + 2*dxt (Fuy + (Et-A)*Duy)
!    BD(tu2,xy1)= -BD(tu1,xy1)/2
!    BD(tu2,xy2)=
!       -Fxtuy - (Ex+Et-A)*Gxtuy + 2*dxt (Fuy + (Ex-A)*Duy)
! where A=EASUM=Sum(w) of (Ew*Dww)


! Loop over superindex symmetry.
      DO ISYM=1,NSYM
        NIN=NINDEP(ISYM,5)
        IF(NIN.EQ.0) CYCLE
        NAS=NTU(ISYM)
        NBD=(2*NAS*(2*NAS+1))/2
        IF(NBD.GT.0) THEN
          CALL mma_allocate(BD,NBD,Label='BD')
!GG.Nov03  Load in SD the diagonal elements of SD matrix:
          NS2=(2*NAS*(2*NAS+1))/2
          NAS2=2*NAS
          CALL mma_allocate(S,NS2,Label='S')
          CALL mma_allocate(SD,NAS2,Label='SD')
          IDS=IDSMAT(ISYM,5)
          CALL DDAFILE(LUSBT,2,S,NS2,IDS)
          IDIAG=0
          DO I=1,NAS2
            IDIAG=IDIAG+I
            SD(I)=S(IDIAG)
          END DO
          CALL mma_deallocate(S)
!GG End
        END IF
        DO ITU=1,NAS
          ITU2=ITU+NAS
          ITUABS=ITU+NTUES(ISYM)
          ITABS=MTU(1,ITUABS)
          IUABS=MTU(2,ITUABS)
          ET=EPSA(ITABS)
          DO IXY=1,ITU
            IXY2=IXY+NAS
            IXYABS=IXY+NTUES(ISYM)
            IXABS=MTU(1,IXYABS)
            IYABS=MTU(2,IXYABS)
            EX=EPSA(IXABS)
            IB11=(ITU*(ITU-1))/2+IXY
            IB21=(ITU2*(ITU2-1))/2+IXY
            IB12=(IXY2*(IXY2-1))/2+ITU
            IB22=(ITU2*(ITU2-1))/2+IXY2
            IUTP=IUABS+NASHT*(ITABS-1)
            IXYP=IXABS+NASHT*(IYABS-1)
            IP1=MAX(IUTP,IXYP)
            IP2=MIN(IUTP,IXYP)
            IP=(IP1*(IP1-1))/2+IP2
            ETX=ET+EX
            B11=Four*(FP(IP)+(ETX-EASUM)*PREF(IP))
            IXTP=IXABS+NASHT*(ITABS-1)
            IUYP=IUABS+NASHT*(IYABS-1)
            IP1=MAX(IXTP,IUYP)
            IP2=MIN(IXTP,IUYP)
            IP=(IP1*(IP1-1))/2+IP2
            B22=-Two*(FP(IP)+(ETX-EASUM)*PREF(IP))
            IF(IXABS.EQ.ITABS) THEN
              ID1=MAX(IUABS,IYABS)
              ID2=MIN(IUABS,IYABS)
              ID=(ID1*(ID1-1))/2+ID2
              FUY=FD(ID)
              DUY=DREF(ID)
              B11=B11+Two*(FUY+(ET-EASUM)*DUY)
              B22=B22+Two*(FUY+(EX-EASUM)*DUY)
            END IF
            BD(IB11)= B11
            BD(IB21)=-Half*B11
            BD(IB12)=-Half*B11
            BD(IB22)= B22
!GG.Nov03
            IF (ITU.eq.IXY) THEN
              IDT=(ITABS*(ITABS+1))/2
              IDU=(IUABS*(IUABS+1))/2
              BD(IB11)=BD(IB11)+ipea_shift*Half*                        &
     &                       (Two-DREF(IDU)+DREF(IDT))*SD(ITU)
              BD(IB22)=BD(IB22)+ipea_shift*Half*                        &
     &                   (Two-DREF(IDU)+DREF(IDT))*SD(ITU+NAS)
            ENDIF
!GG End
          END DO
        END DO

! Write to disk
        IF(NBD.GT.0.and.NINDEP(ISYM,5).GT.0) THEN
          IDISK=IDBMAT(ISYM,5)
          CALL DDAFILE(LUSBT,1,BD,NBD,IDISK)
          CALL mma_deallocate(BD)
!GG.Nov03 DisAlloc SD
          CALL mma_deallocate(SD)
!GG End
        END IF
      END DO


      END SUBROUTINE MKBD
