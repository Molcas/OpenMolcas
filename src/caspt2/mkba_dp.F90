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
      SUBROUTINE MKBA_DP (DREF,NDREF,PREF,NPREF,FD,FP,iSYM,             &
     &                    BA,MBA,iLo,iHi,jLo,jHi,LDA)
      use definitions, only: iwp, wp
      use constants, only: Half, Two, Four
      USE SUPERINDEX, only: MTUV
      use caspt2_global, only:ipea_shift
      use caspt2_module, only: EASUM,NASHT,NTUVES,EPSA,NTUVES
      IMPLICIT None
      INTEGER(KIND=IWP), INTENT(IN):: NDREF, NPREF, iSYM,               &
     &                                MBA,iLo, iHi, jLo, jHi, LDA
      REAL(KIND=WP), INTENT(IN):: DREF(NDREF),PREF(NPREF)
      REAL(KIND=WP), INTENT(IN):: FD(NDREF),FP(NPREF)
      REAL(KIND=WP), INTENT(INOUT):: BA(MBA)

      INTEGER(KIND=IWP) IXYZ,IXYZABS,IXABS,IYABS,IZABS,ITUVABS,ITABS,   &
     &                  IUABS,IVABS,ISADR,IVZ,IXT,IP1,IP2,IP,ID,ID1,    &
     &                  ID2,IDT,IDU,IDV,ITUV,IVT,IVU,IXZ,IYZ
      REAL(KIND=WP)EX,EY,ET,EU,ETU,FACT,VALUE
!SV.20100831: fill in the F2 and F1 corrections for this BA block
! on entry, BA should contain SA!!
      DO IXYZ=jLo,jHi
        IXYZABS=IXYZ+NTUVES(ISYM)
        IXABS=MTUV(1,IXYZABS)
        IYABS=MTUV(2,IXYZABS)
        IZABS=MTUV(3,IXYZABS)
        EX=EPSA(IXABS)
        EY=EPSA(IYABS)
        DO ITUV=iLo,iHi
          ITUVABS=ITUV+NTUVES(ISYM)
          ITABS=MTUV(1,ITUVABS)
          IUABS=MTUV(2,ITUVABS)
          IVABS=MTUV(3,ITUVABS)
          ET=EPSA(ITABS)
          EU=EPSA(IUABS)
          ETU=ET+EU
          FACT=EY+EU+EX+ET-EASUM
          ISADR=1+iTUV-iLo+LDA*(iXYZ-jLo)
          IF (LDA.EQ.0) THEN
            IF (iXYZ.LE.iTUV) THEN
              ISADR=(ITUV*(ITUV-1))/2+IXYZ
            ELSE
              CYCLE
            END IF
          END IF
          VALUE=FACT*BA(ISADR)
          IF(IYABS.EQ.IUABS) THEN
            IVZ=IVABS+NASHT*(IZABS-1)
            IXT=IXABS+NASHT*(ITABS-1)
            IP1=MAX(IVZ,IXT)
            IP2=MIN(IVZ,IXT)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE-Two*(FP(IP)-EU*PREF(IP))
            IF(IXABS.EQ.ITABS) THEN
              ID1=MAX(IVABS,IZABS)
              ID2=MIN(IVABS,IZABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE+Two*(FD(ID)-ETU*DREF(ID))
            END IF
          END IF
! Add  dyt ( -Fvuxz + Et*Gvuxz +dxu (-Fvz+(Et+Eu)*Gvz))
          IF(IYABS.EQ.ITABS) THEN
            IVU=IVABS+NASHT*(IUABS-1)
            IXZ=IXABS+NASHT*(IZABS-1)
            IP1=MAX(IVU,IXZ)
            IP2=MIN(IVU,IXZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE-Two*(FP(IP)-ET*PREF(IP))
            IF(IXABS.EQ.IUABS) THEN
              ID1=MAX(IVABS,IZABS)
              ID2=MIN(IVABS,IZABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE - (FD(ID)-ETU*DREF(ID))
            END IF
          END IF
! Add  dxu ( -Fvtyz + Eu*Gvtyz )
          IF(IXABS.EQ.IUABS) THEN
            IVT=IVABS+NASHT*(ITABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IP1=MAX(IVT,IYZ)
            IP2=MIN(IVT,IYZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE-Two*(FP(IP)-EU*PREF(IP))
          END IF
! Add  2dtx ( Fvuyz-Et*Gvuyz )
          IF(ITABS.EQ.IXABS) THEN
            IVU=IVABS+NASHT*(IUABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IP1=MAX(IVU,IYZ)
            IP2=MIN(IVU,IYZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+Four*(FP(IP)-ET*PREF(IP))
          END IF
!GG.Nov03
          IF (ITUV.eq.IXYZ) THEN
            IDT=(ITABS*(ITABS+1))/2
            IDU=(IUABS*(IUABS+1))/2
            IDV=(IVABS*(IVABS+1))/2
            VALUE=VALUE+ipea_shift*Half*BA(ISADR)*                      &
     &                  (Two-DREF(IDV)+DREF(IDT)+DREF(IDU))
          ENDIF
!GG End
          BA(ISADR)=VALUE
        END DO
      END DO
      END SUBROUTINE MKBA_DP
