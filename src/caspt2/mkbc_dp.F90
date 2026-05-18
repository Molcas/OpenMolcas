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
      SUBROUTINE MKBC_DP (DREF,NDREF,PREF,NPREF,FD,FP,iSYM,             &
     &                    BC,MBC,iLo,iHi,jLo,jHi,LDC)
      use definitions, only: iwp, wp
      use constants, only: Half, Two, Four
      USE SUPERINDEX, only: MTUV
      use caspt2_global, only:ipea_shift
      use caspt2_module, only: EASUM,NASHT,NTUVES,EPSA
      IMPLICIT NONE
      INTEGER(KIND=IWP), intent(in):: NDREF,NPREF, iSYM,MBC,            &
     &                                iLo,iHi,jLo,jHi,LDC
      REAL(KIND=WP), intent(in):: DREF(NDREF),PREF(NPREF)
      REAL(KIND=WP), intent(in):: FD(NDREF),FP(NPREF)
      REAL(KIND=WP), intent(inout):: BC(MBC)

      INTEGER(KIND=IWP) IXYZ,IXYZABS,IXABS,IYABS,IZABS,ITUV,ITUVABS,    &
     &                  ITABS,IUABS,IVABS,ISADR,IVZ,ITX,IVU,ITZ,IP1,    &
     &                  IP2,IP,ID1,ID2,ID,IDT,IDU,IDV,IVX,IYZ
      REAL(KIND=WP) EY,EU, EYU, FACT, VALUE

      DO IXYZ=jLo,jHi
        IXYZABS=IXYZ+NTUVES(ISYM)
        IXABS=MTUV(1,IXYZABS)
        IYABS=MTUV(2,IXYZABS)
        IZABS=MTUV(3,IXYZABS)
        EY=EPSA(IYABS)
        DO ITUV=iLo,iHi
          ITUVABS=ITUV+NTUVES(ISYM)
          ITABS=MTUV(1,ITUVABS)
          IUABS=MTUV(2,ITUVABS)
          IVABS=MTUV(3,ITUVABS)
          EU=EPSA(IUABS)
          EYU=EY + EU
          FACT=EYU-EASUM
          ISADR=1+iTUV-iLo+LDC*(iXYZ-jLo)
          IF (LDC.EQ.0) THEN
            IF (iXYZ.LE.iTUV) THEN
              ISADR=(ITUV*(ITUV-1))/2+IXYZ
            ELSE
              CYCLE
            END IF
          END IF
          VALUE=FACT*BC(ISADR)
! VALUE= Fvutxyz + (EPSA(y)+EPSA(u))*SC(tuv,xyz)
! Add  dyu ( Fvztx - EPSA(u)*Gvztx )
          IF(IYABS.EQ.IUABS) THEN
            IVZ=IVABS+NASHT*(IZABS-1)
            ITX=ITABS+NASHT*(IXABS-1)
            IP1=MAX(IVZ,ITX)
            IP2=MIN(IVZ,ITX)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+Two*(FP(IP)-EU*PREF(IP))
          END IF
! Add  dyx ( Fvutz - EPSA(y)*Gvutz )
          IF(IYABS.EQ.IXABS) THEN
            IVU=IVABS+NASHT*(IUABS-1)
            ITZ=ITABS+NASHT*(IZABS-1)
            IP1=MAX(IVU,ITZ)
            IP2=MIN(IVU,ITZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+Two*(FP(IP)-EY*PREF(IP))
          END IF
! Add  dtu ( Fvxyz - EPSA(u)*Gvxyz + dyx Fvz -
!             (EPSA(u)+EPSA(y)*dyz Gvz)
          IF(ITABS.EQ.IUABS) THEN
            IVX=IVABS+NASHT*(IXABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IP1=MAX(IVX,IYZ)
            IP2=MIN(IVX,IYZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+Two*(FP(IP)-EU*PREF(IP))
            IF(IYABS.EQ.IXABS) THEN
              ID1=MAX(IVABS,IZABS)
              ID2=MIN(IVABS,IZABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE+FD(ID)-EYU*DREF(ID)
            END IF
          END IF
!GG.Nov03
          IF (ITUV.eq.IXYZ) THEN
            IDT=(ITABS*(ITABS+1))/2
            IDU=(IUABS*(IUABS+1))/2
            IDV=(IVABS*(IVABS+1))/2
            VALUE=VALUE+ipea_shift*Half*BC(ISADR)*                      &
     &                    (Four-DREF(IDT)-DREF(IDV)+DREF(IDU))
          ENDIF
!GG End
          BC(ISADR)=VALUE
        END DO
      END DO
      END SUBROUTINE MKBC_DP
