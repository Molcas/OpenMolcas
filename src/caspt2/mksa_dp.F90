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
      SUBROUTINE MKSA_DP (DREF,NDREF,PREF,NPREF,                        &
     &                    iSYM,SA,NSA,iLo,iHi,jLo,jHi,LDA)
! In parallel, this subroutine is called on a local chunk of memory
! and LDA is set. In serial, the whole array is passed but then the
! storage uses a triangular scheme, and the LDA passed is zero.
      use definitions, only: iwp, wp
      use constants, only: Two, Four
      USE SUPERINDEX, only: MTUV
      use caspt2_module, only: NASHT, nTUVES
      IMPLICIT None
      integer(kind=iwp), intent(in):: NDREF,NPREF,iSYM,NSA,             &
     &                                iLo,iHi,jLo,jHi,LDA
      real(kind=wp), intent(in):: DREF(NDREF),PREF(NPREF)
      real(kind=wp), intent(inout):: SA(NSA)

      integer(kind=iwp) ISADR,IXYZ,IXYZABS,IXABS,IYABS,IZABS,ITUV,      &
     &                  ITUVABS,ITABS,IUABS,IVABS,IVU,IYZ,IP1,          &
     &                  IP2,IP,ID1,ID2,ID,IVT,IVZ,IXT,IXZ
      real(kind=wp) VALUE

      ISADR=0
!-SVC20100831: fill in the G2 and G1 corrections for SA
      DO IXYZ=jLo,jHi
        IXYZABS=IXYZ+NTUVES(ISYM)
        IXABS=MTUV(1,IXYZABS)
        IYABS=MTUV(2,IXYZABS)
        IZABS=MTUV(3,IXYZABS)
        DO ITUV=iLo,iHi
          ITUVABS=ITUV+NTUVES(ISYM)
          ITABS=MTUV(1,ITUVABS)
          IUABS=MTUV(2,ITUVABS)
          IVABS=MTUV(3,ITUVABS)
! Add  2 dtx Gvuyz + 2 dtx dyu Gvz
          IF (LDA.NE.0) THEN
            VALUE=SA(1+(iTUV-iLo)+LDA*(iXYZ-jLo))
          ELSE
            IF (IXYZ.LE.ITUV) THEN
              ISADR=(ITUV*(ITUV-1))/2+IXYZ
              VALUE=SA(ISADR)
            ELSE
              CYCLE
            ENDIF
          END IF
          IF(ITABS.EQ.IXABS) THEN
            IVU=IVABS+NASHT*(IUABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IP1=MAX(IVU,IYZ)
            IP2=MIN(IVU,IYZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+Four*PREF(IP)
            IF(IYABS.EQ.IUABS)THEN
              ID1=MAX(IVABS,IZABS)
              ID2=MIN(IVABS,IZABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE+Two*DREF(ID)
            END IF
          END IF
! Add  -dxu Gvtyz -dxu dyt Gvz
          IF(IXABS.EQ.IUABS) THEN
            IVT=IVABS+NASHT*(ITABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IP1=MAX(IVT,IYZ)
            IP2=MIN(IVT,IYZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE - Two*PREF(IP)
            IF(IYABS.EQ.ITABS)THEN
              ID1=MAX(IVABS,IZABS)
              ID2=MIN(IVABS,IZABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE - DREF(ID)
            END IF
          END IF
! Add  -dyt Gvuxz
          IF(IYABS.EQ.ITABS) THEN
            IVU=IVABS+NASHT*(IUABS-1)
            IXZ=IXABS+NASHT*(IZABS-1)
            IP1=MAX(IVU,IXZ)
            IP2=MIN(IVU,IXZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE - Two*PREF(IP)
          END IF
! Add -dyu Gvzxt
          IF(IYABS.EQ.IUABS) THEN
            IVZ=IVABS+NASHT*(IZABS-1)
            IXT=IXABS+NASHT*(ITABS-1)
            IP1=MAX(IVZ,IXT)
            IP2=MIN(IVZ,IXT)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE - Two*PREF(IP)
          END IF
          IF (LDA.NE.0) THEN
            SA(1+(iTUV-iLo)+LDA*(iXYZ-jLo))=VALUE
          ELSE
            SA(ISADR)=VALUE
          END IF
        END DO
      END DO

      END SUBROUTINE MKSA_DP
