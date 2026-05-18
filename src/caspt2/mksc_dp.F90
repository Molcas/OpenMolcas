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
      SUBROUTINE MKSC_DP (DREF,NDREF,PREF,NPREF,                        &
     &                    iSYM,SC,NSC,iLo,iHi,jLo,jHi,LDC)
! In parallel, this subroutine is called on a local chunk of memory
! and LDC is set. In serial, the whole array is passed but then the
! storage uses a triangular scheme, and the LDC passed is zero.
      use definitions, only: iwp, wp
      use constants, only: Two
      USE SUPERINDEX, only: MTUV
      use caspt2_module, only: NASHT, nTUVES
      IMPLICIT None
      integer(kind=iwp), intent(in) :: NDREF,NPREF,iSYM,NSC,            &
     &                                 iLo,iHi,jLo,jHi,LDC
      real(kind=wp), intent(in):: DREF(NDREF),PREF(NPREF)
      real(kind=wp), intent(inout):: SC(NSC)

      integer(kind=iwp) ISADR,IXYZ,IXYZABS,IXABS,IYABS,IZABS,ITUV,      &
     &                  ITUVABS,ITABS,IUABS,IVABS,IVU,IYZ,IP1,          &
     &                  IP2,IP,ID1,ID2,IVZ,ITX,ITZ,IVX
      real(kind=wp) VALUE

      ISADR=0
!-SVC20100831: fill in the G2 and G1 corrections for this SC block
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
          IF (LDC.NE.0) THEN
            VALUE=SC(1+iTUV-iLo+LDC*(iXYZ-jLo))
          ELSE
            IF (IXYZ.LE.ITUV) THEN
              ISADR=(ITUV*(ITUV-1))/2+IXYZ
              VALUE=SC(ISADR)
            ELSE
              CYCLE
            ENDIF
          END IF
! Add  dyu Gvztx
          IF(IYABS.EQ.IUABS) THEN
            IVZ=IVABS+NASHT*(IZABS-1)
            ITX=ITABS+NASHT*(IXABS-1)
            IP1=MAX(IVZ,ITX)
            IP2=MIN(IVZ,ITX)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+Two*PREF(IP)
          END IF
! Add  dyx Gvutz
          IF(IYABS.EQ.IXABS) THEN
            IVU=IVABS+NASHT*(IUABS-1)
            ITZ=ITABS+NASHT*(IZABS-1)
            IP1=MAX(IVU,ITZ)
            IP2=MIN(IVU,ITZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+Two*PREF(IP)
          END IF
! Add  dtu Gvxyz + dtu dyx Gvz
          IF(ITABS.EQ.IUABS) THEN
            IVX=IVABS+NASHT*(IXABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IP1=MAX(IVX,IYZ)
            IP2=MIN(IVX,IYZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+Two*PREF(IP)
            IF(IYABS.EQ.IXABS) THEN
              ID1=MAX(IVABS,IZABS)
              ID2=MIN(IVABS,IZABS)
              VALUE=VALUE+DREF((ID1*(ID1-1))/2+ID2)
            END IF
          END IF
          IF (LDC.NE.0) THEN
            SC(1+iTUV-iLo+LDC*(iXYZ-jLo))=VALUE
          ELSE
            SC(ISADR)=VALUE
          END IF
        END DO
      END DO

      END SUBROUTINE MKSC_DP
