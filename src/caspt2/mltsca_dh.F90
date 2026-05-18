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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
      SUBROUTINE MLTSCA_DH(IMLTOP,LST1,LST2,                            &
     &                     X,NXI,NXA,F,NFI,NFA,                         &
     &                     Y,NAS2,jYLo,jYHi)
      use definitions, only: iwp, wp
      use Sigma_data, only: NLST1, NLST2, VAL1, VAL2
      IMPLICIT None
      integer(kind=iwp), intent(in):: IMLTOP,NXI,NXA,NFI,NFA,NAS2,jYLo, &
     &                                jYHi
      real(kind=wp), intent(inout):: X(NXI,NXA),F(NFI,NFA),             &
     &                               Y(NAS2,jYHi-jYLo+1)
      integer(kind=iwp), intent(in)::  LST1(4,NLST1), LST2(4,NLST2)

      integer(kind=iwp) ILST1, ILST2, JY, L11, L12, L13, L14, L21, L23, &
     &                  L24, L22
      real(kind=wp) V1, V2
! Given two lists with entries LST1(4,ITEM), ITEM=1,NLST1, the
! four entries called L11,L12,L13,L14 for short, for a given
! item, and with V1=VAL1(L14), and similar for the other list,
! compute, for IMLTOP=0 or 1 respectively,
!     X(L11,L21) := Add V1*V2*F(L12,L22)*Y(L13,L23)
!  or Y(L13,L23) := Add V1*V2*F(L12,L22)*X(L11,L21)
! or for IMLTOP=2, compute
!     F(L12,L22) := Add V1*V2*X(L11,L21)*Y(L13,L23)

      IF(IMLTOP.EQ.0) THEN
        DO ILST1=1,NLST1
        L11=LST1(1,ILST1)
        L12=LST1(2,ILST1)
        L13=LST1(3,ILST1)
        L14=LST1(4,ILST1)
        V1=VAL1(L14)
        IF (L13.GE.jYLo .AND. L13.LE.jYHi) THEN
          JY=L13-jYLo+1
          DO ILST2=1,NLST2
          L21=LST2(1,ILST2)
          L22=LST2(2,ILST2)
          L23=LST2(3,ILST2)
          L24=LST2(4,ILST2)
          V2=VAL2(L24)
          X(L11,L21)=X(L11,L21)+V1*V2*F(L12,L22)*Y(L23,JY)
          END DO
        END IF
        END DO
      ELSE IF(IMLTOP.EQ.1) THEN
        DO ILST1=1,NLST1
        L11=LST1(1,ILST1)
        L12=LST1(2,ILST1)
        L13=LST1(3,ILST1)
        L14=LST1(4,ILST1)
        V1=VAL1(L14)
        IF (L13.GE.jYLo .AND. L13.LE.jYHi) THEN
          JY=L13-jYLo+1
          DO ILST2=1,NLST2
          L21=LST2(1,ILST2)
          L22=LST2(2,ILST2)
          L23=LST2(3,ILST2)
          L24=LST2(4,ILST2)
          V2=VAL2(L24)
          Y(L23,JY)=Y(L23,JY)+V1*V2*F(L12,L22)*X(L11,L21)
          END DO
        END IF
        END DO
      ELSE
        DO ILST1=1,NLST1
        L11=LST1(1,ILST1)
        L12=LST1(2,ILST1)
        L13=LST1(3,ILST1)
        L14=LST1(4,ILST1)
        V1=VAL1(L14)
        IF (L13.GE.jYLo .AND. L13.LE.jYHi) THEN
          JY=L13-jYLo+1
          DO ILST2=1,NLST2
          L21=LST2(1,ILST2)
          L22=LST2(2,ILST2)
          L23=LST2(3,ILST2)
          L24=LST2(4,ILST2)
          V2=VAL2(L24)
          F(L12,L22)=F(L12,L22)+V1*V2*X(L11,L21)*Y(L23,JY)
          END DO
        END IF
        END DO
      END IF

!     NFSCA=NFSCA+4*NLST1*NLST2
      END SUBROUTINE MLTSCA_DH
