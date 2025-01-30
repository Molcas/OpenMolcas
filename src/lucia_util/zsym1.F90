!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE ZSYM1(NIRREP,IPRNT)
      use symmetry_info, only: SYMPRO => Mul
      use lucia_data, only: MXPOBS
      use csm_data, only: NSMSX,NSMDX,NSMST,NSMCI,NSMXT,ITSSX,ITSDX,    &
     &                    ITSXT
      use csm_data, only: ADASX,ADSXA,ASXAD,SXDXSX,SXSXDX
      Implicit None
      INTEGER NIRREP,IPRNT
!
      NSMSX = NIRREP
      NSMDX = NIRREP
      NSMST = NIRREP
      NSMCI = NIRREP
      NSMXT = NIRREP
      ITSSX = 1
      ITSDX = 1
      ITSXT = 1

!
      CALL ICPMT2(   SYMPRO,    ADASX,        8,        8,   MXPOBS,    &
     &               MXPOBS,        1)
      CALL ICPMT2(   SYMPRO,    ADSXA,        8,        8,   MXPOBS,    &
     &             2*MXPOBS,        1)
      CALL ICPMT2(   SYMPRO,    ASXAD,        8,        8,   MXPOBS,    &
     &             2*MXPOBS,        1)
      CALL ICPMT2(   SYMPRO,   SXSXDX,        8,        8, 2*MXPOBS,    &
     &             2*MXPOBS,        1)
      CALL ICPMT2(   SYMPRO,   SXDXSX,        8,        8, 2*MXPOBS,    &
     &             4*MXPOBS,        1)
!
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(IPRNT)
      END SUBROUTINE ZSYM1
