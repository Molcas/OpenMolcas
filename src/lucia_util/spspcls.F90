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
! Copyright (C) 1997, Jeppe Olsen                                      *
!***********************************************************************
      SUBROUTINE SPSPCLS(ISPSPCLS,ICLS,NCLS)
!
! Obtain mapping a-supergroup X b-supergroup => class
!
! Classes are specified by ICLS
!
! Jeppe Olsen, Jan 97
!
      use lucia_data, only: NGAS
      use lucia_data, only: IPRDIA
      use lucia_data, only: IBSPGPFTP,ISPGPFTP,NELFGP
      use lucia_data, only: NOCTYP
      use lucia_data, only: MXPNGAS
      IMPLICIT None
!. Specific input
      INTEGER NCLS,ICLS(*)
!. OUtput
      INTEGER ISPSPCLS(*)

      INTEGER IATP,IBTP,NOCTPA,NOCTPB,IOCTPA,IOCTPB
!
      IATP = 1
      IBTP = 2
!
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
!
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)

      CALL SPSPCLS_GAS(  NOCTPA,                                        &
     &                   NOCTPB,                                        &
     &                 ISPGPFTP(1,IOCTPA),                              &
     &                 ISPGPFTP(1,IOCTPB),NELFGP,MXPNGAS,NGAS,ISPSPCLS, &
     &                     ICLS,    NCLS,  IPRDIA)
!
!
      END SUBROUTINE SPSPCLS
!
