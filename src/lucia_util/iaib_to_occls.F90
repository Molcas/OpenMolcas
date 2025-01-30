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
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************
      SUBROUTINE IAIB_TO_OCCLS(IAGRP,IATP,IBGRP,IBTP,IOC)
      use strbas, only: IOCLS
      use lucia_data, only: NGAS,NMXOCCLS
      use lucia_data, only: IBSPGPFTP,NELFSPGP
      use lucia_data, only: MXPNGAS
!
! Find address of the occupation class corresponding to given types
! of alpha and beta strings
!
! Jeppe Olsen, December 2001
!
      Implicit NONE
      INTEGER IAGRP,IATP,IBGRP,IBTP,IOC
!. Local scratch
      INTEGER IABOCC(MXPNGAS)
      INTEGER IATP_ABS,IBTP_ABS,IONE,INUM
!
      IATP_ABS = IATP + IBSPGPFTP(IAGRP) - 1
      IBTP_ABS = IBTP + IBSPGPFTP(IBGRP) - 1
!?    WRITE(6,*) ' IATP, IBTP, IAGRP, IBGRP = ',
!?   &             IATP, IBTP, IAGRP, IBGRP
!?    WRITE(6,*) ' IATP_ABS, IBTP_ABS ', IATP_ABS, IBTP_ABS
!
!  IVCSUM(IA,IB,IC,IFACB,IFACC,NDIM)
      IONE = 1
      CALL IVCSUM(   IABOCC,                                            &
     &            NELFSPGP(1,IATP_ABS),                                 &
     &            NELFSPGP(1,IBTP_ABS),                                 &
     &                 IONE,                                            &
     &                 IONE,                                            &
!
     &                 NGAS)
!. And the address of this occupation class
      CALL CMP_IVEC_ILIST(IABOCC,IOCLS,NGAS,NMXOCCLS,INUM)
!
      IOC = INUM
!
      IF(INUM.EQ.0) THEN
        WRITE(6,*)                                                      &
     &  ' Combination of alpha and beta string not found as occ-class'
        WRITE(6,*) ' Occ of alpha, Occ of beta, Occ of alpha+beta '
        CALL IWRTMA(NELFSPGP(1,IATP_ABS),1,NGAS,1,NGAS)
        CALL IWRTMA(NELFSPGP(1,IBTP_ABS),1,NGAS,1,NGAS)
        CALL IWRTMA(IABOCC,1,NGAS,1,NGAS)
!        STOP
!     &  ' Combination of alpha and beta string not found as occ-class'
         CALL SYSABENDMSG('lucia_util/iaib_to_occls',                   &
     &                    'Internal error',' ')
      END IF
!
      END SUBROUTINE IAIB_TO_OCCLS
