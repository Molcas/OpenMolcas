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
! Copyright (C) 1999, Jeppe Olsen                                      *
!***********************************************************************
      SUBROUTINE CON_BLOCKS(   IATP,   IBTP,   JATP,   JBTP,   IASM,    &
     &                         IBSM,   JASM,   JBSM,                    &
     &                      ICONSPA,ICONSPB, NOCTPA, NOCTPB,  MXEXC,    &
     &                      IH_OCC_CONS,INTERACT)
!
! Does CI blocks IATP IBTP interact with blocks JATP JBTP
!
!. Input
! ======
! IATP IBTP JATP JBTP : Supergroups, relative numbers
! IOCTPA, IOBTPB : Offset for type
! ICONSPA, ICONSPB : Connection matrices giving exciation
!                    level between two string types
! MXEXC : Largest excitation level
! IH_OCC_CONS : = 1 => Use only occupation conserving part of
!                     Hamiltonian
!
!. Output
!. INTERACT : =1 => The two blocks does interact
! Jeppe Olsen, April 99
!
      IMPLICIT NONE
      INTEGER IATP,   IBTP,   JATP,   JBTP,   IASM,                     &
     &        IBSM,   JASM,   JBSM,                                     &
     &        NOCTPA, NOCTPB,  MXEXC,                                   &
     &        IH_OCC_CONS,INTERACT
      INTEGER ICONSPA(NOCTPA,NOCTPA), ICONSPB(NOCTPB,NOCTPB)

      INTEGER IA_EXC,IB_EXC,NTEST
!
      IA_EXC = ICONSPA(IATP,JATP)
      IB_EXC = ICONSPB(IBTP,JBTP)
      IF(IH_OCC_CONS.EQ.0) THEN
!. Usual one- or two- electron operator
        IF(MXEXC.EQ.1) THEN
          IF((IA_EXC.LE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM).OR.        &
     &       (IB_EXC.LE.1.AND.IATP.EQ.JATP.AND.IASM.EQ.JASM)    )       &
     &        INTERACT = 1
        ELSE IF(MXEXC.EQ.2) THEN
          IF((IA_EXC.LE.1.AND.IB_EXC.LE.1)                  .OR.        &
     &       (IA_EXC.EQ.2.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM).OR.        &
     &       (IB_EXC.EQ.2.AND.IATP.EQ.JATP.AND.IASM.EQ.JASM)    )       &
     &        INTERACT = 1
        END IF
!      ELSE
!*. Orbital conserving part of  Hamiltonian
!        IF(IA_EXC.EQ.IB_EXC .AND. IB_EXC.LE.1 ) THEN
!          IATP_ABS = IATP + IOCTPA-1
!          IBTP_ABS = IBTP + IOCTPB-1
!          JATP_ABS = JATP + IOCTPA-1
!          JBTP_ABS = JBTP + IOCTPB-1
!*. Find Orb space where alpha strings differ
!          IPGAS = 0
!          IMGAS = 0
!          DO IGAS = 1, NGAS
!            IAEL = NELFSPGP(IGAS,IATP_ABS)
!            JAEL = NELFSPGP(IGAS,JATP_ABS)
!            IF(IAEL-JAEL.EQ.1) IPGAS = IGAS
!            IF(IAEL-JAEL.EQ.-1)IMGAS = IGAS
!          END DO
!          IF(IPGAS.NE.0) THEN
!            IPDIF = NELFSPGP(IPGAS,IBTP_ABS)-NELFSPGP(IPGAS,JBTP_ABS)
!          ELSE
!            IPDIF = 0
!          END IF
!*. corresponding differences in beta
!          IF(IMGAS.NE.0) THEN
!            IMDIF = NELFSPGP(IMGAS,IBTP_ABS)-NELFSPGP(IMGAS,JBTP_ABS)
!          ELSE
!            IMDIF = 0
!          END IF
!          IF(IPGAS.EQ.0.AND.IMGAS.EQ.0) INTERACT = 1
!          IF(IPGAS.NE.0.AND.IMGAS.NE.0) THEN
!            IF(IPDIF.EQ.-1.AND.IMDIF.EQ.1) INTERACT = 1
!          END IF
!        END IF
      END IF
!
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from CONBLOCKS '
        WRITE(6,*) ' IATP IBTP JATP JBTP ',IATP,IBTP,JATP,JBTP
        WRITE(6,*) ' IH_OCC_CONS, INTERACT = ', IH_OCC_CONS,INTERACT
      END IF
!
      END SUBROUTINE CON_BLOCKS
