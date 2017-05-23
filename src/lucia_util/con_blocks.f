************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1999, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE CON_BLOCKS(   IATP,   IBTP,   JATP,   JBTP,   IASM,
     &                         IBSM,   JASM,   JBSM,
     &                      ICONSPA,ICONSPB, NOCTPA, NOCTPB,  MXEXC,
     &                      IH_OCC_CONS,INTERACT)
*
* Does CI blocks IATP IBTP interact with blocks JATP JBTP
*
*. Input
* ======
* IATP IBTP JATP JBTP : Supergroups, relative numbers
* IOCTPA, IOBTPB : Offset for type
* ICONSPA, ICONSPB : Connection matrices giving exciation
*                    level between two string types
* MXEXC : Largest excitation level
* IH_OCC_CONS : = 1 => Use only occupation conserving part of
*                     Hamiltonian
*
*. Output
*. INTERACT : =1 => The two blocks does interact
* Jeppe Olsen, April 99
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "gasstr.fh"
#include "cgas.fh"
      INTEGER ICONSPA(NOCTPA,NOCTPA), ICONSPB(NOCTPB,NOCTPB)
*
      IA_EXC = ICONSPA(IATP,JATP)
      IB_EXC = ICONSPB(IBTP,JBTP)
      IF(IH_OCC_CONS.EQ.0) THEN
*. Usual one- or two- electron operator
        IF(MXEXC.EQ.1) THEN
          IF((IA_EXC.LE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM).OR.
     &       (IB_EXC.LE.1.AND.IATP.EQ.JATP.AND.IASM.EQ.JASM)    )
     &        INTERACT = 1
        ELSE IF(MXEXC.EQ.2) THEN
          IF((IA_EXC.LE.1.AND.IB_EXC.LE.1)                  .OR.
     &       (IA_EXC.EQ.2.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM).OR.
     &       (IB_EXC.EQ.2.AND.IATP.EQ.JATP.AND.IASM.EQ.JASM)    )
     &        INTERACT = 1
        END IF
c      ELSE
c*. Orbital conserving part of  Hamiltonian
c        IF(IA_EXC.EQ.IB_EXC .AND. IB_EXC.LE.1 ) THEN
c          IATP_ABS = IATP + IOCTPA-1
c          IBTP_ABS = IBTP + IOCTPB-1
c          JATP_ABS = JATP + IOCTPA-1
c          JBTP_ABS = JBTP + IOCTPB-1
c*. Find Orb space where alpha strings differ
c          IPGAS = 0
c          IMGAS = 0
c          DO IGAS = 1, NGAS
c            IAEL = NELFSPGP(IGAS,IATP_ABS)
c            JAEL = NELFSPGP(IGAS,JATP_ABS)
c            IF(IAEL-JAEL.EQ.1) IPGAS = IGAS
c            IF(IAEL-JAEL.EQ.-1)IMGAS = IGAS
c          END DO
c          IF(IPGAS.NE.0) THEN
c            IPDIF = NELFSPGP(IPGAS,IBTP_ABS)-NELFSPGP(IPGAS,JBTP_ABS)
c          ELSE
c            IPDIF = 0
c          END IF
c*. corresponding differences in beta
c          IF(IMGAS.NE.0) THEN
c            IMDIF = NELFSPGP(IMGAS,IBTP_ABS)-NELFSPGP(IMGAS,JBTP_ABS)
c          ELSE
c            IMDIF = 0
c          END IF
c          IF(IPGAS.EQ.0.AND.IMGAS.EQ.0) INTERACT = 1
c          IF(IPGAS.NE.0.AND.IMGAS.NE.0) THEN
c            IF(IPDIF.EQ.-1.AND.IMDIF.EQ.1) INTERACT = 1
c          END IF
c        END IF
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from CONBLOCKS '
        WRITE(6,*) ' IATP IBTP JATP JBTP ',IATP,IBTP,JATP,JBTP
        WRITE(6,*) ' IH_OCC_CONS, INTERACT = ', IH_OCC_CONS,INTERACT
      END IF
*
      RETURN
      END
