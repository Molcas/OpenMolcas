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
* Copyright (C) 1995, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE ZNELFSPGP(NTESTG)
      use lucia_data, only: NGAS
      use lucia_data, only: NSTTP,IBSPGPFTP,ISPGPFTP,NELFGP,
     &                      NELFSPGP,NSPGPFTP
      use lucia_data, only: MXPNGAS
*
* Generate for each supergroup the number of electrons in each active
* orbital space and store in NELFSPGP
*
* Jeppe Olsen, July 1995
*
      IMPLICIT NONE
      INTEGER NTESTG
*. input
*. Input and Output ( NELFSPGP(MXPNGAS,MXPSTT) )
      INTEGER NTESTL,NTEST,ITP,NSPGP,IBSPGP,ISPGP,IGAS
*
      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
*
      DO ITP = 1, NSTTP
        NSPGP = NSPGPFTP(ITP)
        IBSPGP = IBSPGPFTP(ITP)
        DO ISPGP = IBSPGP,IBSPGP + NSPGP - 1
          DO IGAS = 1, NGAS
            NELFSPGP(IGAS,ISPGP) = NELFGP(ISPGPFTP(IGAS,ISPGP))
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Distribution of electrons in Active spaces '
        DO ITP = 1, NSTTP
          WRITE(6,*) ' String type ', ITP
          WRITE(6,*) ' Row : active space, Column: supergroup '
          NSPGP = NSPGPFTP(ITP)
          IBSPGP = IBSPGPFTP(ITP)
          CALL IWRTMA(NELFSPGP(1,IBSPGP),NGAS,NSPGP,MXPNGAS,NSPGP)
        END DO
      END IF
*
      END SUBROUTINE ZNELFSPGP
