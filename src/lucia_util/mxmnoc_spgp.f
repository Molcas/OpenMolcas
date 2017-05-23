************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE MXMNOC_SPGP(  MINEL,  MAXEL, NORBTP,NORBFTP, NELFTP,
     &                           NTESTG)
*
* Construct accumulated MAX and MIN arrays for a GAS supergroup
*
      IMPLICIT REAL*8           ( A-H,O-Z)
*. Output
      DIMENSION  MINEL(*),MAXEL(*)
*. Input
      INTEGER NORBFTP(*),NELFTP(*)
*
* Some dummy initializations
      IORB_START = 1 ! jwk-cleanup
*
      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==========='
        WRITE(6,*) ' MXMNOC_SPGP'
        WRITE(6,*) ' ==========='
        WRITE(6,*)
C?      WRITE(6,*) ' NORBFTP : '
C?      CALL IWRTMA(NORBFTP,1,NORBTP,1,NORBTP)
      END IF
*
      DO IORBTP = 1, NORBTP
*. Max and min at start of this type and at end of this type
        IF(IORBTP.EQ.1) THEN
          IORB_START = 1
          IORB_END = NORBFTP(1)
          NEL_START = 0
          NEL_END   = NELFTP(1)
        ELSE
          IORB_START =  IORB_START + NORBFTP(IORBTP-1)
          IORB_END   =  IORB_START + NORBFTP(IORBTP)-1
          NEL_START = NEL_END
          NEL_END   = NEL_START + NELFTP(IORBTP)
        END IF
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IORBTP,IORB_START-IORB_END,NEL_START,NEL_END '
          WRITE(6,*)   IORBTP,IORB_START-IORB_END,NEL_START,NEL_END
        END IF
*
        DO IORB = IORB_START, IORB_END
          MAXEL(IORB) = MIN(IORB,NEL_END)
          MINEL(IORB) = NEL_START
          IF(NEL_END-MINEL(IORB).GT. IORB_END-IORB)
     &    MINEL(IORB) = NEL_END - ( IORB_END - IORB )
        END DO
      END DO
*
      IF( NTEST .GE. 100 ) THEN
        NORB = IELSUM(NORBFTP,NORBTP)
        WRITE(6,*) ' MINEL : '
        CALL IWRTMA(MINEL,1,NORB,1,NORB)
        WRITE(6,*) ' MAXEL : '
        CALL IWRTMA(MAXEL,1,NORB,1,NORB)
      END IF
*
      RETURN
      END
