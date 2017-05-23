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
      SUBROUTINE PNT2DM(   I12SM,   NSMOB,   NSMSX,    OSXO,    IPSM,
     &                      JPSM,    IJSM,    ISM2,   IPNTR,  MXPOBS)
*
* Pointer to two dimensional array
*
* =====
* Input
* =====
* I12SM  : ne.0 => restrict to lower half
*          eq.0 => complete matrix
* NSMOB : Number of orbital symmetries
* NSMSX : Number of SX      symmetries
* OSXO  : Symmetry of orbital, SX => symmetry of other orbital
* IPSM : Number of orbitals per symmetry for index 1
* JPSM : Number of orbitals per symmetry for index 2
* IJSM  : Symmetry of two index array
*
* =======
* Output
* =======
* IPNTR : Pointer to block with first index of given symmetry
*         = 0 indicates forbidden block
* ISM2  : symmetry of second index for given first index
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER OSXO(MXPOBS,2*MXPOBS),IPSM(*),JPSM(*)
*.Output
      DIMENSION IPNTR(*),ISM2(*)
*
      CALL ISETVC(IPNTR,0,NSMOB)
      CALL ISETVC(ISM2 ,0,NSMOB)
      IOFF = 1
      DO 100 ISM = 1,NSMOB
        JSM = OSXO(ISM,IJSM)
        IF(JSM.EQ.0) GOTO 100
        IF(I12SM.EQ.0.OR.ISM.GE.JSM) THEN
*. Allowed block
          IPNTR(ISM) = IOFF
          ISM2(ISM) = JSM
          IF(I12SM.GT.0.AND.ISM.EQ.JSM) THEN
            IOFF = IOFF + IPSM(ISM)*(IPSM(ISM)+1)/2
          ELSE
            IOFF = IOFF + IPSM(ISM)*JPSM(JSM)
          END IF
        END IF
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' dimension of two-dimensional array ',IOFF-1
      END IF
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Pointer '
        CALL IWRTMA(IPNTR,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' Symmetry of other array '
        CALL IWRTMA(ISM2,1,NSMOB,1,NSMOB)
      END IF
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NSMSX)
      END
