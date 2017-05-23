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
      SUBROUTINE IAIBCM_MCLR(MNRS1,MXRS3,NOCTPA,NOCTPB,NEL1A,NEL3A,
     *                  NEL1B,NEL3B,IOCOC,IPRNT)
*
* RAS allowed combinations of alpha and beta types
*
      Implicit integer (a-z)
*
* =====
*.Input
* =====
*
* NOCTPA : Number of alpha types
* NEL1A  : Number of electrons in RAS 1 for each alpha type
* NEL3A  : Number of electrons in RAS 3 for each alpha type
*
* NOCTPB : Number of beta types
* NEL1B  : Number of electrons in RAS 1 for each beta type
* NEL3B  : Number of electrons in RAS 3 for each beta type
*
* ======
*.Output
* ======
* IOCOC(IATP,IBTP)  = 1 =>      allowed combination
* IOCOC(IATP,IBTP)  = 0 => not allowed combination
*
*.Input
      INTEGER NEL1A(*),NEL3A(*),NEL1B(*),NEL3B(*)
*.Output
      INTEGER IOCOC(NOCTPA,NOCTPB)
*
*     Call qEnter('IAIBCM')
*
      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)
*
      Call iCopy(nOctPA*nOctPB,0,0,iococ,1)
      DO 100 IATP = 1,NOCTPA
         IAEL1 = NEL1A(IATP)
         IAEL3 = NEL3A(IATP)
         DO 90 IBTP = 1,NOCTPB
            IBEL1 = NEL1B(IBTP)
            IBEL3 = NEL3B(IBTP)
            IF( (IAEL1+IBEL1).GE.MNRS1 .AND.
     &          (IAEL3+IBEL3).LE.MXRS3       )
     &          IOCOC(IATP,IBTP) = 1
90       CONTINUE
100   CONTINUE
*
      IF ( NTEST .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Matrix giving allowed combinations of types '
        WRITE(6,*)
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF

*
*     Call qExit('IAIBCM')
*
      RETURN
      END
