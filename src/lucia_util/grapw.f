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
* Copyright (C) Jeppe Olsen                                            *
************************************************************************
      SUBROUTINE GRAPW(W,Y,MINEL,MAXEL,NORB,NEL,IPRNT)
*
* A graph of strings has been defined from
*
*      MINEL(I) is the smallest allowed number of electrons in
*      orbitals 1 through I
*
*      MAXEL(I) is the largest allowed number of electrons in
*      orbitals 1 through I
*
* Set up vertex weights W
* Set up arc weights    Y
*
* Reverse lexical ordering is used with
* weights of unoccupied orbitals set to 0
*
* Jeppe Olsen
*
       IMPLICIT REAL*8(A-H,O-Z)
       INTEGER W(NORB+1,NEL+1)
       INTEGER Y(NORB,NEL)
       INTEGER MAXEL(NORB),MINEL(NORB)
*
       NTEST = 0
       NTEST = MAX(NTEST,IPRNT)
*
      Call iCopy((NEL+1)*(NORB+1),[0],0,W,1)
      Call iCopy(NEL*NORB,[0],0,Y,1)
*
*================
*  Vertex weights
*================
*
*. (Weight for vertex(IEL,IORB) is stored in W(IORB+1,IEL+1) )
      W(1,1) = 1
      DO 300 IEL = 0, NEL
        DO 200 IORB = 1, NORB
          IF(MINEL(IORB).LE.IEL .AND. IEL .LE. MAXEL(IORB) ) THEN
            IF( IEL .GT. 0 ) THEN
              W(IORB+1,IEL+1) = W(IORB,IEL+1)
     *                        + W(IORB,IEL)
            ELSE
              W(IORB+1,1) = W(IORB,1)
            END IF
          END IF
200     CONTINUE
300   CONTINUE
*
*. Weight for arc connecting vertices (IORB-1,IEL-1) and(IORB,IEL)
*. is stored in Y(IORB,IEL)
*. Y(IORB,IEL) = W(IORB-1,IEL)
*
      DO 1300 IEL = 1, NEL
        DO 1200 IORB = 1, NORB
          IF(MINEL(IORB).LE.IEL .AND. IEL .LE. MAXEL(IORB) ) THEN
            Y(IORB,IEL) = W(IORB-1+1,IEL+1)
          END IF
1200    CONTINUE
1300  CONTINUE
*
      IF ( NTEST.GE.100 ) THEN
         WRITE(6,*) ' vertex weights'
         CALL IWRTMA(W,NORB+1,NEL+1,NORB+1,NEL+1)
         WRITE(6,*) ' arc weights'
         CALL IWRTMA(Y,NORB,NEL,NORB,NEL)
      END IF
*
      RETURN
      END
