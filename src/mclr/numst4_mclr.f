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
      SUBROUTINE NUMST4_MCLR(NEL,NORB1,NEL1MN,NEL1MX,NORB2,
     &                  NORB3,NEL3MN,NEL3MX,NSTTP)
*
* Number of strings per type for class of strings with
*
* Between NEL1MN AND NEL1MX electrons in the first NORB1 orbitals
* Between NEL3MN AND NEL3MX electrons in the last  NORB3 orbitals
*
      INTEGER NSTTP(*)
      NTEST = 0
      NSTRIN = 0
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' NUMST4 '
        WRITE(6,*) ' NEL NEL1MN NEL1MX NEL3MN NEL3MX '
        WRITE(6,*) NEL,NEL1MN,NEL1MX,NEL3MN,NEL3MX
      END IF
*
      ITPMAX = 0
      DO 100 IEL1 = NEL1MN,MIN(NEL1MX,NORB1,NEL)
        NSTIN1 = IBION(NORB1,IEL1)
        IEL3MN = MAX ( NEL3MN,NEL-(IEL1+NORB2) )
        IEL3MX = MIN ( NEL3MX,NEL-IEL1)
        DO 80 IEL3 = IEL3MN, IEL3MX
         IEL2 = NEL - IEL1-IEL3
         NSTINT = NSTIN1*IBION(NORB2,IEL2)*IBION(NORB3,IEL3)
         NSTRIN = NSTRIN + NSTINT
         ITP = (NEL1MX-IEL1)
     &         * (NEL3MX-NEL3MN+1)
     &         + IEL3-NEL3MN + 1
         NSTTP(ITP) = NSTINT
         ITPMAX= MAX(ITPMAX,ITP)
  80   CONTINUE
 100  CONTINUE
*
      IF( NTEST .GE.1 ) THEN
      WRITE(6,'(/A,I6)') '  Number of strings generated ... ', NSTRIN
         WRITE(6,*)
         WRITE(6,*) ' Largest type number ', ITPMAX
         WRITE(6,*) ' Number of strings per type '
         CALL IWRTMA(NSTTP,1,ITPMAX,1,ITPMAX)
      END IF
*
      RETURN
      END
