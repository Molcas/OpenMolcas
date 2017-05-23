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
      FUNCTION NUMST3(NEL,NORB1,NEL1MN,NEL1MX,NORB2,
     &                NORB3,NEL3MN,NEL3MX)
*
* Number of strings with NEL electrons that fullfills
*
* Between NEL1MN AND NEL1MX electrons in the first NORB1 orbitals
* Between NEL3MN AND NEL3MX electrons in the last  NORB3 orbitals
*
*
*
      NTEST = 0
      NSTRIN = 0
*
      DO 100 IEL1 = NEL1MN,MIN(NEL1MX,NORB1,NEL)
        NSTIN1 = IBION_LUCIA(NORB1,IEL1)
        IEL3MN = MAX ( NEL3MN,NEL-(IEL1+NORB2) )
        IEL3MX = MIN ( NEL3MX,NEL-IEL1)
        DO 80 IEL3 = IEL3MN, IEL3MX
         IEL2 = NEL - IEL1-IEL3
         NSTINT = NSTIN1*IBION_LUCIA(NORB2,IEL2)*IBION_LUCIA(NORB3,IEL3)
         NSTRIN = NSTRIN + NSTINT
  80   CONTINUE
 100  CONTINUE
      NUMST3 = NSTRIN
*
      IF( NTEST .GE.1 )
     &WRITE(6,'(/A,I6)') '  Number of strings generated ... ', NSTRIN
*
      RETURN
      END
