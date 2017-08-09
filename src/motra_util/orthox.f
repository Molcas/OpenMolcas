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
      SUBROUTINE ORTHOX_MOTRA(S,C,NORB,NBAS)
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(NBAS,NORB),S(NORB,NORB)
*
      Call qEnter('ORTHOX')
*
      DO IORB=1,NORB
         F=1.0/SQRT(S(IORB,IORB))
         DO IBAS=1,NBAS
            C(IBAS,IORB)=F*C(IBAS,IORB)
         END DO
         DO JORB=1,NORB
            S(IORB,JORB)=F*S(IORB,JORB)
            S(JORB,IORB)=F*S(JORB,IORB)
         END DO
         DO JORB=IORB+1,NORB
            A=S(IORB,JORB)
            DO IBAS=1,NBAS
               C(IBAS,JORB)=C(IBAS,JORB)-A*C(IBAS,IORB)
            END DO
            DO KORB=1,NORB
               S(JORB,KORB)=S(JORB,KORB)-A*S(IORB,KORB)
            END DO
            DO KORB=1,NORB
               S(KORB,JORB)=S(KORB,JORB)-A*S(KORB,IORB)
            END DO
         END DO
      END DO
*
      Call qExit('ORTHOX')
*
      RETURN
      END
