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
      SUBROUTINE GTJK(       RJ,       RK,    NTOOB,      SCR,   IREOTS,
     &                   IREOST)
*
* Interface routine for obtaining Coulomb (RJ) and
* Exchange integrals (RK)
*
* Ordering of integrals is in the internal order
      IMPLICIT None
*
*.Input
      INTEGER NTOOB
      INTEGER IREOTS(*), IREOST(*)
*.Output
      REAL*8 RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
*.Scratch
      REAL*8 SCR(2*NTOOB ** 2)

      INTEGER NTEST

      CALL GTJK_RASSCF(RJ,RK,NTOOB,IREOST)
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' RJ and RK from GTJK '
        CALL WRTMAT(RJ,NTOOB,NTOOB,NTOOB,NTOOB)
        CALL WRTMAT(RK,NTOOB,NTOOB,NTOOB,NTOOB)
      END IF
*
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_real_array(SCR)
        CALL Unused_integer_array(IREOTS)
      END IF
      END SUBROUTINE GTJK

* Working on EXPHAM
* some known problems :
*     1 : if CSF are used diagonal is not delivered to H0mat
