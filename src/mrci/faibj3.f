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
         subroutine faibj3(NSIJ,IFT,
     &  AIBJ,FSEC,FAC,IN,INS,IPOA,IPOF)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION IPOF(9),IPOA(9),
     & AIBJ(*),FSEC(*)
         CALL IPO(IPOA,NVIR,MUL,NSYM,NSIJ,IFT)
C INTEGRAL COMBINATION APPROPRIATE FOR SINGLET-COUPLING:
         DO 170 IASYM=1,NSYM
            IBSYM=MUL(NSIJ,IASYM)
            IF(IBSYM.GT.IASYM)GO TO 170
            IAB=IPOA(IASYM+1)-IPOA(IASYM)
            IF(IAB.EQ.0)GO TO 170
            IF(NSIJ.EQ.1) THEN
              CALL SECEQ(AIBJ(IPOF(IASYM)+1),AIBJ(IPOF(IBSYM)+1),
     *             FSEC(IN+1),NVIR(IASYM),0,FAC)
            ELSE
              CALL SECNE(AIBJ(IPOF(IASYM)+1),AIBJ(IPOF(IBSYM)+1),
     *             FSEC(IN+1),NVIR(IASYM),NVIR(IBSYM),NSIJ,0)
            END IF
            IN=IN+IAB
170      CONTINUE
         INS=IN
C INTEGRAL COMBINATION APPROPRIATE FOR TRIPLET-COUPLING:
         DO 180 IASYM=1,NSYM
            IBSYM=MUL(NSIJ,IASYM)
            IF(IBSYM.GT.IASYM)GO TO 180
            IAB=IPOA(IASYM+1)-IPOA(IASYM)
            IF(IAB.EQ.0)GO TO 180
            IF(NSIJ.EQ.1) THEN

              CALL SECEQ(AIBJ(IPOF(IASYM)+1),AIBJ(IPOF(IBSYM)+1),
     *             FSEC(IN+1),NVIR(IASYM),1,DUMMY)
            ELSE
              CALL SECNE(AIBJ(IPOF(IASYM)+1),AIBJ(IPOF(IBSYM)+1),
     *             FSEC(IN+1),NVIR(IASYM),NVIR(IBSYM),NSIJ,1)
            END IF
            IN=IN+IAB
180      CONTINUE
         return
         end
