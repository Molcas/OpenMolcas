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
      FUNCTION IBASSPC_FOR_CLS(ICLS)
*
*. Obtain base space for occupation class ICLS
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.fh"
#include "cgas.fh"
*. Specific input
      INTEGER ICLS(NGAS)
*
* Some dummy initializations
*
      NEL = 0 ! jwk-cleanup
*
      IBASE = 0
      DO ISPC = 1, NCMBSPC
        DO JJSPC = 1, LCMBSPC(ISPC)
          JSPC = ICMBSPC(JJSPC,ISPC)
*. Test for occupation constraints in CI space JSPC
          I_AM_OKAY = 1
          DO IGAS = 1, NGAS
            IF(IGAS.EQ.1) THEN
              NEL = ICLS(IGAS)
            ELSE
              NEL = NEL + ICLS(IGAS)
            END IF
*
            IF(NEL.LT.IGSOCCX(IGAS,1,JSPC).OR.
     &         NEL.GT.IGSOCCX(IGAS,2,JSPC)    ) THEN
                I_AM_OKAY = 0
            END IF
          END DO
*         ^ End of loop over gasspaces for given cispace
*
          IF(I_AM_OKAY.EQ.1.AND.IBASE.EQ.0) THEN
            IBASE = ISPC
          END IF
*
        END DO
*       ^ End of loop over cisspaces for given combination space
      END DO
*     ^ End of loop over combinations apaces
*
      IBASSPC_FOR_CLS = IBASE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation class and its basespace '
        CALL IWRTMA(ICLS,1,NGAS,1,NGAS)
        WRITE(6,*) IBASE
      END IF
*
      RETURN
      END
CTOBE       SUBROUTINE CLS_TO_SPSP
CTOBE      &(ICLS,NCLS,I_CLS_TO_SPSP,N_CLS_TO_SPSP,IB_CLS_TO_SPSP,
CTOBE      & ISPSPCLS,NOCTPA,NOCTPB)
CTOBE *
CTOBE * Combination of supergroups belonging to given class
CTOBE *
CTOBE * Jeppe Olsen, Nov 99
CTOBE *
CTOBE * Input :
CTOBE *
CTOBE * ISPSPCLS : Class for given pair of supergroups
CTOBE *
CTOBE *. Output
CTOBE *
CTOBE * N_CLS_TO_SPSP  : Number of supergroup combinations for class
CTOBE * IB_CLS_TO_SPSP : Base  for supergroup combinations of given class
CTOBE * I_CLS_TO_SPSP  : supergroup combinations of given class
CTOBE *
CTOBE       INCLUDE 'implicit.fh'
CTOBE       INCLUDE 'mxpdim.fh'
CTOBE       INCLUDE 'gasstr.fh'
CTOBE       INCLUDE 'cgas.fh'
CTOBE       INCLUDE 'WrkSpc.fh'
CTOBE *. Input
CTOBE       INTEGER ISPSPCLS(NOCTPA,NOCTPB)
CTOBE *. Output
CTOBE       INTEGER I_CLS_TO_SPSP(2,*),N_CLS_TO_SPSP(NCLS)
CTOBE       INTEGER IB_CLS_TO_SPSP(NCLS)
CTOBE *
CTOBE       IZERO = 0
CTOBE       CALL ISETVC(N_CLS_TO_SPSP,IZERO,NCLS)
*
CTOBE       DO IOCTPA = 1, NOCTPA
CTOBE        DO IOCTPB = 1, NOCTPB
CTOBE          ICLS = ISPSPCLS(IOCTPA,IOCTPB)
CTOBE          N_CLS_TO_SPSP(ICLS) = N_CLS_TO_SPSP(ICLS) + 1
CTOBE        END DO
CTOBE       END DO
CTOBE *
CTOBE       IB_CLS_TO_SPSP(1) = 1
CTOBE       DO ICLS = 2, NCLS
CTOBE         IB_CLS_TO_SPSP(ICLS) =
CTOBE      &  IB_CLS_TO_SPSP(ICLS-1) + N_CLS_TO_SPSP(ICLS-1)
CTOBE       END DO
CTOBE *
CTOBE       CALL ISETVC(N_CLS_TO_SPSP,IZERO,NCLS)
CTOBE       DO IOCTPA = 1, NOCTPA
CTOBE        DO IOCTPB = 1, NOCTPB
CTOBE          ICLS = ISPSPCLS(IOCTPA,IOCTPB)
CTOBE          I_CLS_TO_SPSP(1,
CTOBE          N_CLS_TO_SPSP(ICLS) = N_CLS_TO_SPSP(ICLS) + 1
CTOBE        END DO
CTOBE       END DO
CTOBE
