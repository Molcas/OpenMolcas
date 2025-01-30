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
      INTEGER FUNCTION IBASSPC_FOR_CLS(ICLS)
*
*. Obtain base space for occupation class ICLS
*
      use lucia_data, only: NGAS,NCMBSPC,ICMBSPC,IGSOCCX,LCMBSPC,
     &                      NCMBSPC
      IMPLICIT NONE
*. General input
*. Specific input
      INTEGER ICLS(NGAS)

      INTEGER NEL,IBASE,ISPC,JJSPC,JSPC,I_AM_OKAY,IGAS,NTEST
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
      END FUNCTION IBASSPC_FOR_CLS
