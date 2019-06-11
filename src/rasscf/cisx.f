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
      SUBROUTINE CISX(IDX,D,DS,PS,PA,SCR)

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION IDX(NAC),D(*),DS(*),PS(*),PA(*),SCR(*)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"

*
* Convert from CI to SX ordering
* Note: A factor of 2 arises because matrices are folded
*
      ITR(I)=I*(I-1)/2

*     one-body density
      IJO=0
      NIJ=ITR(NAC+1)
      NIJKL=ITR(NIJ+1)
      DO I=1,NAC
        DO J=1,I
          INEW=IDX(I)
          JNEW=IDX(J)
          IF (JNEW.GT.INEW)THEN
            JNEW=IDX(I)
            INEW=IDX(J)
          ENDIF
          IJNEW=ITR(INEW)+JNEW
          IJO=IJO+1
          SCR(IJNEW)=D(IJO)
        END DO
      END DO
      CALL DCOPY_(NIJ,SCR,1,D,1)

*     spin density
      IJO=0
      NIJ=ITR(NAC+1)
      NIJKL=ITR(NIJ+1)
      DO I=1,NAC
        DO J=1,I
          INEW=IDX(I)
          JNEW=IDX(J)
          IF (JNEW.GT.INEW)THEN
            JNEW=IDX(I)
            INEW=IDX(J)
          ENDIF
          IJNEW=ITR(INEW)+JNEW
          IJO=IJO+1
          SCR(IJNEW)=DS(IJO)
        END DO
      END DO
      CALL DCOPY_(NIJ,SCR,1,DS,1)

*     symmetrized (iCase=1) and antisymmetrized (iCase=2)
*     two-body density
      DO  ICASE=1,2
      IJKLO=0
      CALL DCOPY_(NACPR2,[0.0D0],0,SCR,1)
      DO I=1,NAC
        DO J=1,I
          INEW=IDX(I)
          JNEW=IDX(J)
          SGN0=1.0D0
          IF (JNEW.GT.INEW)THEN
            JNEW=IDX(I)
            INEW=IDX(J)
            SGN0=-1.0D0
          ENDIF
          IJNEW=ITR(INEW)+JNEW
          DO K=1,I
            LLIM=K
            IF (K.EQ.I)LLIM=J
            DO L=1,LLIM
              IJKLO=IJKLO+1
              KNEW=IDX(K)
              LNEW=IDX(L)
              SGN=SGN0
              IF (LNEW.GT.KNEW)THEN
                KNEW=IDX(L)
                LNEW=IDX(K)
                SGN=-SGN0
              ENDIF
              KLNEW=ITR(KNEW)+LNEW
              IF(KLNEW.GT.IJNEW) THEN
                IJKLN=ITR(KLNEW)+IJNEW
              ELSE
                IJKLN=ITR(IJNEW)+KLNEW
              END IF
              IF(ICASE.EQ.1) THEN
                IF (KLNEW.GT.IJNEW)THEN
                  IF (K.EQ.L.AND.I.NE.J) THEN
                    SCR(IJKLN)=2.0D0*PS(IJKLO)
                  ELSE IF (I.EQ.J.AND.K.NE.L)THEN
                    SCR(IJKLN)=0.5D0*PS(IJKLO)
                  ELSE
                    SCR(IJKLN)=PS(IJKLO)
                  ENDIF
                ELSE
                  SCR(IJKLN)=PS(IJKLO)
                ENDIF
              ELSE
                SCR(IJKLN)=SGN*PA(IJKLO)
              END IF
            END DO
          END DO
        END DO
      END DO
      IF(ICASE.EQ.1) THEN
        CALL DCOPY_(NIJKL,SCR,1,PS,1)
      ELSE
        CALL DCOPY_(NIJKL,SCR,1,PA,1)
      END IF

* End of long loop over ICASE.
      END DO

      RETURN
      END
