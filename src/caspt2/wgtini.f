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
* Copyright (C) 2019, Stefano Battaglia                                *
************************************************************************
      SUBROUTINE WGTINI

      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "intgrl.fh"
#include "eqsolv.fh"
#include "warnings.fh"

      CALL QENTER('WGTINI')

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)' Entered WGTINI.'
      END IF

* Initialize array of weights with all zeros
      CALL DCOPY_(NSTATE**2,[0.0D0],0,WORK(LFWGT),1)
      CALL DCOPY_(NSTATE**2,[0.0D0],0,WORK(LVWGT),1)

      DO ISTATE=1,NSTATE

* Compute weights for constructing H0
        IF (IFFDW) THEN
          EBETA = REFENE(ISTATE)
* Compute normalization factor
          DO JSTATE=1,NSTATE
            EALPHA = REFENE(JSTATE)
            FAC = 0.0D0
            DO KSTATE=1,NSTATE
              EGAMMA = REFENE(KSTATE)
              FAC = FAC + EXP(-ZETAF*(EALPHA - EGAMMA)**2)
            END DO
            IJ = (ISTATE-1) + NSTATE*(JSTATE-1)
            WORK(LFWGT+IJ) = EXP(-ZETAF*(EALPHA - EBETA)**2)/FAC
          END DO
* If it is an XMS-CASPT2 calculation, all the weights are 1/NSTATE
        ELSE IF (IFXMS) THEN
          CALL DCOPY_(NSTATE**2,1.0D0/NSTATE,0,WORK(LFWGT),1)
        ELSE
* It is a normal MS-CASPT2 and the weight vectors are the standard e_1, e_2, ...
          WORK(LFWGT + (NSTATE*(ISTATE-1)) + (ISTATE-1)) = 1.0D0
        END IF

* Compute weights for constructing V
        IF (IFVDW) THEN
          EBETA = REFENE(ISTATE)
* Compute normalization factor
          DO JSTATE=1,NSTATE
            EALPHA = REFENE(JSTATE)
            FAC = 0.0D0
            DO KSTATE=1,NSTATE
              EGAMMA = REFENE(KSTATE)
              FAC = FAC + EXP(-ZETAV*(EALPHA - EGAMMA)**2)
            END DO
            IJ = (ISTATE-1) + NSTATE*(JSTATE-1)
            WORK(LVWGT+IJ) = EXP(-ZETAV*(EALPHA - EBETA)**2)/FAC
          END DO
* It is a normal (X)MS-CASPT2 and the weight vectors are the standard e_1, e_2, ...
        ELSE
          WORK(LVWGT + (ISTATE-1) + NSTATE*(ISTATE-1)) = 1.0D0
        END IF

      END DO

      ! IF (IFFDW) THEN
        WRITE(6,*)
        WRITE(6,*)' Weights used for H0:'
        DO I=1,NSTATE
          WRITE(6,'(1x,10f8.4)')(WORK(LFWGT + (I-1) + NSTATE*(J-1)),
     &    J=1,NSTATE)
        END DO
        WRITE(6,*)
      ! END IF

      ! IF (IFVDW) THEN
        WRITE(6,*)
        WRITE(6,*)' Weights used for V:'
        DO I=1,NSTATE
          WRITE(6,'(1x,10f8.4)')(WORK(LVWGT + (I-1) + NSTATE*(J-1)),
     &    J=1,NSTATE)
        END DO
        WRITE(6,*)
      ! END IF

      CALL QEXIT('WGTINI')
      RETURN
      END
