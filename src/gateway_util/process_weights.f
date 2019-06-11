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
* Copyright (C) 2013, Ignacio Fdez. Galvan                             *
*               2016, Roland Lindh                                     *
************************************************************************
*  Process_Weights
*
*> @brief
*>   Process and store the weights used for alignment
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Process the weights according to the `WEIG` option in `GATEWAY`.
*> These weights are used in alignment and in the "sphere" constraint.
*> The list of weights is stored in the runfile, first the symmetry-unique
*> atoms and then the symmetric images, in the manner of ::Expand_Coor.
*>
*> @param[in] iPrint Print level
************************************************************************
      SUBROUTINE Process_Weights(iPrint)
      IMPLICIT REAL*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "constants2.fh"
#include "real.fh"
#include "stdalloc.fh"
      LOGICAL Small
      PARAMETER ( thr=1.0d-6 )
      REAL*8, DIMENSION(:), ALLOCATABLE :: W
*
*---- Count the total and symmetry-unique number of atoms
      nAt=0
      nSymAt=0
      ndc=0
      DO i=1,nCnttp
        DO j=1,nCntr(i)
          ndc=ndc+1
          IF (.NOT.(pChrg(i).OR.FragCnttp(i).OR.AuxCnttp(i))) THEN
            nAt=nAt+nIrrep/nStab(ndc)
            nSymAt=nSymAt+1
          END IF
        END DO
      END DO
      CALL mma_allocate(W,nAt,label='W')
*
*---- By default, all weights are 1
      call dcopy_(nAt,[One],0,W,1)
*
      IF (Align_Weights(1:4).EQ.'MASS') Then
*---- Set the weights to the mass of each atom
        j=1
        DO i=1,nCnttp
          IF (.NOT.(pChrg(i).OR.FragCnttp(i).OR.AuxCnttp(i))) THEN
            DO iCnt=1,nCntr(i)
              W(j)=CntMass(i)/UTOAU
              j=j+1
            END DO
          END IF
        END DO
      ELSE IF (Align_Weights(1:5).EQ.'HEAVY') THEN
*---- Set the the weight to 1 for heavy atoms, 0 for hydrogens
        j=1
        DO i=1,nCnttp
          IF (.NOT.(pChrg(i).OR.FragCnttp(i).OR.AuxCnttp(i))) THEN
            DO iCnt=1,nCntr(i)
              IF (iAtmNr(i).LE.1) W(j)=Zero
              j=j+1
            END DO
          END IF
        END DO
      ELSE IF (Align_Weights(1:5).EQ.'EQUAL') THEN
*---- EQUAL is already the default: 1 for all
        CONTINUE
      ELSE
*---- Read the weights from the input line
        READ(Align_Weights,*,IOSTAT=iErr) (W(i),i=1,nAt)
        IF (iErr.GT.0) THEN
          CALL WarningMessage(2,'Unable to read data from WEIG')
          CALL Quit_OnUserError()
        END IF
      END IF
*
*---- Unfold the symmetry
      iSymAt=1
      iAt=1+nSymAt
      ndc=0
      DO i=1,nCnttp
        DO j=1,nCntr(i)
          ndc=ndc+1
          IF (.NOT.(pChrg(i).OR.FragCnttp(i).OR.AuxCnttp(i))) THEN
            DO k=1,nIrrep/nStab(ndc)-1
              W(iAt)=W(iSymAt)
              iAt=iAt+1
            END DO
            iSymAt=iSymAt+1
          END IF
        END DO
      END DO
*
*---- Check for zero total weight
      wTot=Zero
      DO i=1,nAt
        wTot=wTot+W(i)
      END DO
      IF (wTot.LT.thr) THEN
        CALL WarningMessage(1,
     &       'Total weight too small. Setting equal weights.')
        DO i=1,nAt
          W(i)=One
        END DO
      END IF
*---- Prevent zero weights, it could break the "sphere" constraint
*     (a value between 1e-6 and 1e-1 can still be used)
      Small=.FALSE.
      DO i=1,nAt
        IF (W(i).LT.thr) THEN
          W(i)=1.0d-1
          Small=.TRUE.
        END IF
      END DO
      IF (iPrint.GE.6) THEN
        IF (Small) THEN
          CALL WarningMessage(1,
     &         'Small weights were increased to avoid problems with'//
     &         ' constraints.')
        END IF
        CALL RecPrt('Weights used for alignment and distance',' ',
     &              W,nAt,1)
        WRITE(6,*)
      END IF
*
*---- Store weights in the runfile too
      CALL Put_dArray('Weights',W,nAt)
      CALL mma_deallocate(W)
*
      END
