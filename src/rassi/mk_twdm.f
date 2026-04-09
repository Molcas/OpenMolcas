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
      SUBROUTINE MK_TWDM(mSym,TDMZZ,WDMZZ,nTDMZZ,SCR,nSCR,iOFF,NBASF,
     &                   ISY12)
      use definitions, only: iwp, wp
      use constants, only: Zero
      use Symmetry_Info, only: MUL
      IMPLICIT NONE
      INTEGER(KIND=IWP), INTENT(IN):: mSYM, nTDMZZ,nSCR,ISY12
      REAL(KIND=WP), INTENT(IN):: TDMZZ(NTDMZZ),WDMZZ(NTDMZZ)
      REAL(KIND=WP), INTENT(OUT):: SCR(nSCR,4)
      INTEGER(KIND=IWP), INTENT(IN):: IOFF(mSYM), NBASF(mSym)

      INTEGER(KIND=IWP) IOF,ITD,ISY,NB,J,I,IJ,ISY1,NB1,ISY2,NB2
      REAL(KIND=WP) TDM,WDM
C CALCULATE THE SYMMETRIC AND ANTISYMMETRIC FOLDED TRANS D MATRICES
C AND SIMILAR WE-REDUCED SPIN DENSITY MATRICES

C SPECIAL CASE: DIAGONAL SYMMETRY BLOCKS.
      SCR(:,:)=Zero
      IF(ISY12.EQ.1) THEN
        IOF=0
        ITD=0
        DO ISY=1,mSym
          NB=NBASF(ISY)
          IF(NB.EQ.0) CYCLE
          DO J=1,NB
            DO I=1,NB
              ITD=ITD+1
              TDM=TDMZZ(ITD)
              WDM=WDMZZ(ITD)
              IF(I.GE.J) THEN
                IJ=IOF+(I*(I-1))/2+J
                IF(I.GT.J) THEN
                  SCR(IJ,2)=SCR(IJ,2)+TDM
                  SCR(IJ,4)=SCR(IJ,4)+WDM
                END IF
              ELSE
                IJ=IOF+(J*(J-1))/2+I
                SCR(IJ,2)=SCR(IJ,2)-TDM
                SCR(IJ,4)=SCR(IJ,4)-WDM
              END IF
              SCR(IJ,1)=SCR(IJ,1)+TDM
              SCR(IJ,3)=SCR(IJ,3)+WDM
            END DO
          END DO
          IOF=IOF+(NB*(NB+1))/2
        END DO
      ELSE
C GENERAL CASE, NON-DIAGONAL SYMMETRY BLOCKS
C THEN LOOP OVER ELEMENTS OF TDMZZ
        ITD=0
        DO ISY1=1,mSym
          NB1=NBASF(ISY1)
          IF(NB1.EQ.0) CYCLE
          ISY2=MUL(ISY1,ISY12)
          NB2=NBASF(ISY2)
          IF(NB2.EQ.0) CYCLE
          IF(ISY1.GT.ISY2) THEN
            DO J=1,NB2
              DO I=1,NB1
                ITD=ITD+1
                TDM=TDMZZ(ITD)
                WDM=WDMZZ(ITD)
                IJ=IOFF(ISY1)+I+NB1*(J-1)
                SCR(IJ,1)=SCR(IJ,1)+TDM
                SCR(IJ,2)=SCR(IJ,2)+TDM
                SCR(IJ,3)=SCR(IJ,3)+WDM
                SCR(IJ,4)=SCR(IJ,4)+WDM
              END DO
            END DO
          ELSE
            DO J=1,NB2
              DO I=1,NB1
                ITD=ITD+1
                TDM=TDMZZ(ITD)
                WDM=WDMZZ(ITD)
                IJ=IOFF(ISY2)+J+NB2*(I-1)
                SCR(IJ,1)=SCR(IJ,1)+TDM
                SCR(IJ,2)=SCR(IJ,2)-TDM
                SCR(IJ,3)=SCR(IJ,3)+WDM
                SCR(IJ,4)=SCR(IJ,4)-WDM
              END DO
            END DO
          END IF
        END DO
      END IF
*
      END SUBROUTINE MK_TWDM
