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
      IMPLICIT REAL*8 (A-H,O-Z)
#include "symmul.fh"
      REAL*8 TDMZZ(NTDMZZ),WDMZZ(NTDMZZ),SCR(nSCR,4)
      DIMENSION IOFF(mSYM), NBASF(mSym)

C CALCULATE THE SYMMETRIC AND ANTISYMMETRIC FOLDED TRANS D MATRICES
C AND SIMILAR WE-REDUCED SPIN DENSITY MATRICES

C SPECIAL CASE: DIAGONAL SYMMETRY BLOCKS.
      SCR(:,:)=0.0D0
      IF(ISY12.EQ.1) THEN
        IOF=0
        ITD=0
        DO 100 ISY=1,mSym
          NB=NBASF(ISY)
          IF(NB.EQ.0) GOTO 100
          DO 90 J=1,NB
            DO 90 I=1,NB
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
90        CONTINUE
          IOF=IOF+(NB*(NB+1))/2
100     CONTINUE
      ELSE
C GENERAL CASE, NON-DIAGONAL SYMMETRY BLOCKS
C THEN LOOP OVER ELEMENTS OF TDMZZ
        ITD=0
        DO 200 ISY1=1,mSym
          NB1=NBASF(ISY1)
          IF(NB1.EQ.0) GOTO 200
          ISY2=MUL(ISY1,ISY12)
          NB2=NBASF(ISY2)
          IF(NB2.EQ.0) GOTO 200
          IF(ISY1.GT.ISY2) THEN
            DO 180 J=1,NB2
              DO 180 I=1,NB1
                ITD=ITD+1
                TDM=TDMZZ(ITD)
                WDM=WDMZZ(ITD)
                IJ=IOFF(ISY1)+I+NB1*(J-1)
                SCR(IJ,1)=SCR(IJ,1)+TDM
                SCR(IJ,2)=SCR(IJ,2)+TDM
                SCR(IJ,3)=SCR(IJ,3)+WDM
                SCR(IJ,4)=SCR(IJ,4)+WDM
180         CONTINUE
          ELSE
            DO 190 J=1,NB2
              DO 190 I=1,NB1
                ITD=ITD+1
                TDM=TDMZZ(ITD)
                WDM=WDMZZ(ITD)
                IJ=IOFF(ISY2)+J+NB2*(I-1)
                SCR(IJ,1)=SCR(IJ,1)+TDM
                SCR(IJ,2)=SCR(IJ,2)-TDM
                SCR(IJ,3)=SCR(IJ,3)+WDM
                SCR(IJ,4)=SCR(IJ,4)-WDM
190         CONTINUE
          END IF
200     CONTINUE
      END IF
*
      RETURN
      END
