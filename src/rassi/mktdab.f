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
* Copyright (C) 1989, Per Ake Malmqvist                                *
************************************************************************
*****************************************************************
*  PROGRAM RASSI        PER-AAKE MALMQVIST
*  SUBROUTINE MKTDAB    IBM-3090 RELEASE 89 01 31
*  PURPOSE: CALCULATE TRANSITION DENSITY MATRIX FOR CI EXPANSIONS IN
*  BIORTHONORMAL ORBITAL BASES A AND B.
*****************************************************************
      SUBROUTINE MKTDAB(OVER,GAMMA1,TDMAB,iRC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TDMAB(NTDMAB)
      DIMENSION GAMMA1(NASHT,NASHT)
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "WrkSpc.fh"
      DIMENSION IOFFA(8)
C IOFFA=NR OF ACTIVE ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
      IOFFA(1)=0
      DO I=1,NSYM-1
        IOFFA(I+1)=IOFFA(I)+NASH(I)
      end do
C  INITIALIZE TRANSITION DENSITY MATRIX:
      CALL FZERO(TDMAB,NTDMAB)
C CONTRIBUTION FROM INACTIVE ORBITALS:
      IF (LSYM1.EQ.LSYM2) THEN
         IF (OVER.NE.0.0D0) THEN
            IOFFTD=0
            DO ISY=1,NSYM
               II=0
               DO I=1,NISH(ISY)
                  II=II+1
                  IPOS=IOFFTD+(II-1)*NOSH(ISY)+II
                  TDMAB(IPOS)=2.0D0*OVER
               end do
               IOFFTD=IOFFTD+NOSH(ISY)**2
            end do
         END IF
      END IF
C THEN ADD CONTRIBUTION FROM ACTIVE SPACE.
      ISY12=MUL(LSYM1,LSYM2)
      IOFFTD=0
      DO ISY1=1,NSYM
         NO1=NOSH(ISY1)
         IF(NO1.EQ.0) cycle
         ISY2=MUL(ISY1,ISY12)
         NO2=NOSH(ISY2)
         IF(NO2.EQ.0) cycle
         NA1=NASH(ISY1)
         if (NA1 /= 0) then
           NA2=NASH(ISY2)
           if (NA2 /= 0) then
             NI1=NISH(ISY1)
             NI2=NISH(ISY2)
             do I=1,NA1
               IA=IOFFA(ISY1)+I
               II=NI1+I
               do J=1,NA2
                 JA=IOFFA(ISY2)+J
                 JJ=NI2+J
                 IPOS=IOFFTD+II+(JJ-1)*NO1
                 TDMAB(IPOS)=GAMMA1(IA,JA)
               end do
             end do
           end if
         end if
      IOFFTD=IOFFTD+NO1*NO2
      end do
*
      iRC=1
      If (DDot_(nTDMAB,TDMAB,1,TDMAB,1).le.0.0D0) iRC=0
      RETURN
      END
