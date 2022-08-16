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
* Copyright (C) 2021, Bruno Tenorio                                    *
************************************************************************
      SUBROUTINE SDCHS(IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,
     &                   PSI1,PSI2,IF20,IF02,NSDCHSM,SDCHSM)

      IMPLICIT NONE
      REAL*8 PSI1(*),PSI2(*),SDCHSM(*)
    
      REAL*8 COEFF,OVERLAP_RASSI,OVLP
      !REAL*8 PRTHR
      !INTEGER MAPORB(*)
      INTEGER IORBTAB(*),NASORB
      INTEGER NSDCHSM,NASHT
      INTEGER ISSTAB(*)
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      INTEGER FSBOP,IMODE
      !INTEGER NDETS1,NDETS2
      INTEGER LFSBANN1,LFSBANN2
      INTEGER ISORB,JSORB,IJ
      INTEGER LANN1,LANN2
      INTEGER ND1,ND2
      INTEGER KOINFO
      INTEGER ISMLAB,ISPLAB,JSMLAB,JSPLAB,LSMLAB,LSPLAB
      LOGICAL IF20,IF02

#include "SysDef.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
      EXTERNAL OVERLAP_RASSI

C Calculates the DCH  matrix elements between two states with
C N and N-1 electrons, defined as:
C IF02 D = < N-2 | anni_right anni_right | N >, or
C IF20 D = < N | anni_left anni_left  | N-2 >
 
C reduced 2-electron tdm in the space of active spin-orbitals
      !Nr of active spin-orbitals
      NASORB=IORBTAB(4)

C Nr of active orbital:
      NASHT=NASORB/2

C IF02 = eliminte one electron to the right: < N-2 | anni_right
C anni_right | N >
      IF(IF02) THEN

       KOINFO=19

       DO ISORB=1,NASORB
C Symmetry properties:
        ISMLAB=IORBTAB(KOINFO+1+8*(ISORB-1))
        ISPLAB=IORBTAB(KOINFO+3+8*(ISORB-1))
C Annihilate a single spin orbital from PSI2, the spin orbital ISORB:
        IMODE=-1
        LFSBANN1=FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB2)
        ND1=IWORK(LFSBANN1+4)
        COEFF=1.0D0
        CALL GETMEM('ANN1','Allo','Real',LANN1,ND1)
        CALL DCOPY_(ND1,[0.0D0],0,WORK(LANN1),1)
        CALL PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,IWORK(LFSBANN1),
     &                IFSBTAB2,COEFF,WORK(LANN1),PSI2)
CTEST       WRITE(*,*)' The ANN1 wave function, with ISORB=',ISORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN1),PRTHR,WORK(LANN1))

        DO JSORB=1,ISORB-1
C Symmetry properties:
         JSMLAB=IORBTAB(KOINFO+1+8*(JSORB-1))
         JSPLAB=IORBTAB(KOINFO+3+8*(JSORB-1))
C Pair index J,L:
         !IJSORB=((ISORB-1)*(ISORB-2))/2+JSORB

         OVLP=0.0D0
C Annihilate another spin orbital from PSI2, LSORB:
         IMODE=-1
         LFSBANN2=FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,IWORK(LFSBANN1))
         ND2=IWORK(LFSBANN2+4)
         COEFF=1.0D0
         CALL GETMEM('ANN2','Allo','Real',LANN2,ND2)
         CALL DCOPY_(ND2,[0.0D0],0,WORK(LANN2),1)
         CALL PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,IWORK(LFSBANN2),
     &                   IWORK(LFSBANN1),COEFF,WORK(LANN2),WORK(LANN1))
          
!          write(6,*) 'BRN print the WF table WORK(LANN1)'
!          PRTHR=0.01D0
!          CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN1),PRTHR,WORK(LANN1))
!          write(6,*) 'BRN print the WF table WORK(LANN3)'
!          CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN3),PRTHR,WORK(LANN3))

C Compute the spin transition density matrix element:
         OVLP=OVERLAP_RASSI(IFSBTAB1,
     &                  IWORK(LFSBANN2),PSI1,WORK(LANN2))

         IJ=((ISORB-1)*(ISORB-2))/2+JSORB
         SDCHSM(IJ)=SDCHSM(IJ)+OVLP

         !Write(6,*)'BRNCAT I,J, SDCHSM(IJ)',ISORB,JSORB,SDCHSM(IJ)
! 100      CONTINUE
         !END DO
   
         CALL GETMEM('ANN2','Free','Real',LANN2,ND2)
         CALL KILLOBJ(LFSBANN2)
        !CALL GETMEM('ANN1','Free','Real',LANN1,ND1)
        !CALL KILLOBJ(LFSBANN1)
        END DO
        CALL GETMEM('ANN1','Free','Real',LANN1,ND1)
        CALL KILLOBJ(LFSBANN1)
       END DO

!################################################################################
C IF02 = Eliminate to the right (state 2)
      ELSE IF(IF20) THEN
       WRITE(6,*) 'Invalid state combination. 
     &      Please, give PSI1=(N-2) and PSI2=(N) '

      ELSE
       WRITE(6,*)'Invalid state combination in DCH states'
       !WRITE(6,*)'(No such Dyson orbital can exist!)'

      END IF 

!################################################################################

      RETURN
      END
