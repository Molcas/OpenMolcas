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
      SUBROUTINE MKDCHS(IFSBTAB1,IFSBTAB2,ISSTAB,
     &                  MAPORB,DET1,DET2,
     &                  IF20,IF02,NDCHSM,DCHSM)

      IMPLICIT NONE
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      INTEGER ISSTAB(*),MAPORB(*),NDCHSM
      REAL*8 DET1(*),DET2(*)
      REAL*8 DCHSM(NDCHSM)
      INTEGER NASHT,NASORB,LORBTB
      REAL*8 GVAL,GAB,GBA
      INTEGER IAJB,IBJA
      INTEGER JORB,IORB
      INTEGER JORBA,JORBB,IORBA,IORBB
      INTEGER ITABS,JTABS,IJTABS
      INTEGER NSDCHSM
      LOGICAL IF20,IF02
#include "symmul.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
      Real*8, Allocatable:: SDCHSM(:)

C Given two CI expansions, using a biorthonormal set of SD''s,
C calculate the matrix elements relevant to DCH state intensities
C in the biorthonormal active orbital basis.
C
C  'I,J,|< N-2 | anni_right anni_right | N >|**2'

      LORBTB=ISSTAB(3)
C Pick out nr of active orbitals from orbital table:
      NASORB=IWORK(LORBTB+3)
      NASHT=NASORB/2
      NSDCHSM= NASORB*(NASORB-1)/2
      Call mma_allocate(SDCHSM,nSDCHSM,Label='SDCHSM')
      SDCHSM(:)=0.0D0

        CALL SDCHS(IWORK(LORBTB),ISSTAB,
     &              IFSBTAB1,IFSBTAB2,DET1,DET2,
     &              IF20,IF02,SDCHSM)

C Mapping from active spin-orbital to active orbital in external order.
C Note that these differ, not just because of the existence of two
C spin-orbitals for each orbital, but also because the active orbitals
C (external order) are grouped by symmetry and then RAS space, but the
C spin orbitals are grouped by subpartition.

      IAJB=0      ! dummy initialize
      IBJA=0      ! dummy initialize

      DO IORB=1,NASHT
       IORBA=2*IORB-1
       IORBB=2*IORB
       ITABS=MAPORB(IORBA)
       DO JORB=1,NASHT
        JORBA=2*JORB-1
        JORBB=2*JORB
        JTABS=MAPORB(JORBA)
        GVAL=0.0D0
        IF(IORB.GT.JORB) THEN
         IAJB=((IORBA-1)*(IORBA-2)/2)+JORBB
         IBJA=((IORBB-1)*(IORBB-2)/2)+JORBA
        ELSE IF(JORB.EQ.IORB) THEN
         IAJB=((IORBA-1)*(IORBA-2)/2)+JORBB
         IBJA=((IORBB-1)*(IORBB-2)/2)+JORBA
         GAB=SDCHSM(IAJB)
         GBA=SDCHSM(IBJA)
         GVAL=GAB+GBA
        ELSE IF(IORB.LT.JORB) THEN
         IBJA=((JORBA-1)*(JORBA-2)/2)+IORBB
         IAJB=((JORBB-1)*(JORBB-2)/2)+IORBA
        END IF
        IJTABS=JTABS+NASHT*(ITABS-1)
        DCHSM(IJTABS)=GVAL**2
       END DO
      END DO

      CALL mma_deallocate(SDCHSM)

      RETURN
      END
