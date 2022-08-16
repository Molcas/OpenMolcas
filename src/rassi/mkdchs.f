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
*    Copyright (C) 2021, Bruno Tenorio                                 *
************************************************************************
      SUBROUTINE MKDCHS(IFSBTAB1,IFSBTAB2,ISSTAB,
     &                  MAPORB,DET1,DET2,
     &                  IF20,IF02,NDCHSM,DCHSM,
     &                  ISTATE,JSTATE)

      IMPLICIT NONE
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      INTEGER ISSTAB(*),MAPORB(*),NDCHSM
      REAL*8 DET1(*),DET2(*)
      REAL*8 DCHSM(NDCHSM)
      INTEGER NASHT,NASORB,LORBTB,IP
!      REAL*8 SGNJL,SGNIK
      REAL*8 GVAL,GAB,GBA, GAA,GBB
!      !INTEGER LSYM1,MSPOJ1,LSYM2,MSPROJ2,ISYOP,MS2OP
!      !INTEGER MPLET1,MPLET2, MPLETD

      INTEGER IAJB,IBJA, IAJA,IBJB

      INTEGER JORB,IORB 
      INTEGER JORBA,JORBB,IORBA,IORBB
      INTEGER ITABS,JTABS,IJTABS
!      INTEGER IBKA,IBKB,IJ,IJIJ,IORBA,IORBB,ITU,ITUVX
!      INTEGER IVABS,IVX,IXABS,JALA,JALB,JBJA,JBLA,JBLB
!      INTEGER JORBA,JORBB,KORB,KORBA,KORBB,LORB,LORBA,LORBB

      INTEGER NASGEM,NSDCHSM,ISTATE,JSTATE,SIGNLJ
      LOGICAL IF20,IF02
      !INTEGER job1,job2,ist,jst
#include "symmul.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
      Real*8, Allocatable:: SDCHSM(:)

C Given two CI expansions, using a biorthonormal set of SD''s,
C calculate the matrix elements relevant to DCH state intensities
C in the biorthonormal active orbital basis.
      Write(6,*)' '
      Write(6,*)'*****************************************'
      Write(6,*) 'Printout Sqr DCH matrix elements' 
      Write(6,*) 'ISTATE,JSTATE',ISTATE,JSTATE
      Write(6,*) 'I,J,|< N-2 | anni_right anni_right | N >|**2'
      Write(6,*)'*****************************************'

      LORBTB=ISSTAB(3)
C Pick out nr of active orbitals from orbital table:
      NASORB=IWORK(LORBTB+3)
      NASHT=NASORB/2
      !NASGEM=(NASORB*(NASORB-1))/2
      NSDCHSM= NASORB*(NASORB-1)/2
      Call mma_allocate(SDCHSM,nSDCHSM,Label='SDCHSM')
      SDCHSM(:)=0.0D0
      !ISYOP   = MUL(LSYM1,LSYM2)
      !MS2OP   = MSPROJ1-MSPROJ2
      !MPLETD =  MPLET1 - MPLET2

        CALL SDCHS(IWORK(LORBTB),ISSTAB,
     &              IFSBTAB1,IFSBTAB2,DET1,DET2,
     &              IF20,IF02,NSDCHSM,SDCHSM)
      
C Mapping from active spin-orbital to active orbital in external order.
C Note that these differ, not just because of the existence of two
C spin-orbitals for each orbital, but also because the active orbitals
C (external order) are grouped by symmetry and then RAS space, but the
C spin orbitals are grouped by subpartition.
      IP=0
      IAJB=0      ! dummy initialize
      IBJA=0      ! dummy initialize
      IAJA=0      ! dummy initialize
      IBJB=0      ! dummy initialize

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
         SIGNLJ=1
         IAJB=((IORBA-1)*(IORBA-2)/2)+JORBB
         IBJA=((IORBB-1)*(IORBB-2)/2)+JORBA
         IAJA=((IORBA-1)*(IORBA-2)/2)+JORBA
         IBJB=((IORBB-1)*(IORBB-2)/2)+JORBB
         !GAB=SDCHSM(IAJB)
         !GBA=SDCHSM(IBJA)
         !GAA=SDCHSM(IAJA)
         !GBB=SDCHSM(IBJB)
         !GVAL=GAB+GBA +GAA+GBB
         !Write(6,*) IORB,JORB,GVAL**2
        ELSE IF(JORB.EQ.IORB) THEN
         IAJB=((IORBA-1)*(IORBA-2)/2)+JORBB
         IBJA=((IORBB-1)*(IORBB-2)/2)+JORBA
         GAB=SDCHSM(IAJB)
         GBA=SDCHSM(IBJA)
         GVAL=GAB+GBA
!
         !I,J,DCHSM**2
         Write(6,*) IORB,JORB,GVAL**2
!
        ELSE IF(IORB.LT.JORB) THEN
         SIGNLJ=-1
         IBJA=((JORBA-1)*(JORBA-2)/2)+IORBB
         IAJB=((JORBB-1)*(JORBB-2)/2)+IORBA
         IBJB=((JORBB-1)*(JORBB-2)/2)+IORBB
         IAJA=((JORBA-1)*(JORBA-2)/2)+IORBA
         !GAB=SDCHSM(IAJB)
         !GBA=SDCHSM(IBJA)
         !GAA=SDCHSM(IAJA)
         !GBB=SDCHSM(IBJB)
         !GVAL=SIGNLJ*(GAB+GBA+ GAA+GBB)
        END IF
        IJTABS=JTABS+NASHT*(ITABS-1)
        DCHSM(IJTABS)=GVAL

       END DO
      END DO
   
      CALL mma_deallocate(SDCHSM)

      RETURN
      END
