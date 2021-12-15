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
      SUBROUTINE TSHop(CI1,CI2)
      use rassi_global_arrays, only: JBNUM, LROOT
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='TSHOP')
#include "Molcas.fh"
#include "cntrl.fh"
#include "rasdef.fh"
#include "rassi.fh"
#include "Files.fh"
#include "WrkSpc.fh"
#include "tshcntrl.fh"
      REAL*8       CI1,CI1pr,CI2,CI2pr,prdct(2,2)
      INTEGER      I,JOB1,JOB2,file,file2,maxHop,IAD3,IADR3
      CHARACTER    filnam*80,filother*80
      DIMENSION    CI1(NCI1),CI1pr(NCI1),CI2(NCI2),CI2pr(NCI2)
      DIMENSION    IADR3(3)
      LOGICAL      lMaxHop,lAllowHop,fexist,lHopped
*
*
C Skip the test if a hop has occurred
C  (this should happen when testing for state n+1
C  right after a hop to state n-1)
      CALL Get_iScalar('Relax CASSCF root',I)
      IF (I.NE.ISTATE1) THEN
         IF (IPGLOB.GE.USUAL) WRITE(6,'(6X,A)')
     &        'A hop has just been detected, skipping this state.'
         RETURN
      ENDIF
*
C Initialization
      DO i=1, NCI1
         CI1pr(I)=0.0d0
      END DO
      DO i=1, NCI2
         CI2pr(I)=0.0d0
      END DO
      file=83
      file2=83
      lHopped=.FALSE.
C
C Get the CI coefficients for current state
C
C Open JOBIPH file:
      JOB1=JBNUM(ISTATE1)
      CALL DANAME(LUIPH,JBNAME(JOB1))
C Read table of contents on this JOBIPH file:
      IAD=0
      CALL IDAFILE(LUIPH,2,ITOC15,15,IAD)
C Read CI coefficients from interface.
      IDISK=ITOC15(4)
      LROOT1=LROOT(ISTATE1)
      DO I=1,LROOT1-1
         CALL DDAFILE(LUIPH,0,CI1,NCI1,IDISK)
      END DO
      CALL DDAFILE(LUIPH,2,CI1,NCI1,IDISK)
      CALL DACLOS(LUIPH)
C
C Get the CI coefficients for state2
C
C Open JOBIPH file:
      JOB2=JBNUM(ISTATE2)
      CALL DANAME(LUIPH,JBNAME(JOB2))
C Read table of contents on this JOBIPH file:
      IAD=0
      CALL IDAFILE(LUIPH,2,ITOC15,15,IAD)
C Read CI coefficients from interface.
      IDISK=ITOC15(4)
      LROOT1=LROOT(ISTATE2)
      DO I=1,LROOT1-1
         CALL DDAFILE(LUIPH,0,CI2,NCI2,IDISK)
      END DO
      CALL DDAFILE(LUIPH,2,CI2,NCI2,IDISK)
      CALL DACLOS(LUIPH)
C
C Check if it is a hop up or hop down
C
      I=ISTATE1-ISTATE2
      IF (I.GT.0) THEN
         IF (IPGLOB.GE.USUAL)
     &      WRITE(6,*) 'Checking for a hop to a root lower in energy.'
         filnam='CIVECTOR'
         filother='CIVECTUP'
      ELSEIF (I.LT.0) THEN
         IF (IPGLOB.GE.USUAL)
     &      WRITE(6,*) 'Checking for a hop to a root higher in energy.'
         filnam='CIVECTUP'
         filother='CIVECTOR'
      ELSE
         WRITE(6,*) 'Unknown problem'
         WRITE(6,*) 'ISTATE1 = ',ISTATE1
         WRITE(6,*) 'ISTATE2 = ',ISTATE2
      ENDIF
C
C Open the file with the previous CI-vectors
C
      CALL f_inquire(filnam,fexist)
      IF (fexist) THEN
         CALL DANAME(file,filnam)
         IF (IPGLOB.GE.VERBOSE)
     &      WRITE(6,*) filnam(:mylen(filnam))//' file exists.'
      ELSE
C If the file does not exist, create a new one with the
C current vectors
         CALL DANAME(file,filnam)
C Dummy table of contents is written
         DO I=1,3
            IADR3(I)=0
         END DO
         IAD3=0
         CALL IDAFILE(file,1,IADR3,3,IAD3)
         IADR3(1)=IAD3
C Current CI coefficients are written
         CALL DDAFILE(file,1,CI1,NCI1,IAD3)
         IADR3(2)=IAD3
         CALL DDAFILE(file,1,CI2,NCI2,IAD3)
         IADR3(3)=IAD3
C Write the real table of contents
         IAD3=0
         CALL IDAFILE(file,1,IADR3,3,IAD3)
         IF (IPGLOB.GE.VERBOSE)
     &      WRITE(6,*) filnam(:mylen(filnam))//' file created.'
      ENDIF
C
C Check for surface hop if the energy difference is smaller than
C the threshold.
C
      IF (ChkHop) THEN
C
C Read table of contents on this file
         IAD3=0
         CALL IDAFILE(file,2,IADR3,3,IAD3)
C Read the CI coefficients of ISTATE1 from previuos step
         IAD3=IADR3(1)
         CALL DDAFILE(file,2,CI1pr,NCI1,IAD3)
C Read the CI coefficients of ISTATE2 from previuos step
         IAD3=IADR3(2)
         CALL DDAFILE(file,2,CI2pr,NCI2,IAD3)
C
C Calculate the scalar product of the CI coefficient vectors.
         prdct(1,1)=0.0d0
         prdct(1,2)=0.0d0
         prdct(2,1)=0.0d0
         prdct(2,2)=0.0d0
         IF (IPGLOB.GE.VERBOSE)
     &      WRITE(6,'(4(A16))') "CI1","CI1pr","CI2","CI2pr"
         DO i=1, NCI1
            IF (IPGLOB.GE.VERBOSE)
     &         WRITE(6,'(4(3X,E13.6))') CI1(i),CI1pr(i),CI2(i),CI2pr(i)
            prdct(1,1) = prdct(1,1) + CI1pr(i) * CI1(i)
            prdct(1,2) = prdct(1,2) + CI1pr(i) * CI2(i)
            prdct(2,1) = prdct(2,1) + CI2pr(i) * CI1(i)
            prdct(2,2) = prdct(2,2) + CI2pr(i) * CI2(i)
         END DO
         IF (IPGLOB.GE.USUAL) THEN
            WRITE(6,'(6X,A)')'The scalar products of the CI-vectors:'
            WRITE(6,3000)'CIpr(state1) * CI(state1) =',prdct(1,1)
            WRITE(6,3000)'CIpr(state1) * CI(state2) =',prdct(1,2)
            WRITE(6,3000)'CIpr(state2) * CI(state1) =',prdct(2,1)
            WRITE(6,3000)'CIpr(state2) * CI(state2) =',prdct(2,2)
            WRITE(6,*)
         END IF
C Check the conditions for a surface hop
         IF (ABS(prdct(1,2)).GE.0.25.AND.ABS(prdct(2,1)).GE.0.25) THEN
            WRITE(6,'(6X,80A1)')'+',('-',i=1,78),'+'
            WRITE(6,'(6X,A1,T86,A1)')'|','|'
            WRITE(6,'(6X,A1,T35,A,T86,A1)')'|',
     &           'A HOP event is detected!','|'
            WRITE(6,'(6X,A1,T86,A1)')'|','|'
            WRITE(6,'(6X,A1,T32,2(A,I3,4X),T86,A1)')'|','From state:',
     &           ISTATE1,'To state:',ISTATE2,'|'
C Check if the number of Hops is limited:
            CALL qpg_iScalar('MaxHops',lMaxHop)
            IF (lMaxHop) THEN
               CALL Get_iScalar('MaxHops',maxHop)
               lHop=.FALSE.
               CALL qpg_iScalar('Number of Hops',lHop)
               IF (lHop) THEN
                  CALL Get_iScalar('Number of Hops',nHop)
               ELSE
                  nHop=0
               END IF
               IF (maxHop.LE.nHop) THEN
                  lAllowHop=.FALSE.
                  WRITE(6,'(6X,A1,T40,A,T86,A1)') '|','maxHop > nHop',
     &            '|'
                  WRITE(6,'(6X,A1,T31,A,T86,A1)') '|',
     &            'This surface HOP is not allowed','|'
                  WRITE(6,'(6X,A1,T24,A,T86,A1)') '|',
     &            'because the number of allowed Hops is exceeded','|'
               ELSE
                  lAllowHop=.TRUE.
               END IF
            ELSE
               lAllowHop=.TRUE.
            END IF
C Bottom of the printed box
            WRITE(6,'(6X,A1,T86,A1)')'|','|'
            WRITE(6,'(6X,80A1,//)')'+',('-',i=1,78),'+'
C Set the numbers of Hops
            IF (lAllowHop) THEN
               CALL Put_iScalar('Relax CASSCF root',ISTATE2)
               CALL Put_iScalar('NumGradRoot',ISTATE2)
               nHop=nHop+1
               CALL Put_iScalar('Number of Hops',nHop)
               lHopped=.TRUE.
            END IF
         END IF
      END IF
      IF (IPGLOB.GE.VERBOSE) THEN
         WRITE(6,'(2(6X,A8,I3))')'ISTATE1=',ISTATE1,'ISTATE2=',ISTATE2
         DO i=1, NCI1
            WRITE(6,'(6X,E12.5,8X,E12.5)') CI1(i),CI2(i)
         END DO
      END IF
C
C Save the CI coefficients to the right file,
C dependending on whether or not a hop has occurred
C
      IF (lHopped) THEN
C
C Delete the file with the CI-vectors
         CALL DaEras(file)
C
C Note that, if a hop occurred, the vectors are written
C in the *other* file, and in reversed order
         CALL f_inquire(filother,fexist)
         IF (fexist) THEN
            CALL DANAME(file2,filother)
            IAD3=0
            CALL IDAFILE(file2,2,IADR3,3,IAD3)
            IAD3=IADR3(1)
            CALL DDAFILE(file2,1,CI2,NCI2,IAD3)
            IAD3=IADR3(2)
            CALL DDAFILE(file2,1,CI1,NCI1,IAD3)
            CALL DACLOS(file2)
         ELSE
C The file does not exist, create a new one
            CALL DANAME(file2,filother)
            IAD3=0
            CALL IDAFILE(file2,1,IADR3,3,IAD3)
            IADR3(1)=IAD3
            CALL DDAFILE(file2,1,CI2,NCI2,IAD3)
            IADR3(2)=IAD3
            CALL DDAFILE(file2,1,CI1,NCI1,IAD3)
            IADR3(3)=IAD3
            IAD3=0
            CALL IDAFILE(file2,1,IADR3,3,IAD3)
            CALL DACLOS(file2)
         END IF
      ELSE
C
C Write the CI-vectors normally if no hop occurred
         IAD3=0
         CALL IDAFILE(file,2,IADR3,3,IAD3)
         IAD3=IADR3(1)
         CALL DDAFILE(file,1,CI1,NCI1,IAD3)
         IAD3=IADR3(2)
         CALL DDAFILE(file,1,CI2,NCI2,IAD3)
         CALL DACLOS(file)
      ENDIF
*
      RETURN
3000  FORMAT(6X,A,F7.2)
*
      END
