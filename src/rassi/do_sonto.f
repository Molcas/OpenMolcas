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
* Copyright (C) 2021, Rulin Feng                                       *
************************************************************************
*       ****************************************************
*                           Do SO-NTO
*       ****************************************************
*        This routine is made to prepare the transition density
*        This routine is modified from do_sonatorb.f.
*        matrices, transition spin density matrices.
*        In the future, this may include antisymmetric, transition
*        densities. Namely the keyword antisin or antitrip
*
*                                                      -RF 8/18,2021
      SUBROUTINE DO_SONTO(NSS, USOR, USOI)
      use rassi_global_arrays, only: JBNUM, EIGVEC
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='DO_SONTO')
      Real*8 USOR(NSS,NSS), USOI(NSS,NSS)
      Real*8 IDENTMAT(3,3)

c Calculates natural orbitals, including spinorbit effects
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '*****************************************'
      WRITE(6,*) '* RUNNING SONTO CODE ********************'
      WRITE(6,*) '*****************************************'
      WRITE(6,*)

      IDENTMAT(:,:)=0.0D0
      FOR ALL (I=1:3) IDENTMAT(I,I)=1.0D0

      CALL GETMEM('UMATR2','ALLO','REAL',LUMATR,NSS**2)
      CALL GETMEM('UMATI2','ALLO','REAL',LUMATI,NSS**2)
      CALL GETMEM('EIGVEC2','ALLO','REAL',LVMAT,NSS**2)

      CALL DCOPY_(NSS**2,[0.0d0],0,WORK(LVMAT),1)

c transform V matrix in SF basis to spin basis
c This was taken from smmat.f and modified slightly
      ISS=0
      DO ISTATE=1,NSTATE
       JOB1=JBNUM(ISTATE)
       MPLET1=MLTPLT(JOB1)
c       S1=0.5D0*DBLE(MPLET1-1)

       DO MSPROJ1=-MPLET1+1,MPLET1-1,2
c        SM1=0.5D0*DBLE(MSPROJ1)
        ISS=ISS+1
        JSS=0

        DO JSTATE=1,NSTATE
         JOB2=JBNUM(JSTATE)
         MPLET2=MLTPLT(JOB2)
c         S2=0.5D0*DBLE(MPLET2-1)

         DO MSPROJ2=-MPLET2+1,MPLET2-1,2
c          SM2=0.5D0*DBLE(MSPROJ2)
          JSS=JSS+1

          IF (MPLET1.EQ.MPLET2 .AND. MSPROJ1.EQ.MSPROJ2) THEN
           IJ=(JSS-1)*NSS+ISS
           WORK(LVMAT-1+IJ)=EIGVEC(JSTATE,ISTATE)
          END IF ! IF (MPLET1.EQ.MPLET2 .AND. MSPROJ1.EQ.MSPROJ2)
         END DO ! DO MSPROJ2=-MPLET2+1,MPLET2-1,2
        END DO ! end DO JSTATE=1,NSTATE
       END DO ! DO MSPROJ1=-MPLET1+1,MPLET1-1,2
      END DO ! DO ISTATE=1,NSTATE

c combine this matrix with the SO eigenvector matrices
      IF(.not.NOSO) THEN
        CALL DGEMM_('N','N',NSS,NSS,NSS,
     &      1.0d0,WORK(LVMAT),NSS,USOR,NSS,0.0d0,
     &      WORK(LUMATR),NSS)
        CALL DGEMM_('N','N',NSS,NSS,NSS,
     &      1.0d0,WORK(LVMAT),NSS,USOI,NSS,0.0d0,
     &      WORK(LUMATI),NSS)
      ELSE
c Spinorbit contributions to this are disabled
        CALL DCOPY_(NSS,WORK(LVMAT),1,WORK(LUMATR),1)
        CALL DCOPY_(NSS,[0.0d0],0,WORK(LUMATI),1)
      END IF

c SONTONSTATE = number of state pairs to calculate.
c These states are stored as pairs beginning in IWORK(LSONTO)
      DO I=1,SONTOSTATES
        INTOSTATE=IWORK(LSONTO-1+I*2-1)
        JNTOSTATE=IWORK(LSONTO-1+I*2)
        WRITE(6,*)
        WRITE(6,*) "CALCULATING SO-NTOs BETWEEM SO STATES: ",
     &              INTOSTATE,JNTOSTATE
        IF((INTOSTATE.GT.NSS.OR.INTOSTATE.LE.0)
     & .or.(JNTOSTATE.GT.NSS.OR.JNTOSTATE.LE.0)) THEN
          WRITE(6,*) "...WHICH DOES NOT EXIST!"
          CALL ABEND()
        END IF
        WRITE(6,*)
        iOpt=0
c Currently only HERMISING TDMs are dealt with here
        call GETMEM('TDMAO','ALLO','REAL',LTDMAO,6*NBST**2)
        call GETMEM('TSDMAO','ALLO','REAL',LTSDMAO,6*NBST**2)
        call GETMEM('ANTSIN','ALLO','REAL',LANTSIN,6*NBST**2)
c Initialization is important
        call DCOPY_(6*NBST**2,[0.0D0],0,WORK(LTDMAO),1)
        call DCOPY_(6*NBST**2,[0.0D0],0,WORK(LTSDMAO),1)
        call DCOPY_(6*NBST**2,[0.0D0],0,WORK(LANTSIN),1)
c
        Call MAKETDMAO('HERMSING',WORK(LUMATR),WORK(LUMATI),
     &                        INTOSTATE,JNTOSTATE,NSS,iOpt,IDENTMAT,
     &                        WORK(LTDMAO))
c Following codes are left for other types of SO-TDMs
c        Call print_matrixt('TDM after MAKETDMAO 1',nbst,nbst**2,1,
c     &                    WORK(LTDMAO))
c        CALL MAKETDMAO('HERMTRIP',WORK(LUMATR),WORK(LUMATI),
c     &                        INTOSTATE,JNTOSTATE,NSS,iOpt,IDENTMAT,
c     &                        WORK(LTSDMAO),NBST)
c        Call print_matrixt('TSDM after MAKETDMAO 1',nbst,nbst**2,1,
c     &                    WORK(LTSDMAO))
c        CALL MAKETDMAO('ANTISING',WORK(LUMATR),WORK(LUMATI),
c     &                        INTOSTATE,JNTOSTATE,NSS,iOpt,IDENTMAT,
c     &                        WORK(LANTSIN),NBST)
c        Call print_matrixt('ANTITDM after MAKETDMAO 1',nbst,nbst**2,1,
c     &                    WORK(LANTSIN))
        Call DO_AOTDMNTO(WORK(LTDMAO),WORK(LTSDMAO),WORK(LANTSIN),
     &                     INTOSTATE,JNTOSTATE,NBST,NBST**2)
        call GETMEM('TDMAO','FREE','REAL',LTDMAO,6*NBST**2)
        call GETMEM('TSDMAO','FREE','REAL',LTSDMAO,6*NBST**2)
        call GETMEM('ANTSIN','FREE','REAL',LANTSIN,6*NBST**2)
      END DO
      CALL GETMEM('UMATR2','FREE','REAL',LUMATR,NSS**2)
      CALL GETMEM('UMATI2','FREE','REAL',LUMATI,NSS**2)
      CALL GETMEM('EIGVEC2','FREE','REAL',LVMAT,NSS**2)
      CALL GETMEM('SONTO','FREE','INTE',LSONTO,2*SONTOSTATES)
      RETURN
      END

