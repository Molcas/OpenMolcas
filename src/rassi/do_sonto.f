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
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer NSS
      Real*8 USOR(NSS,NSS), USOI(NSS,NSS)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
      Real*8 IDENTMAT(3,3)
      Real*8, Allocatable:: UMATR(:), UMATI(:), VMAT(:,:)
      Real*8, Allocatable:: TDMAO(:), TSDMAO(:)

c Calculates natural orbitals, including spinorbit effects
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '*****************************************'
      WRITE(6,*) '* RUNNING SONTO CODE ********************'
      WRITE(6,*) '*****************************************'
      WRITE(6,*)

      IDENTMAT(:,:)=0.0D0
      FOR ALL (I=1:3) IDENTMAT(I,I)=1.0D0

      CALL mma_allocate(UMATR,NSS**2,Label='UMATR')
      CALL mma_allocate(UMATI,NSS**2,Label='UMATI')
      CALL mma_allocate(VMAT,NSS,NSS,Label='VMAT')

      VMAT(:,:)=0.0D0

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
           VMAT(ISS,JSS)=EIGVEC(JSTATE,ISTATE)
          END IF ! IF (MPLET1.EQ.MPLET2 .AND. MSPROJ1.EQ.MSPROJ2)
         END DO ! DO MSPROJ2=-MPLET2+1,MPLET2-1,2
        END DO ! end DO JSTATE=1,NSTATE
       END DO ! DO MSPROJ1=-MPLET1+1,MPLET1-1,2
      END DO ! DO ISTATE=1,NSTATE

c combine this matrix with the SO eigenvector matrices
      IF(.not.NOSO) THEN
        CALL DGEMM_('N','N',NSS,NSS,NSS,
     &      1.0d0,VMAT,NSS,USOR,NSS,0.0d0,
     &      UMATR,NSS)
        CALL DGEMM_('N','N',NSS,NSS,NSS,
     &      1.0d0,VMAT,NSS,USOI,NSS,0.0d0,
     &      UMATI,NSS)
      ELSE
c Spinorbit contributions to this are disabled
        CALL DCOPY_(NSS,VMAT,1,UMATR,1)
        CALL DCOPY_(NSS,[0.0d0],0,UMATI,1)
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
        call mma_allocate(TDMAO,6*NBST**2,Label='TDMAO')
        call mma_allocate(TSDMAO,6*NBST**2,Label='TSDMAO')
        call GETMEM('ANTSIN','ALLO','REAL',LANTSIN,6*NBST**2)
c Initialization is important
        TDMAO(:)=0.0D0
        TSDMAO(:)=0.0D0
        call DCOPY_(6*NBST**2,[0.0D0],0,WORK(LANTSIN),1)
c
        Call MAKETDMAO('HERMSING',UMATR,UMATI,
     &                        INTOSTATE,JNTOSTATE,NSS,iOpt,IDENTMAT,
     &                        TDMAO)
c Following codes are left for other types of SO-TDMs
c        Call print_matrixt('TDM after MAKETDMAO 1',nbst,nbst**2,1,
c     &                    TDMAO)
c        CALL MAKETDMAO('HERMTRIP',UMATR,UMATI,
c     &                        INTOSTATE,JNTOSTATE,NSS,iOpt,IDENTMAT,
c     &                        TSDMAO,NBST)
c        Call print_matrixt('TSDM after MAKETDMAO 1',nbst,nbst**2,1,
c     &                    TSDMAO)
c        CALL MAKETDMAO('ANTISING',UMATR,UMATI,
c     &                        INTOSTATE,JNTOSTATE,NSS,iOpt,IDENTMAT,
c     &                        WORK(LANTSIN),NBST)
c        Call print_matrixt('ANTITDM after MAKETDMAO 1',nbst,nbst**2,1,
c     &                    WORK(LANTSIN))
        Call DO_AOTDMNTO(TDMAO,TSDMAO,WORK(LANTSIN),
     &                     INTOSTATE,JNTOSTATE,NBST,NBST**2)
        Call mma_deallocate(TDMAO)
        Call mma_deallocate(TSDMAO)
        call GETMEM('ANTSIN','FREE','REAL',LANTSIN,6*NBST**2)
      END DO
      Call mma_deallocate(UMATR)
      Call mma_deallocate(UMATI)
      Call mma_deallocate(VMAT)
      CALL GETMEM('SONTO','FREE','INTE',LSONTO,2*SONTOSTATES)

      END SUBROUTINE DO_SONTO

