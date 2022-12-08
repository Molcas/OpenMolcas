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

      SUBROUTINE TRINT(CMO1,CMO2,ECORE,NGAM1,FOCKMO,NGAM2,TUVX)
      USE Fock_util_global, only: Fake_CMO2
#if defined (_MOLCAS_MPP_)
      USE Para_Info, ONLY: nProcs
#endif
      use Data_structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 CMO1(NCMO),CMO2(NCMO),FOCKMO(NGAM1),TUVX(NGAM2)
      Integer KEEP(8),NBSX(8), nAux(8)
      LOGICAL   ISQARX
      Type (DSBA_Type) Ash(2), MO1(2), MO2(2), DLT, FLT(1), KSQ,
     &                 FAO, Temp_SQ, DInAO
#include "real.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "stdalloc.fh"
      Logical IfTest,FoundTwoEls,DoCholesky

      Real*8, Dimension(:), Allocatable:: Prod

#include "chorassi.fh"

*****************************************************************
*  CALCULATE AND RETURN ECORE, FOCKMO, AND TUVX. ECORE IS THE
*  INACTIVE-INACTIVE ENERGY CONTRIBUTION, WHICH IS TO BE MULTI-
*  PLIED WITH AN OVERLAP TO GIVE A CONTRIBUTION TO THE HAMILTONIAN
*  MATRIX ELEMENT. SIMILARLY, THE INACTIVE-ACTIVE CONTRIBUTION IS
*  GIVEN BY THE FOCKMO ARRAY CONTRACTED WITH AN ACTIVE TRANSITION
*  DENSITY MATRIX, AND THE ACTIVE-ACTIVE IS THE TRANSFORMED
*  INTEGRALS IN ARRAY TUVX CONTRACTED WITH THE TRANSITION DENSITY
*  TWO-ELECTRON MATRIX. THEREFORE, THE STORAGE OF THE FOCKMO AND
*  THE TUVX MATRICES ARE IN THE SAME FORMAT AS THE DENSITY MATRICES.
*****************************************************************

      IfTest=.False.
#ifdef _DEBUGPRINT_
      IfTest=.True.
#endif
C THE FOLLOWING PROGRAMS USE THE ORDERED INTEGRAL FILE FOR BOTH
C THE FOCKMO MATRIX AND THE TUVX ARRAY. THEREFORE, IT IS IMPERATIVE
C THAT SYMMETRY BLOCKS HAVE NOT BEEN EXCLUDED IN THE GENERATION OF
C THIS ORDERED INTEGRAL FILE, UNLESS BOTH ACTIVE AND INACTIVE ORBITALS
C ARE MISSING FOR THAT SYMMETRY LABEL. THIS IS CHECKED FIRST:

C OPEN THE ELECTRON REPULSION INTEGRAL FILE
      Call f_Inquire('ORDINT',FoundTwoEls)
c      Call DecideOnDirect(.False.,FoundTwoEls,DoDirect,DoCholesky)
      Call DecideOnCholesky(DoCholesky)

      If (.not.DoCholesky) then

         IOPT=0
         IRC=0
         CALL OPNORD(IRC,IOPT,'ORDINT',LUORD)
          IF(IRC.NE.0) THEN
             WRITE(6,*)' ERROR: RETURN CODE=',IRC
             WRITE(6,*)' RECEIVED WHEN OPENING ORDINT.'
             CALL ABEND()
          END IF

          IRC=0
          CALL GETORD(IRC,ISQARX,NSYMX,NBSX,KEEP)
          IF(IRC.NE.0) THEN
            WRITE(6,*)' ERROR: RETURN CODE=',IRC
            WRITE(6,*)' RECEIVED WHEN CALLING GETORD.'
            CALL ABEND()
          END IF
          IF ( NSYMX.NE.NSYM ) THEN
            WRITE(6,*)' '
            WRITE(6,*)'     *** ERROR IN SUBROUTINE TRINT ***'
            WRITE(6,*)'     INCOMPATIBLE NUMBERS OF IRRED. REP.'
            WRITE(6,*)' '
            CALL ABEND()
          END IF
          DO ISYM=1,NSYM
            NB1=NBASF(ISYM)
            NB2=NBSX(ISYM)
            IF ( NB1.NE.NB2 ) THEN
              WRITE(6,*)' '
              WRITE(6,*)'     *** ERROR IN SUBROUTINE TRINT ***'
              WRITE(6,*)'   INCOMPATIBLE NUMBERS OF BASIS FUNCTION'
              WRITE(6,*)' '
              CALL ABEND()
            END IF
          END DO
          DO 10 I=1,NSYM
            IF(KEEP(I).EQ.0) GOTO 10
            IF(NOSH(I).GT.0) GOTO 901
10        CONTINUE

      EndIf


C CALCULATE AN INACTIVE TRANSITION DENSITY MATRIX IN AO BASIS:
      NDINAO=NBSQ
      Call Allocate_DT(DInAO,nBasF,nBasF,nSym)
      CALL DIMAT(CMO1,CMO2,DINAO%A0)

      NFAO=NBSQ
      Call Allocate_DT(FAO,nBasF,nBasF,nSym)
*                                                                      *
************************************************************************
*                                                                      *
      IF (.not.DoCholesky) THEN     ! Conventional integrals
*                                                                      *
************************************************************************
*                                                                      *

         If ( IfTest ) Call dVcPrt('Done',' ',DINAO%A0,NDINAO)
C GET THE ONE-ELECTRON HAMILTONIAN MATRIX FROM ONE-EL FILE AND
C PUT IT INTO A FOCK MATRIX IN AO BASIS:
C Note: GETH1 also adds the reaction field contribution to the
C 1-electron hamiltonian, and the variable ERFNuc in common /general/,
C which is the RF contribution to the nuclear repulsion
         CALL GETH1_RASSI(FAO%A0)
         If ( IfTest ) Call dVcPrt('h0',' ',FAO%A0,NFAO)
C ONE CONTRIBUTION TO ECORE MUST BE CALCULATED FROM THE NAKED
C ONE-EL HAMILTONIAN:
         ECORE1=DDOT_(NBSQ,FAO%A0,1,DINAO%A0,1)
         If ( IfTest ) Write (6,*) '      ECore1=',ECORE1


C ADD IN THE TWO-ELECTRON CONTRIBUTIONS TO THE FOCKAO MATRIX:
*
         Call Allocate_DT(Temp_SQ,nBasF,nBasF,nSym)
         Temp_SQ%A0(:)=Zero
         CALL FOCK_RASSI(DINAO%A0,Temp_SQ%A0)

c --- FAO already contains the one-electron part
         FAO%A0(:) = FAO%A0(:) + Temp_SQ%A0(:)
         Call Deallocate_DT(Temp_SQ)

#ifdef _DEBUGPRINT_
         Do i=1,nSym
         call CHO_OUTPUT(FAO%SB(i)%A2,1,nBasF(i),1,nBasF(i),
     &                                  nBasF(i),nBasF(i),1,6)
         End Do
#endif

         ECORE2=DDOT_(NBSQ,FAO%A0,1,DINAO%A0,1)
*                                                                      *
************************************************************************
*                                                                      *
      Else       ! RI/CD integrals
*                                                                      *
************************************************************************
*                                                                      *
* ------ Initialize Cholesky information

         CALL CHO_X_INIT(irc,ChFracMem)
         if (irc.ne.0) then
            write(6,*)'TrInt: Cho_X_Init returns error code ',irc
            Call AbEnd()
         endif

         Call Allocate_DT(DLT,nBasF,nBasF,nSym,aCase='TRI')
         CALL Fold_Mat(nSym,nBasF,DINAO%A0,DLT%A0)

#ifdef _DEBUGPRINT_

         do i=1,nSym
            call cho_output(DInAO%SB(i)%A2,1,nBasF(i),1,nBasF(i),
     &                      nBasF(i),nBasF(i),1,6)
            call triprt('DLT','',DLT%SB(i)%A1,nBasF(i))
         end do

#endif

         Call Allocate_DT(FLT(1),nBasF,nBasF,nSym,aCase='TRI')

C GET THE ONE-ELECTRON HAMILTONIAN MATRIX FROM ONE-EL FILE AND
C PUT IT INTO A FOCK MATRIX IN AO BASIS:
C Note: CHO_GETH1 also adds the reaction field contribution to the
C 1-electron hamiltonian, and the variable ERFNuc in common /general/,
C which is the RF contribution to the nuclear repulsion

         CALL CHO_GETH1(nBtri,FLT(1)%A0,RFpert,ERFNuc)

         ECORE1=DDOT_(nBtri,FLT(1)%A0,1,DLT%A0,1)
         If ( IfTest ) Write (6,*) '      ECore1=',ECORE1,ALGO
         If ( IfTest ) Write (6,*) '      FAKE_CMO2=',FAKE_CMO2

#if defined (_MOLCAS_MPP_)
         If (nProcs.gt.1) Then
             scx=1.0/dble(nProcs)
C --- to avoid double counting when using gadsum
             FLT(1)%A0(:) = scx * FLT(1)%A0(:)
         EndIf
#endif

         FAO%A0(:)=Zero ! Used as Exchange F matrix
!                                                                      !
!)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()(!
!                                                                      !
         If (ALGO.eq.1) Then
!                                                                      !
!)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()(!
!                                                                      !
c -------- reorder the MOs to fit Cholesky needs

           Call Allocate_DT(MO1(1),nIsh,nBasF,nSym)
           Call Allocate_DT(MO1(2),nIsh,nBasF,nSym)
           Call Allocate_DT(MO2(1),nAsh,nBasF,nSym)
           Call Allocate_DT(MO2(2),nAsh,nBasF,nSym)

           ioff=0
           Do iSym=1,nSym

            do ikk=1,nIsh(iSym)

               ioff1=ioff+nBasF(iSym)*(ikk-1)

               MO1(1)%SB(iSym)%A2(ikk,:) =
     &            CMO1(ioff1+1:ioff1+nBasF(iSym))

               MO1(2)%SB(iSym)%A2(ikk,:) =
     &            CMO2(ioff1+1:ioff1+nBasF(iSym))
            end do

            ioff2=ioff+nBasF(iSym)*nIsh(iSym)

            do ikk=1,nAsh(iSym)

               ioff3=ioff2+nBasF(iSym)*(ikk-1)

               MO2(1)%SB(iSym)%A2(ikk,:) =
     &            CMO1(ioff3+1:ioff3+nBasF(iSym))

               MO2(2)%SB(iSym)%A2(ikk,:) =
     &            CMO2(ioff3+1:ioff3+nBasF(iSym))
            end do

            ioff=ioff+nBasF(iSym)*(nIsh(iSym)+nAsh(iSym))

           End Do


c --- Add the two-electron contribution to the Fock matrix
c ---     and compute the (tu|vx) integrals

           If (Fake_CMO2) Then
              CALL CHO_FOCK_RASSI(DLT,MO1,MO2,FLT,TUVX)
           Else
              CALL CHO_FOCK_RASSI_X(DLT,MO1,MO2,FLT,FAO,TUVX)
           EndIf

           Call Deallocate_DT(MO2(2))
           Call Deallocate_DT(MO2(1))
           Call Deallocate_DT(MO1(2))
           Call Deallocate_DT(MO1(1))
!                                                                      !
!)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()(!
!                                                                      !
         Else  ! algo=2 (local exchange algorithm)
!                                                                      !
!)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()(!
!                                                                      !

           nAux(:) = nIsh(:) + nAsh(:)
           Call Allocate_DT(MO1(1),nBasF,nAux,nSym,Ref=CMO1)
           Call Allocate_DT(MO1(2),nBasF,nAux,nSym,Ref=CMO2)

C *** Only the active orbitals MO coeff need reordering
           Call Allocate_DT(Ash(1),nAsh,nBasF,nSym)
           Call Allocate_DT(Ash(2),nAsh,nBasF,nSym)

           Do iSym=1,nSym

            do ikk=1,nAsh(iSym)
               jkk = nIsh(iSym) + ikk

               Ash(1)%SB(iSym)%A2(ikk,:) =
     &             MO1(1)%SB(iSym)%A2(:,jkk)

               Ash(2)%SB(iSym)%A2(ikk,:) =
     &             MO1(2)%SB(iSym)%A2(:,jkk)

            end do

           End Do

c --- Add the two-electron contribution to the Fock matrix
c ---     and compute the (tu|vx) integrals

           If (Fake_CMO2) Then

             CALL CHO_LK_RASSI(DLT,MO1,FLT,FAO,TUVX,Ash,nScreen,dmpk)
           Else

             Call Allocate_DT(KSQ,nBasF,nBasF,nSym)
             KSQ%A0(:)=Zero

             CALL CHO_LK_RASSI_X(DLT,MO1,FLT,KSQ,FAO,TUVX,Ash,nScreen,
     6                           dmpk)

             Call Deallocate_DT(KSQ)
           EndIf

           Call Deallocate_DT(Ash(2))
           Call Deallocate_DT(Ash(1))
           Call Deallocate_DT(MO1(2))
           Call Deallocate_DT(MO1(1))

!                                                                      !
!)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()(!
!                                                                      !
         EndIf
!                                                                      !
!)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()(!
!                                                                      !

         If (Fake_CMO2) Then
            Do i=1,nSym
               CALL SQUARE(FLT(1)%SB(i)%A1,FAO%SB(i)%A2,
     &                     1,nBasF(i),nBasF(i))
            End Do
         EndIf

         Call Deallocate_DT(FLT(1))
         Call Deallocate_DT(DLT)

         Call GADSUM(FAO%A0,NBSQ)
         Call GADSUM(TUVX,NGAM2)

         ECORE2=DDOT_(NBSQ,FAO%A0,1,DINAO%A0,1)

* --- Finalize Cholesky information

         CALL CHO_X_FINAL(irc)
         if (irc.ne.0) then
            write(6,*)'TrInt: Cho_X_Final returns error code ',irc
            write(6,*)'Try recovery -- continue.'
         endif
*                                                                      *
************************************************************************
*                                                                      *
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
      If ( IfTest ) Write (6,*) '      Etwo  =',ECORE2
      ECORE=0.5D0*(ECORE1+ECORE2)
      If ( IfTest ) Write (6,*) '      Ecore =',ECORE
C (NOTE COMPENSATION FOR DOUBLE-COUNTING OF TWO-ELECTRON CONTRIBUTION.)
      Call Deallocate_DT(DInAO)


C TRANSFORM THE FOCK MATRIX TO MO BASIS:
      NPROD=0
      DO 20 I=1,NSYM
        NPROD=MAX(NPROD,NASH(I)*NBASF(I))
20    CONTINUE
*
      Call mma_allocate(Prod,nProd,Label='Prod')
      Call FZero(PROD,NPROD)
      ISTFMO=1
      ISTC=1
      CALL FZERO(FOCKMO,NGAM1)
      DO 30 ISYM=1,NSYM
         NA=NASH(ISYM)
         NI=NISH(ISYM)
         NO=NI+NA
         NB=NBASF(ISYM)
         IF (NA.NE.0) THEN
            ISTA=ISTC+NI*NB
C ISTFMO: BEGINNING OF F(MO,MO) BLOCK OF SYMMETRY ISYM.
C ISTC: BEGINNING OF MO-S OF SYMMETRY ISYM.
C ISTA: BEGINNING OF ACTIVE MO-S OF SYMMETRY ISYM.
C MATRIX MULT. PROD(AO,ACTIVE MO)=F(AO,AO)*CMO2(AO,ACTIVE MO)
            CALL DGEMM_('N','N',NB,NA,NB,
     &                  1.0D0,FAO%SB(ISYM)%A2,NB,
     &                        CMO2(ISTA),NB,
     &                  0.0D0,PROD,NB)
C -- MATRIX MULT. F(ACT MO,ACT MO)=CMO1(AO,ACT MO)(TRP)*PROD(AO,ACT MO)
            CALL DGEMM_('T','N',NA,NA,NB,
     &                  1.0D0,CMO1(ISTA),NB,
     &                        PROD,NB,
     &                  0.0D0,FOCKMO(ISTFMO),NASHT)
         END IF
         ISTFMO=ISTFMO+NA*(NASHT+1)
         ISTC=ISTC+NO*NB
30    CONTINUE
      Call mma_deallocate(Prod)
      Call Deallocate_DT(FAO)


      If (.not.DoCholesky) then
C TRANSFORM TWO-ELECTRON INTEGRALS:
         CALL TRAINT(CMO1,CMO2,NGAM2,TUVX)

         IRC=0
         CALL CLSORD(IRC)

      End If
      Call Chk4NaN(nasht*(nasht+1)/2,TUVX,iErr)
      If (iErr.ne.0) Then
         Write (6,*) 'TrInt: TUVX corrupted'
         Call Abend()
      End If
      Call Chk4NaN(nGam1,FOCKMO,iErr)
      If (iErr.ne.0) Then
         Write (6,*) 'TrInt: FOCKMO corrupted'
         Call Abend()
      End If
c     Call triprt('tuvx',' ',TUVX,nasht)

      RETURN
901   CONTINUE
      WRITE(6,*)' ERROR IN KEEP PARAMETERS ON ORDINT FILE.'
      WRITE(6,'(A,8I5)')' KEEP ARRAY:',(KEEP(I),I=1,NSYM)
      WRITE(6,'(A,8I5)')' NASH ARRAY:',(NASH(I),I=1,NSYM)
      WRITE(6,*)' A NON-ZERO KEEP PARAMETER INDICATES THAT BASIS'
      WRITE(6,*)' FUNCTIONS WITH A SPECIFIC SYMMETRY LABEL HAS'
      WRITE(6,*)' BEEN SKIPPED. THIS IS POSSIBLE ONLY IF INACTIVE'
      WRITE(6,*)' AND ACTIVE ORBITALS OF THAT SYMMETRY ARE MISSING.'
      CALL ABEND()
      END
