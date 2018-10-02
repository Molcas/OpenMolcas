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
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CMO1(NCMO),CMO2(NCMO),FOCKMO(NGAM1),TUVX(NGAM2)
      DIMENSION KEEP(8),NBSX(8),ipAsh(2)
      LOGICAL   ISQARX
#include "rassi.fh"
#include "symmul.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "para_info.fh"
      Logical IfTest,FoundTwoEls,DoCholesky

      Integer ALGO,Nscreen
      Real*8  dmpk
      Logical Fake_CMO2
      Real*8, Dimension(:), Allocatable:: Prod, FAO, DInAO

      COMMON / CHO_JOBS / Fake_CMO2
      COMMON /CHORASSI / ALGO,Nscreen,dmpk

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

      Call qEnter('TRINT')
      IfTest=.False.
#ifdef _DEBUG_
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
            CALL ERRTRA
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
              CALL ERRTRA
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
      Call mma_allocate(DInAO,nDInAO,Label='DInAO')
      CALL DIMAT(CMO1,CMO2,DINAO)

      NFAO=NBSQ
      Call mma_allocate(FAO,nFAO,Label='FAO')
      LFAO=ip_of_Work(FAO)

      IF (.not.DoCholesky) THEN

         If ( IfTest ) Call dVcPrt('Done',' ',DINAO,NDINAO)
C GET THE ONE-ELECTRON HAMILTONIAN MATRIX FROM ONE-EL FILE AND
C PUT IT INTO A FOCK MATRIX IN AO BASIS:
C Note: GETH1 also adds the reaction field contribution to the
C 1-electron hamiltonian, and the variable ERFNuc in common /general/,
C which is the RF contribution to the nuclear repulsion
         CALL GETH1_RASSI(FAO)
         If ( IfTest ) Call dVcPrt('h0',' ',FAO,NFAO)
C ONE CONTRIBUTION TO ECORE MUST BE CALCULATED FROM THE NAKED
C ONE-EL HAMILTONIAN:
         ECORE1=DDOT_(NBSQ,FAO,1,DINAO,1)
         If ( IfTest ) Write (6,*) '      ECore1=',ECORE1


C ADD IN THE TWO-ELECTRON CONTRIBUTIONS TO THE FOCKAO MATRIX:
*
         Call Allocate_Work(ipTemp,nFAO)
         Call FZero(Work(ipTemp),nFAO)
         CALL FOCK(DINAO,WORK(ipTemp))

c --- FAO already contains the one-electron part
         Call DaXpY_(nFAO,1.0D+0,Work(ipTemp),1,FAO,1)
         Call Free_Work(ipTemp)

#ifdef _DEBUG_
         ioff=1
         Do i=1,nSym
         call CHO_OUTPUT(fao(ioff),1,nBasF(i),1,nBasF(i),
     &                   nBasF(i),nBasF(i),1,6)
         ioff=ioff+nBasF(i)**2
         End Do
#endif

         ECORE2=DDOT_(NBSQ,FAO,1,DINAO,1)

      Else

* --- Initialize Cholesky information

         CALL CHO_X_INIT(irc,ChFracMem)
         if (irc.ne.0) then
            write(6,*)'TrInt: Cho_X_Init returns error code ',irc
            Call AbEnd()
         endif

         CALL GETMEM('DLT','ALLO','REAL',ipDLT,nBTri)
         CALL Fold_Mat(nSym,nBasF,DINAO,Work(ipDLT))

#ifdef _DEBUG_

         ioff1=0
         ioff2=0
         do i=1,nSym
            ipDL = ipDLT + ioff1
            ipDS = 1 + ioff2
            call cho_output(DInAO(ipDS),1,nBasF(i),1,nBasF(i),
     &                      nBasF(i),nBasF(i),1,6)
            call triprt('DLT','',Work(ipDL),nBasF(i))
            ioff1=ioff1+nBasF(i)*(nBasF(i)+1)/2
            ioff2=ioff2+nBasF(i)**2
         end do

#endif

         CALL GETMEM('FLT','ALLO','REAL',ipFLT,nBTri)

C GET THE ONE-ELECTRON HAMILTONIAN MATRIX FROM ONE-EL FILE AND
C PUT IT INTO A FOCK MATRIX IN AO BASIS:
C Note: CHO_GETH1 also adds the reaction field contribution to the
C 1-electron hamiltonian, and the variable ERFNuc in common /general/,
C which is the RF contribution to the nuclear repulsion

         CALL CHO_GETH1(nBtri,ipFLT,RFpert,ERFNuc)

         ECORE1=DDOT_(nBtri,WORK(ipFLT),1,WORK(ipDLT),1)
         If ( IfTest ) Write (6,*) '      ECore1=',ECORE1,ALGO
         If ( IfTest ) Write (6,*) '      FAKE_CMO2=',FAKE_CMO2

#if defined (_MOLCAS_MPP_)
         If (nProcs.gt.1) Then
             scx=1.0/dble(nProcs)
C --- to avoid double counting when using gadsum
             call dscal_(nBtri,scx,WORK(ipFLT),1)
         EndIf
#endif

         call Fzero(FAO,nFAO) ! Used as Exchange F matrix

         If (ALGO.eq.1) Then
c --- reorder the MOs to fit Cholesky needs

           nKB=0
           Do iSym=1,nSym
            nKB = nKB + (nIsh(iSym)+nAsh(iSym))*nBasF(iSym)
           End Do
           Call GetMem('Cka','Allo','Real',ipMO1,nKB)
           Call GetMem('Ckb','Allo','Real',ipMO2,nKB)

           ioff=0
           Do iSym=1,nSym

            do ikk=1,nIsh(iSym)

               ioff1=ioff+nBasF(iSym)*(ikk-1)

               call dcopy_(nBasF(iSym),CMO1(ioff1+1),1,
     &                    Work(ipMO1+ioff+ikk-1),nIsh(iSym))

               call dcopy_(nBasF(iSym),CMO2(ioff1+1),1,
     &                    Work(ipMO2+ioff+ikk-1),nIsh(iSym))
            end do

            ioff2=ioff+nBasF(iSym)*nIsh(iSym)

            do ikk=1,nAsh(iSym)

               ioff3=ioff2+nBasF(iSym)*(ikk-1)

               call dcopy_(nBasF(iSym),CMO1(ioff3+1),1,
     &                    Work(ipMO1+ioff2+ikk-1),nAsh(iSym))

               call dcopy_(nBasF(iSym),CMO2(ioff3+1),1,
     &                    Work(ipMO2+ioff2+ikk-1),nAsh(iSym))
            end do

            ioff=ioff+nBasF(iSym)*(nIsh(iSym)+nAsh(iSym))

           End Do


c --- Add the two-electron contribution to the Fock matrix
c ---     and compute the (tu|vx) integrals

           If (Fake_CMO2) Then
              CALL CHO_FOCK_RASSI(ipDLT,ipMO1,ipMO2,ipFLT,LTUVX)
           Else
              CALL CHO_FOCK_RASSI_X(ipDLT,ipMO1,ipMO2,ipFLT,LFAO,LTUVX)
           EndIf

           Call GetMem('Ckb','Free','Real',ipMO2,nKB)
           Call GetMem('Cka','Free','Real',ipMO1,nKB)

         Else  ! algo=2 (local exchange algorithm)

           ipMO1 = ip_of_Work(CMO1(1))  ! Cak storage
           ipMO2 = ip_of_Work(CMO2(1))

C *** Only the active orbitals MO coeff need reordering
           nVB=0
           Do iSym=1,nSym
            nVB = nVB + nAsh(iSym)*nBasF(iSym)
           End Do
           Call GetMem('Cva','Allo','Real',ipAsh(1),nVB)
           Call GetMem('Cvb','Allo','Real',ipAsh(2),nVB)

           ioff=0
           ioff1=0
           Do iSym=1,nSym

            ioff2 = ioff + nBasF(iSym)*nIsh(iSym)

            do ikk=1,nAsh(iSym)

               ioff3=ioff2+nBasF(iSym)*(ikk-1)

               call dcopy_(nBasF(iSym),CMO1(ioff3+1),1,
     &                 Work(ipAsh(1)+ioff1+ikk-1),nAsh(iSym))

               call dcopy_(nBasF(iSym),CMO2(ioff3+1),1,
     &                    Work(ipAsh(2)+ioff1+ikk-1),nAsh(iSym))

            end do

            ioff=ioff+(nIsh(iSym)+nAsh(iSym))*nBasF(iSym)
            ioff1=ioff1+nAsh(iSym)*nBasF(iSym)

           End Do

c --- Add the two-electron contribution to the Fock matrix
c ---     and compute the (tu|vx) integrals

           If (Fake_CMO2) Then

             CALL CHO_LK_RASSI(ipDLT,ipMO1,ipMO2,ipFLT,LFAO,LTUVX,
     &                         ipAsh,nScreen,dmpk)
           Else

             CALL GetMem('K-mat','Allo','Real',ipK,NBSQ)
             Call FZero(Work(ipK),NFAO)

             CALL CHO_LK_RASSI_X(ipDLT,ipMO1,ipMO2,ipFLT,ipK,LFAO,LTUVX,
     &                         ipAsh,nScreen,dmpk)

             CALL GetMem('K-mat','Free','Real',ipK,NBSQ)
           EndIf

           Call GetMem('Cva','Free','Real',ipAsh(1),nVB)
           Call GetMem('Cvb','Free','Real',ipAsh(2),nVB)

         EndIf

         If (Fake_CMO2) Then
            ioff1=0
            ioff2=1
            Do i=1,nSym
               CALL SQUARE(Work(ipFLT+ioff1),FAO(ioff2),
     &                     1,nBasF(i),nBasF(i))
               ioff1=ioff1+nBasF(i)*(nBasF(i)+1)/2
               ioff2=ioff2+nBasF(i)**2
            End Do
         EndIf

         CALL GETMEM('FLT','Free','REAL',ipFLT,nBTri)
         CALL GETMEM('DLT','Free','REAL',ipDLT,nBTri)

         Call GADSUM(FAO,NBSQ)
         Call GADSUM(TUVX,NGAM2)

         ECORE2=DDOT_(NBSQ,FAO,1,DINAO,1)

* --- Finalize Cholesky information

         CALL CHO_X_FINAL(irc)
         if (irc.ne.0) then
            write(6,*)'TrInt: Cho_X_Final returns error code ',irc
            write(6,*)'Try recovery -- continue.'
         endif

      EndIf


      If ( IfTest ) Write (6,*) '      Etwo  =',ECORE2
      ECORE=0.5D0*(ECORE1+ECORE2)
      If ( IfTest ) Write (6,*) '      Ecore =',ECORE
C (NOTE COMPENSATION FOR DOUBLE-COUNTING OF TWO-ELECTRON CONTRIBUTION.)
      Call mma_deallocate(DInAO)


C TRANSFORM THE FOCK MATRIX TO MO BASIS:
      NPROD=0
      DO 20 I=1,NSYM
        NPROD=MAX(NPROD,NASH(I)*NBASF(I))
20    CONTINUE
*
      Call mma_allocate(Prod,nProd,Label='Prod')
      Call FZero(PROD,NPROD)
      ISTFMO=1
      ISTFAO=1
      ISTC=1
      CALL FZERO(FOCKMO,NGAM1)
      DO 30 ISYM=1,NSYM
         NA=NASH(ISYM)
         NI=NISH(ISYM)
         NO=NI+NA
         NB=NBASF(ISYM)
         IF (NA.NE.0) THEN
            ISTA=ISTC+NI*NB
C ISTFAO: BEGINNING OF F(AO,AO) BLOCK OF SYMMETRY ISYM.
C ISTFMO: BEGINNING OF F(MO,MO) BLOCK OF SYMMETRY ISYM.
C ISTC: BEGINNING OF MO-S OF SYMMETRY ISYM.
C ISTA: BEGINNING OF ACTIVE MO-S OF SYMMETRY ISYM.
C MATRIX MULT. PROD(AO,ACTIVE MO)=F(AO,AO)*CMO2(AO,ACTIVE MO)
            CALL DGEMM_('N','N',NB,NA,NB,
     &                  1.0D0,FAO(ISTFAO),NB,
     &                        CMO2(ISTA),NB,
     &                  0.0D0,PROD,NB)
C -- MATRIX MULT. F(ACT MO,ACT MO)=CMO1(AO,ACT MO)(TRP)*PROD(AO,ACT MO)
            CALL DGEMM_('T','N',NA,NA,NB,
     &                  1.0D0,CMO1(ISTA),NB,
     &                        PROD,NB,
     &                  0.0D0,FOCKMO(ISTFMO),NASHT)
         END IF
         ISTFAO=ISTFAO+NB**2
         ISTFMO=ISTFMO+NA*(NASHT+1)
         ISTC=ISTC+NO*NB
30    CONTINUE
      Call mma_deallocate(Prod)
      Call mma_deallocate(FAO)

      IOPT=0
      IRC=0

      If (.not.DoCholesky) then
C TRANSFORM TWO-ELECTRON INTEGRALS:
         CALL TRAINT(CMO1,CMO2,NGAM2,TUVX)

         CALL CLSORD(IRC,IOPT)

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


      Call qExit('TRINT')
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
