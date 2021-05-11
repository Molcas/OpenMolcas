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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_ComputeFittingCoefficients(irc)
C
C     Thomas Bondo Pedersen, June/July 2010.
C
C     Driver of the Local Density Fitting.
C
C     Atom and Atom Pair data must have been set up before calling this
C     routine.
C
#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      Implicit None
      Integer irc
#include "localdf.fh"
#include "localdf_print.fh"
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"

      Character*30 SecNam
      Parameter (SecNam='LDF_ComputeFittingCoefficients')

#if defined (_DEBUGPRINT_)
      Logical  LDF_AtomInfoIsSet, LDF_AtomPairInfoIsSet
      External LDF_AtomInfoIsSet, LDF_AtomPairInfoIsSet
#endif
      Logical  Rsv_Tsk
      External Rsv_Tsk

      Integer  LDF_OpenC, LDF_CloseC
      External LDF_OpenC, LDF_CloseC

      Logical Modified
      Logical Timing

      Integer nTask
      Parameter (nTask=7)

      Real*8 Tol2C
      Parameter (Tol2C=1.0d-14)

      Integer ID
      Integer iAtomPair
      Integer ip_C, l_C
      Integer ip_Z, l_Z
      Integer iAddr, LuC
      Integer ip_T, l_T
      Integer iTask
      Integer AB_MAE, AB_MRNrm

      Real*8 tC0, tW0, tC1, tW1
      Real*8 totC0, totW0, totC1, totW1
      Real*8 tC_Replicate, tC_Copy, tC_Constraint
      Real*8 tW_Replicate, tW_Copy, tW_Constraint
      Real*8 MAE, MRNrm

      Integer i, j, k
      Logical isUnique
      Integer iTime
      isUnique(i)=iWork(ip_AP_Unique-1+i).eq.i
      iTime(i,j,k)=ip_T-1+2*nTask*(k-1)+2*(j-1)+i

      ! Init return code
      irc=0

#if defined (_DEBUGPRINT_)
      If (.not.LDF_AtomInfoIsSet() .or.
     &    .not.LDF_AtomPairInfoIsSet()) Then
         irc=-1
         Return
      End If
#endif

      ! Init timing
      Timing=iPrint.ge.Inf_DetailedTiming
      If (Timing) Then
         l_T=2*nTask*NumberOfAtomPairs
         Call GetMem('LDFCFCT','Allo','Real',ip_T,l_T)
         Call Cho_dZero(Work(ip_T),l_T)
      Else
         ip_T=0
         l_T=0
      End If

      ! Open coefficient file
      LuC=LDF_OpenC()
      ! Set unit number in ldf_cio.fh (to allow usage of LDF_CIO_ReadC)
      Call LDF_CIO_SetUnit(LuC)
      ! Init disk address for coefficients
      iAddr=0

      ! Set info for constraints
      If (Timing) Call CWTime(tC0,tW0)
      Call LDF_SetConstraint(LDF_Constraint)
      If (Timing) Then
         Call CWTime(tC1,tW1)
         tC_Constraint=tC1-tC0
         tW_Constraint=tW1-tW0
      End If

      ! Compute fitting coefficients for each unique atom pair
      Call Init_Tsk(ID,NumberOfAtomPairs)
      Do While (Rsv_Tsk(ID,iAtomPair))
         If (isUnique(iAtomPair)) Then

            If (Timing) Call CWTime(totC0,totW0)

            ! Compute CBar and Z vectors (lower triangle):
            ! Cbar[uv,J]={(uv|J)-sum(I=1,J-1)CBar[uv,I]*Z[J,I]}/Z[J,J]
            ! where
            ! (J|K)=sum(I) Z[J,I]*Z[K,I]
            If (Timing) Call CWTime(tC0,tW0)
            Call LDF_ComputeCBar(iAtomPair,ip_C,l_C,ip_Z,l_Z,irc)
            If (irc.ne.0) Then
               Write(6,'(A,A,I8)')
     &         SecNam,': LDF_ComputeCBar returned code',irc
               irc=1
               Return
            End If
            If (iPrint.ge.Inf_AuxBas) Then
               Call Cho_Head('Auxiliary Basis Info after Initial Fit',
     &                       '-',80,6)
               Call LDF_PrintAuxBasInfo(iAtomPair)
            End If
            If (Timing) Then
               Call CWTime(tC1,tW1)
               iTask=1
               Work(iTime(1,iTask,iAtomPair))=tC1-tC0
               Work(iTime(2,iTask,iAtomPair))=tW1-tW0
            End If

            ! Update diagonal integrals
            If (Timing) Call CWTime(tC0,tW0)
            Call LDF_UpdateDiagonal(iAtomPair,l_C,Work(ip_C),irc)
            If (irc.ne.0) Then
               Write(6,'(A,A)')
     &         SecNam,
     &         ': LDF_UpdateDiagonal found too negative diagonals!'
               Write(6,'(A,I10)')
     &         'Number of too negative diagonals:',irc
               Call LDF_PrintAtomPairDiagonal(iAtomPair)
               Call WarningMessage(2,'Too negative diagonals')
               irc=1
               Return
            End If
            If (Timing) Then
               Call CWTime(tC1,tW1)
               iTask=2
               Work(iTime(1,iTask,iAtomPair))=tC1-tC0
               Work(iTime(2,iTask,iAtomPair))=tW1-tW0
            End If
            ! Print diagonal information
            If (iPrint.ge.Inf_AP) Then
               Call Cho_Head('Atom Pair Diagonal Info (1-Center Only)',
     &                       '-',80,6)
               Call LDF_PrintAtomPairDiagonal(iAtomPair)
            End If

            ! Check accuracy and include 2-center functions if needed
            ! (and requested). This is based on the updated diagonal and
            ! new CBar and Z arrays are computed.
            If (LDF2) Then
               If (Timing) Call CWTime(tC0,tW0)
               Modified=.False.
               Call LDF_Add2CenterFunctions(iAtomPair,ip_C,l_C,ip_Z,l_Z,
     &                                      Modified,irc)
               If (irc.ne.0) Then
                  Write(6,'(A,A,I8)')
     &            SecNam,': LDF_Add2CenterFunctions returned code',irc
                  irc=1
                  Return
               End If
               If (Timing) Then
                  Call CWTime(tC1,tW1)
                  iTask=3
                  Work(iTime(1,iTask,iAtomPair))=tC1-tC0
                  Work(iTime(2,iTask,iAtomPair))=tW1-tW0
               End If
               ! Update diagonal info (if aux bas was modified)
               If (Modified) Then
                  ! Restore diagonal
                  If (Timing) Call CWTime(tC0,tW0)
                  Call LDF_RestoreDiagonal(iAtomPair)
                  ! Update diagonal
                  Call LDF_UpdateDiagonal(iAtomPair,l_C,Work(ip_C),irc)
                  If (irc.ne.0) Then
                     Write(6,'(A,A)')
     &               SecNam,
     &              ': LDF_UpdateDiagonal found too negative diagonals!'
                     Write(6,'(A,I10)')
     &               'Number of too negative diagonals:',irc
                     Call LDF_PrintAtomPairDiagonal(iAtomPair)
                     Call WarningMessage(2,'Too negative diagonals')
                     irc=1
                     Return
                  End If
                  If (Timing) Then
                     Call CWTime(tC1,tW1)
                     iTask=4
                     Work(iTime(1,iTask,iAtomPair))=tC1-tC0
                     Work(iTime(2,iTask,iAtomPair))=tW1-tW0
                  End If
                  ! Print diagonal information
                  If (iPrint.ge.Inf_AP) Then
                     Call Cho_Head(
     &                    'Atom Pair Diagonal Info (2-Center Included)',
     &                    '-',80,6)
                     Call LDF_PrintAtomPairDiagonal(iAtomPair)
                  End If
               End If
            End If

            ! Compute C from CBar and Z vectors:
            ! C[uv,J]={CBar[uv,J]-sum(I=J+1,M)C[uv,I]*Z[I,J]}/Z[J,J]
            If (Timing) Call CWTime(tC0,tW0)
            Call LDF_ComputeC(iAtomPair,ip_C,l_C,ip_Z,l_Z,irc)
            If (irc.ne.0) Then
               Write(6,'(A,A,I8)')
     &         SecNam,': LDF_ComputeC returned code',irc
               irc=1
               Return
            End If
            If (iPrint.ge.Inf_AuxBas) Then
               Call Cho_Head('Auxiliary Basis Info after Final Fit',
     &                       '-',80,6)
               Call LDF_PrintAuxBasInfo(iAtomPair)
            End If
            If (Timing) Then
               Call CWTime(tC1,tW1)
               iTask=5
               Work(iTime(1,iTask,iAtomPair))=tC1-tC0
               Work(iTime(2,iTask,iAtomPair))=tW1-tW0
            End If

            ! Debug option: write unconstrained coefficients to disk
            If (WriteUnconstrainedC) Then
               Call LDF_WriteUnconstrainedCoefficients(iAtomPair,l_C,
     &                                                 Work(ip_C),irc)
               If (irc.ne.0) Then
                  Call WarningMessage(2,
     &                           SecNam//': unconstrained write failed')
                  Write(6,'(A,I10)') 'irc=',irc
                  Call LDF_Quit(1)
               End If
            End If

            ! Add constraint correction to coefficients
            If (Timing) Call CWTime(tC0,tW0)
            Call LDF_AddConstraintCorrection(LDF_Constraint,iAtomPair,
     &                                       l_C,Work(ip_C))
            If (Timing) Then
               Call CWTime(tC1,tW1)
               tC_Constraint=tC_Constraint+(tC1-tC0)
               tW_Constraint=tW_Constraint+(tW1-tW0)
            End If

            ! Update diagonal from coefficients
            If (Timing) Call CWTime(tC0,tW0)
            Call LDF_RestoreDiagonal(iAtomPair)
            Call LDF_UpdateDiagonalFromC(LDF_Constraint,iAtomPair,
     &                                   l_C,Work(ip_C),irc)
            If (irc.ne.0) Then
               Write(6,'(A,A)')
     &         SecNam,
     &         ': LDF_UpdateDiagonalFromC found too negative diagonals!'
               Write(6,'(A,I10)')
     &         'Number of too negative diagonals:',irc
               Call LDF_PrintAtomPairDiagonal(iAtomPair)
               Call WarningMessage(2,'Too negative diagonals')
               irc=1
               Return
            End If
            If (Timing) Then
               Call CWTime(tC1,tW1)
               iTask=4
               Work(iTime(1,iTask,iAtomPair))=
     &                          Work(iTime(1,iTask,iAtomPair))+(tC1-tC0)
               Work(iTime(2,iTask,iAtomPair))=
     &                          Work(iTime(2,iTask,iAtomPair))+(tW1-tW0)
            End If
            ! Print diagonal information
            If (iPrint.ge.Inf_AP) Then
               Call Cho_Head('Atom Pair Diagonal Info (Final)','-',80,6)
               Call LDF_PrintAtomPairDiagonal(iAtomPair)
            End If

            ! Debug: check integrals for this atom pair
            If (CheckPairIntegrals) Then
               If (LDF_Constraint.eq.-1) Then
                  ! unconstrained LDF, use nonrobust LDF integrals
                  Call LDF_CheckPairIntegrals(2,iAtomPair,
     &                                        l_C,Work(ip_C),irc)
               Else
                  ! constrained LDF, use robust LDF integrals
                  Call LDF_CheckPairIntegrals(1,iAtomPair,
     &                                        l_C,Work(ip_C),irc)
               End If
               If (irc.ne.0) Then
                  Write(6,'(A,A,I8)')
     &            SecNam,': LDF_CheckPairIntegrals returned code',irc
                  irc=1
                  Return
               End If
            End If

            ! Zero negative diagonals
            Call LDF_CleanDiagonal(iAtomPair)

            ! Debug: verify fit for this atom pair
            If (VerifyFit) Then
               Call LDF_VerifyFit(.True.,iPrint.lt.3,LDF_Constraint,
     &                            1.0d-10,iAtomPair,l_C,Work(ip_C),irc)
               If (irc.ne.0) Then
                  Write(6,'(A,A,I8)')
     &            SecNam,': LDF_VerifyFit returned code',irc
                  irc=1
                  Return
               End If
            End If

            ! Save fitting coefficients on disk, and save disk address
            If (Timing) Call CWTime(tC0,tW0)
            Call LDF_WriteC(iAtomPair,l_C,Work(ip_C),LuC,iAddr)
            If (Timing) Then
               Call CWTime(tC1,tW1)
               iTask=6
               Work(iTime(1,iTask,iAtomPair))=tC1-tC0
               Work(iTime(2,iTask,iAtomPair))=tW1-tW0
            End If

            ! Deallocate C and Z
            Call GetMem('LDFC','Free','Real',ip_C,l_C)
            Call GetMem('ZVec','Free','Real',ip_Z,l_Z)

            If (Timing) Then
               Call CWTime(totC1,totW1)
               iTask=7
               Work(iTime(1,iTask,iAtomPair))=totC1-totC0
               Work(iTime(2,iTask,iAtomPair))=totW1-totW0
            End If

         End If
      End Do
      Call Free_Tsk(ID)

      ! Parallel: replicate data on all nodes.
      If (Timing) Call CWTime(tC0,tW0)
      Call LDF_ReplicateData(LuC,irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': LDF_ReplicateData returned code',irc
         irc=1
         Return
      End If
      If (Timing) Then
         Call CWTime(tC1,tW1)
         tC_Replicate=tC1-tC0
         tW_Replicate=tW1-tW0
      End If

      ! Copy unique atom pair info
      If (Timing) Call CWTime(tC0,tW0)
      Call LDF_CopyUniqueAtomPairs(irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': LDF_CopyUniqueAtomPairs returned code',irc
         irc=1
         Return
      End If
      If (Timing) Then
         Call CWTime(tC1,tW1)
         tC_Copy=tC1-tC0
         tW_Copy=tW1-tW0
      End If

      ! Print atom pair info
      If (iPrint.ge.Inf_AP) Then
         Call LDF_PrintAtomPairInfo()
      End If

      ! Check overlap integrals
      If (CheckOverlapIntegrals) Then
#if defined (_MOLCAS_MPP_)
         If (nProcs.gt.1 .and. Is_Real_Par()) Then
            ! Start coefficient IO utility
            Call LDF_CIO_SetUnit(0)
            Call LDF_CIO_Init(0.0d0,irc)
            If (irc.ne.0) Then
               Call WarningMessage(2,SecNam//' LDF_CIO_Init failed')
               Write(6,'(A,I10)') 'Return code=',irc
               Call LDF_Quit(1)
            End If
         End If
#endif
         Call LDF_CheckAllOverlapIntegrals(iPrint.ge.3,Tol2C,
     &                                     MAE,AB_MAE,MRNrm,AB_MRNrm)
         If (iPrint.lt.3) Then
            Call Cho_Head('Results from overlap integral check','-',80,
     &                    6)
            Write(6,'(/,2X,A,1P,D20.10,2X,A)')
     &      'Tolerance for 2C errors..',Tol2C,'(all 2C passed)'
            Write(6,'(2X,A,1P,D20.10,2X,A,I10)')
     &      'Max abs error............',MAE,'@AB=',AB_MAE
            Write(6,'(2X,A,1P,D20.10,2X,A,I10)')
     &      'Max relative norm error..',MRNrm,'@AB=',AB_MRNrm
            Call xFlush(6)
         End If
         If (LDF_Constraint.eq.0) Then ! charge constraint => ovlp exact
            If (MAE.gt.1.0d-12) Then
               Call WarningMessage(2,
     &          SecNam//': MAE too large with Charge Constraint active')
               Call LDF_Quit(1)
            End If
         End If
#if defined (_MOLCAS_MPP_)
         If (nProcs.gt.1 .and. Is_Real_Par()) Then
            ! Shut down coefficient IO utility
            Call LDF_CIO_Final()
            Call LDF_CIO_SetUnit(LuC)
         End If
#endif
      End If

      ! Unset info for constraints
      If (Timing) Call CWTime(tC0,tW0)
      Call LDF_UnsetConstraint(LDF_Constraint)
      If (Timing) Then
         Call CWTime(tC1,tW1)
         tC_Constraint=tC_Constraint+(tC1-tC0)
         tW_Constraint=tW_Constraint+(tW1-tW0)
      End If

      ! Close coefficient file
      irc=LDF_CloseC(LuC)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': LDF_CloseC returned code',irc
         irc=1
         Return
      End If
      Call LDF_CIO_SetUnit(0)

      ! Print timing
      If (Timing) Then
         Write(6,'(/,A)')
     &   'Detailed Timing of LDF Fitting Coefficient Computation:'
         Write(6,'(/,A,3X,A,3X,A,4X,A,3X,A,7X,A,9X,A,3X,A)')
     &   'Atom Pair','Compute CBar','Upd. Diagonal','Add 2C Func.',
     &   '2C Upd. Diag.','Compute C','Write C','Total AP Time'
         Write(6,'(8X,7A16)') ('     CPU    Wall',i=1,7)
         Write(6,'(120A1)') ('-',i=1,120)
         Do iAtomPair=1,NumberOfAtomPairs
            Write(6,'(I8,14(1X,F7.1))')
     &      iAtomPair,(Work(ip_T+2*nTask*(iAtomPair-1)+i),i=0,2*nTask-1)
         End Do
         Write(6,'(120A1)') ('-',i=1,120)
         Do iAtomPair=2,NumberOfAtomPairs
            Do iTask=1,nTask
               Work(iTime(1,iTask,1))=Work(iTime(1,iTask,1))
     &                               +Work(iTime(1,iTask,iAtomPair))
               Work(iTime(2,iTask,1))=Work(iTime(2,iTask,1))
     &                               +Work(iTime(2,iTask,iAtomPair))
            End Do
         End Do
         Write(6,'(A,14F8.1)') 'Sum     ',(Work(ip_T+i),i=0,2*nTask-1)
         Write(6,'(120A1)') ('-',i=1,120)
         Write(6,'(A,2(1X,F7.1))')
     &   'Constraints (CPU,Wall in s).........',
     &   tC_Constraint,tW_Constraint
         Write(6,'(A,2(1X,F7.1))')
     &   'Data Replication (CPU,Wall in s)....',
     &   tC_Replicate,tW_Replicate
         Write(6,'(A,2(1X,F7.1))')
     &   'Copy Unique Pairs (CPU,Wall in s)...',
     &   tC_Copy,tW_Copy
         Call xFlush(6)
         Call GetMem('LDFCFCT','Free','Real',ip_T,l_T)
      End If

      End
