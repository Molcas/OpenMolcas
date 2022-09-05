!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine  Gen_QVec(nIrrep,nBas_Aux)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Integer nBas_Aux(0:nIrrep-1), Lu_Q(0:7), Lu_A(0:7)
      Logical Out_Of_Core
      Character*6 Name_Q

      Real*8, Allocatable, Target :: Mem(:)
      Real*8, Pointer:: A_l(:)=>Null(), Q_l(:)=>Null()
      Real*8, Pointer:: A_k(:)=>Null(), Q_k(:)=>Null()
      Real*8, Pointer:: Am(:)=>Null(), Qm(:)=>Null()
      Integer, Allocatable:: iDiag(:)
      Real*8, Allocatable:: Scr(:), X(:), Z(:)
!                                                                      *
!***********************************************************************
!                                                                      *
      INTERFACE
      SUBROUTINE SORT_mat(irc,nDim,nVec,iD_A,nSym,lu_A0,mode,lScr,Scr,  &
     &                    Diag)
      Integer irc
      Integer nSym
      Integer nDim(nSym)
      Integer nVec(nSym)
      Integer iD_A(*)
      Integer lu_A0(nSym)
      Character(LEN=7) mode
      Integer lScr
      Real*8  Scr(lScr)
      Real*8, Optional ::  Diag(*)
      END SUBROUTINE SORT_mat
      END INTERFACE
!                                                                      *
!***********************************************************************
!                                                                      *
!
      ThrQ=1.0D-14 ! Threshold for Inv_Cho_Factor
!
      mB=0
      nA_Diag=0
      Do iIrrep = 0, nIrrep-1
         nB = nBas_Aux(iIrrep)
         nA_Diag = nA_Diag + nB
         mB=Max(mB,nB)
!
         iSeed=55+iIrrep
         Lu_Q(iIrrep)=IsFreeUnit(iSeed)
         Write(Name_Q,'(A4,I2.2)') 'QMAT',iIrrep
         Call DaName_MF_WA(Lu_Q(iIrrep),Name_Q)
!
         iSeed=63+iIrrep
         Lu_A(iIrrep)=IsFreeUnit(iSeed)
         Write(Name_Q,'(A4,I2.2)') 'AVEC',iIrrep
         Call DaName_MF_WA(Lu_A(iIrrep),Name_Q)
      End Do
      nBfn2=mB**2
!
      Call mma_allocate(Z,mB,Label='Z')
      Call mma_allocate(X,mB,Label='X')
      lScr=3*mB
      Call mma_allocate(Scr,lScr,Label='Scr')
      Call mma_maxDBLE(Mem_Max)
      Call mma_allocate(Mem,Mem_Max,Label='Mem')
!
      Do iIrrep = 0, nIrrep-1
         nB = nBas_Aux(iIrrep)
         nQm=nB*(nB+1)/2
!
         Out_Of_Core = 2*nQm .gt. Mem_Max
         If (Out_Of_Core) Then
            MaxMem = Mem_Max - 2*nB
            mQm=MaxMem/2
            a=One
            b=-Two*DBLE(mQm)
            mB=INT(-a/Two + Sqrt( (a/Two)**2 - b ))
            kQm=mB*(mB+1)/2
            If (kQm.gt.mQm) Then
               Call WarningMessage(2,'Error in Gen_QVec')
               Write (6,*) 'kQm.gt.mQm!'
               Write (6,*) 'MaxMem=',MaxMem
               Write (6,*) 'nQm,mQm,kQm=',nQm,mQm,kQm
               Write (6,*) 'nB,mB=',nB,mB
               Call Abend()
            End If
            iE = 2*kQm
            iS = iE + 1
            iE = iE + mB
            Q_k(1:mB) => Mem(iS:iE)
            iS = iE + 1
            iE = iE + mB
            A_k(1:mB) => Mem(iS:iE)
         Else
            mB = nB
            kQm = nQm
         End If
!
         iS = 1
         iE = kQm
         Qm(1:kQm) => Mem(iS:iE)
         iS = iE + 1
         iE = iE + kQm
         Am(1:kQm) => Mem(iS:iE)
!
         iAddr=0
         Do kCol = 1, nB
!
            If (kCol.le.mB) Then
               iOff = (kCol-1)*kCol/2
               A_l(1:) => Am(1+iOff:)
            Else
               A_l(1:) => A_k(1:)
            End If
!
            iAddr_ = iAddr
            If (kCol.le.mB.and.kCol.eq.1) Then
               Call dDaFile(Lu_A(iIrrep),2,Am,kQm,iAddr_)
            Else If (kCol.gt.mB) Then
               Call dDaFile(Lu_A(iIrrep),2,A_l,kCol,iAddr_)
            End If
#ifdef _DEBUGPRINT_
            Write (6,*) 'kCol=',kCol
            Call TriPrt('Am',' ',Am,mB)
            Call RecPrt('Al',' ',A_l,1,kCol)
#endif
!
            If (kCol.le.mB) Then
               iOff = (kCol-1)*kCol/2
               Q_l(1:) => Qm(1+iOff:)
            Else
               Q_l(1:) => Q_k(1:)
            End If
!
            LinDep=2
            Call Inv_Cho_Factor(A_l,kCol,                               &
     &                          Am,Qm,mB,                               &
     &                          Lu_A(iIrrep),Lu_Q(iIrrep),              &
     &                          Scr,lScr,                               &
     &                          Z,X,ThrQ,                               &
     &                          Q_l,LinDep)

            If (LinDep.ne.0) Then
               Call WarningMessage(2,'Error in Gen_QVec')
               Write(6,*) 'Inv_Cho_Factor found linear dependence!'
               Call Abend()
            End If
#ifdef _DEBUGPRINT_
            Call TriPrt('Qm',' ',Qm,Min(mB,kCol))
            Call RecPrt('Ql',' ',Q_l,1,kCol)
#endif
!
!           Write the new A/Q-vector to file
!
            iAddr_=iAddr
            If (kCol.eq.mB) Then
               lQm=kCol*(kCol+1)/2
               Call dDaFile(Lu_Q(iIrrep),1,Qm,lQm,iAddr )
               Call dDaFile(Lu_A(iIrrep),1,Am,lQm,iAddr_)
            Else If (kCol.gt.mB) Then
               nQ_k=kCol
               Call dDaFile(Lu_Q(iIrrep),1,Q_l,nQ_k,iAddr )
               Call dDaFile(Lu_A(iIrrep),1,A_l,nQ_k,iAddr_)
            End If
!
         End Do    ! kCol
         Call DaClos(Lu_A(iIrrep))
      End Do       ! iIrrep
!
      A_l=>Null()
      Q_l=>Null()
      A_k=>Null()
      Q_k=>Null()
      Am=>Null()
      Qm=>Null()
      Call mma_deallocate(Mem)
      Call mma_deallocate(Scr)
      Call mma_deallocate(X)
      Call mma_deallocate(Z)
!
!     Sort the Q-matrix to square storage.
!
      Call mma_allocate(iDiag,nA_Diag,Label='iDiag')
      ik = 0
      Do iIrrep=0,nIrrep-1
         Do k=1,nBas_Aux(iIrrep)
            ik = ik + 1
            iDiag(ik) = k  ! dummy assignement
         End Do
      End Do
      Call mma_maxDBLE(MaxMem2)
      lScr=Min(MaxMem2,nBfn2)
      Call mma_allocate(Scr,lScr,Label='Scr')
!
      Call SORT_Mat(irc,nBas_Aux,nBas_Aux,                              &
     &              iDiag,nIrrep,Lu_Q,'Restore',                        &
     &              lScr,Scr)
!
!     Note: after the 'Restore' call to Sort_mat, the Q-matrix is
!           no longer stored as upper-triangular but as squared
!           (zeros have been added).
!
      Call mma_deallocate(Scr)
      Call mma_deallocate(iDiag)
!
      Do iIrrep=0,nIrrep-1
         Call DaClos(Lu_Q(iIrrep))
      End Do
!
      Return
      End
