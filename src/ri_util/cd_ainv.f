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
      Subroutine CD_AInv(A,n,AInV,Thr_CD)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Real*8 A(n,n), AInv(n,n)
*
      Call GetMem('ADiag','Allo','Real',ipADiag,n)
      Call GetMem('iADiag','Allo','Inte',ipiADiag,n)
*
      iSeed=77
      Lu_A=IsFreeUnit(iSeed)
      Call DaName_MF_WA(Lu_A,'AMat09')
*
      iDisk=0
      Call dDaFile(Lu_A,1,A,n**2,iDisk)
*
C     Call RecPrt('A',' ',A,n,n)
*
      iSeed=iSeed+1
      Lu_Q=IsFreeUnit(iSeed)
      Call DaName_MF_WA(Lu_Q,'QMat09')
*
      call dcopy_(n,A,n+1,Work(ipADiag),1)
C     Call RecPrt('ADiag',' ',Work(ipADiag),1,n)
*
      Call CD_AInv_(n,m,Work(ipADiag),iWork(ipiADiag),Lu_A,Lu_Q,Thr_CD)
*
      Call GetMem('iADiag','Free','Inte',ipiADiag,n)
      Call GetMem('ADiag','Free','Real',ipADiag,n)
*
      Call GetMem('QVec','Allo','Real',ipQVec,n*m)
*
      iDisk=0
      Call dDaFile(Lu_Q,2,Work(ipQVec),n*m,iDisk)
*
C     Call RecPrt('QVec','(6G20.10)',Work(ipQVec),n,m)
      Call DGEMM_('N','T',n,n,m,
     &            1.0D0,Work(ipQVec),n,
     &                  Work(ipQVec),n,
     &            0.0D0,AInv,n)
C     Call RecPrt('AInv',' ',AInv,n,n)
      Call DaEras(Lu_Q)
      Call GetMem('QVec','Free','Real',ipQVec,n*n)
*                                                                      *
************************************************************************
*                                                                      *
*     Check the accuracy I-AA^1
*
#ifdef _ACCURACY_
      Call Allocate_Work(ipTmp,n*n)
*---
      Call FZero(Work(ipTmp),n*n)
*     I
      call dcopy_(n,1.0D0,0,Work(ipTmp),n+1)
*     I-AA^-1
      Call DGEMM_('N','N',n,n,n,
     &           -1.0D0,A,n,
     &                  AInv,n,
     &            1.0D0,Work(ipTmp),n)
      Call RecPrt('I-AA^-1','(6G20.12)',Work(ipTmp),n,n)
*
      Call DGEMM_('N','N',n,n,n,
     &            1.0D0,A,n,
     &                  AInv,n,
     &            0.0D0,Work(ipTmp),n)
      Call Allocate_Work(ipTmp2,n*n)
      Call FZero(Work(ipTmp2),n*n)
      call dcopy_(n,1.0D0,0,Work(ipTmp2),n+1)
      Call DGEMM_('N','N',n,n,n,
     &           -1.0D0,Work(ipTmp),n,
     &                  Work(ipTmp),n,
     &            1.0D0,Work(ipTmp2),n)
      Call RecPrt('I-AA^-1AA^-1','(6G20.12)',Work(ipTmp2),n,n)
*---
      Call Free_Work(ipTmp2)
      Call Free_Work(ipTmp)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
      Subroutine CD_AInv_(n,m,ADiag,iADiag,Lu_A,Lu_Q,Thr_CD)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "real.fh"
      Real*8 ADiag(n)
      Integer iADiag(n)
      Logical Out_of_Core
*
      nScr=3*n
      Call GetMem('AMax','Max','Real',iDummy,MaxMem)
      lScr=Min(MaxMem,nScr)
      Call GetMem('AScr','Allo','Real',ipScr,lScr)
*
      nDim=n
*
      Thr=Thr_CD*1.0D-1
      Lu_Z=7
      Call DaName_MF_WA(Lu_Z,'ZMAT09')
      Call Get_Pivot_idx(ADiag,nDim,nVec,Lu_A,Lu_Z,iADiag,Work(ipScr),
     &                   lScr,Thr)
      m=nVec
      If (nDim.ne.nVec) Then
         Write (6,*)
         Write (6,*) 'Detected lin. dep. in the auxiliary basis'
         Write (6,'(A,I6)')' # of aux. bfns before lin. dep. removal: ',
     &                     nDim
         Write (6,'(A,I6)')' # of aux. bfns after  lin. dep. removal: ',
     &                     nVec
      End If
*
      Call Pivot_Mat(nDim,nVec,Lu_A,Lu_Z,iADiag,Work(ipScr),lScr)
*
      Call GetMem('AScr','Free','Real',ipScr,lScr)
*
************************************************************************
*     A-vectors are now on disk. Go ahead and compute the Q-vectors!
************************************************************************
*
      ThrQ=Thr_CD*1.0D-1 ! Threshold for Inv_Cho_Factor
*
      nB=nVec
      If (nB.eq.0) Go To 777
      nQm=nB*(nB+1)/2
*
      nXZ=nB
      nQm_full= nB*(nB+1)/2
*
      Out_of_Core=2*nQm_full+5*nXZ.gt.MaxMem
*
      If (Out_Of_Core) Then
         mQm=(nQm*MaxMem-5*nXZ)/(2*nQm_full)
         a=One
         b=-Two*DBLE(mQm)
         mB=INT(-a/Two + Sqrt( (a/Two)**2 - b ))
         kQm=mB*(mB+1)/2
         If (kQm.gt.mQm) Then
            Call WarningMessage(2,'Error in CD_AInv')
            Write (6,*) 'kQm.gt.mQm!'
            Write (6,*) 'MaxMem=',MaxMem
            Write (6,*) 'nQm,mQm,kQm=',nQm,mQm,kQm
            Write (6,*) 'nB,mB=',nB,mB
            Call Abend()
         End If
      Else
         mB = nB
         kQm = nQm
      End If
*
      lQm=kQm
      lAm=lQm
*
      If (lQm.lt.1) Then
         Call WarningMessage(2,'Error in CD_AInv')
         Write (6,*) 'lQm.lt.1'
         Call Abend()
      End If
*
*     Some of memory for scratch arrays for Inv_Cho_Factor
*     Allocate memory for the A- and Q-vectors and initialize.
*
      Call GetMem('ICF','Allo','Real',ip_ICF,lAm+lQm+5*nXZ)
      ipZ    = ip_ICF
      ipX    = ip_ICF + 1*nXZ
      ip_Q_k = ip_ICF + 2*nXZ
      ip_A_k = ip_ICF + 3*nXZ
      ip_Scr = ip_ICF + 4*nXZ
      lScr=nXZ
      ip_Am  = ip_ICF + 5*nXZ
      ip_Qm  = ip_ICF + 5*nXZ + lAm
*
      Call FZero(Work(ip_Am),lAm)
      Call FZero(Work(ip_Qm),lQm)
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Process the A_ks to generate Q_ks.
*
      iAddr=0
      nMem=mB
      Do kCol = 1, nB
*
         iAddr_=iAddr
         If (kCol.le.nMem) Then
*           Point to A_k in Am
            iOff = (kCol-1)*kCol/2
            ip_A_l=ip_Am+iOff
            If (kCol.eq.1) Then
               nAm=nMem*(nMem+1)/2
               Call dDaFile(Lu_Z,2,Work(ip_Am),nAm,iAddr_)
            End If
*           Point to Q_k in Qm
            iOff = (kCol-1)*kCol/2
            ip_Q_l=ip_Qm+iOff
         Else If (kCol.gt.nMem) Then
*           Use special scratch for A_k
            ip_A_l=ip_A_k
            Call dDaFile(Lu_Z,2,Work(ip_A_l),kCol,iAddr_)
*           Use special scratch for Q_k
            ip_Q_l=ip_Q_k
         End If
*
         LinDep=2
         Call Inv_Cho_Factor(Work(ip_A_l),kCol,
     &                       Work(ip_Am),Work(ip_Qm),nMem,
     &                       Lu_Z,Lu_Q,
     &                       Work(ip_Scr),lScr,
     &                       Work(ipZ),Work(ipX),ThrQ,
     &                       Work(ip_Q_l),LinDep)
*
         If (LinDep.ne.0) Then
            Call WarningMessage(2,'Error in CD_AInv')
            Write(6,*) 'Inv_Cho_Factor found linear dependence!'
            Call Abend()
         End If
*
*        Write the new A/Q-vector to file
*
         iAddr_=iAddr
         If (kCol.eq.nMem) Then
            nQm=kCol*(kCol+1)/2
            Call dDaFile(Lu_Q,1,Work(ip_Qm),nQm,iAddr )
            Call dDaFile(Lu_Z,1,Work(ip_Am),nQm,iAddr_)
         Else If (kCol.gt.nMem) Then
            Call dDaFile(Lu_Q,1,Work(ip_Q_k),kCol,iAddr )
            Call dDaFile(Lu_Z,1,Work(ip_A_l),kCol,iAddr_)
         End If
*
      End Do
*
      Call GetMem('ICF','Free','Real',ip_ICF,lAm+lQm+5*nXZ)
 777  Continue
      Call DaEras(Lu_Z)
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Sort the Q-matrix back to the original order.
*
      Call GetMem('MemMax','Max','Real',iDummy,MaxMem2)
*
      nBfnTot=n
      nBfn2=n**2
      lScr=Min(MaxMem2,Max(nBfn2,2*nBfnTot))
      Call GetMem('Scr','Allo','Real',ip_Scr,lScr)
*
      Call Restore_Mat(nDim,nVec,Lu_Q,Lu_A,iADiag,Work(ip_Scr),lScr,
     &                 .true.)
      Call DaEras(Lu_Q)
      Lu_Q=Lu_A
*
*     Note: after the 'Restore' call to Sort_mat, the Q-matrix is
*           no longer stored as upper-triangular but as squared
*           (zeros have been added because the corresponding argument
*            is set to .true.). The column index is still pivoted.
*
      Call GetMem('Scr','Free','Real',ip_Scr,lScr)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Return
      End
