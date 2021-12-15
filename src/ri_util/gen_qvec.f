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
      Subroutine  Gen_QVec(nIrrep,nBas_Aux)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Integer nBas_Aux(0:nIrrep-1), Lu_Q(0:7), Lu_A(0:7)
      Logical Out_Of_Core
      Character*6 Name_Q
*
      ThrQ=1.0D-14 ! Threshold for Inv_Cho_Factor
*
      mB=0
      nA_Diag=0
      Do iIrrep = 0, nIrrep-1
         nB = nBas_Aux(iIrrep)
         nA_Diag = nA_Diag + nB
         mB=Max(mB,nB)
*
         iSeed=55+iIrrep
         Lu_Q(iIrrep)=IsFreeUnit(iSeed)
         Write(Name_Q,'(A4,I2.2)') 'QMAT',iIrrep
         Call DaName_MF_WA(Lu_Q(iIrrep),Name_Q)
*
         iSeed=63+iIrrep
         Lu_A(iIrrep)=IsFreeUnit(iSeed)
         Write(Name_Q,'(A4,I2.2)') 'AVEC',iIrrep
         Call DaName_MF_WA(Lu_A(iIrrep),Name_Q)
      End Do
      nBfn2=mB**2
*
      Call Allocate_Work(ipZ,mB)
      Call Allocate_Work(ipX,mB)
      lScr=3*mB
      Call Allocate_Work(ip_Scr,lScr)
      Call GetMem('MemX','Max','Real',iDum,Mem_Max)
      Call GetMem('MemX','Allo','Real',ip_Mem,Mem_Max)
*
      Do iIrrep = 0, nIrrep-1
         nB = nBas_Aux(iIrrep)
         nQm=nB*(nB+1)/2
*
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
            ip_Q_k = ip_Mem + 2*kQm
            ip_A_k = ip_Mem + 2*kQm + mB
         Else
            mB = nB
            kQm = nQm
            ip_Q_k = ip_Dummy
            ip_A_k = ip_Dummy
         End If
*
         ip_Qm=ip_Mem
         ip_Am=ip_Mem + kQm
         Call FZero(Work(ip_Qm),kQm)
         Call FZero(Work(ip_Am),kQm)
*
         iAddr=0
         Do kCol = 1, nB
*
            If (kCol.le.mB) Then
               iOff = (kCol-1)*kCol/2
               ip_A_l = ip_Am + iOff
            Else
               ip_A_l = ip_A_k
            End If
*
            iAddr_ = iAddr
            If (kCol.le.mB.and.kCol.eq.1) Then
               Call dDaFile(Lu_A(iIrrep),2,Work(ip_Am),kQm,iAddr_)
            Else If (kCol.gt.mB) Then
               Call dDaFile(Lu_A(iIrrep),2,Work(ip_A_l),kCol,iAddr_)
            End If
#ifdef _DEBUGPRINT_
            Write (6,*) 'kCol=',kCol
            Call TriPrt('Am',' ',Work(ip_Am),mB)
            Call RecPrt('Al',' ',Work(ip_A_l),1,kCol)
#endif
*
            If (kCol.le.mB) Then
               iOff = (kCol-1)*kCol/2
               ip_Q_l = ip_Qm + iOff
            Else
               ip_Q_l = ip_Q_k
            End If
*
            LinDep=2
            Call Inv_Cho_Factor(Work(ip_A_l),kCol,
     &                          Work(ip_Am),Work(ip_Qm),mB,
     &                          Lu_A(iIrrep),Lu_Q(iIrrep),
     &                          Work(ip_Scr),lScr,
     &                          Work(ipZ),Work(ipX),ThrQ,
     &                          Work(ip_Q_l),LinDep)

            If (LinDep.ne.0) Then
               Call WarningMessage(2,'Error in Gen_QVec')
               Write(6,*) 'Inv_Cho_Factor found linear dependence!'
               Call Abend()
            End If
#ifdef _DEBUGPRINT_
            Call TriPrt('Qm',' ',Work(ip_Qm),Min(mB,kCol))
            Call RecPrt('Ql',' ',Work(ip_Q_l),1,kCol)
#endif
*
*           Write the new A/Q-vector to file
*
            iAddr_=iAddr
            If (kCol.eq.mB) Then
               lQm=kCol*(kCol+1)/2
               Call dDaFile(Lu_Q(iIrrep),1,Work(ip_Qm),lQm,iAddr )
               Call dDaFile(Lu_A(iIrrep),1,Work(ip_Am),lQm,iAddr_)
            Else If (kCol.gt.mB) Then
               nQ_k=kCol
               Call dDaFile(Lu_Q(iIrrep),1,Work(ip_Q_l),nQ_k,iAddr )
               Call dDaFile(Lu_A(iIrrep),1,Work(ip_A_l),nQ_k,iAddr_)
            End If
*
         End Do    ! kCol
         Call DaClos(Lu_A(iIrrep))
      End Do       ! iIrrep
*
      Call GetMem('MemX','Free','Real',ip_Mem,Mem_Max)
      Call Free_Work(ip_Scr)
      Call Free_Work(ipX)
      Call Free_Work(ipZ)
*
*     Sort the Q-matrix to square storage.
*
      Call GetMem('iDiag','Allo','Inte',ip_iDiag,nA_Diag)
      ik = ip_iDiag
      Do iIrrep=0,nIrrep-1
         Do k=1,nBas_Aux(iIrrep)
            iWork(ik) = k  ! dummy assignement
            ik = ik + 1
         End Do
      End Do
      Call GetMem('MemMax','Max','Real',iDummy,MaxMem2)
      lScr=Min(MaxMem2,nBfn2)
      Call GetMem('Scr','Allo','Real',ip_Scr,lScr)
*
      Call SORT_Mat(irc,Work(ip_Dummy),nBas_Aux,nBas_Aux,
     &              iWork(ip_iDiag),nIrrep,Lu_Q,'Restore',
     &              lScr,Work(ip_Scr))
*
*     Note: after the 'Restore' call to Sort_mat, the Q-matrix is
*           no longer stored as upper-triangular but as squared
*           (zeros have been added).
*
      Call GetMem('Scr','Free','Real',ip_Scr,lScr)
      Call GetMem('iD_Diag','Free','Inte',ip_iDiag,nA_Diag)
*
      Do iIrrep=0,nIrrep-1
         Call DaClos(Lu_Q(iIrrep))
      End Do
*
      Return
      End
