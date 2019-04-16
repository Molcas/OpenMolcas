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
      SubRoutine Cho_Reorder_RI(Vec,lVec,nVec,iSym)
      Implicit Real*8 (a-h,o-z)
      Real*8 Vec(lVec,nVec)
#include "cholesky.fh"
#include "choorb.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
      iRS2F(i,j)=iWork(ip_iRS2F-1+2*(j-1)+i)

      If (nVec .lt. 1) Return
      If (lVec .lt. 1) Return
      If (lVec.ne.nnBstR(iSym,1) .or. nVec.gt.NumCho(iSym)) Then
         Call SysAbendMsg('Cho_Reorder_RI','Input argument error!',' ')
      End If
      If (nnShl .ne. nnShl_Tot) Then
         Call SysAbendMsg('Cho_Reorder_RI','Screening is not allowed!',
     &                    '(nnShl.ne.nnShl_Tot)')
      End If

C     Set mapping from global address to reduced set.
C     -----------------------------------------------

      liF2RS = nBasT*(nBasT+1)/2
      Call GetMem('CR_RI_F2RS','Allo','Inte',ipiF2RS,liF2RS)
      Call iCopy(liF2RS,[0],0,iWork(ipiF2RS),1)
      Do iRS = 1,nnBstR(iSym,1)
         iRS_tot = iiBstR(iSym,1) + iRS
         na = iRS2F(1,iRS_tot)
         nb = iRS2F(2,iRS_tot)
         nab = iTri(na,nb)
         iWork(ipiF2RS-1+nab) = iRS
      End Do

C     Reorder.
C     --------

      lScr = lVec
      Call GetMem('CR_RI_Scr','Allo','Real',ipScr,lScr)
      Do iVec = 1,nVec

         Call dCopy_(lVec,Vec(1,iVec),1,Work(ipScr),1)
         kFrom = ipScr - 1
         Do iSymb = 1,nSym

            iSyma = MulD2h(iSymb,iSym)

            If (iSyma .gt. iSymb) Then
               Do ib = 1,nBas(iSymb)
                  nb = iBas(iSymb) + ib
                  Do ia = 1,nBas(iSyma)
                     na = iBas(iSyma) + ia
                     nab = iTri(na,nb)
                     iRS = iWork(ipiF2RS-1+nab)
#if defined (_DEBUG_)
                     If (iRS.lt.1 .or. iRS.gt.nnBstR(iSym,1)) Then
                        Call SysAbendMsg('Cho_Reorder_RI',
     &                                   'Index out of bounds',' ')
                     End If
#endif
                     kFrom = kFrom + 1
                     Vec(iRS,iVec) = Work(kFrom)
                  End Do
               End Do
            Else If (iSyma .eq. iSymb) Then
               Do ia = 1,nBas(iSyma)
                  na = iBas(iSyma) + ia
                  Do ib = 1,ia
                     nb = iBas(iSymb) + ib
                     nab = iTri(na,nb)
                     iRS = iWork(ipiF2RS-1+nab)
#if defined (_DEBUG_)
                     If (iRS.lt.1 .or. iRS.gt.nnBstR(iSym,1)) Then
                        Call SysAbendMsg('Cho_Reorder_RI',
     &                                   'Index out of bounds',' ')
                     End If
#endif
                     kFrom = kFrom + 1
                     Vec(iRS,iVec) = Work(kFrom)
                  End Do
               End Do
            End If

         End Do

      End Do
      Call GetMem('CR_RI_Scr','Free','Real',ipScr,lScr)

      Call GetMem('CR_RI_F2RS','Free','Inte',ipiF2RS,liF2RS)

      End
