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
      SubRoutine PickMO(COUT,nOut,
     &                  nAcO,
     &                  ishell,icmp,iBasi,iBasn,jBasj,jBasn,
     &                  kBask,kBasn,lBasl,lBasn,iaoii)
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "etwas.fh"
#include "pso.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 COUT(nOut)
      Integer iCmp(4),iBas(4),nBs(4)
      Integer iTwoj(0:7),ishell(4),iaoii(4)
      Data iTwoj/1,2,4,8,16,32,64,128/
*
      iBas(1)=iBasi
      iBas(2)=jBasj
      iBas(3)=kBask
      iBas(4)=lBasl
      nBs(1)=iBasn
      nBs(2)=jBasn
      nBs(3)=kBasn
      nBs(4)=lBasn
      ip2=1
      nA=0
      Do iIrrep=0,nIrrep-1
          nA=nA+nAsh(iIrrep)
      End Do
      Do iCnt=3,4
         ipC=ipCMO-1
         Do iIrrep=0,nIrrep-1
            iOrb=nIsh(iIrrep)
            Do iAsh=1,nAsh(iIrrep)
               jj=iCmp(iCnt)
               Do i1=1,jj
                  iSO=iAOtSO(iAOii(iCnt)+i1,iIrrep)+iBas(iCnt)-1
                  ip1=ipC+(iOrb+iAsh-1)*nBas(iIrrep)+iso
                  If (iAnd(IrrCmp(IndS(iShell(iCnt))+i1),
     &                iTwoj(iIrrep)).ne.0) Then
                     call dcopy_(nBs(iCnt),Work(ip1),1,COUT(ip2),1)
                  Else
                     call dcopy_(nBs(iCnt),[0.0d0],0,COUT(ip2),1)
                  End If
                  ip2=ip2+nBs(iCnt)
               End Do
            End Do
            ipc=ipc+nBas(iIrrep)**2
         End Do
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nAcO)
      End
