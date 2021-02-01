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
      SubRoutine PickMO(COUT,nOut,nAcO,icmp,iBasi,iBasn,jBasj,jBasn,
     &                  kBask,kBasn,lBasl,lBasn,iaoii)
      use Basis_Info, only: nBas
      use SOAO_Info, only: iAOtSO
      use pso_stuff
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (a-h,o-z)
#include "etwas.fh"
#include "real.fh"
      Real*8 COUT(nOut)
      Integer iCmp(4),iBas(4),nBs(4)
      Integer iAOii(4)
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

      Do iCnt=3,4
         ipC=0
         Do iIrrep=0,nIrrep-1
            iOrb=nIsh(iIrrep)
            Do iAsh=1,nAsh(iIrrep)
               jj=iCmp(iCnt)
               Do i1=1,jj
                  iSO=iAOtSO(iAOii(iCnt)+i1,iIrrep)+iBas(iCnt)-1
                  If (iSO>0) Then
                     ip1=ipC+(iOrb+iAsh-1)*nBas(iIrrep)+iSO
                     call dcopy_(nBs(iCnt),CMO(ip1,1),1,COUT(ip2),1)
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
