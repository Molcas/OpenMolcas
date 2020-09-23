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
      SubRoutine PckMO2(COUT,nAcO,icmpi,iBasi,jcmpj,jBasj,iAOi,jAOj)
      use Basis_Info, only: nBas
      use SOAO_Info, only: iAOtSO
      use pso_stuff
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "etwas.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 COUT(*)
      Integer iCmp(4),nBs(4)
      Integer iTwoj(0:7),iaoii(4)
      Data iTwoj/1,2,4,8,16,32,64,128/
*
      nBs(1)=iBasi
      nBs(2)=jBasj
      iAOii(1)=iAOi
      iAOii(2)=jAOj
      icmp(1)=icmpi
      icmp(2)=jcmpj
      ip2=1
      nA=0
      Do iIrrep=0,nIrrep-1
         nA=nA+nAsh(iIrrep)
      End Do
      Do iCnt=1,2
         ipC=0
         Do iIrrep=0,nIrrep-1
            iOrb=nIsh(iIrrep)
            Do iAsh=1,nAsh(iIrrep)
               jj=iCmp(iCnt)
               Do i1=1,jj
                  iSO=iAOtSO(iAOii(iCnt)+i1,iIrrep)
                  If (iSO>0) Then
                     ip1=ipC+(iOrb+iAsh-1)*nBas(iIrrep)+iso
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
