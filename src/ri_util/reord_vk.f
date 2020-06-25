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
      SUBROUTINE Reord_Vk(ip_V_k,nProcs,myProc,nV_k,nV_t,nA,jSym,Array)
      Implicit None
      Integer nProcs, myProc, nV_k(*), nV_t(*), nA(*), jSym
      Integer ip_V_k(nProcs)
      Real*8 Array(*)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"
*
      Integer ipScr, ik, ifr, ito, nAV_t, jOff, kOff, iSym
*
      Integer N2
      Parameter (N2 = InfVec_N2)
      Integer InfVec, i, j, k
*********************************************************************
      InfVec(i,j,k)=iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
*********************************************************************
*
      nAV_t=0
      Do iSym=1,jSym
         nAV_t = nAV_t + nA(iSym)*nV_t(iSym)
      End Do
      Call GetMem('Vk_scr','Allo','Real',ipScr,nAV_t)
      Call FZero(Work(ipScr),nAV_t)
*
      jOff=0
      kOff=0
      Do iSym=1,jSym
         Do ik=1,nV_k(iSym)
*
            ifr = ip_V_k(myProc) + jOff + nA(iSym)*(ik-1)
            ito = ipScr + kOff + nA(iSym)*(InfVec(ik,5,iSym)-1)
            call dcopy_(nA(iSym),Array(ifr),1,Work(ito),1)
*
         End Do
         jOff=jOff+nA(iSym)*nV_k(iSym)
         kOff=kOff+nA(iSym)*nV_t(iSym)
      End Do
*
      call dcopy_(nAV_t,Work(ipScr),1,Array(ip_V_k(1)),1)
      Call GADGOP(Array(ip_V_k(1)),nAV_t,'+')
*
      Call GetMem('Vk_scr','Free','Real',ipScr,nAV_t)
*
      Return
      End
