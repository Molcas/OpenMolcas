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
      Subroutine FormNumHess(nIter,Grdn,Shift,nInter,Delta,Stop,
     &                       qInt,nAtom,Cubic,iNeg,DipM,mTR,Smmtrc,
     &                       Degen,UserT,UserP,nUserPT,nsRot,lTherm,
     &                       lDoubleIso,Curvilinear)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Grdn(nInter,nIter), Shift(nInter,nIter), UserT(64), UserP,
     &       qInt(nInter,nIter+1), DipM(3,nIter), Degen(3,nAtom)
      Logical Stop, Cubic, Smmtrc(3,nAtom), lTherm, lDoubleIso, Found,
     &        Curvilinear
      Integer nUserPT, nsRot
*                                                                      *
************************************************************************
*                                                                      *
      Call QEnter('FormNumHess')
      iRout = 182
      iPrint = nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
      Call Allocate_Work(ipdDipM,3*(nInter+mTR))
      Call FZero(Work(ipdDipM),3*(nInter+mTR))
*                                                                      *
************************************************************************
*                                                                      *
*-----Form the Hessian matrix via a 2-point numerical differentiation.
*
      Stop = .False.
      Call Allocate_Work(ipH,nInter**2)
      If (Cubic) Then
         Call Allocate_Work(ipFEq,nInter**3)
      Else
         ipFEq=ip_Dummy
      End If
      Call NmHess(Shift,nInter,Grdn,nIter,Work(ipH),Delta,qInt,
     &            Work(ipFEq),Cubic,DipM,Work(ipdDipM))
      Write (6,*)
      Write (6,*) ' Numerical differentiation is finished!'
      If (iPrint.GE.98) Call RecPrt(' Numerical force constant matrix',
     &     ' ',Work(ipH),nInter,nInter)
*
      Call Add_Info('Numerical Hessian',Work(ipH),nInter**2,2)
      Call Put_dArray('Hss_Q',Work(ipH),nInter**2)
      Call Put_dArray('Hss_upd',Work(ip_Dummy),0)
*
*-----That is in internal coordinates, now transform it to Cartesians
*     d^2E/dx^2 = dQ/dx d^2E/dQ^2 dQ/dx + d^2Q/dx^2 dE/dQ
*
      Call Qpg_dArray('KtB',Found,nKB)
      If (.Not.Found) Call Abend()
      nDim=nKB/nInter
      Call Allocate_Work(ipKB,nDim*nInter)
      Call Allocate_Work(ipHB,nDim*nInter)
      Call Allocate_Work(ipHx,nDim**2)
      Call Allocate_Work(ipDegen,nDim)
      Call Get_dArray('KtB',Work(ipKB),nKB)
      Call DGeMM_('N','T',nInter,nDim,nInter,One,Work(ipH),nInter,
     &                    Work(ipKB),nDim,Zero,Work(ipHB),nInter)
      Call DGeMM_('T','T',nDim,nDim,nInter,One,Work(ipHB),nInter,
     &                    Work(ipKB),nDim,Zero,Work(ipHx),nDim)
      i=0
      Do ii=1,nAtom
         Do ij=1,3
            If (Smmtrc(ij,ii)) Then
               Work(ipDegen+i) = Degen(ij,ii)
               i=i+1
            End If
         End Do
      End Do
*
      If (Curvilinear) Then
*        Compute and add the d^2Q/dx^2 dE/dQ part
         Call dBuu(Work(ipDegen),nInter,nDim,Grdn(1,1),Work(ipHx),
     &             .True.)
      End If
*
      Call Put_dArray('Hss_X',Work(ipHx),nDim**2)
      Call Free_Work(ipKB)
      Call Free_Work(ipHB)
      Call Free_Work(ipHx)
      Call Free_Work(ipDegen)
      Call Free_Work(ipH)
*
      If (Cubic) Then
         Call RecPrt(' Numerical cubic force constant matrix',' ',
     &               Work(ipFEq),nInter**2,nInter)
         Call Add_Info('Numerical anharm. cons.',Work(ipFEq),
     &                 nInter**3,2)
         Call Free_Work(ipFEq)
      End If
*
*---- Do an harmonic frequency analysis
*
      Call Allocate_Work(ipIRInt,nInter+mTR)
      Call HrmFrq(nAtom,nInter,iNeg,
     &            Work(ipdDipM),mTR,Smmtrc,DipM,Work(ipIRInt),
     &            UserT, UserP, nUserPT, nsRot, lTherm, lDoubleIso)
      Call Add_Info('Numerical IR Intensities',Work(ipIRInt),
     &              nInter,2)
      Call Free_Work(ipIRInt)
      Write (6,*)
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_Work(ipdDipM)
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('FormNumHess')
      Return
      End
