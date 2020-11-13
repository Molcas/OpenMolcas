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
      Subroutine ValANM(nAtom,nInter,nIter,Bmx,Degen,rInt,Cx,Label,
     &                  nWndw)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 BMx(3*nAtom,3*nAtom), rInt(nInter,nIter),
     &     Degen(3*nAtom), Cx(3*nAtom,nIter)
      Character Label*(*)
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*
*
*-----Values:    q=BuX
*-----Gradients: g=Bu(dE/dX)
*
*                                                                      *
************************************************************************
*                                                                      *
      iSt = nIter
      iEnd = iSt - Min(nIter,nWndw+1) + 1
      If (Label.eq.'Values') Then
*                                                                      *
************************************************************************
*                                                                      *
*-----Values:    q=BX
*
         Call GetMem('ScrC','Allo','Real',ipB,3*nAtom*(iSt-iEnd+1))
*
         Do iIter = iSt, iEnd, -1
            Do j = 1, 3*nAtom
               ij = (iIter-iEnd)*3*nAtom + j-1 + ipB
               Work(ij)=Cx(j,iIter)*Degen(j)
            End Do
         End Do
*
         Call DGEMM_('T','N',
     &               nInter,iSt-iEnd+1,3*nAtom,
     &               1.0d0,BMx,3*nAtom,
     &                     Work(ipB),3*nAtom,
     &               0.0d0,rInt(1,iEnd),nInter)
*
         Call GetMem('ScrC','Free','Real',ipB,3*nAtom*(iSt-iEnd+1))
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
         M = 3*nAtom
         N = nInter
*        NRHS=nIter
         NRHS=iSt-iEnd+1
         Call Eq_Solver('N',M,N,NRHS,BMx,.FALSE.,Degen,
     &                  Cx(1,iEnd),rInt(1,iEnd))
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Write (6,'(A)')' In ValANM: New '
      Call RecPrt(Label,' ',rInt(1,iEnd),nInter,iSt-iEnd+1)
#endif
*
      Return
      End
