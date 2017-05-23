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
      Subroutine BMtrx_User_Defined(
     &                 nLines,nBVec,ipBMx,nAtom,nInter,
     &                 ip_rInt,Lbl,Coor,nDim,rMass,
     &                 Name,nSym,iOper,Smmtrc,
     &                 Degen,BSet,HSet,nIter,ip_drInt,
     &                 Gx,Cx,mTtAtm,iAnr,
     &                 nStab,jStab,Numerical,
     &                 HWRS,Analytic_Hessian,
     &                 iOptC,PrQ,mxdc,iCoSet,lOld,
     &                 nFix,mTR,ip_KtB,nQQ,Redundant,nqInt,MaxItr)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Coor(3,nAtom), rMass(nAtom), Degen(3*nAtom),
     &       Gx(3*nAtom,nIter), Cx(3*nAtom,nIter)
      Character Lbl(nInter)*8, Name(nAtom)*(LENIN)
      Integer   iOper(0:nSym-1), iAnr(nAtom),
     &          nStab(nAtom), jStab(0:7,nAtom), iCoSet(0:7,nAtom)
      Logical Smmtrc(3*nAtom), BSet, HSet, Redundant,
     &        Numerical, HWRS, Analytic_Hessian, PrQ, lOld
*                                                                      *
************************************************************************
*                                                                      *
      iRout=133
      iPrint=nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
*.... Section for user defined internal coordinates
*
      Call Rd_UDIC(nLines,iInt,nFix,nRowH)
      nQQ=iInt+nFix
*
      If (ip_rInt.eq.ip_Dummy) Then
         nqInt=nQQ*MaxItr
         Call GetMem(' qInt','Allo','Real',ip_rInt, nqInt)
         Call GetMem('dqInt','Allo','Real',ip_drInt,nqInt)
      End If
      Call Allocate_Work(ipBmx,3*nAtom*nQQ)
      Call FZero(Work(ipBMx),3*nAtom*nQQ)
      Call GetMem('rMult','Allo','Real',ipMult,nBVec)
      Call GetMem('BVec','Allo','Real',ipBVec,nBVec*3*nAtom)
      Call FZero(Work(ipBVec),nBVec*3*nAtom)
      Call GetMem('Val','Allo','Real',ipVal,nBVec)
      Call FZero(Work(ipVal),nBVec)
      Call GetMem('Lab','Allo','Char',ipLab,nBVec*8)
*
*-----Compute the B matrix in symmetry distinct basis and the
*     internal coordinates.
*
      ip = ip_rInt + (nIter-1)*nQQ
      Call DefInt(Work(ipBVec),nBVec,cWork(ipLab),Work(ipBMx),nQQ,
     &            nAtom,nLines,Work(ipVal),Work(ip),Lbl,Name,
     &            Coor,rMass,nSym,iOper,jStab,nStab,mxdc,Work(ipMult),
     &            nDim-mTR,Redundant)
*
      Call GetMem('Lab','Free','Char',ipLab,nBVec*8)
      Call GetMem('Val','Free','Real',ipVal,nBVec)
      Call GetMem('BVec','Free','Real',ipBVec,nBVec*3*nAtom)
      Call GetMem('rMult','Free','Real',ipMult,nBVec)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the gradient
*
      If (BSet) Call Force(nFix,Gx(1,nIter),nAtom,nQQ,Work(ipBMx),
     &                     Name,nIter,Work(ip_drInt),Lbl,Degen)
*                                                                      *
************************************************************************
*                                                                      *
      If (HSet.and..NOT.lOld.and.BSet) Then
         Call Allocate_Work(ip_KtB,nDim*nQQ)
*
         Call Allocate_Work(ipDegen,nDim)
         i=0
         Do ix = 1, 3*nAtom
            If (Smmtrc(ix)) Then
               Work(ipDegen+i) = Degen(ix)
               i = i + 1
            End If
         End Do
*
         Do j = 1, nQQ
            i = 0
            Do ix = 1, 3*nAtom
               If (Smmtrc(ix)) Then
                  i = i + 1
                  ij = (j-1)*nDim + i - 1 + ip_KtB
                  ixj= (j-1)*3*nAtom + ix - 1 + ipBmx
                  Work(ij) = Work(ixj)
               End If
            End Do
         End Do
*
         Do iInter = 1, nQQ
            Do iDim = 1, nDim
               ij = (iInter-1)*nDim + iDim - 1 + ip_KtB
*              Work(ij) = Work(ij) / Sqrt(Work(ipDegen+iDim-1))
               Work(ij) = Work(ij) / Work(ipDegen+iDim-1)
            End Do
         End Do
         Call Free_Work(ipDegen)
      Else
         ip_KtB = ip_Dummy
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Cx)
         Call Unused_integer(mTtAtm)
         Call Unused_integer_array(iAnr)
         Call Unused_logical(Numerical)
         Call Unused_logical(HWRS)
         Call Unused_logical(Analytic_Hessian)
         Call Unused_integer(iOptC)
         Call Unused_logical(PrQ)
         Call Unused_integer_array(iCoSet)
      End If
      End
