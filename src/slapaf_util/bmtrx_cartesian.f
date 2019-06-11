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
      Subroutine BMtrx_Cartesian(
     &                 nLines,ipBMx,nAtom,nInter,ip_rInt,Coor,nDim,
     &                 dMass,Name,nSym,iOper,Smmtrc,Degen,BSet,HSet,
     &                 nIter,ip_drInt,Gx,Cx,mTtAtm,iAnr,nStab,jStab,
     &                 Numerical,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &                 iCoSet,lOld,nFix,mTR,TRVec,ipEVal,ip_Hss_x,
     &                 ip_KtB,nQQ,Redundant,nqInt,MaxItr,nWndw)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Coor(3,nAtom), dMass(nAtom), Degen(3*nAtom),
     &       Gx(3*nAtom,nIter), Cx(3*nAtom,nIter), TRVec(nDim,mTR)
      Character Name(nAtom)*(LENIN)
      Integer   iOper(0:nSym-1), iAnr(nAtom),
     &          nStab(nAtom), jStab(0:7,nAtom), iCoSet(0:7,nAtom)
      Logical Smmtrc(3*nAtom), BSet, HSet, Redundant,
     &        Numerical, HWRS, Analytic_Hessian, PrQ, lOld
*                                                                      *
************************************************************************
*                                                                      *
      iRout=133
      iPrint=nPrint(iRout)
*
*-----Recompute the B matrix once each macroIteration, this is
*     not done if a numerical Hessian is computed.
*                                                                      *
************************************************************************
*                                                                      *
*     R E D U N D A N T  C A R T E S I A N  C O O R D S
*
      If (Redundant) Then
         nQQ = nDim
         If (ip_rInt.eq.ip_Dummy) Then
            nqInt=nQQ*MaxItr
            Call GetMem(' qInt','Allo','Real',ip_rInt, nqInt)
            Call GetMem('dqInt','Allo','Real',ip_drInt,nqInt)
         End If
         Call Allocate_Work(ipEVec,nDim**2)
         Call FZero(Work(ipEVec),nDim**2)
         call dcopy_(nDim,[One],0,Work(ipEVec),nDim+1)
*                                                                      *
************************************************************************
*                                                                      *
*------- Move over the eigenvectors putting to BMx
*
         Call Allocate_Work(ipBMx,(3*nAtom)*nQQ)
         Call FZero(Work(ipBMx), (3*nAtom)*nQQ)
         ipFrom = ipEVec
         Call BPut(Work(ipFrom),nDim,Work(ipBMx),3*nAtom,Smmtrc,
     &             nQQ,Degen)
         If (iPrint.ge.19) Call RecPrt('In Bmtrx: B',' ',Work(ipBMx),
     &                                 3*nAtom,nQQ)
*                                                                      *
************************************************************************
*                                                                      *
         If (PrQ.and.nAtom.le.5)
     &      Call List2('Cartesian Redundant',
     &                 Name,Work(ipBMx),nAtom,nQQ,Smmtrc)
*                                                                      *
************************************************************************
*                                                                      *
*------- Project the model Hessian with respect to rotations and
*        translations. The eigenvalues are shifted to large positive
*        eigenvalues to effectively remove any displacements in the
*        rotational and translations directions and to make sure that
*        the matrix is not singular.
*
         Call Allocate_iWork(ipInd,nDim)
         iInd=0
         Do i = 1, 3*nAtom
            If (Smmtrc(i)) Then
               iWork(ipInd+iInd)=i
               iInd=iInd+1
            End If
         End Do
*
*        Compute H|i>
*
         Call Allocate_Work(ipHi,mTR*nDim)
         Call Allocate_Work(ipiHi,mTR)
         Call FZero(Work(ipHi),mTR*nDim)
*
         Do j = 1, mTR
            Do i = 1, nDim
               Temp = 0.0D0
               Do k = 1, nDim
                  kx = iWork(ipInd+k-1)
                  ik=(i-1)*nDim+k -1 + ip_Hss_x
                  Temp = Temp
     &                 + Work(ik) * Sqrt(Degen(kx))
     &                 * TRVec(k,j)
               End Do
               Work(ipHi+(j-1)*nDim+i-1) = Temp
            End Do
         End Do
*        Call RecPrt('Hi',' ',Work(ipHi),nDim,mTR)
         Do iTR = 1, mTR
            Work(ipiHi+iTR-1) = DDot_(nDim,TRVec(1,iTR),1,
     &                                    Work(ipHi+(iTR-1)*nDim),1)
         End Do
*        Call RecPrt('iHi',' ',Work(ipiHi),mTR,1)
*
         Do i = 1, nDim
            ix = iWork(ipInd+i-1)
            Do j = 1, i
               jx = iWork(ipInd+j-1)
               ij = (j-1)*nDim + i -1 + ip_Hss_x
               ji = (i-1)*nDim + j -1 + ip_Hss_x
               Temp = Half*(Work(ij)+Work(ji))
*define UNIT_MM
#ifndef UNIT_MM
*
*              Here we shift the eigenvectors corresponding to
*              translations and rotations to a large positive values.
*
               Do iTR = 1, mTR
                  Omega = 1.0D+5
                  Hii   = Work(ipiHi+iTR-1)
                  Temp = Temp
     &                 + Sqrt(Degen(ix)) * (
     &                 - TRVec(i,iTR) * Work(ipHi+(iTR-1)*nDim+j-1)
     &                 - Work(ipHi+(iTR-1)*nDim+i-1) * TRVec(j,iTR)
     &                 + TRVec(i,iTR) * (Omega+Hii) * TRVec(j,iTR)
     &                                     )* Sqrt(Degen(jx))
               End Do
#endif
*
               Work(ij)=Temp
               Work(ji)=Temp
            End Do
         End Do
         Call Free_Work(ipiHi)
         Call Free_Work(ipHi)
         Call Free_iWork(ipInd)
*
*        Clean up the gradient wrt translational and rotational
*        component.
*
*        |g> = |g> - Sum(TR) |i><i|g>
*
         If (BSet) Then
*
*           Call RecPrt('Gx',' ',Gx(1,nIter),1,3*nAtom)
*           Call RecPrt('TRVec',' ',TRVec,nDim,mTR)
*
            Do iTR = 1, mTR
*
*              <i|g>
*
               Temp=0.0D0
               iInd=0
               Do i = 1, 3*nAtom
                  If (Smmtrc(i)) Then
                     iInd=iInd+1
                     Temp = Temp + Degen(i)*Gx(i,nIter)*TRVec(iInd,iTR)
                  End If
               End Do
*
               iInd=0
               Do i = 1, 3*nAtom
                  If (Smmtrc(i)) Then
                     iInd=iInd+1
                     Gx(i,nIter) = Gx(i,nIter) - TRVec(iInd,iTR)*Temp
                  End If
               End Do
*
            End Do
*
*           Call RecPrt('Gx',' ',Gx(1,nIter),1,3*nAtom)
*
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
*     N O N - R E D U N D A N T  C A R T E S I A N  C O O R D S
*
         nQQ=nInter
         If (ip_rInt.eq.ip_Dummy) Then
            nqInt=nQQ*MaxItr
            Call GetMem(' qInt','Allo','Real',ip_rInt, nqInt)
            Call GetMem('dqInt','Allo','Real',ip_drInt,nqInt)
         End If
*
*------- Project the model Hessian with respect to rotations and
*        translations. The eigenvalues are shifted to negative
*        eigenvalues.
*
         Call Allocate_iWork(ipInd,nDim)
         iInd=0
         Do i = 1, 3*nAtom
            If (Smmtrc(i)) Then
               iWork(ipInd+iInd)=i
               iInd=iInd+1
            End If
         End Do
*
*        Compute H|i>
*
         Call Allocate_Work(ipHi,mTR*nDim)
         Call Allocate_Work(ipiHi,mTR)
         Call FZero(Work(ipHi),mTR*nDim)
*
         Do j = 1, mTR
            Do i = 1, nDim
               Temp = 0.0D0
               Do k = 1, nDim
                  kx = iWork(ipInd+k-1)
                  ik=(i-1)*nDim+k -1 + ip_Hss_x
                  Temp = Temp
     &                 + Work(ik) * Sqrt(Degen(kx))
     &                 * TRVec(k,j)
               End Do
               Work(ipHi+(j-1)*nDim+i-1) = Temp
            End Do
         End Do
*        Call RecPrt('Hi',' ',Work(ipHi),nDim,mTR)
         Do iTR = 1, mTR
            Work(ipiHi+iTR-1) = DDot_(nDim,TRVec(1,iTR),1,
     &                                    Work(ipHi+(iTR-1)*nDim),1)
         End Do
*        Call RecPrt('iHi',' ',Work(ipiHi),mTR,1)
*
         Do i = 1, nDim
            ix = iWork(ipInd+i-1)
            Do j = 1, i
               jx = iWork(ipInd+j-1)
               ijTri=i*(i-1)/2 + j -1 + ipEVal
               ij = (j-1)*nDim + i -1 + ip_Hss_x
               ji = (i-1)*nDim + j -1 + ip_Hss_x
               Work(ijTri) = Half*(Work(ij)+Work(ji))
*
*              Here we shift the eigenvectors corresponding to tran-
*              lations and rotations down to negative faked eigen-
*              values.
*
               Do iTR = 1, mTR
                  Omega = -DBLE(iTR)
                  Hii   = Work(ipiHi+iTR-1)
                  Work(ijTri) = Work(ijTri)
     &                 + Sqrt(Degen(ix)) * (
     &                 - TRVec(i,iTR) * Work(ipHi+(iTR-1)*nDim+j-1)
     &                 - Work(ipHi+(iTR-1)*nDim+i-1) * TRVec(j,iTR)
     &                 + TRVec(i,iTR) * (Omega+Hii) * TRVec(j,iTR)
     &                                     )* Sqrt(Degen(jx))
               End Do
*
               Work(ij)=Work(ijTri)
               Work(ji)=Work(ijTri)
            End Do
         End Do
         Call Free_Work(ipiHi)
         Call Free_Work(ipHi)
         Call Free_iWork(ipInd)
#ifdef _DEBUG_
         Call TriPrt(' The Projected Model Hessian','(5G20.10)',
     &               Work(ipEVal),nDim)
         Call RecPrt(' The Projected Model Hessian','(5G20.10)',
     &               Work(ip_Hss_x),nDim,nDim)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*        Compute the eigen vectors for the Cartesian Hessian
*
         Call GetMem('EVec','Allo','Real',ipEVec,(3*mTtAtm)**2)
         Call Hess_Vec(mTtAtm,Work(ipEVal),Work(ipEVec),nAtom,nDim)
*                                                                      *
************************************************************************
*                                                                      *
*------- Move over the eigenvectors putting to BMx
*
         Call Allocate_Work(ipBMx,(3*nAtom)**2)
         Call FZero(Work(ipBMx), (3*nAtom)**2)
         ipFrom = ipEVec + mTR*nDim
         Call BPut(Work(ipFrom),nDim,Work(ipBMx),3*nAtom,Smmtrc,
     &             nQQ,Degen)
         If (iPrint.ge.19) Call RecPrt('In Bmtrx: B',' ',Work(ipBMx),
     &                                 3*nAtom,nQQ)
*                                                                      *
************************************************************************
*                                                                      *
         If (PrQ.and.nAtom.le.5)
     &      Call List2('Cartesian Approximate Normal Modes',
     &                  Name,Work(ipBMx),nAtom,nQQ,Smmtrc)
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (HSet.and..NOT.lOld) Then
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
         call dcopy_(nDim*nQQ,Work(ipFrom),1,Work(ip_KtB),1)
         Do iInter = 1, nQQ
            Do iDim = 1, nDim
               ij = (iInter-1)*nDim + iDim - 1 + ip_KtB
               Work(ij) = Work(ij) / Sqrt(Work(ipDegen+iDim-1))
            End Do
         End Do
         Call Free_Work(ipDegen)
      Else
         ip_KtB = ip_Dummy
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('EVec','Free','Real',ipEVec,(3*mTtAtm)**2)
*                                                                      *
************************************************************************
*                                                                      *
*
*---- Compute the value and gradient vectors in the new basis.
*
      Call ValANM(nAtom,nQQ,nIter,Work(ipBmx),Degen,
     &            Work(ip_rInt),Cx,'Values',iPrint,nWndw)
      If (BSet) Call ValANM(nAtom,nQQ,nIter,Work(ipBMx),Degen,
     &                      Work(ip_drInt),Gx,'Gradients',iPrint,nWndw)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nLines)
         Call Unused_real_array(Coor)
         Call Unused_real_array(dMass)
         Call Unused_integer_array(iOper)
         Call Unused_integer_array(iAnr)
         Call Unused_integer_array(nStab)
         Call Unused_integer_array(jStab)
         Call Unused_logical(Numerical)
         Call Unused_logical(HWRS)
         Call Unused_logical(Analytic_Hessian)
         Call Unused_integer(iOptC)
         Call Unused_integer(mxdc)
         Call Unused_integer_array(iCoSet)
         Call Unused_integer(nFix)
      End If
      End
