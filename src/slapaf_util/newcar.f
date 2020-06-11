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
      Subroutine NewCar(kIter,nBVec,nLines,nAtom,nDim,nInter,
     &                  Coor,ipBMx,dMass,Lbl,Shift,ip_qInt,ip_dqInt,
     &                  DFC,dss,Tmp,Name,iOper,nSym,iSym,Smmtrc,
     &                  Degen,Gx,Cx,mTtAtm,iANr,iOptH,User_Def,nStab,
     &                  jStab,Curvilinear,Numerical,DDV_Schlegel,HWRS,
     &                  Analytic_Hessian,iOptC,PrQ,mxdc,iCoSet,rHidden,
     &                  ipRef,Redundant,nqInt,MaxItr,iRef)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "Molcas.fh"
      Real*8 Coor(3,nAtom), dMass(nAtom),Shift(nInter,kIter),
     &       DFC(3*nAtom), dss(nInter), Tmp(nInter), Degen(3*nAtom),
     &       Gx(3*nAtom,kIter), Cx(3*nAtom,kIter+1)
      Integer   iOper(0:7), iSym(3), iANr(nAtom),
     &          nStab(nAtom), jStab(0:7,nAtom), iCoSet(0:7,nAtom)
      Character Lbl(nInter)*8, Name(nAtom)*(LENIN)
      Logical Smmtrc(3,nAtom), User_Def, Redundant,
     &        Curvilinear, Numerical, DDV_Schlegel, HWRS,
     &        Analytic_Hessian, PrQ
*                                                                      *
************************************************************************
*                                                                      *
      Call QEnter('NewCar')
      iRout=136
      iPrint=nPrint(iRout)
      If (iPrint.ge.99) Then
         Call RecPrt('NewCar: q',' ',Work(ip_qInt),nInter,kIter+1)
         Call RecPrt('NewCar: g',' ',Work(ip_dqInt),nInter,kIter)
         Call RecPrt('NewCar: Shift',' ',Shift,nInter,kIter)
      End If
*
*-----Form the new set of symmetry distinct cartesian coordinates.
*
      rMax = Zero
      ThrR = 0.1d0
      Do i = 1, nInter
        If (Abs(Shift(i,kIter)).gt.rMax) rMax=Abs(Shift(i,kIter))
      End Do
*
      ip = ip_qInt + (kIter-1)*nInter
      call dcopy_(nInter,Work(ip),1,Tmp,1)
*                                                                      *
************************************************************************
*                                                                      *
      Error = 1.0D-12
      call dcopy_(nInter,Shift(1,kIter),1,dss,1)
      Call DaXpY_(nInter,One,dss,1,Tmp,1)
*
      Call Int2Car(dss,Tmp,nInter,ip_qInt,Coor,nAtom,nBVec,ipBMx,dMass,
     &             nLines,DFC,ndim,Lbl,Name,iOper,nSym,iSym,Smmtrc,
     &             Degen,kIter,ip_dqInt,Gx,Cx,mTtAtm,iANr,iOptH,
     &             User_Def,nStab,jStab,Curvilinear,Numerical,
     &             DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &             iCoSet,rHidden,Error,ipRef,Redundant,nqInt,MaxItr,
     &             iRef)
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('NewCar')
      Return
      End
