************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2013, Roland Lindh                                     *
*               2015, Ignacio Fdez. Galvan (split from gencxctl)       *
************************************************************************
      Subroutine get_drdq(drdq,mLambda)
      Implicit None
************************************************************************
*     subroutine to get the dr/dq vectors for the constraints as given *
*     in the 'UDC' file.                                               *
************************************************************************
#include "info_slapaf.fh"
#include "WrkSpc.fh"
#include "real.fh"
      Real*8, Intent(InOut) :: drdq(mInt,nLambda)
      Integer, Intent(Out) :: mLambda
*
      Integer n3,nBV,i,iLambda,iOff,iOff2
      Integer ipBMx,ipdBMx,ipBVc,ipdBVc,ipcInt,ipcInt0,ipValue,ipValue0,
     &        ipMult,ip_iFlip
      Real*8 RR
      Real*8, External :: DDot_
*
      Character*8 Lbl(mInt), Labels(iRow_c-1)
*
      n3=3*nsAtom
*
      If (nLambda.ne.0) Then
         nBV=iRow_c-nLambda-1
         Call GetMem('BVec','Allo','Real',ipBVc,n3*nBV)
         Call GetMem('dBVec','Allo','Real',ipdBVc,nBV*n3**2)
         Call GetMem('BMx','Allo','Real',ipBMx,n3*nLambda)
         Call GetMem('Value','Allo','Real',ipValue,nBV)
         Call GetMem('Value0','Allo','Real',ipValue0,nBV)
         Call FZero(Work(ipValue0),nBV) ! dummy initiate
         Call GetMem('cInt','Allo','Real',ipcInt,nLambda)
         Call GetMem('cInt0','Allo','Real',ipcInt0,nLambda)
         Call GetMem('Mult','Allo','Real',ipMult,nBV**2)
         Call GetMem('dBMx','Allo','Real',ipdBMx,nLambda*n3**2)
         Call GetMem('iFlip','Allo','Inte',ip_iFlip,nBV)
*
         Call DefInt2(Work(ipBVc),Work(ipdBVc),nBV,Labels,
     &                Work(ipBMx),nLambda,nsAtom,iRow_c,
     &                Work(ipValue),Work(ipcInt),Work(ipcInt0),
     &                Lbl,AtomLbl,
     &                lWrite,jStab,nStab,mxdc,
     &                Work(ipMult),Work(ipdBMx),
     &                Work(ipValue0),Iter,iWork(ip_iFlip),
     &                Work(ipCM))
*        Call RecPrt('dr/dx',' ',Work(ipBMx),n3,nLambda)
*
*        Assemble dr/dq: Solve  B dr/dq = dr/dx
*
         Call FZero(drdq,nLambda*mInt)
*
*        Temporary fix of the dC/dx vector which always
*        is propted up with the full degeneracy factor.
*
         If (.not.Curvilinear) Then
            iOff=ipBMx
            Do iLambda=1,nLambda
               Do i=1,n3
                  Work(iOff+i-1)=Work(iOff+i-1)/Degen(i)
               End Do
               iOff=iOff+n3
            End Do
         End If
         Call Eq_Solver('N',n3,mInt,nLambda,Work(ipB),Curvilinear,Degen,
     &                  Work(ipBMx),drdq)
*        Call RecPrt('drdq',' ',drdq,mInt,nLambda)
*
         Call Free_iWork(ip_iFlip)
         Call Free_Work(ipdBMx)
         Call Free_Work(ipMult)
         Call Free_Work(ipcInt0)
         Call Free_Work(ipcInt)
         Call Free_Work(ipValue0)
         Call Free_Work(ipValue)
         Call Free_Work(ipBMx)
         Call Free_Work(ipdBVc)
         Call Free_Work(ipBVc)
      End If
*
*     Double check that we don't have any null vectors
*
      iOff=1
      iOff2=1
      mLambda=nLambda
      Do iLambda=1,nLambda
         RR=Sqrt(DDot_(mInt,drdq(1,iOff),1,drdq(1,iOff),1))
         If (RR.lt.1.0D-12) Then
            Write (6,*) 'Warning: constraint ',iLambda,
     &                  ' has a null vector, I''ll remove it!'
            mLambda=mLambda-1
         Else
            If (iOff.ne.iOff2)
     &         Call dCopy_(mInt,drdq(1,iOff),1,drdq(1,iOff2),1)
            iOff2=iOff2+1
         End If
         iOff=iOff+1
      End Do
      If (mLambda.lt.nLambda)
     &   Call FZero(drdq(1,mLambda+1),mInt*(nLambda-mLambda))
*
      End Subroutine get_drdq
