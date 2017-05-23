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
* Copyright (C) 1994, Roland Lindh                                     *
************************************************************************
      Subroutine Cnthlf(Coeff1,nCntr1,nPrm1,Coeff2,nCntr2,nPrm2,
     &                  lZeta,nVec,First,IncVec,A1,A2,A3,
     &                  Indij)
************************************************************************
*                                                                      *
* Object: to do a half transformation. The loop over the two matrix-   *
*         matrix multiplications is segmented such that the end of the *
*         intermediate matrix will not push the start of the same out  *
*         from the cache.                                              *
*                                                                      *
* Called from: Cntrct                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
* Author:     Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8 Coeff1(nPrm1,nCntr1), Coeff2(nPrm2,nCntr2),
     &       A1(lZeta,nVec), A2(IncVec,nprm2),
     &       A3(nVec,nCntr1,nCntr2)
      Parameter (mxnprm=1000)   ! be aware of aCD(fat) basis sets.
      Integer Indij(lZeta),idone(mxnprm),nnz2(mxnprm),
     &        ifirst(mxnprm),last(mxnprm)
      Logical First
*     Call QEnter('CntHlf')
*
      If (nPrm1.gt.mxnprm .or.
     &    nPrm2.gt.mxnprm) Then
          Call WarningMessage(2,'CntHlf: nPrm.gt.mxnprm')
          Call Abend()
      End If
*
*     Sparsity check
*
      iz2=0
      nz2=0
      minva=max(4,(nPrm2+1)/2)
      Do iCntr2 = 1,nCntr2
         nnz2(icntr2)=0
         ifirst(icntr2)=nPrm2+1
         last(icntr2)=0
         Do iPrm2 = 1,nPrm2
            If (Coeff2(iPrm2,iCntr2).ne.Zero) Then
               ifirst(icntr2)=min(ifirst(icntr2),iprm2)
               last(icntr2)=max(last(icntr2),iprm2)
               nnz2(icntr2)=nnz2(icntr2)+1
             End If
         End Do
         If (nnz2(icntr2).ge.minva.and.nz2.eq.iCntr2-1) nz2=iCntr2
      End Do
*
*-----Loop sectioning
*
      Do iiVec = 1, nVec, IncVec
         mVec = Min(IncVec,nVec-iiVec+1)
*                                                                      *
************************************************************************
*                                                                      *
*--------First quarter transformation
*
         Do iCntr1 = 1, nCntr1
            Do iprm2=1,nprm2
               idone(iprm2)=0
            End Do
*
            Do iZeta = 1, lZeta
               iPrm2 = (Indij(iZeta)-1)/nPrm1 + 1
               iPrm1 = Indij(iZeta) - (iPrm2-1)*nPrm1
               C1=Coeff1(iPrm1,iCntr1)
               If (Abs(C1).gt.Zero) Then
                  If (idone(iprm2).gt.0) Then
                     Call DaXpY_(mVec,C1,A1(iZeta,iiVec),lZeta,
     &                                  A2(1,iPrm2),1)
                  Else
                     Call DYaX(mVec,C1,A1(iZeta,iiVec),lZeta,
     &                                 A2(1,iPrm2),1)
                     idone(iprm2)=1
                  End If
               End If
            End Do ! iZeta
*                                                                      *
************************************************************************
*                                                                      *
*-----------Second quarter transformation
*
            Do iprm2=1,nprm2
               If (idone(iprm2).eq.0) Call FZero(a2(1,iprm2),mvec)
            End Do
*
            ic1=1
            If (nz2.gt.1) Then
               If (first) Then
                  Call DGEMM_('n','n',mVec,nz2,nprm2,
     &                       1.0d0,A2,IncVec,Coeff2,nprm2,
     &                       0.0d0,A3(iivec,iCntr1,1),nvec*ncntr1)
               Else
C                  Call mxmb(A2,1,IncVec, Coeff2,1,nprm2,
C     &                      A3(iivec,iCntr1,1),1,nvec*ncntr1,
C     &                       mVec,nPrm2,nz2)
                  Call DGEMM_('N','N',mVec,nz2,nPrm2,
     &                         1.0d0,A2,IncVec,Coeff2,nprm2,
     &                         1.0d0,A3(iivec,iCntr1,1),nvec*ncntr1)
               End If
               ic1=nz2+1
            End If
*
            Do iCntr2=ic1,nCntr2
               If (first) Then
                  If (nnz2(icntr2).ge.minva) Then
C                     Call mxva(A2,1,IncVec, Coeff2(1,icntr2),1,
C     &                         A3(iivec,iCntr1,icntr2),1,mVec,nPrm2)
                      Call dGeMV_('N',mVec,nPrm2,1.d0,A2,IncVec,
     &                            Coeff2(1,icntr2),1,0.d0,
     &                            A3(iivec,iCntr1,icntr2),1)
                  Else
                     iprm2=ifirst(icntr2)
                     c2=coeff2(iprm2,icntr2)
                     Call DYaX(mVec,C2,A2(1,iPrm2),1,
     &                                 A3(iiVec,iCntr1,iCntr2),1)
*
                     iPrm2=ifirst(icntr2)+1
                     mPrm2=last(iCntr2)-iPrm2+1
                     If (mPrm2.gt.0)
     &                  Call DNaXpY(mPrm2,mVec,Coeff2(iPrm2,iCntr2),1,
     &                              A2(1,iPrm2),1,IncVec,
     &                              A3(iiVec,iCntr1,iCntr2),1,0)
*
                  End If
               Else
                  If (nnz2(icntr2).ge.minva) Then
C                     Call mxvb(A2,1,IncVec, Coeff2(1,icntr2),1,
C     &                         A3(iivec,iCntr1,icntr2),1,mVec,nPrm2)
                      Call dGeMV_('N',mVec,nPrm2,1.d0,A2,IncVec,
     &                            Coeff2(1,icntr2),1,1.d0,
     &                            A3(iivec,iCntr1,icntr2),1)
                  Else
                     iPrm2=ifirst(icntr2)
                     mPrm2=last(icntr2)-iPrm2+1
                     Call DNaXpY(mPrm2,mVec,Coeff2(iPrm2,iCntr2),1,
     &                           A2(1,iPrm2),1,IncVec,
     &                           A3(iiVec,iCntr1,iCntr2),1,0)
                  End If
               End If
            End Do ! iCntr2
         End Do    ! iCntr1
*
*-----End of loop sectioning
*
      End Do    ! iiVec
*
*     Call QExit('CntHlf')
      Return
      End
