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
      Subroutine get_drdq(drdq,mInt_,nLambda_,mLambda) !fix the _ later!
      use Slapaf_Info, only: BMx, Degen
      use Slapaf_Parameters, only: iRow_c, Curvilinear
      Implicit None
************************************************************************
*     subroutine to get the dr/dq vectors for the constraints as given *
*     in the 'UDC' file.                                               *
************************************************************************
#include "info_slapaf.fh"
#include "stdalloc.fh"
#include "real.fh"
      Real*8, Intent(InOut) :: drdq(mInt_,nLambda_)
      Integer, Intent(In) :: mInt_, nLambda_
      Integer, Intent(Out) :: mLambda
*
      Integer n3,nBV,i,iLambda,iOff,iOff2, iAtom, ixyz
      Real*8 RR
      Real*8, External :: DDot_
*
      Character(LEN=8), Allocatable:: Lbl(:)
      Real*8, Allocatable:: BVc(:), dBVc(:), BMx_t(:,:), Value(:),
     &                      Value0(:), cInt(:), cInt0(:), Mult(:),
     &                      dBMx(:)
      Integer, Allocatable:: iFlip(:)
*
      n3=SIZE(Degen)
*
      If (nLambda_.ne.0) Then
         nBV=iRow_c-nLambda_-1
         Call mma_allocate(BMx_t,n3,nLambda_,Label='BMx_t')

         Call mma_allocate(BVc,n3*nBV,Label='BVc')
         Call mma_allocate(dBVc,nBV*n3**2,Label='dBVc')
         Call mma_allocate(Value,nBV,Label='Value')
         Call mma_allocate(Value0,nBV,Label='Value0')
         Value0(:)=Zero
         Call mma_allocate(cInt,nLambda_,Label='cInt')
         Call mma_allocate(cInt0,nLambda_,Label='cInt0')
         Call mma_allocate(Mult,nBV**2,Label='Mult')
         Call mma_allocate(dBMx,nLambda_*n3**2,Label='dBMx')
         Call mma_allocate(iFlip,nBV,Label='iFlip')
         Call mma_allocate(Lbl,mInt_,Label='Lbl')
*
         Call DefInt2(BVc,dBVc,nBV,BMx_t,nLambda_,
     &                SIZE(Degen,2),iRow_c,
     &                Value,cInt,cInt0,Lbl,lWrite,
     &                Mult,dBMx,Value0,Iter,iFlip)

         Call mma_deallocate(Lbl)
         Call mma_deallocate(iFlip)
         Call mma_deallocate(dBMx)
         Call mma_deallocate(Mult)
         Call mma_deallocate(cInt0)
         Call mma_deallocate(cInt)
         Call mma_deallocate(Value0)
         Call mma_deallocate(Value)
         Call mma_deallocate(dBVc)
         Call mma_deallocate(BVc)

#ifdef _DEBUGPRINT_
         Call RecPrt('BMx_t',' ',BMx_t,n3,nLambda_)
#endif
*
*        Assemble dr/dq: Solve  B dr/dq = dr/dx
*
         Call FZero(drdq,nLambda_*mInt_)
*
*        Temporary fix of the dC/dx vector which always
*        is propted up with the full degeneracy factor.
*
         If (.not.Curvilinear) Then
            Do iLambda=1,nLambda_
               Do i=1,n3
                  iAtom = (i+2)/3
                  ixyz  = i - (iAtom-1)*3
                  BMx_t(i,iLambda)=BMx_t(i,iLambda)/Degen(ixyz,iAtom)
               End Do
            End Do
         End If

         Call Eq_Solver('N',n3,mInt_,nLambda_,BMx,Curvilinear,Degen,
     &                  BMx_t,drdq)
#ifdef _DEBUGPRINT_
         Call RecPrt('drdq',' ',drdq,mInt_,nLambda_)
#endif
*
         Call mma_deallocate(BMx_t)
      End If
*
*     Double check that we don't have any null vectors
*
      iOff=1
      iOff2=1
      mLambda=nLambda_
      Do iLambda=1,nLambda_
         RR=Sqrt(DDot_(mInt_,drdq(1,iOff),1,drdq(1,iOff),1))
         If (RR.lt.1.0D-12) Then
            Write (6,*) 'Warning: constraint ',iLambda,
     &                  ' has a null vector, I''ll remove it!'
            mLambda=mLambda-1
         Else
            If (iOff.ne.iOff2)
     &         Call dCopy_(mInt_,drdq(1,iOff),1,drdq(1,iOff2),1)
            iOff2=iOff2+1
         End If
         iOff=iOff+1
      End Do
      If (mLambda.lt.nLambda_)
     &   Call FZero(drdq(1,mLambda+1),mInt_*(nLambda_-mLambda))
*
      End Subroutine get_drdq
