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
************************************************************************
      Subroutine genCxCTL(iStop,Cartesian,rDelta)
      use Slapaf_Info, only: Coor, Shift, qInt, BMx, Free_Slapaf
      use Slapaf_Parameters, only: Curvilinear, HSet, BSet, PrQ,
     &                             Numerical, nLambda, iRef, nDimBC,
     &                             mTROld, mTtAtm, nWndw, iter
      Implicit Real*8 (a-h,o-z)
************************************************************************
*                                                                      *
*     subroutine for automatic generation of coordinates for numerical *
*     differentiation based on the rlxctl.f routine.                   *
*                                                                      *
*     Author: R. Lindh, Uppsala University                             *
*             2013, November                                           *
************************************************************************
#include "info_slapaf.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
#include "nadc.fh"
#include "weighting.fh"
#include "db.fh"
#include "print.fh"
      Logical Cartesian, Found, TSC, Error
      Real*8, Allocatable:: DList(:), CList(:,:), du(:), TMx(:)
*                                                                      *
************************************************************************
*                                                                      *
*
      Lu=6
*
      iRout = 32
      iPrint=nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
*-----Process the input
*
      LuSpool=21
      Call SpoolInp(LuSpool)
*
      Call RdCtl_Slapaf(LuSpool,.True.)
      mInt = nDimBC - mTROld
      Curvilinear=.FALSE.
      Cartesian  =.NOT.Curvilinear
      Numerical = .False. ! Just to define it, value is irrelevant here!
*
      Call Close_LuSpool(LuSpool)
*                                                                      *
************************************************************************
*                                                                      *
      jPrint=nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the Wilson B-matrix, these describe the transformations
*     between internal and Cartesian coordinates. Values of the
*     Internal coordinates are computed too.
*
      BSet=.True.
      HSet=.False.
      PrQ=.False.
*
      nWndw=iter
      iRef=0
      Call BMtrx(SIZE(Coor,2),Coor,iter,mTtAtm,nWndw)
*
      nPrint(30) = nPrint(30)-1
*                                                                      *
************************************************************************
*                                                                      *
      Call Put_dArray('BMtrx',BMx,SIZE(Coor)*mInt)
      Call Put_iScalar('No of Internal coordinates',mInt)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Generate list of coordinates for numerical differentiation
*                                                                      *
************************************************************************
*                                                                      *
*     Make lists for all Cartesian coordinates and the
*     corresponding displacement in internal coordinates
*
      Call mma_allocate(CList,SIZE(Coor), 2*mInt,Label='CList')
      CList(:,:)=Zero
      Call mma_allocate(DList,mInt,Label='DList')
      DList(:)=Zero
      Call Allocate_Work(ip_RefCoor,SIZE(Coor))
*                                                                      *
************************************************************************
*                                                                      *
*     If in TS-search regime and without TSConstraints, remove
*     constraints (in slapaf this is done differently)
*
      Call qpg_iScalar('TS Search',Found)
      If (Found) Call Get_lScalar('TS Search',Found)
      Call f_Inquire('TSC',TSC)
      If (Found.and..not.TSC)
     &   Call Merge_Constraints('','','UDC',nLambda,iRow_c)
*                                                                      *
************************************************************************
*                                                                      *
*     Get the T-matrix
*
      Call mma_allocate(TMx,mInt**2,Label='TMx')
      Call TMatrix(TMx,mInt)
      Call Put_iScalar('nLambda',nLambda)
      Call Put_dArray('T-matrix',TMx,mInt**2)
*     Call RecPrt('T-matrix',' ',TMx,mInt,mInt)
*                                                                      *
************************************************************************
*                                                                      *
*     Now start generating the displaced structure to be used in the
*     numerical differentiation.
*
*     Take a copy of the current structure - the reference
*     coordinates.
*
      call dcopy_(SIZE(Coor),Coor,1,Work(ip_RefCoor),1)
*
*     Loop over all displacements which are in the subspace in
*     which we like to minimize the energy. Hence, this will
*     eliminate naturally translational and rotational degrees
*     (3N-6) but also eliminate constrained degrees (3N-6-m)
*
      nDisp = mInt-nLambda
      Call mma_allocate(du,mInt,Label='du')
*
*     Loop only over displacement which do not change the constraint.
*     Note that the constraints are placed first.
*
      Do iDisp = 1+2*nLambda, 2*mInt
*
*        Get a virgin copy of the reference structure
*
         call dcopy_(SIZE(Coor),Work(ip_RefCoor),1,Coor,1)
*
*        Compute the effective index where to find the data
*
         Jter=Iter+1
*
*        Update the shift vector, in the space in which we split
*        the constraints and the reduced subspace.
*
         du(:)=Zero
         Shift(:,iter)=Zero
         jInter = (iDisp+1)/2
         If (Mod(iDisp,2).eq.0) Then
            du(jInter) = -rDelta
         Else
            du(jInter) =  rDelta
         End If
*        Call RecPrt('du',' ',du,mInt,1)
*
*        Transform displacement to the internal coordinate
*        space. This is a simple unitary transformation.
*
         Call DGEMM_('N','N',mInt,1,mInt,
     &               One,TMx,mInt,
     &                   du,mInt,
     &               Zero,Shift(:,iter),mInt)
*        Call RecPrt('shf',' ',Shift(:,iter),mInt,1)
*
*        Save the value of the displacement in the list.
*
         DList(jInter) = rDelta
*
*        Take a copy of the current values of the internal
*        coordinates.
*
         call dcopy_(mInt,qInt(:,Iter),1,qInt(:,Jter),1)
*        Call RecPrt('Int_Ref',' ',qInt(:,Jter),1,mInt)
*
*        To the second set of coordinates add the shift.
*        This set of internal coordinates corresponds to
*        the set for which we like to get the Cartesian
*        coordinates.
*
         Call DaXpY_(mInt,One,Shift(:,iter),1,qInt(:,Jter),1)
*        Call RecPrt('Int    ',' ',qInt(:,Jter),1,mInt)
*
*--------Transform the new internal coordinates to Cartesians
*
         PrQ=.False.
         BSet=.False.
         Error=.False.
         nWndw=Iter
         iRef=0
         Call NewCar(Iter,SIZE(Coor,2),Coor,mTtAtm,Error)
*
*        Move the new Cartesian coordinate to the list.
*
         call dcopy_(SIZE(Coor),Coor,1,CList(:,iDisp),1)
      End Do
*
      Call mma_deallocate(du)
      Call mma_deallocate(TMx)
*
*     Call RecPrt('DList',' ',DList,1,mInt)
*     Call RecPrt('CList',' ',CList,SIZE(Coor),2*mInt)
*
*     Save the lists on the runfile. To be used in the
*     numerical gradient module.
*
      Call Put_dArray('DList',DList,SIZE(DList))
      Call Put_dArray('CList',CList,SIZE(CList))
*
*     Deallocate temporary memory.
*
      Call Free_Work(ip_RefCoor)
      Call mma_deallocate(DList)
      Call mma_deallocate(CList)

*     Alaska only
      iStop=3
*
*     Done!
*
      Call Free_Slapaf()
*
*-----Terminate the calculations.
*
      Return
      End
