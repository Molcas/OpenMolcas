!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2013, Roland Lindh                                     *
!***********************************************************************

subroutine genCxCTL(iStop,Cartesian,rDelta)
!***********************************************************************
!                                                                      *
! subroutine for automatic generation of coordinates for numerical     *
! differentiation based on the rlxctl.f routine.                       *
!                                                                      *
! Author: R. Lindh, Uppsala University                                 *
!         2013, November                                               *
!***********************************************************************

use Slapaf_Info, only: Coor, Shift, qInt, BMx, Free_Slapaf
use Slapaf_Parameters, only: Curvilinear, HSet, BSet, PrQ, Numerical, nLambda, iRef, nDimBC, mTROld, mTtAtm, nWndw, iter
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: iStop
logical(kind=iwp), intent(out) :: Cartesian
real(kind=wp), intent(in) :: rDelta
#include "print.fh"
integer(kind=iwp) :: iDisp, iRow_c, jInter, Jter, LuSpool, mInt
logical(kind=iwp) :: Found, TSC, Error
real(kind=wp), allocatable :: DList(:), CList(:,:), du(:), TMx(:), RefCoor(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
!-----Process the input

LuSpool = 21
call SpoolInp(LuSpool)

call RdCtl_Slapaf(LuSpool,.true.)
mInt = nDimBC-mTROld
Curvilinear = .false.
Cartesian = .not. Curvilinear
Numerical = .false. ! Just to define it, value is irrelevant here!

call Close_LuSpool(LuSpool)
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Compute the Wilson B-matrix, these describe the transformations
!     between internal and Cartesian coordinates. Values of the
!     Internal coordinates are computed too.

BSet = .true.
HSet = .false.
PrQ = .false.

nWndw = iter
iRef = 0
call BMtrx(size(Coor,2),Coor,iter,mTtAtm,nWndw)

nPrint(30) = nPrint(30)-1
!                                                                      *
!***********************************************************************
!                                                                      *
call Put_dArray('BMtrx',BMx,size(Coor)*mInt)
call Put_iScalar('No of Internal coordinates',mInt)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Generate list of coordinates for numerical differentiation
!                                                                      *
!***********************************************************************
!                                                                      *
! Make lists for all Cartesian coordinates and the
! corresponding displacement in internal coordinates

call mma_allocate(CList,size(Coor),2*mInt,Label='CList')
CList(:,:) = Zero
call mma_allocate(DList,mInt,Label='DList')
DList(:) = Zero
call mma_allocate(RefCoor,3,size(Coor,2),Label='RefCoor')
!                                                                      *
!***********************************************************************
!                                                                      *
! If in TS-search regime and without TSConstraints, remove
! constraints (in slapaf this is done differently)

call qpg_iScalar('TS Search',Found)
if (Found) call Get_lScalar('TS Search',Found)
call f_Inquire('TSC',TSC)
if (Found .and. (.not. TSC)) call Merge_Constraints('','','UDC',nLambda,iRow_c)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Get the T-matrix

call mma_allocate(TMx,mInt**2,Label='TMx')
call TMatrix(TMx,mInt)
call Put_iScalar('nLambda',nLambda)
call Put_dArray('T-matrix',TMx,mInt**2)
!call RecPrt('T-matrix',' ',TMx,mInt,mInt)
!                                                                      *
!***********************************************************************
!                                                                      *
! Now start generating the displaced structure to be used in the
! numerical differentiation.
!
! Take a copy of the current structure - the reference
! coordinates.

RefCoor(:,:) = Coor(:,:)

! Loop over all displacements which are in the subspace in
! which we like to minimize the energy. Hence, this will
! eliminate naturally translational and rotational degrees
! (3N-6) but also eliminate constrained degrees (3N-6-m)

call mma_allocate(du,mInt,Label='du')

! Loop only over displacement which do not change the constraint.
! Note that the constraints are placed first.

do iDisp=1+2*nLambda,2*mInt

  ! Get a virgin copy of the reference structure

  Coor(:,:) = RefCoor(:,:)

  ! Compute the effective index where to find the data

  Jter = Iter+1

  ! Update the shift vector, in the space in which we split
  ! the constraints and the reduced subspace.

  du(:) = Zero
  Shift(:,iter) = Zero
  jInter = (iDisp+1)/2
  if (mod(iDisp,2) == 0) then
    du(jInter) = -rDelta
  else
    du(jInter) = rDelta
  end if
  !call RecPrt('du',' ',du,mInt,1)

  ! Transform displacement to the internal coordinate
  ! space. This is a simple unitary transformation.

  call DGEMM_('N','N',mInt,1,mInt,One,TMx,mInt,du,mInt,Zero,Shift(:,iter),mInt)
  !call RecPrt('shf',' ',Shift(:,iter),mInt,1)

  ! Save the value of the displacement in the list.

  DList(jInter) = rDelta

  ! Take a copy of the current values of the internal
  ! coordinates.

  call dcopy_(mInt,qInt(:,Iter),1,qInt(:,Jter),1)
  !call RecPrt('Int_Ref',' ',qInt(:,Jter),1,mInt)

  ! To the second set of coordinates add the shift.
  ! This set of internal coordinates corresponds to
  ! the set for which we like to get the Cartesian
  ! coordinates.

  call DaXpY_(mInt,One,Shift(:,iter),1,qInt(:,Jter),1)
  !call RecPrt('Int    ',' ',qInt(:,Jter),1,mInt)

  !------Transform the new internal coordinates to Cartesians

  PrQ = .false.
  BSet = .false.
  Error = .false.
  nWndw = Iter
  iRef = 0
  call NewCar(Iter,size(Coor,2),Coor,mTtAtm,Error)

  ! Move the new Cartesian coordinate to the list.

  call dcopy_(size(Coor),Coor,1,CList(:,iDisp),1)
end do

call mma_deallocate(du)
call mma_deallocate(TMx)

!call RecPrt('DList',' ',DList,1,mInt)
!call RecPrt('CList',' ',CList,SIZE(Coor),2*mInt)

! Save the lists on the runfile. To be used in the
! numerical gradient module.

call Put_dArray('DList',DList,size(DList))
call Put_dArray('CList',CList,size(CList))

! Deallocate temporary memory.

call mma_deallocate(RefCoor)
call mma_deallocate(DList)
call mma_deallocate(CList)

! Alaska only
iStop = 3

! Done!

call Free_Slapaf()

!-----Terminate the calculations.

return

end subroutine GenCxCTL
