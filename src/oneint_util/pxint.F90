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
! Copyright (C) 2006, Roland Lindh                                     *
!***********************************************************************

subroutine PXInt( &
#                define _CALLING_
#                include "int_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of  pX integrals          *
!                                                                      *
! Called from: OneEl                                                   *
!                                                                      *
! Author: Roland Lindh, Dept. Chem. Phys., Lund University,            *
!         June 2006                                                    *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChBas
use Index_Functions, only: nTri_Elem1
use Integral_interfaces, only: int_kernel
use Oneint_interfaces, only: PVInt
use Property_Label, only: PLabel
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
integer(kind=iwp), parameter :: mComp = 200
integer(kind=iwp) :: iComp, ipar_p1, ipar_p2, ipar_p3, iSym_p1, iSym_p2, iSym_p3, iSym_px, iSym_X, jComp1, jComp2, jComp3, &
                     jpar_p1, jpar_p2, jpar_p3, jTemp1, jTemp2, jTemp3, kComp, kIC, kOrdOp, nRys
integer(kind=iwp), allocatable :: kChO(:), kOper(:)
integer(kind=iwp), external :: IrrFnc
procedure(int_kernel) :: CntInt, EFInt, MltInt, NAInt

!                                                                      *
!***********************************************************************
!                                                                      *
! nIC: number of symmetry adapted blocks in total for the nComp
!      elements of the compund operator, pX.
! nComp: is the number of elements of the compund operator
!
! kIC: number of symmetry adapted blocks in total for the kComp
!      elements of the operator X
! kComp: is the number of elements of the operator X.
!                                                                      *
!***********************************************************************
!                                                                      *
! Note that the p operator's each element is only a basis function
! of a single irreducible representation. Hence, 3*kIC=nIC.
!
! In addition if the operator X has kComp elements pX has 3*kComp
! elements.
!                                                                      *
!***********************************************************************
!                                                                      *
nRys = nHer
kIC = nIC/3
kComp = nComp/3
kOrdOp = nOrdOp-1
!                                                                      *
!***********************************************************************
!                                                                      *
! Now produce the kOper array with kComp elements from the lOper
! array. Ditto kChO/iChO.
!
! lOper is an integer which bit pattern indicate to which irreps
! the component of the operator is a basis function. Note that the
! operator is not symmetry adapted, i.e. it can be a basis function
! in more than one irrep.
!
! iChO is an integer which describe the parity character of the
! operator with respect to X, Y, and Z coordinates. For example,
! if the first bit is set this means that the operator change sign
! under a reflection  in the yz-plane, etc.

call mma_allocate(kChO,kComp,label='kChO')
call mma_allocate(kOper,kComp,label='kOper')

! As we remove the p operator (three of them) X should be the same
! regardless of if we remove d/dx, d/dy, or d/dz.

iSym_p1 = IrrFnc(1)   ! d/dx
iSym_p2 = IrrFnc(2)   ! d/dy
iSym_p3 = IrrFnc(4)   ! d/dz
!write(u6,*)
!write(u6,*) 'pXInt******'
!write(u6,*) 'iSym_p=',iSym_p1,iSym_p2,iSym_p3
ipar_p1 = iChBas(2)
ipar_p2 = iChBas(3)
ipar_p3 = iChBas(4)
!write(u6,*) 'ipar_p=',ipar_p1,ipar_p2,ipar_p3
do iComp=1,kComp
  jComp1 = (iComp-1)*3+1
  jComp2 = (iComp-1)*3+2
  jComp3 = (iComp-1)*3+3
  jpar_p1 = iChO(jComp1)
  jpar_p2 = iChO(jComp2)
  jpar_p3 = iChO(jComp3)

  ! Look thru all irreps and check if pX is a basis function in
  ! irrep iSym_pX. If so find the symmetry to which X is a basis function

  jTemp1 = 0
  jTemp2 = 0
  jTemp3 = 0
  !write(u6,*) 'lOper=',lOper(jComp1),lOper(jComp2),lOper(jComp3)
  do iSym_pX=0,nIrrep-1
    if (btest(lOper(jComp1),iSym_pX)) then
      iSym_X = ieor(iSym_pX,iSym_p1)
      !write(u6,*) 'iSym_pX,iSym_X=',iSym_pX,iSym_X
      jTemp1 = ibset(jTemp1,iSym_X)
      !write(u6,*) 'jTemp1=',jTemp1
    end if
    if (btest(lOper(jComp2),iSym_pX)) then
      iSym_X = ieor(iSym_pX,iSym_p2)
      !write(u6,*) 'iSym_pX,iSym_X=',iSym_pX,iSym_X
      jTemp2 = ibset(jTemp2,iSym_X)
      !write(u6,*) 'jTemp2=',jTemp2
    end if
    if (btest(lOper(jComp3),iSym_pX)) then
      iSym_X = ieor(iSym_pX,iSym_p3)
      !write(u6,*) 'iSym_pX,iSym_X=',iSym_pX,iSym_X
      jTemp3 = ibset(jTemp3,iSym_X)
      !write(u6,*) 'jTemp3=',jTemp3
    end if
  end do

  ! Check for consistency!

  if ((jTemp1 /= jTemp2) .or. (jTemp1 /= jTemp3)) then
    call WarningMessage(2,'PXInt: corrupted jTemps!')
    write(u6,*) 'jTemp1,jTemp2,jTemp3=',jTemp1,jTemp2,jTemp3
    call Abend()
  end if

  ! Compute the parity of X

  jpar_p1 = ieor(jpar_p1,ipar_p1)
  jpar_p2 = ieor(jpar_p2,ipar_p2)
  jpar_p3 = ieor(jpar_p3,ipar_p3)

  if ((jpar_p1 /= jpar_p2) .or. (jpar_p1 /= jpar_p3)) then
    call WarningMessage(2,'PXInt: corrupted jpars!')
    call Abend()
  end if

  ! Store the data

  kOper(iComp) = jTemp1
  kChO(iComp) = jpar_p1
end do

!write(u6,*) 'pXpInt'
!do iComp=1,nComp
!  write(u6,*) lOper(iComp),iChO(iComp)
!end do
!write(u6,*)
!do iComp=1,kComp
!  write(u6,*) kOper(iComp),kChO(iComp)
!end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute now the integrals

if (PLabel == 'NAInt ') then
  call PVInt(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,rFinal,nZeta,kIC,kComp,la,lb,A,RB,nRys,Array,nArr,CoorO,kOrdOp,kOper,kChO, &
             iStabM,nStabM,PtChrg,nGrid,iAddPot,NAInt)
else if (PLabel == 'MltInt') then
  call PVInt(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,rFinal,nZeta,kIC,kComp,la,lb,A,RB,nRys,Array,nArr,CoorO,kOrdOp,kOper,kChO, &
             iStabM,nStabM,PtChrg,nGrid,iAddPot,MltInt)
else if (PLabel == 'EFInt ') then
  call PVInt(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,rFinal,nZeta,kIC,kComp,la,lb,A,RB,nRys,Array,nArr,CoorO,kOrdOp,kOper,kChO, &
             iStabM,nStabM,PtChrg,nGrid,iAddPot,EFInt)
else if (PLabel == 'CntInt') then
  call PVInt(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,rFinal,nZeta,kIC,kComp,la,lb,A,RB,nRys,Array,nArr,CoorO,kOrdOp,kOper,kChO, &
             iStabM,nStabM,PtChrg,nGrid,iAddPot,CntInt)
else
  call WarningMessage(2,'PXInt: Illegal type!')
  write(u6,*) '       PLabel=',PLabel
  call Abend()
end if

call mma_deallocate(kChO)
call mma_deallocate(kOper)

return

end subroutine PXInt
