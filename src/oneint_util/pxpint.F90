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
! Copyright (C) 1993, Bernd Artur Hess                                 *
!               1999, Roland Lindh                                     *
!***********************************************************************

subroutine pXpInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the comutation of pXp integrals           *
!                                                                      *
!     Author: Bernd Hess, Institut fuer Physikalische und Theoretische *
!             Chemie, University of Bonn, Germany, April 1993          *
!             R. Lindh, modified to molcas 4.1 form, Oct 1999          *
!***********************************************************************

use Symmetry_Info, only: iChBas, nIrrep
use Index_Functions, only: nTri_Elem1
use Integral_interfaces, only: int_kernel
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: iBeta, iComp, iDum, ipar, ipar_p1, ipar_p2, ipar_p3, ipArr, ipB, ipOff, iPrint, ipS1, ipS2, iRout, iSym_p1, &
                     iSym_p2, iSym_p3, iSym_pX, iSym_pXp, iTemp, jTemp1, jTemp2, jTemp3, kComp, kIC, kOrdOp, mArr, nip
integer(kind=iwp), allocatable :: kChO(:,:), kOper(:,:)
integer(kind=iwp), external :: IrrFnc
procedure(int_kernel) :: pXint

#include "macros.fh"
unused_var(nHer)

iRout = 220
iPrint = nPrint(iRout)

rFinal(:,:,:,:) = Zero
Array(:) = Zero
nip = 1
ipB = nip
nip = nip+nZeta
ipS1 = nip
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb+1)*3*nIC
if (lb > 0) then
  ipS2 = nip
  nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb-1)*3*nIC
else
  ipS2 = ipS1
end if
ipArr = nip
mArr = nArr-(nip-1)/nZeta
if (mArr < 0) then
  call WarningMessage(2,'pXpInt: mArr<0!')
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! nIC: the number of blocks of the symmetry adapted operator pXp
! nComp: number of components of the operator pXp
!
! pXp = d/dx X d/dx + d/dy X d/dy + d/dz X d/dz

kIC = nIC*3
kComp = nComp*3
kOrdOp = nOrdOp-1
call mma_allocate(kChO,3,nComp,label='kChO')
call mma_allocate(kOper,3,nComp,label='kOper')
!write(u6,*)
!write(u6,*) 'pXpInt:**********'

iSym_p1 = IrrFnc(1)
iSym_p2 = IrrFnc(2)
iSym_p3 = IrrFnc(4)
!write(u6,*) 'iSym_p=',iSym_p1,iSym_p2,iSym_p3
ipar_p1 = iChBas(2)
ipar_p2 = iChBas(3)
ipar_p3 = iChBas(4)
!write(u6,*) 'ipar_p=',ipar_p1,ipar_p2,ipar_p3
do iComp=1,nComp
  iTemp = lOper(iComp)
  ipar = iChO(iComp)

  jTemp1 = 0
  jTemp2 = 0
  jTemp3 = 0
  do iSym_pXp=0,nIrrep-1
    if (btest(iTemp,iSym_pXp)) then
      iSym_pX = ieor(iSym_pXp,iSym_p1)
      !write(u6,*) 'iSym_pXp,iSym_pX=',iSym_pXp,iSym_pX
      jTemp1 = ibset(jTemp1,iSym_pX)
      iSym_pX = ieor(iSym_pXp,iSym_p2)
      !write(u6,*) 'iSym_pXp,iSym_pX=',iSym_pXp,iSym_pX
      jTemp2 = ibset(jTemp2,iSym_pX)
      iSym_pX = ieor(iSym_pXp,iSym_p3)
      !write(u6,*) 'iSym_pXp,iSym_pX=',iSym_pXp,iSym_pX
      jTemp3 = ibset(jTemp3,iSym_pX)
    end if
  end do
  kOper(1,iComp) = jTemp1
  kOper(2,iComp) = jTemp2
  kOper(3,iComp) = jTemp3

  kChO(1,iComp) = ieor(ipar,ipar_p1)
  kChO(2,iComp) = ieor(ipar,ipar_p2)
  kChO(3,iComp) = ieor(ipar,ipar_p3)

end do
!                                                                      *
!***********************************************************************
!                                                                      *
call pXint(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS1),nZeta,kIC,kComp,la,lb+1,A,RB,iDum,Array(ipArr),mArr,CoorO, &
           kOrdOp,kOper,kChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)
!                                                                      *
!***********************************************************************
!                                                                      *
if (lb > 0) then
  call pXint(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS2),nZeta,kIC,kComp,la,lb-1,A,RB,iDum,Array(ipArr),mArr,CoorO, &
             kOrdOp,kOper,kChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)
end if
call mma_deallocate(kChO)
call mma_deallocate(kOper)
!                                                                      *
!***********************************************************************
!                                                                      *
ipOff = ipB-1
do iBeta=1,nBeta
  Array(ipOff+1:ipOff+nAlpha) = Beta(iBeta)
  ipOff = ipOff+nAlpha
end do

if (iPrint >= 99) then
  call RecPrt(' In pXpint: Beta (expanded)','(5D20.13)',Array(ipB),nZeta,1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Combine pX integrals to generate the pXp integrals.
!
! Note that the pX integrals have 3*nComp components.

call Ass_pXp(Array(ipB),nZeta,rFinal,la,lb,Array(ipS1),Array(ipS2),nComp)
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 49) call RecPrt('pXpInt: rFinal',' ',rFinal(:,:,:,1),nZeta,nTri_Elem1(la)*nTri_Elem1(lb))

return

end subroutine pXpInt
