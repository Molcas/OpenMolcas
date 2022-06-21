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

use Symmetry_Info, only: nIrrep, iChBas

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Local variables
parameter(mComp=200)
integer kOper(mComp), kChO(mComp)
! Statement function for Cartesian index
nElem(ixyz) = ((ixyz+1)*(ixyz+2))/2

iRout = 220
iPrint = nPrint(iRout)

iSize = nZeta*nElem(la)*nElem(lb)*nComp
call dcopy_(iSize,[Zero],0,final,1)
call dcopy_(nZeta*nArr,[Zero],0,Array,1)
nip = 1
ipB = nip
nip = nip+nZeta
ipS1 = nip
nip = nip+nZeta*nElem(la)*nElem(lb+1)*3*nIC
if (lb > 0) then
  ipS2 = nip
  nip = nip+nZeta*nElem(la)*nElem(lb-1)*3*nIC
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
if (kComp > mComp) then
  write(6,*) 'pxpint: kComp > mComp'
  call Abend()
end if
!write(6,*)
!write(6,*) 'pXpInt:**********'

iSym_p1 = IrrFnc(1)
iSym_p2 = IrrFnc(2)
iSym_p3 = IrrFnc(4)
!write(6,*) 'iSym_p=',iSym_p1,iSym_p2,iSym_p3
ipar_p1 = iChBas(2)
ipar_p2 = iChBas(3)
ipar_p3 = iChBas(4)
!write(6,*) 'ipar_p=',ipar_p1,ipar_p2,ipar_p3
do iComp=1,nComp
  jComp1 = (iComp-1)*3+1
  jComp2 = (iComp-1)*3+2
  jComp3 = (iComp-1)*3+3
  iTemp = lOper(iComp)
  ipar = iChO(iComp)

  jTemp1 = 0
  jTemp2 = 0
  jTemp3 = 0
  do iSym_pXp=0,nIrrep-1
    if (iand(2**iSym_pXp,iTemp) /= 0) then
      iSym_pX = ieor(iSym_pXp,iSym_p1)
      !write(6,*) 'iSym_pXp,iSym_pX=',iSym_pXp,iSym_pX
      jTemp1 = ior(jTemp1,2**iSym_pX)
    end if
    if (iand(2**iSym_pXp,iTemp) /= 0) then
      iSym_pX = ieor(iSym_pXp,iSym_p2)
      !write(6,*) 'iSym_pXp,iSym_pX=',iSym_pXp,iSym_pX
      jTemp2 = ior(jTemp2,2**iSym_pX)
    end if
    if (iand(2**iSym_pXp,iTemp) /= 0) then
      iSym_pX = ieor(iSym_pXp,iSym_p3)
      !write(6,*) 'iSym_pXp,iSym_pX=',iSym_pXp,iSym_pX
      jTemp3 = ior(jTemp3,2**iSym_pX)
    end if
  end do
  kOper(jComp1) = jTemp1
  kOper(jComp2) = jTemp2
  kOper(jComp3) = jTemp3

  kChO(jComp1) = ieor(ipar,ipar_p1)
  kChO(jComp2) = ieor(ipar,ipar_p2)
  kChO(jComp3) = ieor(ipar,ipar_p3)

end do
!                                                                      *
!***********************************************************************
!                                                                      *
call pXint(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS1),nZeta,kIC,kComp,la,lb+1,A,RB,iDum,Array(ipArr),mArr,CCoor, &
           kOrdOp,kOper,kChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)
!                                                                      *
!***********************************************************************
!                                                                      *
if (lb > 0) then
  call pXint(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS2),nZeta,kIC,kComp,la,lb-1,A,RB,iDum,Array(ipArr),mArr,CCoor, &
             kOrdOp,kOper,kChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
ipOff = ipB
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(ipOff),nAlpha)
  ipOff = ipOff+1
end do
!
if (iPrint >= 99) then
  call RecPrt(' In pXpint: Beta (expanded)','(5D20.13)',Array(ipB),nZeta,1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Combine pX integrals to generate the pXp integrals.
!
! Note that the pX integrals have 3*nComp components.

call Ass_pXp(Array(ipB),nZeta,final,la,lb,Array(ipS1),Array(ipS2),nComp)
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 49) call RecPrt('pXpInt: Final',' ',final,nZeta,nElem(la)*nElem(lb))

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nHer)

end subroutine pXpInt
