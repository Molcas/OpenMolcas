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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine DMSInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of diamagnetic shielding  *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, February '91                            *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Local variables
real*8 TC(3,2)
integer iDCRT(0:7), iStabO(0:7)
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

iRout = 230
iPrint = nPrint(iRout)

nRys = nHer

if (iPrint >= 99) then
  call RecPrt(' In DMSInt: Alpha',' ',Alpha,nAlpha,1)
  call RecPrt(' In DMSInt: Beta',' ',Beta,nBeta,1)
end if

nip = 1
ipS1 = nip
nip = nip+nZeta*nElem(la)*nElem(lb+1)*3
ipS2 = nip
nip = nip+nZeta*nElem(la)*nElem(lb)*3
ipRes = nip
nip = nip+nZeta*nElem(la)*nElem(lb)*nComp
if (nip-1 > nZeta*nArr) then
  call WarningMessage(2,'DMSInt: nip-1 > nZeta*nArr')
  write(6,*) 'nip=',nip
  write(6,*) 'nZeta,nArr=',nZeta,nArr
  call Abend()
end if
ipArr = nip
mArr = nZeta*nArr-nip+1

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)

iComp = 1
llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

do lDCRT=0,nDCRT-1
  call OA(iDCRT(lDCRT),Ccoor(1:3,1),TC(1:3,1))
  call OA(iDCRT(lDCRT),Ccoor(1:3,2),TC(1:3,2))

  ! Compute contribution from a,b+1

  nComp_ = 1
  call EFPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS1),nZeta,nComp_,la,lb+1,A,RB,nRys,Array(ipArr),mArr,TC,nOrdOp-1)

  ! Compute contribution from a,b

  call EFPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS2),nZeta,nComp_,la,lb,A,RB,nRys,Array(ipArr),mArr,TC,nOrdOp-1)

  ! Assemble final integral from the derivative integrals

  call Util4(nZeta,Array(ipRes),la,lb,Array(ipS1),Array(ipS2),RB,TC(1,2))

  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array(ipRes),nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,One)

end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(PtChrg)
  call Unused_integer(iAddPot)
end if

end subroutine DMSInt
