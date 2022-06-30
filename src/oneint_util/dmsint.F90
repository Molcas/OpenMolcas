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

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: iComp, iDCRT(0:7), ipArr, ipRes, iPrint, ipS1, ipS2, iRout, iStabO(0:7), lDCRT, llOper, LmbdT, mArr, nComp_, &
                     nDCRT, nip, nOp, nRys, nStabO
real(kind=wp) :: TC(3,2)
integer(kind=iwp), external :: NrOpr

#include "macros.fh"
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 230
iPrint = nPrint(iRout)

nRys = nHer

if (iPrint >= 99) then
  call RecPrt(' In DMSInt: Alpha',' ',Alpha,nAlpha,1)
  call RecPrt(' In DMSInt: Beta',' ',Beta,nBeta,1)
end if

nip = 1
ipS1 = nip
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb+1)*3
ipS2 = nip
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb)*3
ipRes = nip
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb)*nComp
if (nip-1 > nZeta*nArr) then
  call WarningMessage(2,'DMSInt: nip-1 > nZeta*nArr')
  write(u6,*) 'nip=',nip
  write(u6,*) 'nZeta,nArr=',nZeta,nArr
  call Abend()
end if
ipArr = nip
mArr = nZeta*nArr-nip+1

rFinal(:,:,:,:) = Zero

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
  call SymAdO(Array(ipRes),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,One)

end do

return

end subroutine DMSInt
