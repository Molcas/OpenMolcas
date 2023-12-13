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
! Copyright (C) 2002, Roland Lindh                                     *
!***********************************************************************

subroutine dTdmu_int( &
#                    define _CALLING_
#                    include "int_interface.fh"
                    )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of diamagnetic shielding  *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemical Physics, University      *
!             of Lund, Sweden, September 2002.                         *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: iBeta, iComp, iDCRT(0:7), ipArr, ipB, ipOff, ipRes, iPrint, ipS1, ipS2, iRout, iStabO(0:7), lDCRT, llOper, &
                     LmbdT, mArr, nDCRT, nip, nOp, nRys, nStabO
real(kind=wp) :: TC(3,2)
integer(kind=iwp), external :: NrOpr

#include "macros.fh"
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 230
iPrint = nPrint(iRout)

nRys = nHer

if (iPrint >= 99) then
  call RecPrt(' In dTdmu_int: Alpha',' ',Alpha,nAlpha,1)
  call RecPrt(' In dTdmu_int: Beta',' ',Beta,nBeta,1)
end if

nip = 1
ipS1 = nip
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb+1)*3
ipS2 = nip
if (lb >= 1) nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb-1)*3
ipRes = nip
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb)*nComp
ipB = nip
nip = nip+nZeta
if (nip-1 > nZeta*nArr) then
  call WarningMessage(2,'dTdmu_int: nip-1 > nZeta*nArr')
  write(u6,*) 'nip=',nip
  write(u6,*) 'nZeta,nArr=',nZeta,nArr
  call Abend()
end if
ipArr = nip
mArr = nZeta*nArr-nip+1

rFinal(:,:,:,:) = Zero

ipOff = ipB-1
do iBeta=1,nBeta
  Array(ipOff+1:ipOff+nAlpha) = Beta(iBeta)
  ipOff = ipOff+nAlpha
end do

iComp = 1
llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

do lDCRT=0,nDCRT-1
  call OA(iDCRT(lDCRT),CoorO(1:3,1),TC(1:3,1))
  call OA(iDCRT(lDCRT),CoorO(1:3,2),TC(1:3,2))

  ! Compute contribution from a,b+1

  call EFPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS1),nZeta,nComp,la,lb+1,A,RB,nRys,Array(ipArr),mArr,TC,nOrdOp)

  ! Compute contribution from a,b-1

  if (lb >= 1) &
    call EFPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS2),nZeta,nComp,la,lb-1,A,RB,nRys,Array(ipArr),mArr,TC,nOrdOp)

  ! Assemble final integral from the derivative integrals

  call Assemble_dTdmu(nZeta,Array(ipRes),la,lb,Array(ipS1),Array(ipS2),Array(ipB))

  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array(ipRes),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,One)

end do

return

end subroutine dTdmu_int
