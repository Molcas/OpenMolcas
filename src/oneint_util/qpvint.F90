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

subroutine QpVInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of velocity quadrupole    *
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
real*8 TC(3)
integer iStabO(0:7), iDCRT(0:7)
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

#include "macros.fh"
unused_var(PtChrg)
unused_var(iAddPot)

nip = 1
ipB = nip
nip = nip+nZeta
ipS1 = nip
nip = nip+nZeta*nElem(la)*nElem(lb+1)*3
ipS2 = 1
if (lb > 0) then
  ipS2 = nip
  nip = nip+nZeta*nElem(la)*nElem(lb-1)*3
end if
ipRes = nip
nip = nip+nZeta*nElem(la)*nElem(lb)*nComp
if (nip-1 > nZeta*nArr) then
  call WarningMessage(2,' QpVInt: nip-1 > nZeta*nArr')
  call Abend()
end if
ipArr = nip
mArr = (nArr*nZeta-(nip-1))/nZeta

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)

iComp = 3
llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

ipOff = ipB
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(ipOff),nAlpha)
  ipOff = ipOff+1
end do

do lDCRT=0,nDCRT-1
  call OA(iDCRT(lDCRT),Ccoor,TC)

  nHer = (la+(lb+1)+(nOrdOp-1)+2)/2
  call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS1),nZeta,iComp,la,lb+1,A,RB,nHer,Array(ipArr),mArr,TC,nOrdOp-1)

  if (lb > 0) then
    nHer = (la+(lb-1)+(nOrdOp-1)+2)/2
    call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS2),nZeta,iComp,la,lb-1,A,RB,nHer,Array(ipArr),mArr,TC,nOrdOp-1)
  end if

  ! Combine derivatives and dipole integrals to generate the
  ! velocity quadrupole integrals.

  call Util5(Array(ipB),nZeta,Array(ipRes),la,lb,Array(ipS1),Array(ipS2))

  ! Accumulate contributions

  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array(ipRes),nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,One)

end do

return

end subroutine QpVInt
