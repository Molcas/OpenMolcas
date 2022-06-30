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
! Copyright (C) 2015, Lasse Kragh Soerensen                            *
!               2015, Roland Lindh                                     *
!***********************************************************************

subroutine OMQInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of orbital magnetic       *
!         quadrupole integrals => OMQInt                               *
!                                                                      *
!     Author: Lasse Kragh Soerensen and Roland Lindh  2015             *
!             Based on OAMInt                                          *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "int_interface.fh"
integer(kind=iwp) :: iBeta, iComp, iDCRT(0:7), ipArr, ipB, ipOff, ipRes, ipS1, ipS2, ipS3, iStabO(0:7), lDCRT, llOper, LmbdT, &
                     mArr, nDCRT, nip, nOp, nStabO
real(kind=wp) :: TC(3)
integer(kind=iwp), external :: NrOpr

#include "macros.fh"
unused_var(PtChrg)
unused_var(iAddPot)

nip = 1
ipB = nip
nip = nip+nZeta

! L + 1 component
ipS1 = nip
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb+1)*6 ! not ncomp

! L - 1 component
ipS2 = 1
if (lb > 0) then
  ipS2 = nip
  nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb-1)*6 ! not ncomp
end if

! L + 0 component
ipS3 = nip
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb)*3

ipRes = nip
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb)*nComp
if (nip-1 > nZeta*nArr) then
  call WarningMessage(2,' OMQInt: nip-1 > nZeta*nArr')
  call Abend()
end if
ipArr = nip
mArr = (nArr*nZeta-(nip-1))/nZeta

rFinal(:,:,:,:) = Zero

llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

ipOff = ipB-1
do iBeta=1,nBeta
  Array(ipOff+1:ipOff+nAlpha) = Beta(iBeta)
  ipOff = ipOff+nAlpha
end do

do lDCRT=0,nDCRT-1
  call OA(iDCRT(lDCRT),CCoor,TC)

  iComp = 6 ! Why are these here ncomp is passed down?

  nHer = (la+(lb+1)+(nOrdOp-1)+2)/2
  call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS1),nZeta,iComp,la,lb+1,A,RB,nHer,Array(ipArr),mArr,TC,nOrdOp-1)

  if (lb > 0) then
    nHer = (la+(lb-1)+(nOrdOp-1)+2)/2
    call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS2),nZeta,iComp,la,lb-1,A,RB,nHer,Array(ipArr),mArr,TC,nOrdOp-1)
  end if

  iComp = 3 ! Why are these here ncomp is passed down?

  nHer = (la+lb+(nOrdOp-2)+2)/2
  ! check to see dipole integral sure looks a lot like dipole integrals
  call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS3),nZeta,iComp,la,lb,A,RB,nHer,Array(ipArr),mArr,TC,nOrdOp-2)

  ! Combine derivatives of dipole integrals to generate the
  ! orbital angular momentum integrals.

  call Util3(Array(ipB),nZeta,Array(ipRes),la,lb,Array(ipS1),Array(ipS3),Array(ipS2))

  ! Accumulate contributions

  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array(ipRes),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,One)

end do

return

end subroutine OMQInt
