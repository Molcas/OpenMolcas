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
! Copyright (C) 1996, Per Ake Malmqvist                                *
!               1996, Roland Lindh                                     *
!***********************************************************************

subroutine AMPInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for computing matrix elements of the          *
!         six hermitized products of two angular momentum ops          *
!                                                                      *
!     Author: Per-AAke Malmqvist, Dept. of Theoretical Chemistry,      *
!             University of Lund, SWEDEN                               *
!             November '96                                             *
!     After pattern of other SEWARD soubroutines by R. Lindh.          *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Local Variables
real*8 TC(3)
integer iStabO(0:7), iDCRT(0:7)
! Statement function for Cartesian index
nElem(ixyz) = ((ixyz+1)*(ixyz+2))/2

iRout = 220
iPrint = nPrint(iRout)

nip = 1
ipB = nip
nip = nip+nZeta
ipTpp = nip
nip = nip+nZeta*nElem(la)*nElem(lb+2)*6
ipTp = nip
nip = nip+nZeta*nElem(la)*nElem(lb+1)*3
ipT = nip
nip = nip+nZeta*nElem(la)*nElem(lb)*6
ipTm = 1
ipTmm = 1
if (lb > 0) then
  ipTm = nip
  nip = nip+nZeta*nElem(la)*nElem(lb-1)*3
  if (lb > 1) then
    ipTmm = nip
    nip = nip+nZeta*nElem(la)*nElem(lb-2)*6
  end if
end if
ipRes = nip
nip = nip+nZeta*nElem(la)*nElem(lb)*nComp
if (nip-1 > nZeta*nArr) then
  call WarningMessage(2,' AMPInt: nip-1 > nZeta*nArr')
  call Abend()
end if
ipArr = nip
mArr = (nArr*nZeta-(nip-1))/nZeta

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)

ipOff = ipB
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(ipOff),nAlpha)
  ipOff = ipOff+1
end do

llOper = lOper(1)
do iComp=2,nComp
  iDum = lOper(iComp)
  llOper = ior(llOper,iDum)
end do

! Compute stabilizer, and then the double coset representation:
call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

! Loop over the cosets of the stabilizer group:
do lDCRT=0,nDCRT-1
  call OA(iDCRT(lDCRT),Ccoor,TC)

  ! Generate the quadrupole integral tables:
  iComp = 6
  iOrdOp = 2
  nHer = (la+(lb+2)+2+2)/2
  call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipTpp),nZeta,iComp,la,lb+2,A,RB,nHer,Array(ipArr),mArr,TC,iOrdOp)
  nHer = (la+lb+2+2)/2
  call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipT),nZeta,iComp,la,lb,A,RB,nHer,Array(ipArr),mArr,TC,iOrdOp)
  if (lb >= 2) then
    nHer = (la+(lb-2)+2+2)/2
    call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipTmm),nZeta,iComp,la,lb-2,A,RB,nHer,Array(ipArr),mArr,TC,iOrdOp)
  end if
  ! Generate the dipole integral tables:
  iComp = 3
  iOrdOp = 1
  nHer = (la+(lb+1)+1+2)/2
  call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipTp),nZeta,iComp,la,lb+1,A,RB,nHer,Array(ipArr),mArr,TC,iOrdOp)
  if (lb >= 1) then
    nHer = (la+(lb-1)+1+2)/2
    call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipTm),nZeta,iComp,la,lb-1,A,RB,nHer,Array(ipArr),mArr,TC,iOrdOp)
  end if

  if (iprint > 49) write(6,*) ' AMPInt calling AMPr.'
  call AMPr(Array(ipB),nZeta,Array(ipRes),la,lb,Array(ipTpp),Array(ipTp),Array(ipT),Array(ipTm),Array(ipTmm))

  ! Symmetry adaption:
  if (iprint > 49) write(6,*) ' AMPInt calling SymAdO'
  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array(ipRes),nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,One)
  if (iprint > 49) write(6,*) ' Back to AMPInt.'
end do

if (iprint > 49) write(6,*) ' Leaving AMPInt.'

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(nOrdOp)
  call Unused_real_array(PtChrg)
  call Unused_integer(iAddPot)
end if

end subroutine AMPInt
