!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Acore(iang,la,ishll,nordop,TC,A,Array,narr,Alpha,nalpha,fa1,fa2,jfgrad,jfhess,ld,debug)
! Calculates <A'|core> and <A"|core>
!
! @parameter iang Angular momenta for core
! @parameter la Angular momenta for bra
! @parameter ishll identification for core shell
! @parameter nordop order for operator
! @parameter TC Cartesian coordinates for core
! @parameter A Cartesian coordinates for bra
! @parameter Array Scratch
! @parameter narr size for scratch
! @parameter Alpha Bra exponents
! @parameter nalpha number of exponents
! @parameter FA1 First derivatives (out)
! @parameter FA2 2nd derivatives (out)
! @parameter jfgrad true for all 1-derivatives that are needed
! @parameter jfhess true for all 2-derivatives that are needed
! @parameter ld Order of derivatives
! @parameter debug guess

use Index_Functions, only: nTri_Elem1
use Basis_Info, only: Shells
use Her_RW, only: HerR, HerW, iHerR, iHerW
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iang, la, ishll, nordop, narr, nalpha, ld
real(kind=wp), intent(in) :: TC(3), A(3), Alpha(nAlpha)
real(kind=wp), intent(_OUT_) :: Array(*), fa1(*)
real(kind=wp), intent(inout) :: fa2(*)
logical(kind=iwp), intent(in) :: jfgrad(3), jfhess(4,3,4,3), debug
integer(kind=iwp) :: i, iGamma, ip, ipA, ipAxyz, ipCxyz, ipK1, ipP1, ipQ1, ipRxyz, ipV, ipZ1, ipZI1, iStrt, n, nExpi, nHer, nVecAC
logical(kind=iwp) :: ABeq(3)
real(kind=wp), external :: DNrm2_

nExpi = Shells(iShll)%nExp
ip = 1
ipP1 = ip
ip = ip+3*nAlpha*nExpi
ipZ1 = ip
ip = ip+nAlpha*nExpi
ipK1 = ip
ip = ip+nAlpha*nExpi
ipZI1 = ip
ip = ip+nAlpha*nExpi
if (ip-1 > nArr) then
  write(u6,*) ' ip-1 > nArr in acore  (',ip,',',narr,')'
  call Abend()
end if

! Calculate Effective center and exponent for <A|core>

call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,nExpi,Alpha,Shells(iShll)%Exp)
call SetUp1(Alpha,nAlpha,Shells(iShll)%Exp,nExpi,A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))

! Calculate Overlap <A|core> and derivative <A'|core>

nHer = (la+1+iAng+1+ld)/2
ipAxyz = ip
ip = ip+nAlpha*nExpi*3*nHer*(la+1+ld)
ipCxyz = ip
ip = ip+nAlpha*nExpi*3*nHer*(iAng+1)
ipRxyz = ip
ip = ip+nAlpha*nExpi*3*nHer*(nOrdOp+1)
ipQ1 = ip
ip = ip+nAlpha*nExpi*3*(la+1+ld)*(iAng+1)*(nOrdOp+1)
ipA = ip
ip = ip+nAlpha*nExpi
if (ip-1 > nArr) then
  write(u6,*) '  ip-1 > nArr (1b) in acore (',ip,',',narr,')','Order',ld,Shells(ishll)%nExp,nalpha
  call Abend()
end if
ABeq(1) = A(1) == TC(1)
ABeq(2) = A(2) == TC(2)
ABeq(3) = A(3) == TC(3)
call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*nExpi,A,Array(ipAxyz),la+ld,HerR(iHerR(nHer)),nHer,ABeq)
call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*nExpi,TC,Array(ipCxyz),iAng,HerR(iHerR(nHer)),nHer,ABeq)
ABeq(1) = .false.
ABeq(2) = .false.
ABeq(3) = .false.
call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*nExpi,A,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)
if (debug) then
  write(u6,*) ' nAlpha = ',nAlpha,' nExp(',ishll,')=',nExpi,' nHer=',nHer,' la=',la,' iAng=',iAng,' nOrdOp=',nOrdOp

  write(u6,*) ' Array(ipAxyz)=',DNrm2_(nAlpha*nExpi*3*nHer*(la+ld+1),Array(ipAxyz),1)
  write(u6,*) ' Array(ipCxyz)=',DNrm2_(nAlpha*nExpi*3*nHer*(iAng+1),Array(ipCxyz),1)
  write(u6,*) ' Array(ipRxyz)=',DNrm2_(nAlpha*nExpi*3*nHer*(nOrdOp+1),Array(ipRxyz),1)
end if

call Assmbl(Array(ipQ1),Array(ipAxyz),la+ld,Array(ipRxyz),nOrdOp,Array(ipCxyz),iAng,nAlpha*nExpi,HerW(iHerW(nHer)),nHer)
iStrt = ipA
do iGamma=1,nExpi
  Array(iStrt:iStrt+nAlpha-1) = Alpha
  iStrt = iStrt+nAlpha
end do
if (debug) then
  write(u6,*) ' Array(ipA)=',DNrm2_(nAlpha*nExpi,Array(ipA),1)
end if

call rKappa_Zeta(Array(ipK1),Array(ipZ1),nExpi*nAlpha)
call CmbnAC(Array(ipQ1),nAlpha*nExpi,la,iAng,Array(ipK1),FA1,Array(ipA),JfGrad,ld,nVecAC)
if (debug) then
  write(u6,*) 'nVecAC',nvecac
  write(u6,*) ' Array(ipQ1)=',DNrm2_(nAlpha*nExpi*3*(la+ld+1)*(iAng+1)*(nOrdOp+1),Array(ipQ1),1)
  write(u6,*) ' Array(ipA)=',DNrm2_(nAlpha*nExpi,Array(ipA),1)
  do i=1,nvecac
    ipV = 1
    n = nAlpha*nExpi*nTri_Elem1(la)*nTri_Elem1(iAng)
    write(u6,*) 'Cmbn(',i,')=',DNrm2_(n,FA1(ipV),1)
    ipV = ipV+n
  end do
end if

if (ld >= 2) then
  call CmbnS2a(Array(ipQ1),nAlpha*nExpi,la,iAng,Array(ipK1),FA2,Array(ipA),jfHess,ld)
  if (debug) then
    do i=1,6
      ipV = 1
      n = nAlpha*nExpi*nTri_Elem1(la)*nTri_Elem1(iAng)
      write(u6,*) 'Cmbn2(',i,')=',DNrm2_(n,FA2(ipV),1)
      ipV = ipV+n
    end do
  end if
end if

return

end subroutine Acore
