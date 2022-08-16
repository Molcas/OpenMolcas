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

subroutine coreB(iang,lb,ishll,nordop,TC,RB,Array,narr,Beta,nbeta,fb1,fb2,jfgrad,jfhess,ld,debug)
!  Calculates <core|B'> and <core|B">
!
! @parameter iang Angular momenta for core
! @parameter lb Angular momenta for ket
! @parameter ishll identification for core shell
! @parameter nordop order for operator
! @parameter TC Cartesian coordinates for core
! @parameter RB Cartesian coordinates for ket
! @parameter Array Scratch
! @parameter narr size for scratch
! @parameter Beta Ket exponents
! @parameter nbeta number of exponents
! @parameter FB1 First derivatives (out)
! @parameter FB2 2nd derivatives (out)
! @parameter jfgrad true for all 1-deriavtives that are needed
! @parameter jfhess true for all 2-deriavtives that are needed
! @parameter ld Order of derivatives
! @parameter debug guess

use Basis_Info
use Her_RW, only: HerR, HerW, iHerR, iHerW
use Real_Spherical

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "disp.fh"
logical ABeq(3), jfgrad(3), jfhess(4,3,4,3), debug
real*8 TC(3), RB(3), Array(*), fb1(*), fb2(*), beta(*)
! Statement function
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

nExpi = Shells(iShll)%nExp
if (debug) then
  write(6,*) 'Shell: ',ishll,' nBeta:',nbeta,' nExp:',nExpi,'Angular',lb,iang
end if

ip = 1
ipP2 = ip
ip = ip+3*nExpi*nBeta
ipZ2 = ip
ip = ip+nExpi*nBeta
ipK2 = ip
ip = ip+nExpi*nBeta
ipZI2 = ip
ip = ip+nExpi*nBeta
if (ip-1 > nArr) then
  write(6,*) '  ip-1 > nArr*nZeta(2) in bcore (',ip,',',narr,')'
  call Abend()
end if

! Calculate Effective center and exponent for <core|B>

call ZXia(Array(ipZ2),Array(ipZI2),nExpi,nBeta,Shells(iShll)%Exp,Beta)
call SetUp1(Shells(iShll)%Exp,nExpi,Beta,nBeta,TC,RB,Array(ipK2),Array(ipP2),Array(ipZI2))

! Calculate Overlap <core|B> and <core|B'>

nHer = (lb+1+iAng+1+ld)/2
ipCxyz = ip
ip = ip+nBeta*nExpi*3*nHer*(iAng+1)
ipBxyz = ip
ip = ip+nBeta*nExpi*3*nHer*(lb+1+ld)
ipRxyz = ip
ip = ip+nBeta*nExpi*3*nHer*(nOrdOp+1)
ipQ1 = ip
ip = ip+nBeta*nExpi*3*(iAng+1)*(lb+1+ld)*(nOrdOp+1)
ipB = ip
ip = ip+nBeta*nExpi
if (ip-1 > nArr) then
  write(6,*) '  ip-1 > nArr*nZeta(2b) in PrjGrd'
  call Abend()
end if
ABeq(1) = TC(1) == RB(1)
ABeq(2) = TC(2) == RB(2)
ABeq(3) = TC(3) == RB(3)
if (debug) write(6,*) 'shll=',ishll,' nExp=',nExpi,' nBeta=',nBeta
call CrtCmp(Array(ipZ2),Array(ipP2),nExpi*nBeta,TC,Array(ipCxyz),iAng,HerR(iHerR(nHer)),nHer,ABeq)
call CrtCmp(Array(ipZ2),Array(ipP2),nExpi*nBeta,RB,Array(ipBxyz),lb+ld,HerR(iHerR(nHer)),nHer,ABeq)
ABeq(1) = .false.
ABeq(2) = .false.
ABeq(3) = .false.
call CrtCmp(Array(ipZ2),Array(ipP2),nExpi*nBeta,TC,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)
if (debug) then
  write(6,*) ' nbeta  = ',nbeta,' nExp(',ishll,')=',nExpi,' nHer=',nHer,' lb=',lb,' iAng=',iAng,' nOrdOp=',nOrdOp

  write(6,*) ' Array(ipCxyz)=',DNrm2_(nBeta*nExpi*3*nHer*(iAng+1),Array(ipCxyz),1)
  write(6,*) ' Array(ipBxyz)=',DNrm2_(nBeta*nExpi*3*nHer*(lb+2),Array(ipBxyz),1)
  write(6,*) ' Array(ipRxyz)=',DNrm2_(nBeta*nExpi*3*nHer*(nOrdOp+1),Array(ipRxyz),1)
end if

call Assmbl(Array(ipQ1),Array(ipCxyz),iAng,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb+ld,nExpi*nBeta,HerW(iHerW(nHer)),nHer)
iStrt = ipB
do iGamma=1,nExpi
  call dcopy_(nBeta,Beta,1,Array(iStrt),nExpi)
  iStrt = iStrt+1
end do
if (debug) then
  write(6,*) ' Array(ipB)=',DNrm2_(nExpi*nBeta,Array(ipB),1)
end if

call rKappa_Zeta(Array(ipK2),Array(ipZ2),nExpi*nBeta)
call CmbnCB(Array(ipQ1),nExpi*nBeta,iAng,lb,Array(ipK2),FB1,Array(ipB),JfGrad,ld,nVecCB)
if (debug) then
  write(6,*) ' Array(ipQ1)=',DNrm2_(nExpi*nBeta*3*(lb+1+ld+2)*(iAng+1)*(nOrdOp+1),Array(ipQ1),1)
  write(6,*) ' Array(ipB)=',DNrm2_(nExpi*nBeta,Array(ipB),1)
end if
if (ld >= 2) then
  call CmbnS2b(Array(ipQ1),nBeta*nExpi,iang,lb,Array(ipK2),FB2,Array(ipB),jfHess,ld)
  if (debug) then
    do i=1,6
      ipV = 1
      n = nBeta*nExpi*nElem(lb)*nElem(iAng)
      write(6,*) n,nBeta,nExpi,nElem(lb),nElem(iAng)

      write(6,*) 'CmbnB2(',n,')=',DNrm2_(n,FB2(ipV),1)
      ipV = ipV+n
    end do
  end if
end if

return

end subroutine coreB
