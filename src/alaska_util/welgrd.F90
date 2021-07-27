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
! Copyright (C) 1992,1995, Roland Lindh                                *
!***********************************************************************

subroutine WelGrd( &
#                 define _CALLING_
#                 include "grd_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: to compute the Pauli repulsion integrals with the            *
!         Gauss-Hermite quadrature.                                    *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden. October '92.                            *
!                                                                      *
!             Modified to gradients, April '95. R. Lindh               *
!***********************************************************************

use Center_Info

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "wldata.fh"
#include "print.fh"
#include "disp.fh"
#include "grd_interface.fh"
! Statement function for Cartesian index
nElem(i) = (i+1)*(i+2)/2

iRout = 122
iPrint = nPrint(iRout)
!iQ = 1
if (iPrint >= 59) then
  write(6,*) ' In WelGrd'
  write(6,*) ' r0, ExpB=',r0,ExpB
  write(6,*) ' la,lb=',la,lb
  write(6,*) '  A=',A
  write(6,*) ' RB=',RB
end if

k = la+lb+1
jsump = 1
do i=1,k
  jsump = jsump+3**i
end do
k0 = max(la+lb-1,0)
jsumm = 1
do i=1,k0
  jsumm = jsumm+3**i
end do

ip = 1
ipGri = ip
ip = ip+nZeta*jsump
ipTGri = ip
ip = ip+nZeta*jsump
ipGrin = ip
ip = ip+nZeta*(k+1)*(k/2+1)*(k/4+1)
iPxyz = ip
ip = ip+nZeta
if (ip-1 > nZeta*nArr) then
  write(6,*) ' ip-1 > nZeta*nArr(pos.1)'
  write(6,*) ip-1,'>',nZeta*nArr
  call ErrTra()
  call Abend()
end if

call Rowel(nZeta,r0,expB,k,Zeta,P,Array(iPxyz),Array(ipGri),Array(ipGrin),jsump)
ip = ip-nZeta
ip = ip-nZeta*(k+1)*(k/2+1)*(k/4+1)
if (iPrint >= 99) call RecPrt(' In WelInt: Array(ipGri)l',' ',Array(ipGri),nZeta,jSumP)

ipAMx = ip
ip = ip+nZeta*9
ipScr = ip
ip = ip+nZeta*3**k
if (ip-1 > nZeta*nArr) then
  write(6,*) ' ip-1 > nZeta*nArr(pos.2)'
  write(6,*) ip-1,'>',nZeta*nArr
  call ErrTra()
  call Abend()
end if

! Transform each block to the global coordinate system

iOff = ipgri+nZeta
do ik=1,k
  if (ik == 1) call SetUpA(nZeta,Array(ipAMx),P)
  call Traxyz(nZeta,ik,Array(iOff),Array(ipScr),Array(ipAMx))
  iOff = iOff+nZeta*3**ik
end do
if (iPrint >= 99) call RecPrt(' In WelInt: Array(ipGri)g',' ',Array(ipGri),nZeta,jSumP)
call dcopy_(nZeta*jsump,Array(ipGri),1,Array(ipTGri),1)
ip = ip-nZeta*3**k
ip = ip-nZeta*9

ip1 = ip
ip = ip+nZeta
ip2 = ip
ip = ip+nZeta
ip3 = ip
ip = ip+nZeta
ip4 = ip
ip = ip+nZeta
ip5 = ip
ip = ip+nZeta
if (ip-1 > nZeta*nArr) then
  write(6,*) ' ip-1 > nZeta*nArr(pos.3)'
  write(6,*) ip-1,'>',nZeta*nArr
  call ErrTra()
  call Abend()
end if

! Compute <a|O|b+1>

ip0p = ip
ip = ip+nZeta*nElem(la)*nElem(lb+1)
call TraPAB(nZeta,la,lb+1,Array(ip0p),Array(ipgri),jSump,rKappa,Array(ip1),Array(ip2),Array(ip3),Array(ip4),Array(ip5),A,RB,P)

! Compute <a|O|b-1>

if (lb >= 1) then
  ip0m = ip
  ip = ip+nZeta*nElem(la)*nElem(lb-1)
  call dcopy_(nZeta*jsumm,Array(ipTGri),1,Array(ipGri),1)
  call TraPAB(nZeta,la,lb-1,Array(ip0m),Array(ipgri),jSumm,rKappa,Array(ip1),Array(ip2),Array(ip3),Array(ip4),Array(ip5),A,RB,P)
else
  ip0m = 1
end if

! Compute <a+1|O|b>

ipp0 = ip
ip = ip+nZeta*nElem(la+1)*nElem(lb)
call dcopy_(nZeta*jsump,Array(ipTGri),1,Array(ipGri),1)
call TraPAB(nZeta,la+1,lb,Array(ipp0),Array(ipgri),jSump,rKappa,Array(ip1),Array(ip2),Array(ip3),Array(ip4),Array(ip5),A,RB,P)

! Compute <a-1|O|b>

if (la >= 1) then
  ipm0 = ip
  ip = ip+nZeta*nElem(la-1)*nElem(lb)
  call dcopy_(nZeta*jsumm,Array(ipTGri),1,Array(ipGri),1)
  call TraPAB(nZeta,la-1,lb,Array(ipm0),Array(ipgri),jSumm,rKappa,Array(ip1),Array(ip2),Array(ip3),Array(ip4),Array(ip5),A,RB,P)
else
  ipm0 = 1
end if

ipAlph = ip
ip = ip+nZeta
ipBeta = ip
ip = ip+nZeta

jp = ipAlph
do iBeta=1,nBeta
  call dcopy_(nAlpha,Alpha,1,Array(jp),1)
  jp = jp+nAlpha
end do
jp = ipBeta
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(jp),nAlpha)
  jp = jp+1
end do

! Assemble the derivative integrals and distribute contributions.

call CmbnW1(Array(ipp0),Array(ipm0),Array(ip0p),Array(ip0m),nZeta,la,lb,Zeta,rKappa,Final,Array(ipAlph),Array(ipBeta),Grad,nGrad, &
            DAO,IfGrad,IndGrd,dc(mdc)%nStab,dc(ndc)%nStab,kOp)

ip = ip-5*nZeta
ip = ip-2*nZeta
if (la >= 1) ip = ip-nZeta*nElem(la-1)*nElem(lb)
ip = ip-nZeta*nElem(la)*nElem(lb+1)
if (lb >= 1) ip = ip-nZeta*nElem(la)*nElem(lb-1)
ip = ip-nZeta*nElem(la+1)*nElem(lb)
ip = ip-2*nZeta*jsump

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(ZInv)
  call Unused_integer(nHer)
  call Unused_real_array(Ccoor)
  call Unused_integer(nOrdOp)
  call Unused_integer_array(lOper)
  call Unused_integer_array(iStabM)
end if

end subroutine WelGrd
