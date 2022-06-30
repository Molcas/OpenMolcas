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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

subroutine WelInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: to compute the Pauli repulsion integrals with the            *
!         Gauss-Hermite quadrature.                                    *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden. October '92.                            *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "wldata.fh"
#include "print.fh"
integer(kind=iwp) :: i, ik, iOff, ip, ip1, ip2, ip3, ip4, ip5, ipA, ipGri, ipGrin, iPrint, ipScr, iPxyz, iRout, jsum, k

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(ZInv)
unused_var(nHer)
unused_var(Ccoor)
unused_var(nOrdOp)
unused_var(lOper)
unused_var(iChO)
unused_var(iStabM)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 122
iPrint = nPrint(iRout)
!iQ = 1
if (iPrint >= 59) then
  write(u6,*) ' In WelInt'
  write(u6,*) ' r0, ExpB=',r0,ExpB
  write(u6,*) ' la,lb=',la,lb
end if

k = la+lb
jsum = 1
do i=1,k
  jsum = jsum+3**i
end do

ip = 1
ipGri = ip
ip = ip+nZeta*jsum
ipGrin = ip
ip = ip+nZeta*(k+1)*(k/2+1)*(k/4+1)
iPxyz = ip
ip = ip+nZeta
if (ip-1 > nZeta*nArr) then
  call WarningMessage(2,'WelInt:  ip-1 > nZeta*nArr(pos.1)')
  write(u6,*) ip-1,'>',nZeta*nArr
  call Abend()
end if

call Rowel(nZeta,r0,expB,k,Zeta,P,Array(iPxyz),Array(ipGri),Array(ipGrin),jsum)
ip = ip-nZeta
ip = ip-nZeta*(k+1)*(k/2+1)*(k/4+1)

ipA = ip
ip = ip+nZeta*9
ipScr = ip
ip = ip+nZeta*3**k
if (ip-1 > nZeta*nArr) then
  call WarningMessage(2,'WelInt:  ip-1 > nZeta*nArr(pos.2)')
  write(u6,*) ip-1,'>',nZeta*nArr
  call Abend()
end if

! Transform each block to the global coordinate system

iOff = ipgri+nZeta
do ik=1,k
  if (ik == 1) call SetUpA(nZeta,Array(ipA),P)
  call Traxyz(nZeta,ik,Array(iOff),Array(ipScr),Array(ipA))
  iOff = iOff+nZeta*3**ik
end do
if (iPrint >= 99) call RecPrt(' In WelInt: Array(ipGri)',' ',Array(ipGri),nZeta,jSum)
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
  call WarningMessage(2,'WelInt:  ip-1 > nZeta*nArr(pos.3)')
  write(u6,*) ip-1,'>',nZeta*nArr
  call Abend()
end if
call TraPAB(nZeta,la,lb,rFinal,Array(ipgri),jSum,rKappa,Array(ip1),Array(ip2),Array(ip3),Array(ip4),Array(ip5),A,RB,P)
ip = ip-nZeta*5
ip = ip-nZeta*jsum

return

end subroutine WelInt
