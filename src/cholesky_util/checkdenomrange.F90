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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

function CheckDenomRange(xmin,xmax,nSym,EOcc,Evir,iOcc,nOcc,iVir,nVir)

use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: CheckDenomRange
integer(kind=iwp), intent(in) :: nSym, iOcc(nSym), nOcc(nSym), iVir(nSym), nVir(nSym)
real(kind=wp), intent(in) :: xmin, xmax, EOcc(nSym), EVir(nSym)
integer(kind=iwp) :: a, aSym, i, irC, iSym
real(kind=wp) :: e, emax, emin
real(kind=wp), parameter :: Tol = 1.0e-12_wp

emin = 9.9e15_wp
emax = -9.9e15_wp
do iSym=1,nSym
  do i=iOcc(iSym)+1,iOcc(iSym)+nOcc(iSym)
    do aSym=1,nSym
      do a=iVir(aSym)+1,iVir(aSym)+nVir(aSym)
        e = EVir(a)-EOcc(i)
        emin = min(emin,e)
        emax = max(emax,e)
      end do
    end do
  end do
end do
emin = Two*emin
emax = Two*emax

irc = 0
if (abs(emin-xmin) > Tol) irc = irc+1
if (abs(emax-xmax) > Tol) irc = irc+2

!-tbp:
if (irc /= 0) then
  write(u6,'(A,2ES25.16)') 'xmin,xmax=',xmin,xmax
  write(u6,'(A,2ES25.16)') 'emin,emax=',emin,emax
  write(u6,'(A,2ES25.16)') 'diff=     ',xmin-emin,xmax-emax
end if

CheckDenomRange = irc

end function CheckDenomRange
