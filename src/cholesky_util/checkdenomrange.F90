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

integer function CheckDenomRange(xmin,xmax,nSym,EOcc,Evir,iOcc,nOcc,iVir,nVir)

implicit none
real*8 xmin
real*8 xmax
integer nSym
real*8 EOcc(nSym)
real*8 EVir(nSym)
integer iOcc(nSym)
integer nOcc(nSym)
integer iVir(nSym)
integer nVir(nSym)
real*8 Tol
parameter(Tol=1.0d-12)
real*8 e, emin, emax
integer iSym, i
integer aSym, a
integer irc

emin = 9.9d15
emax = -9.9d15
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
emin = 2.0d0*emin
emax = 2.0d0*emax

irc = 0
if (abs(emin-xmin) > Tol) irc = irc+1
if (abs(emax-xmax) > Tol) irc = irc+2

!-tbp:
if (irc /= 0) then
  write(6,'(A,1P,2D25.16)') 'xmin,xmax=',xmin,xmax
  write(6,'(A,1P,2D25.16)') 'emin,emax=',emin,emax
  write(6,'(A,1P,2D25.16)') 'diff=     ',xmin-emin,xmax-emax
end if

CheckDenomRange = irc

end function CheckDenomRange
