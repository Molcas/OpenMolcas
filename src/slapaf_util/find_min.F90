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

subroutine Find_Min(nOrder,XStart,A,XMin,RC,XLow,XHi,ENew)

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nOrder
real(kind=wp), intent(in) :: XStart, A(0:nOrder), XLow, XHi
real(kind=wp), intent(out) :: XMin, ENew
logical(kind=iwp), intent(out) :: RC
#include "print.fh"
integer(kind=iwp) :: i, iPrint, iRout, j, MaxIter
real(kind=wp) :: ddfnc, dfnc, fnc, tmp, X, XInc, XValue, XX
real(kind=wp), parameter :: Thr = 1.0e-12_wp

iRout = 117
iPrint = nPrint(iRout)
if (iPrint >= 99) call RecPrt('Find_Min: A',' ',A,1,nOrder+1)
XValue = XStart
RC = .false.
MaxIter = 100
do i=1,MaxIter
  X = XValue
  fnc = Zero
  XX = One
  do j=0,nOrder
    fnc = fnc+A(j)*XX
    XX = XX*X
  end do
  dfnc = Zero
  XX = One
  do j=1,nOrder
    tmp = real(j,kind=wp)
    dfnc = dfnc+A(j)*tmp*XX
    XX = XX*X
  end do
  ddfnc = Zero
  XX = One
  do j=2,nOrder
    tmp = real(j*j-j,kind=wp)
    ddfnc = ddfnc+A(j)*tmp*XX
    XX = XX*X
  end do
  XInc = dfnc/ddfnc
  XValue = XValue-XInc
  if (iPrint == 99) write(u6,*) 'Fnc,dFnc,ddFnc=',Fnc,dFnc,ddFnc
  if (abs(XInc) < Thr) then
    ENew = fnc
    XMin = XValue
    RC = .true.
    exit
  end if
  if (XValue > XHi) XValue = XHi
  if (XValue < XLow) XValue = XLow
end do
if ((.not. RC) .and. (iPrint >= 6)) write(u6,*) '-- Too many iterations in Find_Min'

return

end subroutine Find_Min
