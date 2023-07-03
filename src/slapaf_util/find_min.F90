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

implicit real*8(a-h,o-z)
#include "real.fh"
#include "print.fh"
real*8 A(0:nOrder)
logical RC

iRout = 117
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt('Find_Min: A',' ',A,1,nOrder+1)
end if
Thr = 1.0D-12
XValue = XStart
RC = .true.
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
    tmp = dble(j)
    dfnc = dfnc+A(j)*tmp*XX
    XX = XX*X
  end do
  ddfnc = Zero
  XX = One
  do j=2,nOrder
    tmp = dble(j*j-j)
    ddfnc = ddfnc+A(j)*tmp*XX
    XX = XX*X
  end do
  XInc = dfnc/ddfnc
  XValue = XValue-XInc
  if (iPrint == 99) then
    write(6,*) 'Fnc,dFnc,ddFnc=',Fnc,dFnc,ddFnc
  end if
  if (abs(XInc) < Thr) then
    ENew = fnc
    XMin = XValue
    return
  end if
  if (XValue > XHi) XValue = XHi
  if (XValue < XLow) XValue = XLow
end do
if (iPrint >= 6) write(6,*) '-- Too many iterations in Find_Min'
RC = .false.

return

end subroutine Find_Min
