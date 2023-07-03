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

subroutine LnSrch(Energy,q,dq,g,nInter,nIter,dqHdq)

implicit real*8(a-h,o-z)
#include "real.fh"
#include "print.fh"
real*8 Energy(nIter), A(0:4), Projg(2), q(nInter,nIter), g(nInter,nIter), dq(nInter,nIter), B(4,0:3), FVal(4)
logical RC

iRout = 200
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  write(6,*) ' Enter LnSrch'
  write(6,*) 'dqHdq=',dqHdq
  call RecPrt('LnSrch: Energy',' ',Energy,1,nIter)
  call RecPrt('LnSrch: q',' ',q,nInter,nIter)
  call RecPrt('LnSrch:dq',' ',dq,nInter,nIter)
  call RecPrt('LnSrch: g',' ',g,nInter,nIter)
end if

Thr = 1.0d-12
A(0) = Zero

iNew = nIter
iOld = nIter-1

! 1) Quartic search, -0.5 <= Xmin <= 0.5
! 2) Cubic search, -1.5 <= Xmin <= 1.5
! 3) E(iNew) is the best, no action
! 4) E(iNew) is not the best,
!    XMin=0.0 between iNew and the Best point

! Quartic line search

Projg(1) = -DDot_(nInter,q(1,iNew),1,g(1,iOld),1)+DDot_(nInter,q(1,iOld),1,g(1,iOld),1)
Projg(2) = -DDot_(nInter,q(1,iNew),1,g(1,iNew),1)+DDot_(nInter,q(1,iOld),1,g(1,iNew),1)
A(3) = -(Two*(Energy(iNew)-Energy(iOld))+(Projg(2)+Projg(1)))
A(1) = (Energy(iNew)-Energy(iOld))-(One/Four)*A(3)
!write(6,*) ' Projg(1),Projg(2),A(3),A(1)=',Projg(1),Projg(2),A(3),A(1)

! Test the gradient difference

Test1 = Projg(2)-Projg(1)
!write(6,*) ' Test1=',Test1
if (Test1 < Thr) then
  if (iPrint >= 6) write(6,*) '-- Line search failed, negative curvature'
  return
end if

! Now check if data is consistent with the condition that
! the second derivative is zero at the minimum.

Test2 = Test1**2-Three*A(3)**2
!write(6,*) ' Test2=',Test2
if (Test2 < Zero) then
  if (iPrint >= 6) write(6,*) '-- Quartic line search failed, nonzero 2nd derivative'
  Go To 100
end if

A(2) = (sqrt(Test2)+Test1)/Four
if (A(2) < Thr) then
  if (iPrint >= 6) write(6,*) '-- Quartic line search failed, A(2) too small'
  Go To 100
end if

A(4) = Half*(Test1-Test2)
nOrder = 4
XStart = Zero
XLow = -Half
XHi = Half
Go To 500

! Cubic line search

! Range: -0.5 - 0.5

100 continue
if (iPrint >= 6) write(6,*) '-- Cubic line search'
FVal(1) = Energy(iOld)
FVal(2) = Energy(iNew)
FVal(3) = Projg(1)
FVal(4) = Projg(2)
!B = Zero
x0 = One
x1 = Half
x2 = x1*Half
x3 = x2*Half
x0g = Zero
x1g = One
x2g = Two*x1
x3g = Three*x2
Fact = One
ii = 1
jj = 3
do i=1,2
  Fact = -One*Fact
  B(ii,0) = x0
  B(ii,1) = Fact*x1
  B(ii,2) = x2
  B(ii,3) = Fact*x3
  B(jj,0) = x0g
  B(jj,1) = x1g
  B(jj,2) = Fact*x2g
  B(jj,3) = x3g
  ii = ii+1
  jj = jj+1
end do
nOrder = 3
!call RecPrt('B-Matrix',' ',B,4,4)
!call RecPrt('FVal',' ',FVal,4,1)
call Gauss(nOrder+1,nOrder+1,B,A,FVal)
!call RecPrt('A-Coefficients',' ',A,4,1)
XStart = Zero
XLow = -One
XHi = Two+Half

500 continue
call Find_Min(nOrder,XStart,A,XMin,RC,XLow,XHi,ENew)
!write(6,*) 'ENew=',ENew
if (.not. RC) then
  !write(6,*) 'RC=.False.'
  return
end if
if (Test2 < Zero) then
  dqHdq = ENew-Energy(nIter)
else
  dqHdq = ENew
end if
XMin = XMin+Half
if (iPrint >= 6) write(6,*) 'Minimum found at -->',XMin,'<--'
!                                                                      *
!***********************************************************************
!                                                                      *
! Update vectors in accordance with the result

! Compute the displacement iOld -> iNew

call dcopy_(nInter,q(1,iNew),1,dq(1,iNew-1),1)
call DaXpY_(nInter,-One,q(1,iOld),1,dq(1,iNew-1),1)
call DScal_(nInter,XMin,dq(1,iNew-1),1)
dqdq = DDot_(nInter,dq(1,iNew-1),1,dq(1,iNew-1),1)

! Compute the new q

call dcopy_(nInter,q(1,iOld),1,q(1,iNew),1)
call DaXpY_(nInter,One,dq(1,iNew-1),1,q(1,iNew),1)

! Update the gradient

call DScal_(nInter,XMin,g(1,iNew),1)
call DaXpY_(nInter,One-XMin,g(1,iOld),1,g(1,iNew),1)
gdq = -DDot_(nInter,g(1,iNew),1,dq(1,iNew-1),1)
call DaXpY_(nInter,gdq/dqdq,dq(1,iNew-1),1,g(1,iNew),1)

! Update the displacement nIter-1 -> nIter

call dcopy_(nInter,q(1,iNew),1,dq(1,iNew-1),1)
call DaXpY_(nInter,-One,q(1,iNew-1),1,dq(1,iNew-1),1)
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 99) then
  call RecPrt('LnSrch: q',' ',q,nInter,nIter)
  call RecPrt('LnSrch:dq',' ',dq,nInter,nIter)
  call RecPrt('LnSrch: g',' ',g,nInter,nIter)
  write(6,*) ' Exit LnSrch'
end if

return

end subroutine LnSrch
