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

subroutine FREQ(nX,H,nDeg,nrvec,Tmp3,EVec,EVal,RedM,iNeg)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, auTocm, uToau
use Definitions, only: wp, u6

implicit real*8(a-h,o-z)
real*8 H(*), Tmp3(nX,nX), EVec(2*nX,nX), EVal(2*nX), RedM(nX)
integer nrvec(*), ndeg(*)
logical Found
real*8, allocatable :: Mass(:)
! Statement function
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

! read masses from runfile

call Qpg_dArray('Isotopes',Found,nIsot)
if (.not. Found) then
  write(u6,*) 'No masses found on RunFile'
  call AbEnd()
end if
call mma_allocate(Mass,nIsot)
call Get_dArray('Isotopes',Mass,nIsot)

! form the Mass Weighted Cartesian force constant matrix.

iprint = 0
do i=1,nX
  rm = Mass(nrvec(i))
  if (rm == Zero) rm = 1.0e7_wp
  do j=1,nX
    Tmp3(i,j) = sqrt(real(nDeg(i)*nDeg(j),kind=wp))*H(itri(i,j))/rm
  end do
end do

! Compute eigenvectors and eigenfunctions

iOpt = 1
if (nX > 0) call Not_DGeEv(iOpt,Tmp3,nX,EVal,EVec,nX,nX)

! Compute the harmonic frequencies

iNeg = 0
do iHarm=1,2*nX,2
  jHarm = (iHarm+1)/2
  temp = EVal(iHarm)
  if (temp >= Zero) then
    EVal(jHarm) = sqrt(temp)*autocm
  else
    iNeg = iNeg+1
    EVal(jHarm) = -sqrt(abs(temp))*autocm
  end if
end do
if (iPrint >= 99) call RecPrt('Converted EVal',' ',EVal,1,nX)
if (iPrint >= 99) call RecPrt('Normalized EVec',' ',EVec,nX*2,nX)

! Normalize

do iHarm=1,nX
  r2 = Zero
  do i=1,nx
    rm = Mass(nrvec(i))/UTOAU
    r2 = r2+Evec(2*(i-1)+1,iHarm)*Evec(2*(i-1)+1,iHarm)*rm
  end do
  RedM(iHarm) = r2
  r2 = One/sqrt(r2)
  call DScal_(nX,r2,EVec(1,iHarm),2)
end do

! Order, from low to high.

do iHarm=1,nX-1
  do jHarm=iHarm+1,nX
    if (EVal(jHarm) < EVal(iHarm)) then
      rlow = EVal(iHarm)
      EVal(iHarm) = EVal(jHarm)
      EVal(jHarm) = rLow
      rlow = RedM(iHarm)
      RedM(iHarm) = RedM(jHarm)
      RedM(jHarm) = rLow
      call DSwap_(nX,EVec(1,iHarm),2,EVec(1,jHarm),2)
    end if
  end do
end do

call mma_deallocate(Mass)

return

end subroutine FREQ
