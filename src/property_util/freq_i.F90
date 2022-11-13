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

subroutine FREQ_i(nX,H,mass,EVec,EVal,iNeg)

use Constants, only: Zero, One, autocm
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nX
real(kind=wp), intent(inout) :: H(nX,nX)
real(kind=wp), intent(in) :: mass(*)
real(kind=wp), intent(out) :: EVec(2*nX,nX), EVal(2*nX)
integer(kind=iwp), intent(out) :: iNeg
integer(kind=iwp) :: i, iHarm, ii, iOpt, iprint, j, jHarm
real(kind=wp) :: r2, rlow, temp
real(kind=wp), external :: DDot_

iprint = 0
do i=1,nX
  ii = (i-1)/3+1
  do j=1,nX
    H(i,j) = H(i,j)/mass(ii)
  end do
end do

! Compute eigenvectors and eigenfunctions

iOpt = 1
if (nX > 0) then
  call Not_DGeEv(iOpt,H,nX,EVal,EVec,nX,nX)
end if

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

! Normalize

do iHarm=1,nX
  r2 = DDot_(nX,EVec(1,iHarm),2,EVec(1,iHarm),2)
  r2 = One/sqrt(r2)
  call DScal_(nX,r2,EVec(1,iHarm),2)
end do
if (iPrint >= 99) call RecPrt('Normalized EVec',' ',EVec,nX*2,nX)

! Order, from low to high.

do iHarm=1,nX-1
  do jHarm=iHarm+1,nX
    if (EVal(jHarm) < EVal(iHarm)) then
      rlow = EVal(iHarm)
      EVal(iHarm) = EVal(jHarm)
      EVal(jHarm) = rLow
      call DSwap_(nX,EVec(1,iHarm),2,EVec(1,jHarm),2)
    end if
  end do
end do

return

end subroutine FREQ_i
