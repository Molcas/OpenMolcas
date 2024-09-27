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

subroutine Int_Setup(iSD,nSkal,iS,jS,kS,lS,Coor,Shijij,iAngV,iCmpV,iShelV,iShllV,iAOV,iStabs)

use Basis_Info, only: dbsc
use Gateway_Info, only: DoFMM, RPQMin
use Gateway_global, only: FMM_shortrange
use iSD_data, only: nSD
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSkal, iSD(0:nSD,nSkal), iS, jS, kS, lS
real(kind=wp), intent(out) :: Coor(3,4)
logical(kind=iwp), intent(out) :: Shijij
integer(kind=iwp), intent(out) :: iAngV(4), iCmpV(4), iShelV(4), iShllV(4), iAOV(4), iStabs(4)
integer(kind=iwp) :: i, iCnt, iCnttp, iQuad, iSkal, jCnt, jCnttp, jQuad(4), kCnt, kCnttp, lCnt, lCnttp
real(kind=wp) :: D, P, Q

iCnttp = iSD(13,iS)
iCnt = iSD(14,iS)
jCnttp = iSD(13,jS)
jCnt = iSD(14,jS)
kCnttp = iSD(13,kS)
kCnt = iSD(14,kS)
lCnttp = iSD(13,lS)
lCnt = iSD(14,lS)

if (dbsc(iCnttp)%Aux) then
  Coor(:,1) = dbsc(jCnttp)%Coor(:,jCnt)
else
  Coor(:,1) = dbsc(iCnttp)%Coor(:,iCnt)
end if
Coor(:,2) = dbsc(jCnttp)%Coor(:,jCnt)

if (dbsc(kCnttp)%Aux) then
  Coor(:,3) = dbsc(lCnttp)%Coor(:,lCnt)
else
  Coor(:,3) = dbsc(kCnttp)%Coor(:,kCnt)
end if
Coor(:,4) = dbsc(lCnttp)%Coor(:,lCnt)

Shijij = ((iSD(0,iS) == iSD(0,kS)) .and. (iSD(10,iS) == iSD(10,kS)) .and. (iSD(0,jS) == iSD(0,lS)) .and. (iSD(10,jS) == iSD(10,lS)))

jQuad(1) = iS
jQuad(2) = jS
jQuad(3) = kS
jQuad(4) = lS
do iQuad=1,4
  iSkal = jQuad(iQuad)
  iAngV(iQuad) = iSD(1,iSkal)
  iCmpV(iQuad) = iSD(2,iSkal)
  iAOV(iQuad) = iSD(7,iSkal)
  iStabs(iQuad) = iSD(10,iSkal)
  iShelV(iQuad) = iSD(11,iSkal)
  iShllV(iQuad) = iSD(0,iSkal)
end do
!MAW start

! For the FMM coulomb integrals <AB(r1)|1/r12|CD(r2)>
! Here we flag the integral routines that we only want to compute
! the short-range non-multipole component of integrals over this
! shell quartet if midpoint(A,B) is sufficiently far from
! midpoint(C,D) for numerical stability.
! Note that midpoint(X,Y) corresponds to the multipole expansion
! centre of an XY AO-pair, regardless of exponents.

FMM_shortrange = .false.
if (DoFMM) then
  D = Zero
  do i=1,3
    P = Half*(Coor(i,1)+Coor(i,2))    ! AB shell-pair
    Q = Half*(Coor(i,3)+Coor(i,4))    ! CD shell-pair
    D = D+(P-Q)*(P-Q)
  end do
  if (D > RPQMIN*RPQMIN) FMM_shortrange = .true.
end if
!MAW end

return

end subroutine Int_Setup
