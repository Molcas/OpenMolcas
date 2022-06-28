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

! All points. A bunch of special cases, see ABOne and ABBoth for more details.
subroutine ABNone(iLA,iLB,dMulA,Rinv,Colle)

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, Two, Three
use Definitions, only: wp, iwp

! Maximum multipole implemented
#define _MxM_ 2

implicit none
integer(kind=iwp), intent(in) :: iLA, iLB
real(kind=wp), intent(in) :: dMulA(nTri_Elem1(_MxM_)), Rinv
real(kind=wp), intent(out) :: Colle(3)
real(kind=wp) :: d3, Pi1, Pi2, Sigma

Colle(:) = Zero

if (iLA == 0) then
  if (iLB == 0) then
    Sigma = dMulA(1)
    Colle(1) = Sigma*Rinv
  else if (iLB == 1) then
    Sigma = dMulA(1)
    Colle(1) = Sigma*Rinv**2
  else if (iLB == 2) then
    Sigma = dMulA(1)
    Colle(1) = Sigma*Rinv**3
  end if
else if (iLA == 1) then
  if (iLB == 0) then
    Sigma = dMulA(3)
    Colle(1) = Sigma*Rinv**2
  else if (iLB == 1) then
    Sigma = dMulA(3)
    Pi1 = dMulA(1)
    Pi2 = dMulA(2)
    Colle(1) = Two*Sigma*(Rinv**3)
    Colle(2) = Pi1*(Rinv**3)
    Colle(3) = Pi2*(Rinv**3)
  else if (iLB == 2) then
    d3 = sqrt(Three)
    Sigma = dMulA(3)
    Pi1 = dMulA(1)
    Pi2 = dMulA(2)
    Colle(1) = Three*Sigma*(Rinv**4)
    Colle(2) = d3*Pi1*(Rinv**4)
    Colle(3) = d3*Pi2*(Rinv**4)
  end if

! Jose* This is for Quadrupoles in Classical. We do not use in QmStat.
!--------------------------------------
!else if (iLA == 2) then
!  if (iLB == 0) then
!    Sigma = dMulA(3)
!    Colle(1) = Sigma*Rinv**3
!  else if (iLB == 1) then
!    d3 = sqrt(Three)
!    Sigma = dMulA(3)
!    Pi1 = dMulA(2)
!    Pi2 = dMulA(4)
!    Colle(1) = Three*Sigma*(Rinv**4)
!    Colle(2) = d3*Pi1*(Rinv**4)
!    Colle(3) = d3*Pi2*(Rinv**4)
!  else if (iLB == 2) then
!    Sigma = dMulA(3)
!    Pi1 = dMulA(2)
!    Pi2 = dMulA(4)
!    ! Jose. Remember dMulB(1)=sqrt(3)*xy
!    !            and dMulB(5)=Half*sqrt(3)*(x2-y2)
!    Del1 = dMulA(1)
!    Del2 = dMulA(5)
!    Colle(1) = Six*Sigma*(Rinv**5)
!    Colle(2) = Four*Pi1*(Rinv**5)
!    Colle(3) = Four*Pi2*(Rinv**5)
!    Colle(4) = Del1*(Rinv**5)
!    Colle(5) = Del2*(Rinv**5)
!  end if
!--------------------------------------
end if

return

end subroutine ABNone
