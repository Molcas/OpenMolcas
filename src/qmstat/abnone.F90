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

implicit real*8(a-h,o-z)
parameter(MxMltp=2)
dimension dMulA((MxMltp+1)*(MxMltp+2)/2)
dimension Colle(3)

do i=1,3
  Colle(i) = 0.0d0
end do

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
    Colle(1) = 2.0d0*Sigma*(Rinv**3)
    Colle(2) = Pi1*(Rinv**3)
    Colle(3) = Pi2*(Rinv**3)
  else if (iLB == 2) then
    d3 = sqrt(3.0d0)
    Sigma = dMulA(3)
    Pi1 = dMulA(1)
    Pi2 = dMulA(2)
    Colle(1) = 3.0d0*Sigma*(Rinv**4)
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
!    d3 = sqrt(3.0d0)
!    Sigma = dMulA(3)
!    Pi1 = dMulA(2)
!    Pi2 = dMulA(4)
!    Colle(1) = 3.0d0*Sigma*(Rinv**4)
!    Colle(2) = d3*Pi1*(Rinv**4)
!    Colle(3) = d3*Pi2*(Rinv**4)
!  else if (iLB == 2) then
!    Sigma = dMulA(3)
!    Pi1 = dMulA(2)
!    Pi2 = dMulA(4)
!    ! Jose. Remember dMulB(1)=sqrt(3)*xy
!    !            and dMulB(5)=0.5*sqrt(3)*(x2-y2)
!    Del1 = dMulA(1)
!    Del2 = dMulA(5)
!    Colle(1) = 6.0d0*Sigma*(Rinv**5)
!    Colle(2) = 4.0d0*Pi1*(Rinv**5)
!    Colle(3) = 4.0d0*Pi2*(Rinv**5)
!    Colle(4) = Del1*(Rinv**5)
!    Colle(5) = Del2*(Rinv**5)
!  end if
!--------------------------------------
end if

return

end subroutine ABNone
