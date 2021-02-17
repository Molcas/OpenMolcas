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

subroutine TransNow(V,S)

use FFPT_Global, only: nBas, TranCoo, ComStk, ComVal
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: V(nBas(1)*(nBas(1)+1)/2)
real(kind=wp), intent(in) :: S(nBas(1)*(nBas(1)+1)/2)
integer(kind=iwp) :: i, iXYZ, j, kaunter
real(kind=wp) :: Trana

!-- Be informative.

if (ComStk(2,1,1,0)) then
  write(u6,*)
  write(u6,*) '    The DIPO perturbation is translated.'
else
  write(u6,*)
  write(u6,*) '    No translation of the perturbation (only implemented for DIPO.'
  write(u6,*) 'OBSERVE! Your result can be origo dependent!'
end if

!-- Translate the perturbation. Why, oh why, is the sign of Trana
!   the opposite of what it should be according to the formula?
!   Answer: Our charge is positive in the MLTPL integrals, hence
!   the field strength should be of opposite sign, see ptdipo.f
!   for example. Hence the minus in the translation formula becomes
!   a plus. Oh yeah!

if ((.not. ComStk(2,1,1,1)) .and. (.not. ComStk(2,1,1,2)) .and. (.not. ComStk(2,1,1,3))) then
  write(u6,*)
  write(u6,*) 'A strange error has occured. ComStk modified?'
  call Abend()
end if
kaunter = 1
!write(u6,*) 'TranCoo',(TranCoo(i),i=1,3)
do i=1,nBas(1)
  do j=1,i
    do iXYZ=1,3
      if (ComStk(2,1,1,iXYZ)) then
        Trana = ComVal(2,1,1,iXYZ)*TranCoo(iXYZ)*S(kaunter)
        !-----Sign of Trana? See source code comment above.
        V(kaunter) = V(kaunter)+Trana
      end if
    end do
    kaunter = kaunter+1
  end do
end do

return

end subroutine TransNow
