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

subroutine G_Nrm(nInter,GNrm,Iter,Grad,mIntEff)

use Slapaf_Info, only: Degen, Gx
use Constants, only: Zero
use Definitions, only: wp, iwp
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nInter, Iter
real(kind=wp), intent(inout) :: GNrm(Iter)
real(kind=wp), intent(in) :: Grad(nInter,Iter)
integer(kind=iwp), intent(out) :: mIntEff
integer(kind=iwp) :: i, j
real(kind=wp) :: Fabs

! Compute the norm of the cartesian force vector.
!
! |dE/dx|=Sqrt(dE/dx|u|dE/dx)

Fabs = Zero
do i=1,size(Gx,2)
  do j=1,3
    Fabs = Fabs+Degen(j,i)*Gx(j,i,Iter)**2
  end do
end do
Fabs = sqrt(Fabs)
#ifdef _DEBUGPRINT_
write(u6,42) Fabs
42 format(/,' Norm of the force vector',F20.15)
#endif
GNrm(iter) = Fabs

! Write out the internal force vector.

mIntEff = 0
do i=1,nInter
  if (abs(Grad(i,Iter)) > 1.0e-6_wp) mIntEff = mIntEff+1
end do
if (mIntEff == 0) mIntEff = 1
#ifdef _DEBUGPRINT_
write(u6,*) ' mIntEff=',mIntEff
#endif

return

end subroutine G_Nrm
