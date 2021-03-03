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

subroutine Init_ave(Title,iPrint,Wset,Wsum,PrOcc,PrEne,DensityBased,ThrOcc,Dummy,iDummy)
!-- Initializations and defaults.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "mxdm.fh"
#include "mxave.fh"
character(len=72), intent(out) :: Title
integer(kind=iwp), intent(out) :: iPrint, iDummy
real(kind=wp), intent(out) :: Wset(MxSets), Wsum, ThrOcc, Dummy
logical(kind=iwp), intent(out) :: PrOcc, PrEne, DensityBased
integer(kind=iwp) :: i

iPrint = 2
do i=1,MxSets
  Wset(i) = Zero
end do
Wsum = Zero
PrOcc = .true.
PrEne = .false.
DensityBased = .true.
ThrOcc = 1.0e-5_wp
Dummy = One
iDummy = 1
Title = ' Untitled job. '

return

end subroutine Init_ave
