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

function nbfshl(iSkal,irp)
!***********************************************************************
!                                                                      *
!  Integer Function NbfShl    returns # of bf for shell,symmetry       *
!                                                                      *
!***********************************************************************

use iSD_data, only: iSD
use SOAO_Info, only: iAOtSO
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nbfshl
integer(kind=iwp), intent(in) :: iSkal, irp
integer(kind=iwp) :: i, iAO, iCmp

! returns number of basis functions for given shell and symmetry

nbfshl = 0
iAO = iSD(7,iSkal)
iCmp = iSD(2,iSkal)
! loop over components of shell...
do i=1,iCmp
  if (iAOtSO(iAO+i,irp) > 0) nbfshl = nbfshl+iSD(3,iSkal)
end do

return

end function nbfshl
