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

subroutine Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ipTmp,nr_of_Densities)

use k2_arrays, only: ipOffD
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ijS, nr_of_Densities
integer(kind=iwp), intent(out) :: ipDij, ipDSij, mDCRij, ipDDij
integer(kind=iwp), intent(inout) :: ipTmp
integer(kind=iwp) :: nDij

ipDij = ipOffD(1,ijS)
mDCRij = ipOffD(2,ijS)
nDij = ipOffD(3,ijS)
if (nr_of_Densities == 2) then
  ipDSij = ipOffD(4,ijS)
else
  ipDSij = 1
end if

if (mDCRij*nDij /= 0) then
  ipDDij = ipTmp
  ipTmp = ipTmp+nDij*mDCRij
else
  ipDDij = 1
end if

return

end subroutine Dens_Info
