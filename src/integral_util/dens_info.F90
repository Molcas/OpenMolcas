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

subroutine Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDij2,ipDDij2)

use k2_arrays, only: ipOffD, ipOffDA
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ijS, nr_of_Densities, nMethod
integer(kind=iwp), intent(out) :: ipDij, ipDSij, mDCRij, ipDDij, ipDij2, ipDDij2
integer(kind=iwp), intent(inout) :: ipTmp, ipTmp2
integer(kind=iwp) :: nDij
integer(kind=iwp), parameter :: SCF = 1, RASSCF = 2

ipDij = ipOffD(1,ijS)
mDCRij = ipOffD(2,ijS)
nDij = ipOffD(3,ijS)
if (nMethod == RASSCF) ipDij2 = ipOffDA(1,ijS)

if (nr_of_Densities == 2) then
  ipDSij = ipOffD(4,ijS)
else
  ipDSij = 1
end if

if (mDCRij*nDij /= 0) then
  ipDDij = ipTmp
  ipTmp = ipTmp+nDij*mDCRij
  if (nMethod == RASSCF) then
    ipDDij2 = ipTmp2
    ipTmp2 = ipTmp2+nDij*mDCRij
  end if

else
  ipDDij = 1
end if

end subroutine Dens_Info
