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

subroutine fokupdate2(foka,symp,i,vint,ndimv1,ndimv2,ndimv3)
! this routine realizes update
! foka(p,q) = foka(p,q) - <ip|qi>
!
! N.B. integrals are of type <symi, symp| symp, symi>
!
! foka    - packed Fokaa matrix (I,O)
! symp    - irrep or p (and also q) index (I)
! i       - value of i, (I)
! vint    - array of integrals <ip|iq> for given i (I)
! ndimv1  - first dimension (norb(symp)) (I)
! ndimv2  - second dimension (norb(symi)) (I)
! ndimv3  - third dimension (norb(symp)) (I)

use ccsort_global, only: NORB
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: foka(*)
integer(kind=iwp), intent(in) :: symp, i, ndimv1, ndimv2, ndimv3
real(kind=wp), intent(in) :: vint(ndimv1,ndimv2,ndimv3)
integer(kind=iwp) :: nhelp1, nhelp2, p, pq, q

! calculate shift

nhelp1 = 0
do nhelp2=1,symp-1
  nhelp1 = nhelp1+(norb(nhelp2)**2+norb(nhelp2))/2
end do

! add integral

pq = nhelp1
do p=1,norb(symp)
  do q=1,p
    pq = pq+1
    foka(pq) = foka(pq)-vint(p,q,i)
  end do
end do

return

end subroutine fokupdate2
