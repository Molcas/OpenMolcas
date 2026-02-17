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

subroutine DBLOCK(D)
! RASSCF program version IBM-3090: SX section
!
! Purpose: To symmetry-block a full matrix d over the active orbital
!          the result is overlaid on the input matrix.
!
! ********** IBM-3090 release 88 10 10 **********

use Index_Functions, only: iTri, nTri_Elem
use general_data, only: NASH, NSYM
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: D(*)
integer(kind=iwp) :: IA, ISYM, ITU, NA, NAT, NAU, NTU

IA = NASH(1)
NTU = nTri_Elem(IA)
do ISYM=2,NSYM
  NA = NASH(ISYM)
  if (NA == 0) cycle
  do NAT=1,NA
    do NAU=1,NAT
      NTU = NTU+1
      ITU = iTri(NAT+IA,NAU+IA)
      D(NTU) = D(ITU)
    end do
  end do
  IA = IA+NASH(ISYM)
end do

end subroutine DBLOCK
