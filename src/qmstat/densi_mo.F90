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

! Here we construct the density matrix given the orbital coefficients.
subroutine DENSI_MO(DENS,ORBCO,IS,IA,NBAS,IDM)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IS, IA, NBAS, IDM
real(kind=wp), intent(out) :: DENS(nTri_Elem(NBAS))
real(kind=wp), intent(in) :: ORBCO(IDM,*)
integer(kind=iwp) :: I, IJ, J, K

DENS(:) = Zero
do I=IS,IS+IA-1
  IJ = 0
  do J=1,NBAS
    do K=1,J
      IJ = IJ+1
      DENS(IJ) = DENS(IJ)+Four*ORBCO(J,I)*ORBCO(K,I)
    end do
    DENS(IJ) = DENS(IJ)-ORBCO(J,I)*ORBCO(J,I)*Two
  end do
end do

return

end subroutine DENSI_MO
