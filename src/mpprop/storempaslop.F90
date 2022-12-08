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

subroutine StoreMpAsLop(nAtoms,ANr,nB,T,Ti,MP,lMax,EC)

use MPProp_globals, only: AtBoMltPl, Cor
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nB, lMax
integer(kind=iwp), intent(out) :: ANr(nAtoms)
real(kind=wp), intent(out) :: T(nB,nB), Ti(nB,nB), MP(nAtoms*(nAtoms+1)/2,(lMax+1)*(lMax+2)*(lMax+3)/6), EC(3,nAtoms*(nAtoms+1)/2)
integer(kind=iwp) :: i, iAt1, iAt2, iAtK, iMu, ix, iy, j, kaunter, kompost, l

!-- Let's fix the ANr.

call Get_iArray('LP_A',ANr,nAtoms)

!-- Let's fix the uber-simple T and T(-1).

call unitmat(T,nB)
Ti(:,:) = T(:,:)

!-- Let's fix the expansion centres.

kaunter = 0
do i=1,nAtoms
  do j=1,i
    kaunter = kaunter+1
    EC(:,kaunter) = Cor(:,i,j)
  end do
end do

! Let's fix the multipole moments. Unlike LoProp, MpProp has here
! included the nuclei contribution, which we have to remove pronto
! to be compatible.

iMu = 0
do l=0,lMax
  kompost = 0
  do ix=l,0,-1
    do iy=l-ix,0,-1
      kompost = kompost+1
      iMu = iMu+1
      iAtK = 0
      do iAt1=1,nAtoms
        do iAt2=1,iAt1
          iAtK = iAtK+1
          MP(iAtK,iMu) = AtBoMltPl(l)%A(kompost,iAtK)
        end do
        if (l == 0) then
          MP(iAtK,iMu) = MP(iAtK,iMu)-real(ANr(iAt1),kind=wp)
        end if
      end do
    end do
  end do
end do

return

end subroutine StoreMpAsLop
