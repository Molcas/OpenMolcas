!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!***********************************************************************

subroutine P2_contraction(D1MO,P2MO)

use Index_Functions, only: iTri
use rasscf_global, only: NAC, NACPAR, NACPR2
use Constants, only: One, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: d1mo(NACPAR)
real(kind=wp), intent(out) :: p2mo(NACPR2)
integer(kind=iwp) :: i, ij, ijkl, j, k, kl, l, lmax
real(kind=wp) :: fact

ijkl = 0
do i=1,nac
  do j=1,i
    ij = iTri(i,j)
    do k=1,i
      if (i == k) then
        lmax = j
      else
        lmax = k
      end if
      do l=1,lmax
        kl = iTri(k,l)
        ijkl = ijkl+1
        fact = One
        if (k == l) fact = Half
        p2MO(ijkl) = fact*D1MO(ij)*D1MO(kl)
      end do
    end do
  end do
end do

end subroutine P2_contraction
