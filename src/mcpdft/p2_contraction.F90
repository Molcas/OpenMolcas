!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!***********************************************************************

Subroutine P2_contraction(D1MO,P2MO)
  use definitions,only:iwp,wp
  use constants,only:one,half
  use rasscf_global,only:NAC

  implicit none

  real(kind=wp),intent(in) :: d1mo(*)
  real(kind=wp),intent(out) :: p2mo(*)

  integer(kind=iwp) :: i,j,k,l,ij,kl,ijkl,lmax
  real(kind=wp) :: fact

  ijkl = 0
  do i = 1,nac
    do j = 1,i
      ij = iTrii(i,j)
      do k = 1,i
        if(i == k) then
          lmax = j
        else
          lmax = k
        endif
        do l = 1,lmax
          kl = iTrii(k,l)
          ijkl = ijkl+1
          fact = one
          if(k == l) fact = half
          p2MO(ijkl) = fact*D1MO(ij)*D1MO(kl)
        enddo
      enddo
    enddo
  enddo
contains
  integer function iTrii(i,j)
    integer,intent(in) :: i,j
    itrii = Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
  endfunction
endSubroutine P2_contraction
