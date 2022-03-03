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

subroutine unpackk_ic_3(i,vint,ndimvi,Vic)
! this routine vint(j,k,l) = <i,j|k,l>
! for given i from incore (reduced) expanded block Vic
! ie. symp=symq=symr=syms
!
! i      - value of pivot index (I)
! vint   - array of integrals (O)
! ndimvi - (norb(symi)) (I)
! Vic    - incore expanded block of integrals (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: i, ndimvi
real(kind=wp), intent(out) :: vint(ndimvi,ndimvi,ndimvi)
real(kind=wp), intent(in) :: Vic((ndimvi*(ndimvi+1)/2)*(ndimvi*(ndimvi+1)/2+1)/2)
integer(kind=iwp) :: ik, ikjl, j, jl, k, l

do k=1,ndimvi

  ! def ik
  if (i >= k) then
    ik = i*(i-1)/2+k
  else
    ik = k*(k-1)/2+i
  end if

  jl = 0
  do j=1,ndimvi
    do l=1,j

      ! def jl
      jl = jl+1

      ! def ikjl
      if (ik >= jl) then
        ikjl = ik*(ik-1)/2+jl
      else
        ikjl = jl*(jl-1)/2+ik
      end if

      vint(j,k,l) = Vic(ikjl)
      vint(l,k,j) = Vic(ikjl)

    end do
  end do
end do

return

end subroutine unpackk_ic_3
