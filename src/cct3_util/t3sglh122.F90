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

subroutine t3sglh122(w,dimab,dimc,s3,d3,ns)
! this routine adds following contribution to W
! for syma=symb > symc
!
! W(ab,c) <-  + S3 _i(c) . D3 _jk(ab)
!
! w     - W matrix (I/O)
! dimab - dimension of ab index (I)
! dimc  - dimension of c index (I)
! s3    - S3 matrix (I)
! d3    - D3 matrix (I)
! ns    - signum of the contribution (+-1) (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimab, dimc, ns
real(kind=wp), intent(inout) :: w(dimab,dimc)
real(kind=wp), intent(in) :: s3(dimc), d3(dimab)
integer(kind=iwp) :: c

if (ns == 1) then
  ! phase +1

  do c=1,dimc
    w(:,c) = w(:,c)+d3(:)*s3(c)
  end do

else
  ! phase -1

  do c=1,dimc
    w(:,c) = w(:,c)-d3(:)*s3(c)
  end do

end if

return

end subroutine t3sglh122
