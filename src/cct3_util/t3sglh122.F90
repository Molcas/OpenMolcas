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

subroutine t3sglh122(w,dima,dimab,dimc,s3,d3,ns)
! this routine adds following contribution to W
! for syma=symb > symc
!
! W(ab,c) <-  + S3 _i(c) . D3 _jk(ab)
!
! w     - W matrix (I/O)
! dima  - dimension of a (b) index (I)
! dimab - dimension of ab index (I)
! dimc  - dimension of c index (I)
! s3    - S3 matrix (I)
! d3    - D3 matrix (I)
! ns    - signum of the contribution (+-1) (I)

integer dima, dimab, dimc, ns
real*8 w(1:dimab,1:dimc)
real*8 s3(1:dimc)
real*8 d3(1:dimab)
! help variables
integer c, ab
real*8 s

if (ns == 1) then
  ! phase +1

  do c=1,dimc
    s = s3(c)
    do ab=1,dimab
      w(ab,c) = w(ab,c)+d3(ab)*s
    end do
  end do

else
  ! phase -1

  do c=1,dimc
    s = s3(c)
    do ab=1,dimab
      w(ab,c) = w(ab,c)-d3(ab)*s
    end do
  end do

end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(dima)

end subroutine t3sglh122
