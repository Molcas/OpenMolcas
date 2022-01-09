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

subroutine Reorder_GW(A,B,k,l,n,m)

implicit real*8(a-h,o-z)
real*8 A(k,l,n,m), B(k,n,l,m)

do ik=1,k
  do il=1,l
    do in=1,n
      do im=1,m

        B(ik,in,il,im) = A(ik,il,in,im)

      end do
    end do
  end do
end do

return

end subroutine Reorder_GW
