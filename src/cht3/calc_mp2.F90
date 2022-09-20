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

subroutine calc_MP2(w,e,no,nv)
! this is primitive checking routine to calculate 2nd order energy

implicit none
integer i, j, a, b, no, nv
real*8 e(1:(no+nv)), w(1:nv,1:no,1:nv,1:no)
real*8 e2, integral, denom

e2 = 0.0d0

do j=1,no
  do i=1,no
    do b=1,nv
      do a=1,nv

        denom = e(no+a)+e(no+b)-e(i)-e(j)
        !mp write(6,'(4(i3,2x),A,3(f17.10,2x))') a,i,b,j,'w1, w2, denom ',w(a,i,b,j),w(a,j,b,i),denom

        integral = (-1.0d0)*w(a,i,b,j)*(2.0d0*w(a,i,b,j)+(-1.0d0)*w(a,j,b,i))

        !write(6,*) integral

        e2 = e2+(integral/denom)

      end do
    end do
  end do
end do

write(6,*) 'Druhy rad je asi = ',e2

return

end subroutine calc_MP2
