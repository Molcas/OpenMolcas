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
subroutine grow_t2_fblocked1(t2,tmp,dima,dimb,nv,no,lasta,lastb,length1,length2,grpa,grpb)

implicit none
integer a, b, dima, dimb, nv, no, i, j
!integer ij
integer lasta, lastb
integer grpa, grpb
integer length1, length2
!mp real*8 t2(1:nv,1:nv,1:no,1:no)
real*8 t2(1:length1,1:length2,1:no,1:no)
real*8 tmp(1:dima,1:dimb,1:no,1:no)

!mp write(6,*) 'lasta+dima, length1 ',lasta+dima,length1
!mp write(6,*) 'lastb+dimb, length2 ',lastb+dimb,length2

!mp write(6,'(A,2(i10,x),i3)') 'length1, length2, no ',length1,length2,no
!mp write(6,'(A,4(i3,x))') 'dima, dimb, lasta, lastb',dima,dimb,lasta,lastb

!mp ij = 0
!mp do i=2,no
!mp   do j=1,i-1
!mp     ij = ij+1
do i=1,no
  do j=1,no
    do b=1,dimb
      do a=1,dima

        !mp t2(lasta+a,lastb+b,ij) = tmp(a,b,i,j)+(-1.0d0*tmp(a,b,j,i))
        t2(lasta+a,lastb+b,i,j) = tmp(a,b,i,j)

      end do
    end do
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(nv)
  call Unused_integer(grpa)
  call Unused_integer(grpb)
end if

end subroutine grow_t2_fblocked1
