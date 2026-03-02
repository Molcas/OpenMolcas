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
! Copyright (C) 2026, Lila Zapp                                        *
!***********************************************************************

subroutine upper_triag2vec(squaremat,matdim,vec,vecdim)
use Definitions, only: u6,wp,iwp

implicit none
integer(kind=iwp),intent(in) :: matdim,vecdim
real(kind=wp),intent(in) :: squaremat(matdim,matdim)
real(kind=wp),intent(out) :: vec(vecdim)
integer(kind=iwp) :: i,j,listindex

! putting the upper triagonal elements & diagonal elements into a list
listindex=0
do i=1,matdim-1
    do j=i+1,matdim
        listindex=listindex+1
        if (.false.) then
            write(u6,"(A,I5,A,I5,A,I5,A,F8.3)") "i=",i ,"j= ",j,"listindex=",listindex,"mat(i,j)=",squaremat(i,j)
        end if
        vec(listindex)=squaremat(i,j)
    end do
end do

if (.true.) then
    write(u6,*) "In upper_triag2vec:"
    call RecPrt("NxN Matrix",' ',squaremat,matdim,matdim)
    call RecPrt("matrix as vector of upper triagonal values:",' ',vec,listindex,1)
end if

end subroutine upper_triag2vec
