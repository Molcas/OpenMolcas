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

subroutine vec2upper_triag(squaremat,matdim,vec,vecdim,antisymmetric)
use Definitions, only: u6,wp,iwp

implicit none
integer(kind=iwp),intent(in) :: matdim,vecdim
real(kind=wp),intent(out) :: squaremat(matdim,matdim) !antisymmetric or symmetric
real(kind=wp),intent(in) :: vec(vecdim)
integer(kind=iwp) :: i,j,listindex
logical, intent(in) :: antisymmetric

! putting data stored as vector back into anti/symmetric matrix format
listindex=0
do i=1,matdim
    do j=i,matdim
        listindex=listindex+1
        if (.false.) then
            write(u6,"(A,I5,A,I5,A,I5,A,F8.3)") "i=",i ,"j= ",j,"listindex=",listindex,"mat(i,j)=",squaremat(i,j)
        end if

        squaremat(i,j)=vec(listindex)

        if (antisymmetric) then
            squaremat(j,i)=-vec(listindex)
        else
            squaremat(j,i)=vec(listindex)
        end if
    end do
end do

if (.false.) then
    write(u6,*) "In vec2upper_triag:  antisymmetric = ",antisymmetric
    call RecPrt("matrix as vector of upper triagonal values:",' ',vec,listindex,1)
    call RecPrt("NxN Matrix",' ',squaremat,matdim,matdim)
end if

end subroutine vec2upper_triag
