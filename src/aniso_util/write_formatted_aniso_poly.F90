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

subroutine write_formatted_aniso_poly(filename,nss,eso,MM,MS)

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: nss
real(kind=8), intent(in) :: eso(nss)
complex(kind=8), intent(in) :: MM(3,nss,nss)
complex(kind=8), intent(in) :: MS(3,nss,nss)
character(len=180), intent(in) :: filename
! local stuff
integer :: l, i, j, LuAniso, IsFreeUnit
integer :: nstate, multiplicity
external :: IsFreeUnit

LuAniso = IsFreeUnit(81)
call molcas_open(LuAniso,filename)
nstate = 1
multiplicity = 1
write(LuAniso,'(2i10)') nstate,nss
write(LuAniso,'(5ES24.14)') (eso(i),i=1,nss)
write(LuAniso,'(30i4)') multiplicity
do l=1,3
  do i=1,nss
    write(LuAniso,'(5ES24.14)') (MM(l,i,j),j=1,nss)
  end do
end do
do l=1,3
  do i=1,nss
    write(LuAniso,'(5ES24.14)') (MS(l,i,j),j=1,nss)
  end do
end do
close(LuAniso)

return

end subroutine write_formatted_aniso_poly
