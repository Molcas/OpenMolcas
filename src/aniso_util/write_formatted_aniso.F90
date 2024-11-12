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

subroutine write_formatted_aniso(nss,nstate,multiplicity,eso,esfs,U,MM,MS,DM,angmom,edmom,amfi,HSO)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nss, nstate, multiplicity(nstate)
real(kind=wp), intent(in) :: eso(nss), esfs(nstate), angmom(3,nstate,nstate), edmom(3,nstate,nstate), amfi(3,nstate,nstate)
complex(kind=wp), intent(in) :: U(nss,nss), MM(3,nss,nss), MS(3,nss,nss), DM(3,nss,nss), HSO(nss,nss)
integer(kind=iwp) :: i, j, l, LuAniso
integer(kind=iwp), external :: IsFreeUnit

LuAniso = IsFreeUnit(81)
call molcas_open(LuAniso,'ANISOINPUT')
write(LuAniso,'(2i10)') nstate,nss
write(LuAniso,'(5ES24.14)') (eso(i),i=1,nss)
write(LuAniso,'(30i4)') (multiplicity(i),i=1,nstate)
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
! add data at the end, so that we do not break the functionality
! with the present format:
write(LuAniso,'(5ES24.14)') (esfs(i),i=1,nstate)
do i=1,nss
  write(LuAniso,'(5ES24.14)') (U(i,j),j=1,nss)
end do
! angmom
do l=1,3
  do i=1,nstate
    write(LuAniso,'(5ES24.14)') (angmom(l,i,j),j=1,nstate)
  end do
end do

! DMmom
do l=1,3
  do i=1,nss
    write(LuAniso,'(5ES24.14)') (DM(l,i,j),j=1,nss)
  end do
end do

! edmom
do l=1,3
  do i=1,nstate
    write(LuAniso,'(5ES24.14)') (edmom(l,i,j),j=1,nstate)
  end do
end do

! amfi
do l=1,3
  do i=1,nstate
    write(LuAniso,'(5ES24.14)') (amfi(l,i,j),j=1,nstate)
  end do
end do

! HSO matrix
do i=1,nss
  write(LuAniso,'(5ES24.14)') (HSO(i,j),j=1,nss)
end do

close(LuAniso)

return

end subroutine write_formatted_aniso
