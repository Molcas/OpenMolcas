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

subroutine write_binary_aniso(nss,nstate,multiplicity,eso,esfs,U,MM,MS,DM,ANGMOM,EDMOM,AMFI,HSO)

implicit none
integer, parameter :: wp = kind(0.d0)
integer :: nss, nstate
integer :: multiplicity(nstate)
real(kind=8) :: eso(nss), esfs(nstate)
real(kind=8) :: angmom(3,nstate,nstate)
real(kind=8) :: edmom(3,nstate,nstate)
real(kind=8) :: amfi(3,nstate,nstate)
complex(kind=8) :: U(nss,nss), HSO(nss,nss)
complex(kind=8) :: MM(3,nss,nss), MS(3,nss,nss), DM(3,nss,nss)
! local variables:
#include "stdalloc.fh"
integer :: i, j, l
integer :: luaniso, idisk, idum(1)
real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:)

LUANISO = 8
! open the binary $Project.aniso file
call DANAME(LUANISO,'POLYFILE')
! writing data to it:
idisk = 0
idum(1) = nstate
call idafile(luaniso,1,idum,1,idisk)
idum(1) = nss
call idafile(luaniso,1,idum,1,idisk)
! spin multiplicity ofthe spin free states:
call idafile(luaniso,1,multiplicity,nstate,idisk)
! spin-orbit energies:
call ddafile(luaniso,1,eso,nss,idisk)
! spin-free energies:
call ddafile(luaniso,1,esfs,nstate,idisk)

call mma_allocate(tmpR,nss,nss,'tmpR')
call mma_allocate(tmpI,nss,nss,'tmpI')
! spin-orbit mixing coefficients:
tmpR = 0.0_wp
tmpI = 0.0_wp
do i=1,nss
  do j=1,nss
    tmpR(i,j) = dble(U(i,j))
    tmpI(i,j) = aimag(U(i,j))
  end do
end do
call ddafile(luaniso,1,tmpR,nss*nss,idisk)
call ddafile(luaniso,1,tmpI,nss*nss,idisk)

! spin-orbit Hamiltonian
tmpR = 0.0_wp
tmpI = 0.0_wp
do i=1,nss
  do j=1,nss
    tmpR(i,j) = dble(HSO(i,j))
    tmpI(i,j) = aimag(HSO(i,j))
  end do
end do
call ddafile(luaniso,1,tmpR,nss*nss,idisk)
call ddafile(luaniso,1,tmpI,nss*nss,idisk)

! angmom
call ddafile(luaniso,1,angmom,3*nstate*nstate,idisk)
! edmom
call ddafile(luaniso,1,edmom,3*nstate*nstate,idisk)
! amfi
call ddafile(luaniso,1,amfi,3*nstate*nstate,idisk)

! magnetic moment:
do l=1,3
  tmpR = 0.0_wp
  tmpI = 0.0_wp
  do i=1,nss
    do j=1,nss
      tmpR(i,j) = dble(MM(l,i,j))
      tmpI(i,j) = aimag(MM(l,i,j))
    end do
  end do
  call ddafile(luaniso,1,tmpR,nss*nss,idisk)
  call ddafile(luaniso,1,tmpI,nss*nss,idisk)
end do

! spin moment:
do l=1,3
  tmpR = 0.0_wp
  tmpI = 0.0_wp
  do i=1,nss
    do j=1,nss
      tmpR(i,j) = dble(MS(l,i,j))
      tmpI(i,j) = aimag(MS(l,i,j))
    end do
  end do
  call ddafile(luaniso,1,tmpR,nss*nss,idisk)
  call ddafile(luaniso,1,tmpI,nss*nss,idisk)
end do

! electric dipole moment:
do l=1,3
  tmpR = 0.0_wp
  tmpI = 0.0_wp
  do i=1,nss
    do j=1,nss
      tmpR(i,j) = dble(DM(l,i,j))
      tmpI(i,j) = aimag(DM(l,i,j))
    end do
  end do
  call ddafile(luaniso,1,tmpR,nss*nss,idisk)
  call ddafile(luaniso,1,tmpI,nss*nss,idisk)
end do

call mma_deallocate(tmpR)
call mma_deallocate(tmpI)
! close the binary $Project.aniso file
call daclos(luaniso)

return

end subroutine write_binary_aniso
