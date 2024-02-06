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

subroutine read_binary_aniso(nss,nstate,multiplicity,eso,esfs,U,MM,MS,ML,DM,angmom,edmom,amfi,HSO)

use Constants, only: Zero, cZero, gElectron
use Definitions, only: wp

implicit none
integer :: nss, nstate
integer :: multiplicity(nstate)
real(kind=8) :: eso(nss), esfs(nstate)
real(kind=8) :: angmom(3,nstate,nstate)
real(kind=8) :: edmom(3,nstate,nstate)
real(kind=8) :: amfi(3,nstate,nstate)
complex(kind=8) :: MM(3,nss,nss), MS(3,nss,nss), ML(3,nss,nss)
complex(kind=8) :: DM(3,nss,nss)
complex(kind=8) :: U(nss,nss)
complex(kind=8) :: HSO(nss,nss)
! local variables:
#include "stdalloc.fh"
integer :: i, j, l
integer :: luaniso, idisk, idum(1)
real(kind=8), parameter :: g_e = 2.0023193043718_wp !IFG -gElectro
real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:)

! initialize:
multiplicity = 0
eso(:) = Zero
esfs(:) = Zero
angmom(:,:,:) = Zero
edmom(:,:,:) = Zero
amfi(:,:,:) = Zero
U(:,:) = cZero
MS(:,:,:) = cZero
ML(:,:,:) = cZero
MM(:,:,:) = cZero
DM(:,:,:) = cZero
HSO(:,:) = cZero
luaniso = 8
idisk = 0
! get the information from binary "$project.aniso" file:
call daname(luaniso,'POLYFILE')
call idafile(luaniso,2,idum,1,idisk)
nstate = idum(1)
call idafile(luaniso,2,idum,1,idisk)
nss = idum(1)

call mma_allocate(tmpR,nss,nss,'tmpR')
call mma_allocate(tmpI,nss,nss,'tmpI')

call idafile(luaniso,2,multiplicity,nstate,idisk)
call ddafile(luaniso,2,eso,nss,idisk)
call ddafile(luaniso,2,esfs,nstate,idisk)

! spin-orbit mixing coefficients:
tmpR(:,:) = Zero
tmpI(:,:) = Zero
call ddafile(luaniso,2,tmpR,nss**2,idisk)
call ddafile(luaniso,2,tmpI,nss**2,idisk)
do i=1,nss
  do j=1,nss
    U(i,j) = cmplx(tmpR(i,j),tmpI(i,j),wp)
  end do
end do

! spin-orbit Hamiltonian
tmpR(:,:) = Zero
tmpI(:,:) = Zero
call ddafile(luaniso,2,tmpR,nss**2,idisk)
call ddafile(luaniso,2,tmpI,nss**2,idisk)
do i=1,nss
  do j=1,nss
    HSO(i,j) = cmplx(tmpR(i,j),tmpI(i,j),wp)
  end do
end do

! angmom
call ddafile(luaniso,3,angmom,3*nstate*nstate,idisk)
! edmom
call ddafile(luaniso,3,edmom,3*nstate*nstate,idisk)
! amfi
call ddafile(luaniso,3,amfi,3*nstate*nstate,idisk)

! magnetic moment
do l=1,3
  tmpR(:,:) = Zero
  tmpI(:,:) = Zero
  call ddafile(luaniso,2,tmpR,nss**2,idisk)
  call ddafile(luaniso,2,tmpI,nss**2,idisk)
  do i=1,nss
    do j=1,nss
      MM(l,i,j) = cmplx(tmpR(i,j),tmpI(i,j),wp)
    end do
  end do
end do

! spin moment
do l=1,3
  tmpR(:,:) = Zero
  tmpI(:,:) = Zero
  call ddafile(luaniso,2,tmpR,nss**2,idisk)
  call ddafile(luaniso,2,tmpI,nss**2,idisk)
  do i=1,nss
    do j=1,nss
      MS(l,i,j) = cmplx(tmpR(i,j),tmpI(i,j),wp)
    end do
  end do
end do
! generate magnetic moment on the fly:
do l=1,3
  do i=1,nss
    do j=1,nss
      ML(l,i,j) = -MM(l,i,j)-g_e*MS(l,i,j)
    end do
  end do
end do

! electric dipole moment
do l=1,3
  tmpR(:,:) = Zero
  tmpI(:,:) = Zero
  call ddafile(luaniso,2,tmpR,nss**2,idisk)
  call ddafile(luaniso,2,tmpI,nss**2,idisk)
  do i=1,nss
    do j=1,nss
      DM(l,i,j) = cmplx(tmpR(i,j),tmpI(i,j),wp)
    end do
  end do
end do

call mma_deallocate(tmpR)
call mma_deallocate(tmpI)
call daclos(luaniso)

return

end subroutine read_binary_aniso
