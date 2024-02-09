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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, gElectron
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: nss, nstate
integer(kind=iwp), intent(out) :: multiplicity(nstate)
real(kind=wp), intent(out) :: eso(nss), esfs(nstate), angmom(3,nstate,nstate), edmom(3,nstate,nstate), amfi(3,nstate,nstate)
complex(kind=wp), intent(out) :: U(nss,nss), MM(3,nss,nss), MS(3,nss,nss), ML(3,nss,nss), DM(3,nss,nss), HSO(nss,nss)
integer(kind=iwp) :: idisk, idum(1), l, luaniso
real(kind=wp), allocatable :: tmp(:,:,:), tmpI(:,:), tmpR(:,:)
real(kind=wp), parameter :: g_e = -gElectron

! initialize:
multiplicity(:) = 0
eso(:) = Zero
esfs(:) = Zero
angmom(:,:,:) = Zero
edmom(:,:,:) = Zero
amfi(:,:,:) = Zero
U(:,:) = cZero
MM(:,:,:) = cZero
MS(:,:,:) = cZero
ML(:,:,:) = cZero
DM(:,:,:) = cZero
HSO(:,:) = cZero
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
call ddafile(luaniso,2,tmpR,nss**2,idisk)
call ddafile(luaniso,2,tmpI,nss**2,idisk)
U(1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)

! spin-orbit Hamiltonian
call ddafile(luaniso,2,tmpR,nss**2,idisk)
call ddafile(luaniso,2,tmpI,nss**2,idisk)
HSO(1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)

call mma_allocate(tmp,3,nstate,nstate,label='tmp')
! angmom
call ddafile(luaniso,3,tmp,3*nstate*nstate,idisk)
angmom(:,1:nstate,1:nstate) = tmp(:,:,:)
! edmom
call ddafile(luaniso,3,tmp,3*nstate*nstate,idisk)
edmom(:,1:nstate,1:nstate) = tmp(:,:,:)
! amfi
call ddafile(luaniso,3,tmp,3*nstate*nstate,idisk)
amfi(:,1:nstate,1:nstate) = tmp(:,:,:)

! magnetic moment
do l=1,3
  call ddafile(luaniso,2,tmpR,nss**2,idisk)
  call ddafile(luaniso,2,tmpI,nss**2,idisk)
  MM(l,1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
end do

! spin moment
do l=1,3
  call ddafile(luaniso,2,tmpR,nss**2,idisk)
  call ddafile(luaniso,2,tmpI,nss**2,idisk)
  MS(l,1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
end do
! generate magnetic moment on the fly:
ML(:,1:nss,1:nss) = -MM(:,1:nss,1:nss)-g_e*MS(:,1:nss,1:nss)

! electric dipole moment
do l=1,3
  call ddafile(luaniso,2,tmpR,nss**2,idisk)
  call ddafile(luaniso,2,tmpI,nss**2,idisk)
  DM(l,1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
end do

call mma_deallocate(tmpR)
call mma_deallocate(tmpI)
call daclos(luaniso)

return

end subroutine read_binary_aniso
