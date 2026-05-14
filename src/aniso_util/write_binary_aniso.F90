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

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nss, nstate
integer(kind=iwp), intent(_IN_) :: multiplicity(nstate)
real(kind=wp), intent(_IN_) :: eso(nss), esfs(nstate), angmom(3,nstate,nstate), edmom(3,nstate,nstate), amfi(3,nstate,nstate)
complex(kind=wp), intent(in) :: U(nss,nss), MM(3,nss,nss), MS(3,nss,nss), DM(3,nss,nss), HSO(nss,nss)
integer(kind=iwp) :: idisk, idum(1), l, luaniso
real(kind=wp), allocatable :: tmpI(:,:), tmpR(:,:)

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
tmpR(:,:) = real(U(:,:))
tmpI(:,:) = aimag(U(:,:))
call ddafile(luaniso,1,tmpR,nss*nss,idisk)
call ddafile(luaniso,1,tmpI,nss*nss,idisk)

! spin-orbit Hamiltonian
tmpR(:,:) = real(HSO(:,:))
tmpI(:,:) = aimag(HSO(:,:))
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
  tmpR(:,:) = real(MM(l,:,:))
  tmpI(:,:) = aimag(MM(l,:,:))
  call ddafile(luaniso,1,tmpR,nss*nss,idisk)
  call ddafile(luaniso,1,tmpI,nss*nss,idisk)
end do

! spin moment:
do l=1,3
  tmpR(:,:) = real(MS(l,:,:))
  tmpI(:,:) = aimag(MS(l,:,:))
  call ddafile(luaniso,1,tmpR,nss*nss,idisk)
  call ddafile(luaniso,1,tmpI,nss*nss,idisk)
end do

! electric dipole moment:
do l=1,3
  tmpR(:,:) = real(DM(l,:,:))
  tmpI(:,:) = aimag(DM(l,:,:))
  call ddafile(luaniso,1,tmpR,nss*nss,idisk)
  call ddafile(luaniso,1,tmpI,nss*nss,idisk)
end do

call mma_deallocate(tmpR)
call mma_deallocate(tmpI)
! close the binary $Project.aniso file
call daclos(luaniso)

return

end subroutine write_binary_aniso
