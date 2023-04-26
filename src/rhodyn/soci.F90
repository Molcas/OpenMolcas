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
! Copyright (C) 2021-2023, Vladislav Kochetov                          *
!***********************************************************************

subroutine soci()
! Purpose: Read in the SO Coefficient from MOLCAS output, The prepare
! for the density matrix which used in propagation
!
!***********************************************************************
!  VR       : eigenvector matrix of total Hamiltonian HTOT_CSF
!  SO_CI    : SO coefficients
!  CH_SO    : H(RASSCF)_diag+H(SO), eigenvalue of the spin-free states
!             plus spin-orbit coupling
!  SO_CI2   : SO_CI^T*SO_CI
!  SO_eig   : SO_CI^T*CH_SO*SO_CI, diagonalize the Hamiltonian
!  Hfull    : HTOT_CSF hamiltonian diagonalized

use rhodyn_data, only: CSF2SO, E_SO, HTOT_CSF, ipglob, lrootstot, nconftot, prep_csfsoi, prep_csfsor, sint, threshold, SO_CI, &
                       U_CI, flag_so
use linalg_mod, only: mult
use mh5, only: mh5_put_dset
use rhodyn_utils, only: dashes, transform
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, INFO, j, LWORK
real(kind=wp), allocatable :: RWORK(:), W(:)
complex(kind=wp), allocatable :: Hdiag(:,:), Hfull(:,:), Hfull2(:,:), SO_CI2(:,:), WORK(:)

call mma_allocate(Hfull,nconftot,nconftot)
call mma_allocate(Hfull2,nconftot,nconftot)
call mma_allocate(Hdiag,nconftot,nconftot)
call mma_allocate(SO_CI2,lrootstot,lrootstot)
call mma_allocate(Work,1)
call mma_allocate(RWORK,3*nconftot-2,label='RWORK')
call mma_allocate(W,nconftot,label='W')

Hfull(:,:) = Zero
Hfull(:,:) = HTOT_CSF
if (ipglob > 4) then
  call dashes()
  write(u6,*) 'Printout Hfull matrix'
  call dashes()
  write(u6,sint) 'Number of the CSFs:',nconftot
  call dashes()
  do i=1,10
    write(u6,*) (Hfull(i,j),j=1,10)
  end do
end if

if (ipglob > 2) then
  call dashes()
  write(u6,*) 'diagonalize the full Hamiltonian HTOT_CSF and ascending sorted eigenvalues'
  call dashes()
end if
LWORK = -1

call zheev('V','L',nconftot,Hfull,nconftot,W,WORK,LWORK,RWORK,INFO)

if (INFO == 0) then
  LWORK = max(1,int(WORK(1))+1)
  call mma_deallocate(WORK)
  call mma_allocate(WORK,LWORK)
else
  call abend()
end if

call zheev('V','L',nconftot,Hfull,nconftot,W,WORK,LWORK,RWORK,INFO)

if (ipglob > 4) then
  call dashes()
  write(u6,*) 'printout the eigenvectors after diagnolize'
  call dashes()
  do i=1,nconftot
    write(u6,*) (Hfull(i,j),j=1,nconftot)
  end do
  call dashes()
  write(u6,*) 'printout the eigenvalues'
  call dashes()
  do i=1,nconftot
    write(u6,*) W(i)
  end do
end if

! eigenvalue of SO states E_SO
E_SO(:) = W(1:lrootstot)

call mult(Hfull,Hfull,Hfull2,.true.,.false.)

if (ipglob > 4) then
  call dashes()
  write(u6,*) 'Printout the overlap of eigenvector'
  call dashes()
  do i=1,nconftot
    write(u6,*) (Hfull2(i,j),j=1,nconftot)
  end do
end if

if (ipglob > 3) then
  call dashes()
  write(u6,*) 'Check whether the eigenvectors of the full Hamiltonian are orthonormalized'
  call dashes()
  do i=1,nconftot
    do j=1,nconftot
      if ((i == j) .and. ((real(Hfull2(I,J))-One >= threshold) .or. (aimag(Hfull2(I,J)) >= threshold))) then
        write(u6,*) 'ERROR INFO!!!: Hfull is not orthonomalized:',I,J,Hfull2(I,J)
      else if ((i /= j) .and. ((real(Hfull2(i,j)) >= threshold) .or. (aimag(Hfull2(i,j)) >= threshold))) then
        write(u6,*) 'ERROR INFO!!!: Hfull is not orhtonomalized:',i,j,Hfull2(i,j)
      end if
    end do
  end do
  call dashes()
  write(u6,*) 'If there is no error info printout'
  write(u6,*) 'Hfull coeffs are orthonormalized'
  call dashes()
end if

Hdiag(:,:) = Zero

call transform(HTOT_CSF,Hfull,Hdiag)

if (ipglob > 4) then
  call dashes()
  write(u6,*) 'printout the Hdiag'
  call dashes()
  do i=1,nconftot
    write(u6,*) (Hdiag(i,j),j=1,nconftot)
  end do
  call dashes()
end if

if (ipglob > 4) then
  call dashes(72)
  write(u6,*) 'Printout the SO Coefficient matrix'
  call dashes(72)
  do i=1,lrootstot
    write(u6,*) (SO_CI(i,j),j=1,lrootstot)
  end do
end if

call mult(SO_CI,SO_CI,SO_CI2,.true.,.false.)

! check whether SO_CI2=(SO_CI*T)*SO_CI equal to unity matric I, overlap
if (ipglob > 4) then
  call dashes()
  write(u6,*) 'check whether SO_CI2=(SO_CI*T)*SO_CI equal to unity'
  call dashes()
  do i=1,lrootstot
    do j=1,lrootstot
      if (i /= j) then
        if ((abs(real(SO_CI2(i,j))) >= threshold) .or. (abs(aimag(SO_CI2(i,j))) >= threshold)) then
          write(u6,*) 'ERROR! SO_CI is not orthonormal:',i,j,SO_CI2(i,j)
        end if
      else if (i == j) then
        if ((abs(real(SO_CI2(i,j))-One) >= threshold) .or. (abs(aimag(SO_CI2(i,j))) >= threshold)) then
          write(u6,*) 'ERROR! SO_CI is not orthonormal',i,j,SO_CI2(i,j)
        end if
      end if
    end do
  end do
  write(u6,*) 'IF THERE IS NOT ANY ERROR INFO PRINTOUT, SO_CI is orthonomalized'
  call dashes()
end if

! construct the transformation matrix CSF2SO from SO to CSF,
! CSF2SO=U_CI*SO_CI and check whether CSF2SO coincide with
! the eigenvector VR of the full Hamiltonian HTOT_CSF which
! was obtained at the begining of this subroutine
call mult(cmplx(U_CI,kind=wp),SO_CI,CSF2SO)

if (ipglob > 4) then
  call dashes()
  write(u6,*) 'Printout the transformation matrix CSF2SO'
  call dashes()
  do i=1,10
    write(u6,*) (CSF2SO(i,j),j=1,10)
  end do
  ! Check
  call dashes()
  write(u6,*) 'Check if CSF2SO coincides with eigenvector Hfull'
  ! vk: attention that Hfull is of nconftot size
  !     probably code is correct as nconftot>=lrootstot
  call dashes()
  do i=1,nconftot
    do j=1,lrootstot
      if ((abs(abs(real(CSF2SO(i,j)))-abs(real(Hfull(i,j)))) >= threshold) .or. &
          (abs(abs(aimag(CSF2SO(i,j)))-abs(aimag(Hfull(i,j)))) >= threshold)) then
        if (((real(CSF2SO(i,j))**2)+(aimag(CSF2SO(i,j))**2))-((real(Hfull(i,j)))**2+(aimag(Hfull(i,j)))**2) >= threshold) then
          write(u6,*) 'WARNING! CSF2SO does not concide with Hfull:',i,j,real(CSF2SO(i,j)),aimag(CSF2SO(i,j)),real(Hfull(i,j)), &
                      aimag(Hfull(i,j)),real(CSF2SO(i,j))**2+aimag(CSF2SO(i,j))**2,real(Hfull(i,j))**2+(aimag(Hfull(i,j)))**2
        end if
      end if
    end do
  end do
  write(u6,*) 'IF THERE IS NO ANY WARNING INFO, CSF2SO IS COINCIDED WITH EIGENVECTOR sortVR!'
  call dashes()
end if

if (flag_so) then
  call mh5_put_dset(prep_csfsor,real(CSF2SO))
  call mh5_put_dset(prep_csfsoi,aimag(CSF2SO))
end if

if (allocated(Hfull)) call mma_deallocate(Hfull)
if (allocated(Hfull2)) call mma_deallocate(Hfull2)
if (allocated(Hdiag)) call mma_deallocate(Hdiag)
if (allocated(SO_CI2)) call mma_deallocate(SO_CI2)
if (allocated(WORK)) call mma_deallocate(WORK)
if (allocated(RWORK)) call mma_deallocate(RWORK)
if (allocated(W)) call mma_deallocate(W)

end subroutine soci
