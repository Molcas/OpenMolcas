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

#include "compiler_features.h"
#ifdef _HDF5_

subroutine read_hdf5_all(file_h5,nss,nstate,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,amfi,HSO)

use mh5, only: mh5_open_file_r, mh5_fetch_attr, mh5_exists_dset, mh5_fetch_dset, mh5_close_file
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, cOne, Onei, auTocm, gElectron
use Definitions, only: wp, iwp, u6

implicit none
character(len=180), intent(in) :: file_h5
integer(kind=iwp), intent(in) :: nss, nstate
integer(kind=iwp), intent(out) :: multiplicity(nstate)
real(kind=wp), intent(out) :: eso(nss), esfs(nstate), angmom(3,nstate,nstate), edmom(3,nstate,nstate), amfi(3,nstate,nstate)
complex(kind=wp), intent(out) :: U(nss,nss), MM(3,nss,nss), MS(3,nss,nss), ML(3,nss,nss), DM(3,nss,nss), HSO(nss,nss)
integer(kind=iwp) :: fileid, i, i1, INRM, ipar, iss, ist, j, j1, jend, jst, l, mult, multI, multJ
real(kind=wp) :: RNRM
logical(kind=iwp) :: found_edmom !, found_amfi, found_angmom, found_esfs, found_eso, found_hso, found_mult, found_sos_coeff
integer(kind=iwp), allocatable :: ibas(:,:)
real(kind=wp), allocatable :: AL(:,:,:), etmp(:), RI(:,:), RR(:,:)
complex(kind=wp), allocatable :: tmp(:,:)
real(kind=wp), parameter :: g_e = -gElectron
logical(kind=iwp), parameter :: DBG = .false.
real(kind=wp), external :: dnrm2_, dznrm2_
complex(kind=wp), external :: spin

found_edmom = .false.
!found_angmom = .false.
!found_hso = .false.
!found_amfi = .false.
!found_sos_coeff = .false.
!found_eso = .false.
!found_esfs = .false.
!found_mult = .false.

write(u6,'(A,A)') 'Read data from rassi.h5 file ',trim(file_h5)

! the presence of the RASSI-HDF5 has already been made in
! read_hdf5_init subroutine, executed earlier
! nss, nstate are already known
!----------------------------------------------------------------------|
! open the file
fileid = mh5_open_file_r(trim(file_h5))
if (DBG) write(u6,'(A,I24)') 'read_hdf5_all:: fileid=',fileid
!----------------------------------------------------------------------|
! read spin multiplicity of each state:
!if (mh5_exists_dset(fileid,'STATE_SPINMULT')) then
!  found_mult = .true.
call mh5_fetch_attr(fileid,'STATE_SPINMULT',multiplicity)
if (DBG) then
  write(u6,'(A)') 'read_hdf5_all:: multiplicity'
  write(u6,'(20I4)') (multiplicity(i),i=1,nstate)
end if
INRM = sum(multiplicity(:))
if (INRM == 0) then
  call WarningMessage(2,'STATE_SPINMULT array read from HDF5 file has norm = zero')
  write(u6,*) 'Norm=',INRM
end if
!else
!  call WarningMessage(2,'State multiplicity array was not found on HDF5 file')
!end if

!----------------------------------------------------------------------|
! read the spin free energies (cm-1)
call mma_allocate(etmp,nstate,'etmp)')
if (mh5_exists_dset(fileid,'SFS_ENERGIES')) then
  !found_esfs = .true.
  call mh5_fetch_dset(fileid,'SFS_ENERGIES',etmp)
  RNRM = dnrm2_(nstate,etmp,1)
  if (RNRM < 1.0e-50_wp) then
    call WarningMessage(2,'ESFS read from HDF5 file  has norm = zero')
    write(u6,*) 'Norm=',RNRM
  end if
  ! compute the energeis in cm-1:
  esfs(:) = (etmp(:)-etmp(1))*auTocm
  if (DBG) then
    write(u6,'(A)') 'read_hdf5_all:: esfs'
    do i=1,nstate,4
      jEND = min(nstate,i+3)
      write(u6,'(4ES24.14)') (esfs(j),j=i,jEnd)
    end do
  end if
else
  call WarningMessage(2,'Spin-free energies were not found on HDF5 file')
end if
call mma_deallocate(etmp)

!----------------------------------------------------------------------|
! read the spin-orbit energies (cm-1)
call mma_allocate(etmp,nss,'tmp)')
etmp(:) = Zero
if (mh5_exists_dset(fileid,'SOS_ENERGIES')) then
  !found_eso = .true.
  call mh5_fetch_dset(fileid,'SOS_ENERGIES',etmp)
  RNRM = dnrm2_(nss,etmp,1)
  if (RNRM < 1.0e-50_wp) then
    call WarningMessage(2,'ESO read from HDF5 file  has norm = zero')
    write(u6,*) 'Norm=',RNRM
  end if
  ! compute the energies in cm-1:
  eso(:) = (etmp(:)-etmp(1))*auTocm
  if (DBG) then
    write(u6,'(A)') 'read_hdf5_all:: eso'
    do i=1,nss,4
      jEnd = min(nss,i+3)
      write(u6,'(4ES24.14)') (eso(j),j=i,jEnd)
    end do
  end if
else
  call WarningMessage(2,'Spin-orbit energies were not found on HDF5 file')
end if
call mma_deallocate(etmp)

!----------------------------------------------------------------------|
! read the spin-orbit mixing coefficient matrix:
call mma_allocate(RR,nss,nss,'RR')
call mma_allocate(RI,nss,nss,'RI')
U(:,:) = Zero
if (mh5_exists_dset(fileid,'SOS_COEFFICIENTS_REAL') .and. mh5_exists_dset(fileid,'SOS_COEFFICIENTS_IMAG')) then
  !found_sos_coeff = .true.
  call mh5_fetch_dset(fileid,'SOS_COEFFICIENTS_REAL',RR)
  call mh5_fetch_dset(fileid,'SOS_COEFFICIENTS_IMAG',RI)
  ! assemble the complex matrix U:
  U(:,:) = cmplx(RR(:,:),RI(:,:),kind=wp)
  RNRM = dznrm2_(nss*nss,U,1)
  if (RNRM < 1.0e-50_wp) then
    call WarningMessage(2,'SOS-U matrix read from HDF5 file  has norm = zero')
    write(u6,*) 'Norm=',RNRM
  end if
  if (DBG) then
    write(u6,'(A)') 'read_hdf5_all:: U'
    do i=1,nss
      do j=1,nss
        write(u6,'(2i4,A,2ES24.14)') i,j,' |',U(i,j)
      end do
    end do
  end if
else
  call WarningMessage(2,'SO mixing coefficeints were not found on HDF5 file')
end if
call mma_deallocate(RR)
call mma_deallocate(RI)

!----------------------------------------------------------------------|
! read the angular momentum integrals (L):
call mma_allocate(AL,nstate,nstate,3,'AL')
if (mh5_exists_dset(fileid,'SFS_ANGMOM')) then
  !found_angmom = .true.
  call mh5_fetch_dset(fileid,'SFS_ANGMOM',AL)
  do i=1,nstate
    do j=1,nstate
      ANGMOM(:,i,j) = AL(i,j,:)
    end do
  end do
  RNRM = dnrm2_(nstate*nstate*3,AL,1)
  if (RNRM < 1.0e-50_wp) then
    call WarningMessage(2,'SFS_ANGMOM read from HDF5 file  has norm = zero')
    write(u6,*) 'Norm=',RNRM
  end if
  if (DBG) then
    write(u6,'(A)') 'read_hdf5_all:: ANGMOM (x,y,z)'
    do i=1,nstate
      do j=1,nstate
        write(u6,'(2i4,A,3ES24.14)') i,j,' |',(ANGMOM(l,i,j),l=1,3)
      end do
    end do
  end if
else
  call WarningMessage(2,'ANGMOM integrals were not found on HDF5 file')
end if
call mma_deallocate(AL)

!----------------------------------------------------------------------|
! read the electric dipole momentum integrals (EDMOM):
call mma_allocate(AL,nstate,nstate,3,'AL')
if (mh5_exists_dset(fileid,'SFS_EDIPMOM')) then
  found_edmom = .true.
  call mh5_fetch_dset(fileid,'SFS_EDIPMOM',AL)
  do i=1,nstate
    do j=1,nstate
      EDMOM(:,i,j) = AL(i,j,:)
    end do
  end do
  RNRM = dnrm2_(nstate*nstate*3,AL,1)
  if (RNRM < 1.0e-50_wp) then
    call WarningMessage(2,'SFS_EDIPMOM read from HDF5 file  has norm = zero')
    write(u6,*) 'Norm=',RNRM
  end if
  if (DBG) then
    write(u6,'(A)') 'read_hdf5_all:: SFS_EDIPMOM(x,y,z)'
    do i=1,nstate
      do j=1,nstate
        write(u6,'(2i4,A,3ES24.14)') i,j,' |',(EDMOM(l,i,j),l=1,3)
      end do
    end do
  end if
else
  call WarningMessage(2,'EDMOM integrals were not found on HDF5 file')
end if
call mma_deallocate(AL)

!----------------------------------------------------------------------|
! read the spin-orbit integrals (AMFI):
call mma_allocate(AL,nstate,nstate,3,'AL')
if (mh5_exists_dset(fileid,'SFS_AMFIINT')) then
  !found_amfi = .true.
  call mh5_fetch_dset(fileid,'SFS_AMFIINT',AL)
  do i=1,nstate
    do j=1,nstate
      AMFI(:,i,j) = AL(i,j,:)
    end do
  end do
  RNRM = dnrm2_(nstate*nstate*3,AL,1)
  if (RNRM < 1.0e-50_wp) then
    call WarningMessage(2,'SFS_AMFIINT read from HDF5 file  has norm = zero')
    write(u6,*) 'Norm=',RNRM
  end if
  if (DBG) then
    write(u6,'(A)') 'read_hdf5_all:: SFS_AMFIINT(x,y,z)'
    do i=1,nstate
      do j=1,nstate
        write(u6,'(2i4,A,3ES24.14)') i,j,' |',(AMFI(l,i,j),l=1,3)
      end do
    end do
  end if
else
  call WarningMessage(2,'AMFI integrals were not found on HDF5 file')
end if
call mma_deallocate(AL)

!----------------------------------------------------------------------|
! read the RASSI SFS Hamiltonian (SFS_HAM):

!----------------------------------------------------------------------|
! read the RASSI SOS Hamiltonian (SOS_HAM):
call mma_allocate(RR,nss,nss,'RR')
call mma_allocate(RI,nss,nss,'RI')
if (mh5_exists_dset(fileid,'HSO_MATRIX_REAL') .and. mh5_exists_dset(fileid,'HSO_MATRIX_IMAG')) then
  !found_hso = .true.
  call mh5_fetch_dset(fileid,'HSO_MATRIX_REAL',RR)
  call mh5_fetch_dset(fileid,'HSO_MATRIX_IMAG',RI)
  ! assemble the complex matrix U:
  HSO(:,:) = cmplx(RR(:,:),RI(:,:),kind=wp)
  RNRM = dznrm2_(nss*nss,HSO,1)
  if (RNRM < 1.0e-50_wp) then
    call WarningMessage(2,'HSO matrix read from HDF5 file  has norm = zero')
    write(u6,*) 'Norm=',RNRM
  end if
  if (DBG) then
    write(u6,'(A)') 'read_hdf5_all:: HSO'
    do i=1,nss
      do j=1,nss
        write(u6,'(2i4,A,2ES24.14)') i,j,' |',HSO(i,j)
      end do
    end do
  end if
else
  call WarningMessage(2,'HSO matrix was not found on HDF5 file')
end if
call mma_deallocate(RR)
call mma_deallocate(RI)

!----------------------------------------------------------------------|
! All info has been read
! transform the SFS data to the SO basis
!----------------------------------------------------------------------|
! generate a local indexing table:
call mma_allocate(Ibas,[1,nstate],[-50,50],label='Ibas')
iss = 0
ibas = 0
ipar = mod(multiplicity(1),2)
do Ist=1,nstate
  Mult = Multiplicity(Ist)
  do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
    if ((Ipar == 0) .and. (I == 0)) cycle
    Iss = Iss+1
    Ibas(Ist,I) = Iss
  end do ! i
end do ! ist

! expand the spin free basis to the spin-orbit basis:
MM(:,:,:) = cZero
ML(:,:,:) = cZero
MS(:,:,:) = cZero
do Ist=1,nstate
  Mult = Multiplicity(Ist)
  do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
    if ((Ipar == 0) .and. (I == 0)) cycle
    do J=-(Mult-Ipar)/2,(Mult-Ipar)/2
      if ((Ipar == 0) .and. (J == 0)) cycle
      do l=1,3
        i1 = Ibas(Ist,I)
        j1 = Ibas(Ist,J)
        MM(l,i1,j1) = -Spin(l,Mult,I,J)*g_e
        MS(l,i1,j1) = Spin(l,Mult,I,J)
      end do ! l
    end do ! J
  end do ! I
end do ! Ist

do Ist=1,nstate
  MultI = Multiplicity(Ist)
  do Jst=1,nstate
    MultJ = Multiplicity(Jst)
    if (MultI == MultJ) then
      do I=-(MultI-Ipar)/2,(MultI-Ipar)/2
        if ((Ipar == 0) .and. (I == 0)) cycle
        i1 = Ibas(Ist,I)
        j1 = Ibas(Jst,I)
        MM(:,i1,j1) = MM(:,i1,j1)-Angmom(:,Ist,Jst)*Onei
        ML(:,i1,j1) = ML(:,i1,j1)+Angmom(:,Ist,Jst)*Onei
        DM(:,i1,j1) = DM(:,i1,j1)+eDmom(:,Ist,Jst)*cOne
      end do   ! I
    end if
  end do   ! Jst
end do   ! Ist

call mma_deallocate(Ibas)

! calculate the matrix elements of the spin and magnetic moment
! in the spin-orbit basis:
call mma_allocate(tmp,nss,nss,'tmp')
do L=1,3
  ! spin moment
  call ZGEMM_('C','N',nss,nss,nss,cOne,U,nss,MS(L,:,:),nss,cZero,TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,cOne,TMP,nss,U,nss,cZero,MS(L,:,:),nss)
  ! orbital moment
  call ZGEMM_('C','N',nss,nss,nss,cOne,U,nss,ML(L,:,:),nss,cZero,TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,cOne,TMP,nss,U,nss,cZero,ML(L,:,:),nss)
  ! magnetic moment
  call ZGEMM_('C','N',nss,nss,nss,cOne,U,nss,MM(L,:,:),nss,cZero,TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,cOne,TMP,nss,U,nss,cZero,MM(L,:,:),nss)

  if (found_EDMOM) then
    ! electric dipole moment
    call ZGEMM_('C','N',nss,nss,nss,cOne,U,nss,DM(L,:,:),nss,cZero,TMP,nss)
    call ZGEMM_('N','N',nss,nss,nss,cOne,TMP,nss,U,nss,cZero,DM(L,:,:),nss)
  end if
end do !L
call mma_deallocate(tmp)

! close the file
call mh5_close_file(fileid)

return

end subroutine read_hdf5_all

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(read_hdf5_all)

#endif
