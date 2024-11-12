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

subroutine read_hdf5_poly(file_h5,nss,nstate,eso,MM,MS,iReturn)

use mh5, only: mh5_open_file_r, mh5_fetch_attr, mh5_exists_dset, mh5_fetch_dset, mh5_close_file
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, cOne, Onei, auTocm, gElectron
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nss, nstate
real(kind=wp), intent(out) :: eso(nss)
complex(kind=wp), intent(out) :: MM(3,nss,nss), MS(3,nss,nss)
integer(kind=iwp), intent(out) :: iReturn
integer(kind=iwp) :: fileid, i, i1, INRM, ipar, iss, ist, j, j1, jst, l, mult, multI, multJ
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: jend
#endif
real(kind=wp) :: RNRM
!logical(kind=iwp) :: found_amfi, found_angmom, found_edmom, found_esfs, found_eso, found_hso, found_mult, found_sos_coeff
character(len=180) :: file_h5
integer(kind=iwp), allocatable :: ibas(:,:), multiplicity(:)
real(kind=wp), allocatable :: AL(:,:,:), angmom(:,:,:), etmp(:), RR(:,:), RI(:,:) !, amfi(:,:,:), edmom(:,:,:), esfs(:)
complex(kind=wp), allocatable :: tmp(:,:), tmp2(:,:), U(:,:) !, ML(:,:,:), DM(:,:,:), HSO(:,:)
real(kind=wp), parameter :: g_e = -gElectron
real(kind=wp), external :: dnrm2_, dznrm2_
complex(kind=wp), external :: spin

!found_edmom = .false.
!found_angmom = .false.
!found_hso = .false.
!found_amfi = .false.
!found_sos_coeff = .false.
!found_eso = .false.
!found_esfs = .false.
!found_mult = .false.
iReturn = 0

call mma_allocate(multiplicity,nstate,label='multiplicity')
call mma_allocate(angmom,3,nstate,nstate,label='angmom')

write(u6,'(A,A)') 'Read data from rassi.h5 file ',trim(file_h5)

! the presence of the RASSI-HDF5 has already been made in
! read_hdf5_init subroutine, executed earlier
! nss, nstate are already known
!----------------------------------------------------------------------|
! open the file
fileid = mh5_open_file_r(trim(file_h5))
#ifdef _DEBUGPRINT_
write(u6,'(A,I24)') 'read_hdf5_all:: fileid=',fileid
#endif
!----------------------------------------------------------------------|
! read spin multiplicity of each state:
!if (mh5_exists_dset(fileid,'STATE_SPINMULT')) then
!  found_mult = .true.
call mh5_fetch_attr(fileid,'STATE_SPINMULT',multiplicity)
#ifdef _DEBUGPRINT_
write(u6,'(A)') 'read_hdf5_all:: multiplicity'
write(u6,'(20I4)') (multiplicity(i),i=1,nstate)
#endif
INRM = 0
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
!call mma_allocate(etmp,nstate,'etmp)')
!etmp(:) = Zero
!if (mh5_exists_dset(fileid,'SFS_ENERGIES')) then
!  found_esfs = .true.
!  call mh5_fetch_dset(fileid,'SFS_ENERGIES',etmp)
!  RNRM = dnrm2_(nstate,etmp,1)
!  if (RNRM < 1.0e-50_wp) then
!    call WarningMessage(2,'ESFS read from HDF5 file has norm = zero')
!    write(u6,*) 'Norm=',RNRM
!  end if
!  ! compute the energeis in cm-1:
!  esfs(:) = (etmp(:)-etmp(1))*auTocm
!# ifdef _DEBUGPRINT_
!  write(u6,'(A)') 'read_hdf5_all:: esfs'
!  do i=1,nstate,4
!    jEND = MIN(nstate,i+3)
!    write(u6,'(4ES24.14)') (esfs(j),j=i,jEnd)
!  end do
!# endif
!else
!  call WarningMessage(2,'Spin-free energies were not found on HDF5 file')
!end if
!call mma_deallocate(etmp)

!----------------------------------------------------------------------|
! read the spin-orbit energies (cm-1)
call mma_allocate(etmp,nss,'tmp)')
etmp(:) = Zero
if (mh5_exists_dset(fileid,'SOS_ENERGIES')) then
  !found_eso = .true.
  call mh5_fetch_dset(fileid,'SOS_ENERGIES',etmp)
  RNRM = dnrm2_(nss,etmp,1)
  if (RNRM < 1.0e-50_wp) then
    call WarningMessage(2,'ESO read from HDF5 file has norm = zero')
    write(u6,*) 'Norm=',RNRM
  end if
  ! compute the energeis in cm-1:
  eso(:) = (etmp(:)-etmp(1))*auTocm
# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'read_hdf5_all:: eso'
  do i=1,nss,4
    jEnd = min(nss,i+3)
    write(u6,'(4ES24.14)') (eso(j),j=i,jEnd)
  end do
# endif
else
  call WarningMessage(2,'Spin-orbit energies were not found on HDF5 file')
end if
call mma_deallocate(etmp)

!----------------------------------------------------------------------|
! read the spin-orbit mixing coefficient matrix:
call mma_allocate(RR,nss,nss,'RR')
call mma_allocate(RI,nss,nss,'RI')
call mma_allocate(U,nss,nss,'U')
U(:,:) = cZero
if (mh5_exists_dset(fileid,'SOS_COEFFICIENTS_REAL') .and. mh5_exists_dset(fileid,'SOS_COEFFICIENTS_IMAG')) then
  !found_sos_coeff = .true.
  call mh5_fetch_dset(fileid,'SOS_COEFFICIENTS_REAL',RR)
  call mh5_fetch_dset(fileid,'SOS_COEFFICIENTS_IMAG',RI)
  ! assemble the complex matrix U:
  U(:,:) = cmplx(RR(:,:),RI(:,:),kind=wp)
  RNRM = dznrm2_(nss*nss,U,1)
  if (RNRM < 1.0e-50_wp) then
    call WarningMessage(2,'SOS-U matrix read from HDF5 file has norm = zero')
    write(u6,*) 'Norm=',RNRM
  end if
# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'read_hdf5_all:: U'
  do i=1,nss
    do j=1,nss
      write(u6,'(2i4,A,2ES24.14)') i,j,' |',U(i,j)
    end do
  end do
# endif
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
      ANGMOM(:,i,j) = AL(:,j,l)
    end do
  end do
  RNRM = dnrm2_(nstate*nstate*3,AL,1)
  if (RNRM < 1.0e-50_wp) then
    call WarningMessage(2,'SFS_ANGMOM read from HDF5 file has norm = zero')
    write(u6,*) 'Norm=',RNRM
  end if
# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'read_hdf5_all:: ANGMOM (x,y,z)'
  do i=1,nstate
    do j=1,nstate
      write(u6,'(2i4,A,3ES24.14)') i,j,' |',(ANGMOM(l,i,j),l=1,3)
    end do
  end do
# endif
else
  call WarningMessage(2,'ANGMOM integrals were not found on HDF5 file')
end if
call mma_deallocate(AL)

!----------------------------------------------------------------------|
! read the electric dipole momentum integrals (EDMOM):
!call mma_allocate(AL,nstate,nstate,3,'AL')
!if (mh5_exists_dset(fileid,'SFS_EDIPMOM')) then
!  found_edmom = .true.
!  call mh5_fetch_dset(fileid,'SFS_EDIPMOM',AL)
!  do i=1,nstate
!    do j=1,nstate
!      EDMOM(:,i,j) = AL(i,j,:)
!    end do
!  end do
!  RNRM = dnrm2_(nstate*nstate*3,AL,1)
!  if (RNRM < 1.0e-50_wp) then
!    call WarningMessage(2,'SFS_EDIPMOM read from HDF5 file has norm = zero')
!    write(u6,*) 'Norm=',RNRM
!  end if
!# ifdef _DEBUGPRINT_
!  write(u6,'(A)') 'read_hdf5_all:: SFS_EDIPMOM(x,y,z)'
!  do i=1,nstate
!    do j=1,nstate
!      write(u6,'(2i4,A,3ES24.14)') i,j,' |',(EDMOM(l,i,j),l=1,3)
!    end do
!  end do
!# endif
!else
!  call WarningMessage(2,'EDMOM integrals were not found on HDF5 file')
!end if
!call mma_deallocate(AL)

!----------------------------------------------------------------------|
! read the spin-orbit integrals (AMFI):
!call mma_allocate(AL,nstate,nstate,3,'AL')
!if (mh5_exists_dset(fileid,'SFS_AMFIINT')) then
!  found_amfi=.true.
!  call mh5_fetch_dset(fileid,'SFS_AMFIINT',AL)
!  do i=1,nstate
!    do j=1,nstate
!      AMFI(:,i,j) = AL(i,j,:)
!    end do
!  end do
!  RNRM = dnrm2_(nstate*nstate*3,AL,1)
!  if (RNRM < 1.0e-50_wp) then
!    call WarningMessage(2,'SFS_AMFIINT read from HDF5 file has norm = zero')
!    write(u6,*) 'Norm=',RNRM
!  end if
!# ifdef _DEBUGPRINT_
!  write(u6,'(A)') 'read_hdf5_all:: SFS_AMFIINT(x,y,z)'
!  do i=1,nstate
!    do j=1,nstate
!      write(u6,'(2i4,A,3ES24.14)') i,j,' |',(AMFI(l,i,j),l=1,3)
!    end do
!  end do
!# endif
!else
!  call WarningMessage(2,'AMFI integrals were not found on HDF5 file')
!end if
!call mma_deallocate(AL)
!
!----------------------------------------------------------------------|
! read the RASSI SFS Hamiltonian (SFS_HAM):

!----------------------------------------------------------------------|
! read the RASSI SOS Hamiltonian (SOS_HAM):
!call mma_allocate(RR,nss,nss,'RR')
!call mma_allocate(RI,nss,nss,'RI')
!if ( mh5_exists_dset(fileid,'HSO_MATRIX_REAL') .and. mh5_exists_dset(fileid,'HSO_MATRIX_IMAG') ) then
!  found_hso = .true.
!  call mh5_fetch_dset(fileid,'HSO_MATRIX_REAL',RR )
!  call mh5_fetch_dset(fileid,'HSO_MATRIX_IMAG',RI )
!  ! assemble the complex matrix U:
!  HSO(:,:) = cmplx(RR(:,:),RI(:,:),kind=wp)
!  RNRM = dznrm2_(nss*nss,HSO,1)
!  if (RNRM < 1.0e-50_wp) then
!    call WarningMessage(2,'HSO matrix read from HDF5 file has norm = zero')
!    write(u6,*) 'Norm=',RNRM
!  end if
!# ifdef _DEBUGPRINT_
!  write(u6,'(A)') 'read_hdf5_all:: HSO'
!  do i=1,nss
!    do j=1,nss
!      write(u6,'(2i4,A,2ES24.14)') i,j,' |',HSO(i,j)
!    end do
!  end do
!# endif
!else
!  call WarningMessage(2,'HSO matrix was not found on HDF5 file')
!end if
!call mma_deallocate(RR)
!call mma_deallocate(RI)

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
!ML(:,:,:) = cZero
MS(:,:,:) = cZero
do Ist=1,nstate
  Mult = Multiplicity(Ist)
  do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
    if ((Ipar == 0) .and. (I == 0)) cycle
    do J=-(Mult-Ipar)/2,(Mult-Ipar)/2
      if ((Ipar == 0) .and. (J == 0)) cycle
      i1 = Ibas(Ist,I)
      j1 = Ibas(Ist,J)
      do l=1,3
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
        !ML(:,i1,j1) = ML(:,i1,j1)+Angmom(:,Ist,Jst)*Onei
        !DM(:,i1,j1) = DM(:,i1,j1)+eDmom(:,Ist,Jst)*cOne
      end do   ! I
    end if
  end do   ! Jst
end do   ! Ist

call mma_deallocate(Ibas)

! calculate the matrix elements of the spin and magnetic moment
! in the spin-orbit basis:
call mma_allocate(tmp,nss,nss,'tmp')
call mma_allocate(tmp2,nss,nss,'tmp2')
do L=1,3
  ! spin moment
  tmp2(:,:) = MS(L,:,:)
  call ZGEMM_('C','N',nss,nss,nss,cOne,U,nss,tmp2,nss,cZero,TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,cOne,TMP,nss,U,nss,cZero,tmp2,nss)
  MS(L,:,:) = tmp2(:,:)
  ! orbital moment
  !tmp2(:,:) = ML(L,:,:)
  !call ZGEMM_('C','N',nss,nss,nss,cOne,U,nss,tmp2,nss,cZero,TMP,nss)
  !call ZGEMM_('N','N',nss,nss,nss,cOne,TMP,nss,U,nss,cZero,tmp2,nss)
  !ML(L,:,:) = tmp2(:,:)
  ! magnetic moment
  tmp2(:,:) = MM(L,:,:)
  call ZGEMM_('C','N',nss,nss,nss,cOne,U,nss,tmp2,nss,cZero,TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,cOne,TMP,nss,U,nss,cZero,tmp2,nss)
  MM(L,:,:) = tmp2(:,:)

  !if (found_EDMOM) then
  !  ! electric dipole moment
  !  tmp2(:,:) = DM(L,:,:)
  !  call ZGEMM_('C','N',nss,nss,nss,cOne,U,nss,tmp2,nss,cZero,TMP,nss)
  !  call ZGEMM_('N','N',nss,nss,nss,cOne,TMP,nss,U,nss,cZero,tmp2,nss)
  !  DM(L,:,:) = tmp2(:,:)
  !end if
end do !L
call mma_deallocate(tmp)
call mma_deallocate(tmp2)
call mma_deallocate(U)
call mma_deallocate(multiplicity)
call mma_deallocate(angmom)

! close the file
call mh5_close_file(fileid)

return

end subroutine read_hdf5_poly

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(read_hdf5_poly)

#endif
