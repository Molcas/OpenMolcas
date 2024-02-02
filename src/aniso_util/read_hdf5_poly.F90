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

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer, intent(in) :: nstate, nss
integer :: multiplicity(nstate)
integer :: iReturn
!real(kind=8), intent(out) :: esfs(nstate)
real(kind=8), intent(out) :: eso(nss)
!real(kind=8), intent(out) ::  edmom(3,nstate,nstate)
!real(kind=8), intent(out) ::   amfi(3,nstate,nstate)
real(kind=8) :: angmom(3,nstate,nstate)
complex(kind=8) :: MM(3,nss,nss)
complex(kind=8) :: MS(3,nss,nss)
!complex(kind=8), intent(out) :: ML(3,nss,nss)
!complex(kind=8), intent(out) :: DM(3,nss,nss) ! electric dipole moment
!complex(kind=8), intent(out) :: HSO(nss,nss)
!complex(kind=8) :: U(nLoc,nLoc)
real(kind=8) :: AU2CM
real(kind=8), allocatable :: etmp(:)
real(kind=8), allocatable :: RR(:,:), RI(:,:)
real(kind=8), allocatable :: AL(:,:,:)
complex(kind=8), allocatable :: U(:,:)
integer :: fileid, jend, INRM
character(Len=180) :: file_h5
real(kind=8) :: RNRM
real(kind=8), external :: dnrm2_, dznrm2_
complex(kind=8), external :: spin
! local variables:
integer :: iss, ibas(nstate,-50:50)
integer :: i, j, i1, j1, ist, jst, mult, multI, multJ
integer :: l, ipar
real(kind=8) :: g_e
complex(kind=8), allocatable :: tmp(:,:)
!logical :: Exist
!logical :: found_edmom, found_angmom, found_hso, found_amfi, found_sos_coeff, found_eso, found_esfs, found_mult
logical :: DBG

DBG = .false.
AU2CM = 219474.6313702_wp
g_e = 2.00231930437180_wp
!found_edmom = .false.
!found_angmom = .false.
!found_hso = .false.
!found_amfi = .false.
!found_sos_coeff = .false.
!found_eso = .false.
!found_esfs = .false.
!found_mult = .false.
iReturn = 0

write(6,'(A,A)') 'Read data from rassi.h5 file ',trim(file_h5)

! the presence of the RASSI-HDF5 has already been made in
! read_hdf5_init subroutine, executed earlier
! nss, nstate are already known
!----------------------------------------------------------------------|
! open the file
fileid = mh5_open_file_r(trim(file_h5))
if (DBG) write(6,'(A,I24)') 'read_hdf5_all:: fileid=',fileid
!----------------------------------------------------------------------|
! read spin multiplicity of each state:
!if (mh5_exists_dset(fileid,'STATE_SPINMULT')) then
!  found_mult = .true.
call mh5_fetch_attr(fileid,'STATE_SPINMULT',multiplicity(1:nstate))
if (DBG) then
  write(6,'(A)') 'read_hdf5_all:: multiplicity'
  write(6,'(20I4)') (multiplicity(i),i=1,nstate)
end if
INRM = 0
INRM = sum(multiplicity(:))
if (INRM == 0) then
  call WarningMessage(2,'STATE_SPINMULT array read from HDF5 file has norm = zero')
  write(6,*) 'Norm=',INRM
end if
!else
!  call WarningMessage(2,'State multiplicity array was not found on HDF5 file')
!end if

!----------------------------------------------------------------------|
! read the spin free energies (cm-1)
!call mma_allocate(etmp,nstate,'etmp)')
!call dcopy_(nstate,[0.0_wp],0,etmp,1)
!if (mh5_exists_dset(fileid,'SFS_ENERGIES')) then
!  found_esfs = .true.
!  call mh5_fetch_dset(fileid,'SFS_ENERGIES',etmp)
!  RNRM = 0.0_wp
!  RNRM = dnrm2_(nstate,etmp,1)
!  if (RNRM < 1.0D-50) then
!    call WarningMessage(2,'ESFS read from HDF5 file has norm = zero')
!    write(6,*) 'Norm=',RNRM
!  end if
!  ! compute the energeis in cm-1:
!  do i=1,nstate
!    esfs(i) = (etmp(i)-etmp(1))*AU2CM
!  end do
!  if (DBG) then
!    write(6,'(A)') 'read_hdf5_all:: esfs'
!    do i=1,nstate,4
!      jEND = MIN(nstate,i+3)
!      write(6,'(4ES24.14)') (esfs(j),j=i,jEnd)
!    end do
!  end if
!else
!  call WarningMessage(2,'Spin-free energies were not found on HDF5 file')
!end if
!call mma_deallocate(etmp)

!----------------------------------------------------------------------|
! read the spin-orbit energies (cm-1)
call mma_allocate(etmp,nss,'tmp)')
call dcopy_(nss,[0.0_wp],0,etmp,1)
if (mh5_exists_dset(fileid,'SOS_ENERGIES')) then
  !found_eso = .true.
  call mh5_fetch_dset(fileid,'SOS_ENERGIES',etmp)
  RNRM = 0.0_wp
  RNRM = dnrm2_(nss,etmp,1)
  if (RNRM < 1.0D-50) then
    call WarningMessage(2,'ESO read from HDF5 file has norm = zero')
    write(6,*) 'Norm=',RNRM
  end if
  ! compute the energeis in cm-1:
  do i=1,nss
    eso(i) = (etmp(i)-etmp(1))*AU2CM
  end do
  if (DBG) then
    write(6,'(A)') 'read_hdf5_all:: eso'
    do i=1,nss,4
      jEnd = min(nss,i+3)
      write(6,'(4ES24.14)') (eso(j),j=i,jEnd)
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
call mma_allocate(U,nss,nss,'U')
call dcopy_(nss*nss,[0.0_wp],0,RR,1)
call dcopy_(nss*nss,[0.0_wp],0,RI,1)
call zcopy_(nss*nss,[(0.0_wp,0.0_wp)],0,U,1)

if (mh5_exists_dset(fileid,'SOS_COEFFICIENTS_REAL') .and. mh5_exists_dset(fileid,'SOS_COEFFICIENTS_IMAG')) then
  !found_sos_coeff = .true.
  call mh5_fetch_dset(fileid,'SOS_COEFFICIENTS_REAL',RR)
  call mh5_fetch_dset(fileid,'SOS_COEFFICIENTS_IMAG',RI)
  ! assemble the complex matrix U:
  do i=1,nss
    do j=1,nss
      U(i,j) = cmplx(RR(i,j),RI(i,j),wp)
    end do
  end do
  RNRM = 0.0_wp
  RNRM = dznrm2_(nss*nss,U,1)
  if (RNRM < 1.0D-50) then
    call WarningMessage(2,'SOS-U matrix read from HDF5 file has norm = zero')
    write(6,*) 'Norm=',RNRM
  end if
  if (DBG) then
    write(6,'(A)') 'read_hdf5_all:: U'
    do i=1,nss
      do j=1,nss
        write(6,'(2i4,A,2ES24.14)') i,j,' |',U(i,j)
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
call dcopy_(nstate*nstate*3,[0.0_wp],0,AL,1)
if (mh5_exists_dset(fileid,'SFS_ANGMOM')) then
  !found_angmom = .true.
  call mh5_fetch_dset(fileid,'SFS_ANGMOM',AL)
  do i=1,nstate
    do j=1,nstate
      do l=1,3
        ANGMOM(l,i,j) = AL(i,j,l)
      end do
    end do
  end do
  RNRM = 0.0_wp
  RNRM = dnrm2_(nstate*nstate*3,AL,1)
  if (RNRM < 1.0D-50) then
    call WarningMessage(2,'SFS_ANGMOM read from HDF5 file has norm = zero')
    write(6,*) 'Norm=',RNRM
  end if
  if (DBG) then
    write(6,'(A)') 'read_hdf5_all:: ANGMOM (x,y,z)'
    do i=1,nstate
      do j=1,nstate
        write(6,'(2i4,A,3ES24.14)') i,j,' |',(ANGMOM(l,i,j),l=1,3)
      end do
    end do
  end if
else
  call WarningMessage(2,'ANGMOM integrals were not found on HDF5 file')
end if
call mma_deallocate(AL)

!----------------------------------------------------------------------|
! read the electric dipole momentum integrals (EDMOM):
!call mma_allocate(AL,nstate,nstate,3,'AL')
!call dcopy_(nstate*nstate*3,[0.0_wp],0,AL,1)
!if (mh5_exists_dset(fileid,'SFS_EDIPMOM')) then
!  found_edmom = .true.
!  call mh5_fetch_dset(fileid,'SFS_EDIPMOM',AL)
!  do i=1,nstate
!    do j=1,nstate
!      do l=1,3
!        EDMOM(l,i,j) = AL(i,j,l)
!      end do
!    end do
!  end do
!  RNRM = 0.0_wp
!  RNRM = dnrm2_(nstate*nstate*3,AL,1)
!  if (RNRM < 1.0D-50) then
!    call WarningMessage(2,'SFS_EDIPMOM read from HDF5 file has norm = zero')
!    write(6,*) 'Norm=',RNRM
!  end if
!  if (DBG) then
!    write(6,'(A)') 'read_hdf5_all:: SFS_EDIPMOM(x,y,z)'
!    do i=1,nstate
!      do j=1,nstate
!        write(6,'(2i4,A,3ES24.14)') i,j,' |',(EDMOM(l,i,j),l=1,3)
!      end do
!    end do
!  end if
!else
!  call WarningMessage(2,'EDMOM integrals were not found on HDF5 file')
!end if
!call mma_deallocate(AL)

!----------------------------------------------------------------------|
! read the spin-orbit integrals (AMFI):
!call mma_allocate(AL,nstate,nstate,3,'AL')
!call dcopy_(nstate*nstate*3,[0.0_wp],0,AL,1)
!if (mh5_exists_dset(fileid,'SFS_AMFIINT')) then
!  found_amfi=.true.
!  call mh5_fetch_dset(fileid,'SFS_AMFIINT',AL)
!  do i=1,nstate
!    do j=1,nstate
!      do l=1,3
!        AMFI(l,i,j) = AL(i,j,l)
!      end do
!    end do
!  end do
!  RNRM = 0.0_wp
!  RNRM = dnrm2_(nstate*nstate*3,AL,1)
!  if (RNRM < 1.0D-50) then
!    call WarningMessage(2,'SFS_AMFIINT read from HDF5 file has norm = zero')
!    write(6,*) 'Norm=',RNRM
!  end if
!  if (DBG) then
!    write(6,'(A)') 'read_hdf5_all:: SFS_AMFIINT(x,y,z)'
!    do i=1,nstate
!      do j=1,nstate
!        write(6,'(2i4,A,3ES24.14)') i,j,' |',(AMFI(l,i,j),l=1,3)
!      end do
!    end do
!  end if
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
!call dcopy_(nss*nss,[0.0_wp],0,RR,1)
!call dcopy_(nss*nss,[0.0_wp],0,RI,1)
!
!if ( mh5_exists_dset(fileid,'HSO_MATRIX_REAL') .and. mh5_exists_dset(fileid,'HSO_MATRIX_IMAG') ) then
!  found_hso = .true.
!  call mh5_fetch_dset(fileid,'HSO_MATRIX_REAL',RR )
!  call mh5_fetch_dset(fileid,'HSO_MATRIX_IMAG',RI )
!  ! assemble the complex matrix U:
!  do i=1,nss
!    do j=1,nss
!      HSO(i,j) = cmplx(RR(i,j),RI(i,j),wp)
!    end do
!  end do
!  RNRM = 0.0_wp
!  RNRM = dznrm2_(nss*nss,HSO,1)
!  if (RNRM < 1.0D-50) then
!    call WarningMessage(2,'HSO matrix read from HDF5 file has norm = zero')
!    write(6,*) 'Norm=',RNRM
!  end if
!  if (DBG) then
!    write(6,'(A)') 'read_hdf5_all:: HSO'
!    do i=1,nss
!      do j=1,nss
!        write(6,'(2i4,A,2ES24.14)') i,j,' |',HSO(i,j)
!      end do
!    end do
!  end if
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
iss = 0
ibas = 0
ipar = mod(multiplicity(1),2)
do Ist=1,nstate
  Mult = Multiplicity(Ist)
  do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
    if ((Ipar == 0) .and. (I == 0)) go to 310
    Iss = Iss+1
    Ibas(Ist,I) = Iss
310 continue
  end do ! i
end do ! ist

! expand the spin free basis to the spin-orbit basis:
MM = (0.0_wp,0.0_wp)
!ML = (0.0_wp,0.0_wp)
MS = (0.0_wp,0.0_wp)
do Ist=1,nstate
  Mult = Multiplicity(Ist)
  do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
    if ((Ipar == 0) .and. (I == 0)) go to 301
    do J=-(Mult-Ipar)/2,(Mult-Ipar)/2
      if ((Ipar == 0) .and. (J == 0)) go to 302
      do l=1,3
        i1 = Ibas(Ist,I)
        j1 = Ibas(Ist,J)
        MM(l,i1,j1) = -Spin(l,Mult,I,J)*g_e
        MS(l,i1,j1) = Spin(l,Mult,I,J)
      end do ! l
302   continue
    end do ! J
301 continue
  end do ! I
end do ! Ist

do Ist=1,nstate
  MultI = Multiplicity(Ist)
  do Jst=1,nstate
    MultJ = Multiplicity(Jst)
    if (MultI == MultJ) then
      do I=-(MultI-Ipar)/2,(MultI-Ipar)/2
        if ((Ipar == 0) .and. (I == 0)) go to 303
        do l=1,3
          i1 = Ibas(Ist,I)
          j1 = Ibas(Jst,I)
          MM(l,i1,j1) = MM(l,i1,j1)-cmplx(0.0_wp,Angmom(l,Ist,Jst),wp)
          !ML(l,i1,j1) = ML(l,i1,j1)+cmplx(0.0_wp,Angmom(l,Ist,Jst),wp)
          !DM(l,i1,j1) = DM(l,i1,j1)+cmplx(eDmom(l,Ist,Jst),0.0_wp,wp)
        end do   ! l
303     continue
      end do   ! I
    end if
  end do   ! Jst
end do   ! Ist

! calculate the matrix elements of the spin and magnetic moment
! in the spin-orbit basis:
call mma_allocate(tmp,nss,nss,'tmp')
do L=1,3
  TMP = (0.0_wp,0.0_wp)
  ! spin moment
  call ZGEMM_('C','N',nss,nss,nss,(1.0_wp,0.0_wp),U,nss,MS(L,:,:),nss,(0.0_wp,0.0_wp),TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,(1.0_wp,0.0_wp),TMP,nss,U,nss,(0.0_wp,0.0_wp),MS(L,:,:),nss)
  ! orbital moment
  ! TMP = (0.0_wp,0.0_wp)
  ! call ZGEMM_('C','N',nss,nss,nss,(1.0_wp,0.0_wp),U,nss,ML(L,:,:),nss,(0.0_wp,0.0_wp),TMP,nss)
  ! call ZGEMM_('N','N',nss,nss,nss,(1.0_wp,0.0_wp),TMP,nss,U,nss,(0.0_wp,0.0_wp),ML(L,:,:),nss)
  ! magnetic moment
  TMP = (0.0_wp,0.0_wp)
  call ZGEMM_('C','N',nss,nss,nss,(1.0_wp,0.0_wp),U,nss,MM(L,:,:),nss,(0.0_wp,0.0_wp),TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,(1.0_wp,0.0_wp),TMP,nss,U,nss,(0.0_wp,0.0_wp),MM(L,:,:),nss)

  !if (found_EDMOM) then
  !  ! electric dipole moment
  !  TMP = (0.0_wp,0.0_wp)
  !  call ZGEMM_('C','N',nss,nss,nss,(1.0_wp,0.0_wp),U,nss,DM(L,:,:),nss,(0.0_wp,0.0_wp),TMP,nss)
  !  call ZGEMM_('N','N',nss,nss,nss,(1.0_wp,0.0_wp),TMP,nss,U,nss,(0.0_wp,0.0_wp),DM(L,:,:),nss)
  !end if
end do !L
call mma_deallocate(tmp)
call mma_deallocate(U)

! close the file
call mh5_close_file(fileid)

return

end subroutine read_hdf5_poly

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(read_hdf5_poly)

#endif
