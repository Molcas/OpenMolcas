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
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************
subroutine get_dipole
!
! Purpose :  Read in the dipole matrix from the MOLCAS output (SO)
!
  use rhodyn_data
  use rhodyn_utils, only: transform, dashes
  use definitions, only: wp, iwp, u6
  use mh5, only: mh5_open_file_r, mh5_exists_dset, mh5_fetch_dset, &
                 mh5_close_file, mh5_put_dset
  use stdalloc, only: mma_allocate, mma_deallocate
  implicit none
  integer(kind=iwp) :: fileid
  real(kind=wp),allocatable,dimension(:,:,:) :: DIPR, DIPI

  call dashes()
  write(u6,*) 'Begin get_dipole'
  call dashes()

  call mma_allocate (DIPR,lrootstot,lrootstot,3)
  call mma_allocate (DIPI,lrootstot,lrootstot,3)

! Read in the components dipole matrix (X,Y,Z) from the MOLCAS output
! probably better move that to read_rassisd.f90
  fileid = mh5_open_file_r('RASSISD')
  write(u6,*) 'dipole real read'
  if (flag_so) then
    if (mh5_exists_dset(fileid,'SOS_EDIPMOM_REAL').and. &
      mh5_exists_dset(fileid,'SOS_EDIPMOM_IMAG')) then
      call mh5_fetch_dset(fileid,'SOS_EDIPMOM_REAL',DIPR)
      call mh5_fetch_dset(fileid,'SOS_EDIPMOM_IMAG',DIPI)
    else
      write(u6,*) 'Error in reading RASSISD file, no dipole matrix in SO basis'
      call abend()
    endif
  else ! to read SFS_EDIPMOM (flag_so = off)
    if (mh5_exists_dset(fileid,'SFS_EDIPMOM')) then
      call mh5_fetch_dset(fileid,'SFS_EDIPMOM',DIPR)
      DIPI=0d0
    else
      write(u6,*) 'Error in reading RASSISD file, no dipole matrix in SF basis'
      call abend()
    endif
  endif
  !write(u6,*) 'dysorb read'
  if (mh5_exists_dset(fileid,'DYSORB').and.flag_dyson) then
    call mh5_fetch_dset(fileid,'DYSORB',dysamp)
  else if (mh5_exists_dset(fileid,'DYSAMP').and.flag_dyson) then
    call mh5_fetch_dset(fileid,'DYSAMP',dysamp)
  else
    write(u6,*) 'Ionization is not taken into account (set flag DYSO) and/or'
    write(u6,*) 'RASSI file does not contain Dyson amplitudes'
    flag_dyson=.False.
  endif
  !write(u6,*) 'dysorb has been read'
  call mh5_close_file(fileid)
  dipole = dcmplx(DIPR,DIPI)

! To put the imaginary part of diagonal elements to 0 just in case
  do j=1,lrootstot
    dipole(j,j,:)=dble(dipole(j,j,:))
  enddo
! process Dyson amplitudes
  if (flag_dyson.and.N>2) then
! Transformation of Dyson amplitudes matrix to SF basis
    call transform(dcmplx(dysamp,0d0),SO_CI,dysamp_bas,.False.)
! nullify non-neighbouring SF blocks, ii,jj - block indices
    ii=0
    jj=0
    do k=1,N
      do l=1,k
        if ((k-l)>1) then
          do i=ii,(ii+lroots(k)*ispin(k))
            do j=jj,(jj+lroots(l)*ispin(l))
              dysamp_bas(i,j)=zero
              dysamp_bas(j,i)=zero
            enddo
          enddo
        end if
        jj=jj+lroots(l)*ispin(l)
      enddo
      ii=ii+lroots(k)*ispin(k)
      jj=0
    enddo
  else
    dysamp_bas = dysamp
  endif
  if (flag_dyson) dysamp_bas=abs(dysamp_bas**2)
  write(u6,*)'dysorb processing has been successfully finished'

! calculate matrix of Einstein coefficient A if emission spectrum needed
!      if ((DM_basis/='SO').and.(DM_basis/='CSF_SO').and.
!  &    (DM_basis/='SF_SO')) then
!     flag_emiss=.False.
!   write(u6,*) 'Emission spectra can be calculated only if SOC'//
!  &             'density matrix available. Set DMBasis keyword' //
!  &             'to SO, CSF_SO, or SF_SO'
!  endif
  if (flag_emiss) then
    a_einstein = 0d0
    emiss = 0d0
    ii = 1
    do j=1,(lrootstot-1)
      do i=(j+1),lrootstot
        do k=1,3
          a_einstein(i,j) = a_einstein(i,j) + abs(dipole(i,j,k))**2
        enddo
        ! write down frequencies
        emiss(ii) = abs(E_SO(i)-E_SO(j))
        a_einstein(i,j) = a_einstein(i,j) * emiss(ii)**3
        ii = ii + 1
      enddo
    enddo
  endif

  if (preparation/=4) then
    ! not CM case
    call mh5_put_dset(prep_dipoler, dble(dipole))
    call mh5_put_dset(prep_dipolei, aimag(dipole))
    if (flag_dyson) then
      call mh5_put_dset(prep_do, dble(dysamp_bas))
    endif
  endif

  call mma_deallocate (DIPR)
  call mma_deallocate (DIPI)

  write(u6,*) 'End get_dipole'

end
