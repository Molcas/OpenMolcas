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
  use definitions, only: u6
  use mh5, only: mh5_put_dset
  use stdalloc, only: mma_allocate, mma_deallocate
  implicit none

  if (ipglob>2) then
    call dashes()
    write(u6,*) 'Begin get_dipole'
    call dashes()
  endif

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
  else if (flag_dyson) then
    dysamp_bas(:,:) = dysamp
  endif
  if (flag_dyson) then
    dysamp_bas = abs(dysamp_bas**2)
    if (ipglob>2) write(u6,*)'dysamp processing successfully finished'
  endif

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

  if (runmode/=4) then
    ! not CM case
    call mh5_put_dset(prep_dipoler, dble(dipole))
    call mh5_put_dset(prep_dipolei, aimag(dipole))
    if (flag_dyson) then
      call mh5_put_dset(prep_do, dble(dysamp_bas))
    endif
  endif

  if (ipglob>2) write(u6,*) 'End get_dipole'

end
