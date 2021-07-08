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
subroutine get_hcsf()
!***********************************************************************
!
! Purpose: Construct the H(RASSCF) matrix in the total CSF basis
! of different spin manifolds. H(RASSCF) has been already read
! in the matrix H_CSF(:,:,N) for different spin manifolds
!
!***********************************************************************
  use rhodyn_data
  use rhodyn_utils, only: dashes
  use mh5
  implicit none

  call dashes()
  write(*,*) 'Begin get_Hcsf'

! Construct the Hamiltonian matrix H(RASSCF) in CSF basis include all
! spin manifolds to ! HTOTRE_CSF
  HTOTRE_CSF=0d0
  ii=0
  jj=0
  kk=0
  if (flag_so) then
    do i=1,N
      if (i/=1) ii=ii+nconf(i-1)*ispin(i-1)
      do j=1,nconf(i)
        do k=1,nconf(i)
          do l=1,ispin(i)
            jj=ii+(j-1)*ispin(i)+l
            kk=ii+(k-1)*ispin(i)+l
            HTOTRE_CSF(jj,kk)=H_CSF(j,k,i)
          enddo
        enddo
      enddo
    enddo
  else
! without accounting for spin-degeneracy
    do i=1,N
      if (i/=1) ii=ii+nconf(i-1)
      do j=1,nconf(i)
        do k=1,nconf(i)
          jj=ii+j
          kk=ii+k
          HTOTRE_CSF(jj,kk)=H_CSF(j,k,i)
        enddo
      enddo
    enddo
  endif

  write(*,*) 'Construct the full Hamiltonian!'
  call dashes()
  write(*,sint)'number of NCSF:',nconftot
  call dashes()
  HTOT_CSF=(0d0,0d0)

! If consider the spin-orbit coupling
! Computing the total Hamiltonian H(RASSCF)+H(SO) in CSF basis
! to HTOT_CSF
  if (.not.flag_so) then
    HTOT_CSF = HTOTRE_CSF
  else
    HTOT_CSF = HTOTRE_CSF + V_CSF
  endif
  write(*,*) 'end contructing full Hamiltonian'

! Check whether total Hamiltonian is Hermitian
  if (ipglob>3) then
    call dashes()
    print*,'Check whether total Hamiltonian HTOT_CSF is Hermitian'
    call dashes()
    do  i=1,nconftot
      do j=1,nconftot
        if (abs(dble(HTOT_CSF(i,j)-HTOT_CSF(j,i)))>=threshold)then
        write(*,int2real) 'WARNING!!!: HTOT_CSF matrix is not'// &
                          ' hermitian in real part:',i,j, &
                          dble(HTOT_CSF(i,j)),dble(HTOT_CSF(j,i))
        endif
        if(abs(aimag(HTOT_CSF(i,j)+HTOT_CSF(j,i)))>=threshold)then
        write(*,int2real) 'WARNING!!!: HTOT_CSF matrix is not'// &
                          ' hermitian in imag part:',i,j, &
                          aimag(HTOT_CSF(i,j)),aimag(HTOT_CSF(j,i))
        endif
      enddo
    enddo

    write(*,*) 'If there is no warning info, total'// &
              ' Hamiltonian matrix HTOT_CSF is hermitian!'
  endif

  call mh5_put_dset_array_real(prep_fhr, dble(HTOT_CSF))
  call mh5_put_dset_array_real(prep_fhi, aimag(HTOT_CSF))

end
