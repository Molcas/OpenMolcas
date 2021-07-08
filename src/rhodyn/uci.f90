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
subroutine uci
  use rhodyn_data
  use rhodyn_utils, only: mult, dashes
  use stdalloc, only: mma_allocate, mma_deallocate
  use mh5
  implicit none
!
!***********************************************************************
! Purpose: construct the transformation matrix U_CI from the basis
!          of CSFs to the basis of spin-free states (accounting for
!          spin-degeneracy)
!***********************************************************************
!
!  UTU    : overlap matrix UTU=U_CI**T*U_CI for the spin-free states
!
  real(8),dimension(:,:),allocatable :: UTU

  call dashes()
  write(*,*) 'Begin of uci'

  write(*,*)'Dimensions of transformation matrix accounting'// &
                                      ' for spin-degeneracy'
  call dashes()
  write(*,sint)'Number of total CSFs:', nconftot
  write(*,sint)'Number of total states (roots):', lrootstot
  call dashes()

  U_CI=0d0
  ii=0
  jj=0
  kk=0
  if (.not.flag_so) then
    do k=1,N
      if (k/=1) then
        kk=kk+lroots(k-1)
      endif
      do i=1,nconf(k)
        ii=ii+1
        do j=1,lroots(k)
          U_CI(ii,j+kk) = CI(i,j,k)
        enddo
      enddo
    enddo
  else
    write(*,*) 'Construct transformation matrix U_CI'
    ll=0
    do l=1,N
      if (l/=1) then
        ll=ll+lroots(l-1)*ispin(l-1)
        kk=kk+nconf(l-1)*ispin(l-1)
      endif
      do i=1,nconf(l)
        do j=1,lroots(l)
          do k=1,ispin(l)
            ii=kk+(i-1)*ispin(l)+k
            jj=ll+(j-1)*ispin(l)+k
            U_CI(ii,jj)=CI(i,j,l)
          enddo
        enddo
      enddo
    enddo
  endif

  if (preparation/=2.and.preparation/=4) &
     call mh5_put_dset_array_real(prep_uci, U_CI)

! Check whether the trafo matrix U_CI is orthonormalized
  if (ipglob>2) then
    call mma_allocate(UTU,lrootstot,lrootstot)
    call mult(U_CI,U_CI,UTU,.True.,.False.)
    call dashes()
    write(*,*)'Internal check for CI coefficients'
    do i=1,lrootstot
      do j=1,lrootstot
        if (i/=j) then
          if (abs(UTU(i,j))>=threshold) then
            write(*,*)'ERROR INFO!!! CI coeffs are not'// &
                      ' orthonormalized',i,j,UTU(i,j)
          endif
        elseif (i==j) then
          if (abs(UTU(i,j)-1d0)>=threshold) then
            write(*,*)'ERROR INFO!!! CI coeffs are not'// &
                      ' orthonormalized',i,j,UTU(i,j)
          endif
        endif
      enddo
    enddo
    call dashes()
    write(*,*)'If there is no error info printout'
    write(*,*)'CI coeffs are orthonormalized'

    if (preparation/=2.and.preparation/=4) then
      prep_utu = mh5_create_dset_real (prep_id, &
          'UTU', 2, [lrootstot,lrootstot])
      call mh5_init_attr(prep_utu, 'description', &
          'UTU=U_CI**T*U_CI overlap matrix')
      call mh5_put_dset_array_real(prep_utu, UTU)
    endif

    call dashes()
    write(*,*) 'Overlap of transformation matrix is saved in file'
    call dashes()
    call mma_deallocate(UTU)
  endif

  write(*,*) 'End of uci'

end
