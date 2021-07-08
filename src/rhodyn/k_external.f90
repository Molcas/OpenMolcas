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
subroutine k_external
!***********************************************************************
!
! Purpose :  calculate dissipation rates k_ab
!
!***********************************************************************
  use rhodyn_data
  use rhodyn_utils, only: transform, dashes
  use stdalloc, only: mma_allocate, mma_deallocate
  implicit none

  integer :: max_i, max_j, iii, jjj, n_sf
  real(8) ::  max_k
  real(8), dimension(:,:), allocatable :: omega_ab,kab_real
  complex(8), dimension(:,:),allocatable :: k_ab
  character(len=256):: format1 = '(2(I8),2(G15.8,X))'

  n_sf = sum(lroots)

  write(*,*) 'N_SF=', n_sf
  write(*,*) 'Nconftot', nconftot
  write(*,*) 'Nstate', Nstate

  call mma_allocate(kab_real,n_sf,n_sf)
  call mma_allocate(k_ab,Nstate,Nstate)
  call mma_allocate(omega_ab,Nstate,Nstate)

  call dashes()
  write(*,*) ' Begin reading k-matrix data '
  call dashes()

  open(11,file='HRFACT',status='old',iostat=error)
  if (error/=0) write(*,*) 'reading file HRFACT failed!'
  do i=1,n_sf
    read(11,*) (kab_real(i,j),j=1,n_sf)
  enddo
  close(11)
	  
  call dashes()
  write(*,*)' End read data '
  call dashes()
	  
! expand Kab to pseudo SF
	k_ab=zero
	ii=0
    jj=0
    kk=0
    ll=0
	iii=0
	jjj=0
  do k=1,n ! manifolds
    if (k/=1) then
      kk=kk+lroots(k-1)*ispin(k-1)
      ll=ll+lroots(k-1)*ispin(k-1)
	    iii=iii+lroots(k-1)
	    jjj=jjj+lroots(k-1)
    endif
    do i=1,lroots(k)
      do j=1,lroots(k)
        do l=1,ispin(k)
            ii=kk+(i-1)*ispin(k)+l
            jj=ll+(j-1)*ispin(k)+l
            k_ab(ii,jj)=kab_real(iii+i,jjj+j)
        enddo
      enddo
    enddo
  enddo

  open(20,file='kab_out.dat',status='replace')
  do i=1,Nstate
    do j=1,Nstate
      omega_ab(i,j)=E_SO(i)-E_SO(j)
      if (real(k_ab(i,j))>=(0.01/autoeV)) then
          write (20,format1) i,j,dble(k_ab(i,j))*autoev, &
             omega_ab(i,j)
      endif
    enddo
  enddo
  close(20)

  if (ipglob>3) then
    call dashes()
    write(*,*)' Print matrix k_ab '
    call dashes()

    open(13,file='Kab_matrix_eV.dat',status='replace')
!!vk!! write procedure for printing matrices
    do i=1,Nstate
      write(13,*)(dble(k_ab(i,j))*autoev,j=1,Nstate)
    enddo
    max_k=0d0
    do i=1,Nstate
      do j=1,Nstate
        if (real(k_ab(i,j))>=max_k) then
          max_k=k_ab(i,j)
          max_i=i
          max_j=j
        endif
      enddo
    enddo
    write(*,*) Max_I,Max_J,Max_K*autoev,' eV', &
    omega_ab(max_i,max_j)
    close(13)
  endif

! transform the k_ab matrix to the required basis

  select case (basis)
  case ('CSF')
    call transform(k_ab,dcmplx(U_CI,0d0),kab_basis,.False.)
  case ('SO')
    call transform(k_ab,SO_CI,kab_basis)
  case ('SF')
      kab_basis=k_ab
  end select

! print out the bigger kab_basis

  open (22,file='max_kab_basis.dat',status='replace')
  max_k=0d0
  do i=1,Nstate
    do j=1,Nstate
      if (abs(kab_basis(i,j))>=max_k) then
        max_k=abs(kab_basis(i,j))
        max_i=i
        max_j=j
      endif
    enddo
  enddo
  write(22,*)'the maximum of Kab in ', basis
  write(22,'(2(I8),G15.8,A)')Max_I,Max_J,Max_K*autoev,' eV'
  do i=1,Nstate
    do j=1,Nstate
      if(abs(Kab_basis(i,j))>=(0.01/autoeV))then
        write(22,'(2(I8),3(G15.8,X))')i,j,abs(kab_basis(i,j)), &
      dble(Kab_basis(i,j))*autoev,aimag(Kab_basis(i,j))*autoev
      endif
    enddo
  enddo
  close(22)
	  
! contruct the matrix (k_bar)_ij=0.5*sum_k[(kab_basis)_ik+(kab_basis)_jk]
  do i=1,Nstate
    do j=1,Nstate
      do k=1,Nstate
        k_bar_basis(j,i) = abs(k_bar_basis(j,i)+0.5d0* &
                       (kab_basis(j,k)+kab_basis(i,k)))
      enddo
    enddo
  enddo

  write(*,*) 'End k_external'

	if (allocated(kab_real)) call mma_deallocate(kab_real)
  if (allocated(k_ab)) call mma_deallocate(k_ab)
	if (allocated(omega_ab)) call mma_deallocate(omega_ab)

end
