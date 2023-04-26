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

subroutine k_external()
!***********************************************************************
! Purpose :  calculate dissipation rates k_ab
!***********************************************************************

use rhodyn_data, only: basis, E_SO, ipglob, ispin, k_bar_basis, kab_basis, lroots, n, nconftot, Nstate, SO_CI, U_CI
use rhodyn_utils, only: dashes, transform
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half, cZero, auToeV
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, j, k, l, ii, jj, kk, ll, iii, jjj, lu, max_i, max_j, n_sf
real(kind=wp) :: max_k
real(kind=wp), allocatable :: omega_ab(:,:), kab_real(:,:)
complex(kind=wp), allocatable :: k_ab(:,:)
character(len=256), parameter :: format1 = '(2(i8),2(g15.8,1x))'
integer(kind=iwp), external :: isFreeUnit

n_sf = sum(lroots)

write(u6,*) 'N_SF=',n_sf
write(u6,*) 'Nconftot',nconftot
write(u6,*) 'Nstate',Nstate

call mma_allocate(kab_real,n_sf,n_sf)
call mma_allocate(k_ab,Nstate,Nstate)
call mma_allocate(omega_ab,Nstate,Nstate)

call dashes()
write(u6,*) ' Begin reading k-matrix data '
call dashes()

lu = isFreeUnit(11)
call molcas_open(lu,'HRFACT')
do i=1,n_sf
  read(lu,*) (kab_real(i,j),j=1,n_sf)
end do
close(lu)

call dashes()
write(u6,*) ' End read k-matrix data '
call dashes()

! expand Kab to pseudo SF
k_ab = cZero
ii = 0
jj = 0
kk = 0
ll = 0
iii = 0
jjj = 0
do k=1,n ! manifolds
  if (k /= 1) then
    kk = kk+lroots(k-1)*ispin(k-1)
    ll = ll+lroots(k-1)*ispin(k-1)
    iii = iii+lroots(k-1)
    jjj = jjj+lroots(k-1)
  end if
  do i=1,lroots(k)
    do j=1,lroots(k)
      do l=1,ispin(k)
        ii = kk+(i-1)*ispin(k)+l
        jj = ll+(j-1)*ispin(k)+l
        k_ab(ii,jj) = kab_real(iii+i,jjj+j)
      end do
    end do
  end do
end do

lu = isFreeUnit(20)
call molcas_open(lu,'kab_out.dat')
do i=1,Nstate
  do j=1,Nstate
    omega_ab(i,j) = real(E_SO(i)-E_SO(j))
    if (real(k_ab(i,j)) >= (0.01_wp/autoeV)) then
      write(lu,format1) i,j,real(k_ab(i,j))*autoev,omega_ab(i,j)
    end if
  end do
end do
close(lu)

if (ipglob > 3) then
  call dashes()
  write(u6,*) ' Print matrix k_ab '
  call dashes()

  lu = isFreeUnit(13)
  call molcas_open(lu,'Kab_matrix_eV.dat')
  !!vk!! write procedure for printing matrices
  do i=1,Nstate
    write(lu,*) (real(k_ab(i,j))*autoev,j=1,Nstate)
  end do
  close(lu)
  max_k = Zero
  do i=1,Nstate
    do j=1,Nstate
      if (real(k_ab(i,j)) >= max_k) then
        max_k = real(k_ab(i,j))
        max_i = i
        max_j = j
      end if
    end do
  end do
  write(u6,*) Max_I,Max_J,Max_K*autoev,' eV',omega_ab(max_i,max_j)
end if

! transform the k_ab matrix to the required basis

select case (basis)
  case ('CSF')
    call transform(k_ab,cmplx(U_CI,kind=wp),kab_basis,.false.)
  case ('SO')
    call transform(k_ab,SO_CI,kab_basis)
  case ('SF')
    kab_basis(:,:) = k_ab
end select

! print out the bigger kab_basis

lu = isFreeUnit(22)
call molcas_open(lu,'max_kab_basis.dat')
max_k = Zero
do i=1,Nstate
  do j=1,Nstate
    if (abs(kab_basis(i,j)) >= max_k) then
      max_k = abs(kab_basis(i,j))
      max_i = i
      max_j = j
    end if
  end do
end do
write(lu,*) 'the maximum of Kab in ',basis
write(lu,'(2(i8),g15.8,a)') Max_I,Max_J,Max_K*autoev,' eV'
do i=1,Nstate
  do j=1,Nstate
    if (abs(Kab_basis(i,j)) >= (0.01_wp/autoeV)) then
      write(lu,'(2(i8),3(g15.8,1x))') i,j,abs(kab_basis(i,j)),real(Kab_basis(i,j))*autoev,aimag(Kab_basis(i,j))*autoev
    end if
  end do
end do
close(lu)

! construct the matrix (k_bar)_ij=0.5*sum_k[(kab_basis)_ik+(kab_basis)_jk]
do i=1,Nstate
  do j=1,Nstate
    do k=1,Nstate
      k_bar_basis(j,i) = abs(k_bar_basis(j,i)+Half*(kab_basis(j,k)+kab_basis(i,k)))
    end do
  end do
end do

if (allocated(kab_real)) call mma_deallocate(kab_real)
if (allocated(k_ab)) call mma_deallocate(k_ab)
if (allocated(omega_ab)) call mma_deallocate(omega_ab)

end subroutine k_external
