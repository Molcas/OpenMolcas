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

subroutine kab()
!***********************************************************************
! Purpose :  calculate dissipation rates k_ab
!***********************************************************************

use rhodyn_data, only: basis, CSF2SO, d, E_SO, cgamma, HRSO, ipglob, ispin, k_b, K_bar_basis, kab_basis, n, nconf, nconftot, &
                       Nmode, Nstate, SO_CI, T
use rhodyn_utils, only: dashes, transform
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, auToCm, auToeV
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, j, k, ii, jj, iii, jjj, lu, max_i, max_j, n_sf
real(kind=wp) :: max_k
real(kind=wp), allocatable :: freq(:), G(:,:), G_SF(:,:), J_w(:,:), n_w(:,:), omega_ab(:,:), r_ab_SO(:,:), temp_gk(:)
complex(kind=wp), allocatable :: G_SO(:,:,:), gamma_pd(:,:), gamma_pd_basis(:,:), k_ab(:,:)
character(len=256), parameter :: format1 = '(2(i8),4(g15.8,1x))'
integer(kind=iwp), external :: isFreeUnit

n_sf = sum(nconf)

write(u6,*) 'N_SF=',n_sf
write(u6,*) 'Nmode=',Nmode
write(u6,*) 'Nconftot=',nconftot
write(u6,*) 'Nstate=',Nstate

call mma_allocate(G,Nmode,n_sf)
call mma_allocate(G_SF,Nmode,Nstate)
call mma_allocate(G_SO,Nmode,Nstate,Nstate)
call mma_allocate(freq,Nmode)
call mma_allocate(temp_gk,Nstate)
call mma_allocate(r_ab_SO,Nstate,Nstate)

! reading quintets and triplets

call dashes()
write(u6,*) ' Begin reading HR-factors data '
call dashes()

jj = 0
iii = 0

lu = isFreeUnit(11)
call molcas_open(lu,'HRFACT')
read(lu,*)
do k=1,Nmode
  read(lu,*) freq(k)
end do
read(lu,*)
if (HRSO) then
  do i=1,NState
    read(lu,*)
    do k=1,Nmode
      write(u6,*) 'k=',k
      do j=1,Nstate
        read(lu,'(ES16.8)',advance='no') temp_gk(j)
        write(u6,*) temp_gk(j)
        G_SO(k,i,j) = cmplx(temp_gk(j),kind=wp)
      end do
      read(lu,*)
    end do
  end do
else
  do k=1,Nmode
    read(lu,*) (G(k,j),j=2,n_sf)
    G(k,1) = Zero
    write(u6,*) k,(G(k,j),j=2,n_sf)
  end do
  !!!!
  ! the gradient for each SF states is calculated using RASSCF while
  ! the other data are calculated in DFT (Gaussian), so there is a shift
  ! from the ground to ground state G(I,1)!=0, here we shift
  ! all the SF states with respect to ground states
  do k=1,n
    if (k == 1) then
      iii = 0
    else
      iii = iii+nconf(k-1)*ispin(k-1)
    end if
    jj = jj+nconf(k)
    do i=1,Nmode
      jjj = iii
      do j=jj-nconf(k)+1,jj
        do ii=1,ispin(k)
          jjj = jjj+1
          G_SF(i,jjj) = G(i,j)
        end do
      end do
    end do
  end do
  call dashes(72)
  do i=1,Nmode
    write(u6,*) (G_SF(i,j),j=1,nconftot)
  end do
  ! construct the G_SO matrix
  do i=1,NState
    do j=1,NState
      do k=1,Nmode
        do ii=1,Nstate
          G_SO(k,i,j) = G_SO(k,i,j)+conjg(SO_CI(ii,i))*G_SF(k,ii)*SO_CI(ii,j)
        end do
      end do
    end do
  end do
end if !HRSO
close(lu) ! close HRFACT file

! change the unit to au
Freq(:) = Freq/auToCm

call dashes()
write(u6,*) ' End read data '
call dashes()

! Spectrum density

call mma_allocate(k_ab,Nstate,Nstate)
call mma_allocate(n_w,Nstate,Nstate)
call mma_allocate(J_w,Nstate,Nstate)
call mma_allocate(omega_ab,Nstate,Nstate)

write(u6,*) 'Gamma=',cgamma,'Hartree',cgamma*auToCm,'cm-1'

lu = isFreeUnit(20)
call molcas_open(lu,'kab_out.dat')
J_w = Zero
do i=1,Nstate
  ! let us cut the spectral density (no), if yes then j=1,i
  do j=1,Nstate
    omega_ab(i,j) = real(E_SO(i)-E_SO(j))
    do k=1,Nmode
      J_w(i,j) = J_w(i,j)+abs(G_SO(k,i,j))**2*Freq(k)**2* &
                 (abs(omega_ab(i,j))*Freq(k)*cgamma/((abs(omega_ab(i,j))**2-Freq(k)**2)**2+omega_ab(i,j)**2*cgamma**2))
    end do
    if (i == j) then
      n_w(i,j) = Zero
    else
      n_w(i,j) = One/(exp(abs(omega_ab(i,j))/(k_b*T))-One)
    end if
    ! note that the K_ab matrix should be real matrix, but here we define
    ! it to be complex for the convernient of the following transform
    if (omega_ab(i,j) >= Zero) then
      k_ab(i,j) = 4*(1+n_w(i,j))*J_w(i,j)
    else
      k_ab(i,j) = 4*n_w(i,j)*J_w(i,j)
    end if
    if (real(k_ab(i,j)) >= (0.01_wp/autoeV)) then
      if (omega_ab(i,j) >= Zero) then
        write(lu,format1) i,j,real(k_ab(i,j))*autoev,omega_ab(i,j)**2,(One+n_w(i,j)),J_w(i,j)
      else
        write(lu,format1) i,j,real(k_ab(i,j))*autoev,omega_ab(i,j)**2,n_w(i,j),J_w(i,j)
      end if
    end if

  end do
end do

close(lu) ! close kab_out.dat

! for this  spectral density the pure decay matrix gamma_pd equal
! to 0.0 always because when the frequency equal
! to 0.0, J_w(0.0)=0.0 so gamma_pd=0.0
call mma_allocate(gamma_pd,Nstate,Nstate)
gamma_pd = Zero

if (ipglob > 3) then
  call dashes()
  write(u6,*) ' Print matrix K_ab '
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
  write(u6,*) Max_I,Max_J,Max_K*autoev,' eV',omega_ab(max_i,max_j),1+n_w(Max_I,Max_J),J_w(Max_I,Max_J)
end if

! transform the k_ab and gamma_pd matrix to the required basis

call mma_allocate(gamma_pd_basis,d,d)

select case (basis)
  case ('CSF')
    ! k_ab SO->CSF
    call transform(k_ab,CSF2SO,kab_basis,.false.)
    ! gamma_pd SO->CSF
    call transform(gamma_pd,CSF2SO,gamma_pd_basis,.false.)
  case ('SF')
    ! K_ab SO->SF
    call transform(k_ab,SO_CI,kab_basis,.false.)
    ! gamma_pd SO->SF
    call transform(gamma_pd,SO_CI,gamma_pd_basis,.false.)
    ! Kab_basis = real(Kab_basis)
  case ('SO')
    kab_basis(:,:) = k_ab
    gamma_pd_basis(:,:) = gamma_pd
end select

! print out the bigger kab_basis

lu = isFreeUnit(22)
call molcas_open(lu,'max_Kab_basis.dat')

max_k = Zero
do i=1,d
  do j=1,d
    if (abs(kab_basis(i,k)) >= max_k) then
      max_k = abs(kab_basis(i,j))
      max_i = i
      max_j = j
    end if
  end do
end do

write(lu,*) 'the maximum of Kab in ',basis
write(lu,'(2(i8),g15.8,a)') max_i,max_j,max_k*autoev,' eV'

do i=1,d
  do j=1,d
    if (abs(kab_basis(i,j)) >= (0.01_wp/autoeV)) then
      write(lu,'(2(i8),3(g15.8,1x))') i,j,abs(kab_basis(i,j)),real(kab_basis(I,J))*autoev,aimag(kab_basis(I,J))*autoev
    end if
  end do
end do

!if (basis == 'SO') then
!  Max_K = Zero
!  do i=161,Nstate
!    do j=161,Nstate
!      if (abs(Kab_basis(i,j)) >= Max_K) then
!        Max_K = abs(Kab_basis(i,j))
!        Max_I = i
!        Max_J = j
!      end if
!    end do
!  end do

!  write(lu,*) 'Maximum of Kab in core states in basis of ',Basis
!  write(lu,'(2(I8),G15.8,A)') Max_I,Max_J,Max_K*autoev,' eV'

!  do i=1,Nstate
!    do j=1,Nstate
!      if ((abs(Kab_basis(i,j)) <= 0.01_wp/autoeV) .and. (abs(Kab_basis(i,j)) >= 0.005_wp/autoeV)) then
!        write(lu,'(2(I8),3(G15.8,X))') I,J,abs(Kab_basis(i,j)),real(Kab_basis(i,j))*autoev,aimag(Kab_basis(i,j))*autoev
!      end if
!    end do
!  end do

!else if (basis == 'SF') then

!  Max_k = Zero
!  do i=26,175
!    do j=26,175
!      if (abs(Kab_basis(I,J)) >= Max_K) then
!        Max_K = abs(Kab_basis(i,j))
!        Max_i = i
!        Max_j = j
!      end if
!    end do
!  end do

!  write(lu,*) 'Maximum of Kab in quintet core states in basis of ',Basis
!  write(lu,'(2(I8),G15.8,A)') Max_I,Max_J,Max_K*autoev,' eV'

!  Max_k = Zero
!  do i=311,Nstate
!    do j=311,Nstate
!      if (abs(Kab_basis(I,J)) >= Max_K) then
!        Max_K = abs(Kab_basis(I,J))
!        Max_i = i
!        Max_j = j
!      end if
!    end do
!  end do

!  write(lu,*) 'Maximum of Kab in triplets core states in basis of ',Basis
!  write(lu,'(2(I8),G15.8,A)') Max_i,Max_j,Max_k*autoev,' eV'

!  do i=1,Nstate
!    do j=1,Nstate
!      if ((abs(Kab_basis(i,J)) <= 0.01_wp/autoeV) .and. (abs(Kab_basis(I,J)) >= 0.005_wp/autoeV)) then
!        write(lu,'(2(I8),3(G15.8,X))') I,J,abs(Kab_basis(I,J)),real(Kab_basis(I,J))*autoev,aimag(Kab_basis(I,J))*autoev
!      end if
!    end do
!  end do
!end if

close(lu) ! close file max_Kab_basis.dat

! construct the matrix (k_bar)_ij=0.5*sum_k[(kab_basis)_ik+(kab_basis)_jk]

do i=1,d
  do j=1,d
    do k=1,d
      k_bar_basis(j,i) = abs(k_bar_basis(j,i)+Half*(kab_basis(j,k)+kab_basis(i,k)))
    end do
  end do
end do

K_bar_basis(:,:) = K_bar_basis+gamma_pd_basis
lu = isFreeUnit(23)
call molcas_open(lu,'r_ab_SO.dat')

! check whether coherence dephasing r_ab=0.5*[sum_e(K_ae+K_be)]>0 in SO basis
do i=1,Nstate
  do j=1,Nstate
    do k=1,Nstate
      r_ab_SO(j,i) = r_ab_SO(j,i)+Half*(real(K_ab(j,k))+real(K_ab(i,k)))
    end do
  end do
end do

do i=1,Nstate
  do j=1,Nstate
    if (r_ab_SO(i,j) <= Zero) then
      write(lu,*) 'R_ab_SO<0',i,j,r_ab_SO(i,j)
    end if
  end do
end do
close(lu)

if (allocated(G)) call mma_deallocate(G)
if (allocated(G_SF)) call mma_deallocate(G_SF)
if (allocated(G_SO)) call mma_deallocate(G_SO)
if (allocated(Freq)) call mma_deallocate(Freq)
if (allocated(temp_gk)) call mma_deallocate(temp_gk)
if (allocated(r_ab_SO)) call mma_deallocate(r_ab_SO)
if (allocated(k_ab)) call mma_deallocate(k_ab)
if (allocated(n_w)) call mma_deallocate(n_w)
if (allocated(J_w)) call mma_deallocate(J_w)
if (allocated(omega_ab)) call mma_deallocate(omega_ab)
if (allocated(gamma_pd)) call mma_deallocate(gamma_pd)
if (allocated(gamma_pd_basis)) call mma_deallocate(gamma_pd_basis)

end subroutine kab
