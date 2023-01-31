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

#include "macros.fh"

module rhodyn_utils
!***********************************************************************
! module contains some auxiliary routines
!***********************************************************************

use Definitions, only: wp, iwp, u6
use Constants, only: Zero, One, cZero
use linalg_mod, only: mult
use stdalloc, only: mma_allocate, mma_deallocate
#ifdef _ADDITIONAL_RUNTIME_CHECK_
use linalg_mod, only: abort_
#endif

implicit none
private

public :: dashes, removeColumn, removeLineAndColumn, sortci, transform, print_c_matrix, &
          check_hermicity, compare_matrices, WERDM, WERDM_back, WERSO, WERSO_back, W3J, W6J, get_kq_order

interface removeLineAndColumn
  module procedure removeLineAndColumnR, removeLineAndColumnZ
end interface
interface removeColumn
  module procedure removeColumnR, removeColumnZ
end interface
interface transform
  module procedure transformR, transformZ
end interface

contains

subroutine dashes(length)
  ! print the line of dashes of length 'length' to the standard output

  integer(kind=iwp), intent(in), optional :: length
  integer(kind=iwp) :: i, l

  if (present(length)) then
    l = length
  else
    l = 72
  end if
  do i=1,l
    write(u6,'(A)',advance='no') '-'
  end do
  write(u6,*)

end subroutine dashes

subroutine removeLineAndColumnZ(a,remLCarray)

  complex(kind=wp), allocatable, intent(inout) :: a(:,:)
  integer(kind=iwp), intent(in) :: remLCarray(:)
  integer(kind=iwp) :: i, j, l, m, n, k
  complex(kind=wp), allocatable :: b(:)
  logical(kind=iwp), allocatable :: mask(:,:)

  ! sizes
  n = size(a,1)
  m = size(a,2)
  k = size(remLCarray)
  call mma_allocate(mask,n,m,label='mask')
  ! mask creation
  l = 0
  mask(:,:) = .false.
  do i=1,n
    do j=1,m
      if (any(remLCarray == i) .and. any(remLCarray == j)) then
        mask(i,j) = .true.
        l = l+1
      end if
    end do
  end do
  call mma_allocate(b,l,label='b')
  ! copy remaining elements into temp b
  b(:) = pack(a,mask)
  call mma_deallocate(a)
  call mma_allocate(a,k,k,label='a')
  ! copy back
  a(:,:) = reshape(b,[k,k])
  call mma_deallocate(b)
  call mma_deallocate(mask)

end subroutine removeLineAndColumnZ

subroutine removeLineAndColumnR(a,remLCarray)

  real(kind=wp), allocatable, intent(inout) :: a(:,:)
  integer(kind=iwp), intent(in) :: remLCarray(:)
  integer(kind=iwp) :: i, j, l, m, n, k
  real(kind=wp), allocatable :: b(:)
  logical(kind=iwp), allocatable :: mask(:,:)

  ! sizes
  n = size(a,1)
  m = size(a,2)
  k = size(remLCarray)
  call mma_allocate(mask,n,m,label='a')
  ! mask creation
  l = 0
  mask(:,:) = .false.
  do i=1,n
    do j=1,m
      if (any(remLCarray == i) .and. any(remLCarray == j)) then
        mask(i,j) = .true.
        l = l+1
      end if
    end do
  end do
  call mma_allocate(b,l,label='b')
  ! copy remaining elements into temp b
  b(:) = pack(a,mask)
  call mma_deallocate(a)
  call mma_allocate(a,k,k,label='a')
  ! copy back
  a(:,:) = reshape(b,[k,k])
  call mma_deallocate(b)
  call mma_deallocate(mask)

end subroutine removeLineAndColumnR

subroutine removeColumnZ(a,remCarray)

  complex(kind=wp), allocatable, intent(inout) :: a(:,:)
  integer(kind=iwp), intent(in) :: remCarray(:)
  integer(kind=iwp) :: j, l, m, n, k
  complex(kind=wp), allocatable :: b(:)
  logical(kind=iwp), allocatable :: mask(:,:)

  ! sizes
  n = size(a,1)
  m = size(a,2)
  k = size(remCarray)
  call mma_allocate(mask,n,m,label='mask')
  ! mask creation
  l = 0
  mask(:,:) = .false.
  do j=1,m
    if (any(remCarray == j)) then
      mask(:,j) = .true.
      l = l+n
    end if
  end do
  call mma_allocate(b,l,label='b')
  ! copy remaining elements into temp b
  b(:) = pack(a,mask)
  call mma_deallocate(a)
  call mma_allocate(a,n,k,label='a')
  ! copy back
  a(:,:) = reshape(b,[n,k])
  call mma_deallocate(b)
  call mma_deallocate(mask)

end subroutine removeColumnZ

subroutine removeColumnR(a,remCarray)

  real(kind=wp), allocatable, intent(inout) :: a(:,:)
  integer(kind=iwp), intent(in) :: remCarray(:)
  integer(kind=iwp) :: j, l, m, n, k
  real(kind=wp), allocatable :: b(:)
  logical(kind=iwp), allocatable :: mask(:,:)

  ! sizes
  n = size(a,1)
  m = size(a,2)
  k = size(remCarray)
  call mma_allocate(mask,n,m,label='mask')
  ! mask creation
  l = 0
  mask(:,:) = .false.
  do j=1,m
    if (any(remCarray == j)) then
      mask(:,j) = .true.
      l = l+n
    end if
  end do
  call mma_allocate(b,l,label='b')
  ! copy remaining elements into temp b
  b(:) = pack(a,mask)
  call mma_deallocate(a)
  call mma_allocate(a,n,k,label='a')
  ! copy back
  a(:,:) = reshape(b,[n,k])
  call mma_deallocate(b)
  call mma_deallocate(mask)

end subroutine removeColumnR

subroutine transformZ(a,u,b,order)
  ! perform orthogonal transformation with transformation matrix
  ! order = True (default): direct transform : b = u^C * a * u
  ! order = False         : inverse transform: b = u   * a * u^C

  complex(kind=wp), intent(in) :: a(:,:), u(:,:)
  complex(kind=wp), intent(out) :: b(:,:)
  logical(kind=iwp), optional, intent(in) :: order
  integer(kind=iwp) :: m, n
  complex(kind=wp), allocatable :: temp(:,:)
  logical(kind=iwp) :: order_

  if (present(order)) then
    order_ = order
  else
    order_ = .true.
  end if
  m = size(u,1)
  n = size(u,2)
  if (order_) then
    ASSERT(m == size(a,1))
    call mma_allocate(temp,n,m,label='temp')
    call mult(u,a,temp,.true.,.false.)
    call mult(temp,u,b,.false.,.false.)
  else
    ASSERT(n == size(a,1))
    call mma_allocate(temp,m,n,label='temp')
    call mult(u,a,temp,.false.,.false.)
    call mult(temp,u,b,.false.,.true.)
  end if
  call mma_deallocate(temp)

end subroutine transformZ

subroutine transformR(a,u,b,order)
  ! perform orthogonal transformation with transformation matrix
  ! order = True (default): direct transform : b = u^T * a * u
  ! order = False       : inverse transform: b = u   * a * u^T

  real(kind=wp), intent(in) :: a(:,:), u(:,:)
  real(kind=wp), intent(out) :: b(:,:)
  logical(kind=iwp), optional, intent(in) :: order
  integer(kind=iwp) :: m, n
  real(kind=wp), allocatable :: temp(:,:)
  logical(kind=iwp) :: order_

  if (present(order)) then
    order_ = order
  else
    order_ = .true.
  end if
  m = size(u,1)
  n = size(u,2)
  if (order_) then
    ASSERT(m == size(a,1))
    call mma_allocate(temp,n,m,label='temp')
    call mult(u,a,temp,.true.,.false.)
    call mult(temp,u,b,.false.,.false.)
  else
    ASSERT(n == size(a,1))
    call mma_allocate(temp,m,n,label='temp')
    call mult(u,a,temp,.false.,.false.)
    call mult(temp,u,b,.false.,.true.)
  end if
  call mma_deallocate(temp)

end subroutine transformR

subroutine sortci(N1,A,WR,C,print_level)

  integer(kind=iwp), intent(in) :: N1, print_level
  real(kind=wp), intent(inout) :: A(N1,N1)
  real(kind=wp), intent(out) :: WR(N1), C(N1,N1)
  integer(kind=iwp) :: INFO, k, l, LWORK
  real(kind=wp), allocatable :: diag(:,:), B(:,:), WORK(:)

  if (print_level > 3) then
    call mma_allocate(B,N1,N1,label='B')
    call mma_allocate(diag,N1,N1,label='diag')
    B(:,:) = A
  end if
  LWORK = 2*N1
  call mma_allocate(WORK,LWORK,label='WORK')
  call dsyev_('V','U',N1,A,N1,WR,WORK,LWORK,INFO)
  if (INFO /= 0) then
    write(u6,*) 'ERROR in sortci'
    call Abend()
  end if
  call dsyev_('V','U',N1,A,N1,WR,WORK,LWORK,INFO)
  call mma_deallocate(WORK)
  C = A
  if (print_level > 3) then
    call transform(B,C,diag)
    call dashes()
    write(u6,*) 'Printout the diagonalized matrix:'
    call dashes()
    do k=1,10
      write(u6,*) (diag(k,l),l=1,10)
    end do
    call mma_deallocate(B)
    call mma_deallocate(diag)
  end if

end subroutine sortci

subroutine compare_matrices(A, B, n, header, thrs)
  integer(kind=iwp), intent(in) :: n
  complex(kind=wp), dimension(n,n), intent(in) :: A, B
  character(len=*), intent(in) :: header
  real(kind=wp), intent(in) :: thrs
  integer(kind=iwp) :: i
  logical :: AB_equal
  call dashes()
  write(u6,*) header
  AB_equal = .true.
  do i=1,n
    if (all(abs(A(:,i)-B(:,i)) < thrs)) then
      cycle
    else
      AB_equal = .false.
      exit
    end if
  end do
  if (AB_equal) write(u6,*) "matrices are equal"
  call dashes()
end subroutine compare_matrices

! routines for spherical tensor basis !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function DCLEBS(XJ1,XJ2,XJ3,XM1,XM2,XM3)
  real(kind=wp) :: DCLEBS
  real(kind=wp), intent(in) :: XJ1,XJ2,XJ3,XM1,XM2,XM3
  integer, parameter :: MAXJ=10, MAXF=3*MAXJ+1
  integer(kind=iwp), save :: icall = 0
  real(kind=wp), save :: DFACT(0:MAXF)
  real(kind=wp) :: DF, den, PRE, PRE2, SUMMA, TERM, XJSUM
  integer(kind=iwp) :: i, IA1, IA2, IA3, IB1, IB2, IB3, IX, IX1, IX2, IY, IY0, JSUM
! DCLEBS: real Clebsch-Gordan coefficients
! From a modification of Racah''s formula. Coded: Malmqvist 1998
! Note carefully: The input values XJ1..XM3 are REAL, not integers. Half-integer spins are allowed
! Half-integers are assumed exactly represented
  if (icall == 0) then
    icall = icall+1
    DF = One
    DFACT(0) = DF
    do i=1,MAXF
      DF = DBLE(i)*DF
      DFACT(i) = DF
    end do
  end if
  DCLEBS = Zero
  XJSUM = XJ1+XJ2+XJ3
  JSUM = NINT(XJSUM)
  if (XJSUM /= DBLE(JSUM)) return
  if (XM1+XM2 /= XM3) return
  IA1 = nint(XJ1+XM1)
  if (IA1 < 0) return
  IB1 = nint(XJ1-XM1)
  if (IB1 < 0) return
  IA2 = nint(XJ2+XM2)
  if (IA2 < 0) return
  IB2 = nint(XJ2-XM2)
  if (IB2 < 0) return
  IA3 = nint(XJ3-XM3)
  if (IA3 < 0) return
  IB3 = nint(XJ3+XM3)
  if (IB3 < 0) return
  if (JSUM-IA1-IB1 < 0) return
  if (JSUM-IA2-IB2 < 0) return
  if (JSUM-IA3-IB3 < 0) return
  PRE2 = real(1+IA3+IB3)*DFACT(JSUM-IA1-IB1) &
         *DFACT(JSUM-IA2-IB2)*DFACT(JSUM-IA3-IB3) &
         *DFACT(IA1)*DFACT(IA2)*DFACT(IA3) &
         *DFACT(IB1)*DFACT(IB2)*DFACT(IB3) &
         /DFACT(JSUM+1)
  PRE = sqrt(PRE2)
  IY0 = (JSUM-IA3-IB3)
  IX1 = (IA2+IB1-JSUM)+IB2
  IX2 = (IA2+IB1-JSUM)+IA1
  IX = MAX(0,IX1,IX2)
  IY = MIN(IY0,IB1,IA2)
  SUMMA = Zero
  do I=IX,IY
    DEN = DFACT(I)*DFACT(I-IX1)*DFACT(I-IX2)*DFACT(IY0-I)*DFACT(IB1-I)*DFACT(IA2-I)
    TERM = One/DEN
    SUMMA = SUMMA+real((-1)**I)*TERM
  end do
  DCLEBS = PRE*SUMMA
end function DCLEBS

function W3J(j1,j2,j3,m1,m2,m3)
! Calculates a Wigner 3-j symbol in the form
! { j1 j2 j3 }
! { m1 m2 m3 }
  real(kind=wp) :: W3j
  real(kind=wp), intent(in) :: j1, j2, j3, m1, m2, m3

  W3J = Zero
  !coeffCG=Zero
  !Call Clebsh_Gordan(j1, m1, j2, m2, j3,-m3, coeffCG)
  !If(coeffCG==Zero) Return
  W3J=DBLE((-1)**(nint(j1-j2-m3)))*DCLEBS(j1, j2, j3, m1, m2,-m3)/SQRT(DBLE(2*j3+1))
end function W3J

subroutine ITO(n,k,q,spins,projs,T)
! calculates the matrix <SM|T^K_Q|S'M'> of irreducible tensor operator
  integer(kind=iwp), intent(in) :: n, k, q
  real(kind=wp), dimension(n), intent(in) :: spins, projs
  real(kind=wp), intent(out) :: T(n,n)
  integer(kind=iwp) :: i, j
  real(kind=wp) :: s1, s2, m1, m2, fact

  do i=1,n
    do j=1,n
      s1 = spins(i)
      s2 = spins(j)
      m1 = projs(i)
      m2 = projs(j)
      fact = sqrt(dble(2*k+1))
      if (mod(int(s1-m1),2)==1) fact = -fact
      T(i,j) = fact * W3J(s1,dble(k),s2,-m1,dble(q),m2)
    enddo
  enddo
  return
end subroutine ITO

subroutine WERDM(rho,n_so,n_sf,k,q,spins,projs,so_sf,RED)
  ! calculates the elements of Wigner-Eckart reduced density matrix for fixed values
  ! of K, Q
  integer(kind=iwp), intent(in) :: n_so, n_sf, k, q
  integer(kind=iwp), dimension(n_so), intent(in) :: so_sf
  complex(kind=wp), intent(in) :: rho(n_so,n_so)
  real(kind=wp), dimension(n_so), intent(in) :: spins, projs
  complex(kind=wp), intent(out) :: RED(n_sf,n_sf)
  integer(kind=iwp) :: i, j, ii, jj
  real(kind=wp) :: T(n_so,n_so)

  RED = cZero

  ! calculate matrix of ITO
  call ITO(n_so,k,q,spins,projs,T)

  do i=1,n_so
    do j=1,n_so
      ! determine sf indices
      ii = so_sf(i)
      jj = so_sf(j)
      RED(ii,jj) = RED(ii,jj) + rho(i,j)*T(i,j)
!TEST
!      write(u6,*) 'i,j,ii,jj:',i,j,ii,jj
!      write(u6,*) 'rho, T, rho_red',rho(i,j),T(i,j),RED(ii,jj)
!TEST
    enddo
  enddo
  return
end subroutine WERDM

subroutine WERDM_back(RED,n_so,n_sf,len_sph,k_ranks,q_proj,spins,projs,so_sf,rho_back)
  ! calculates the density matrix in SO basis from the elements of
  ! Wigner-Eckart reduced density matrix, which are stored in 3d-matrix
  ! (len_sph,n_sf,n_sf)
  integer(kind=iwp), intent(in) :: n_so, n_sf, len_sph
  integer(kind=iwp), dimension(n_so), intent(in) :: so_sf
  integer(kind=iwp), dimension(len_sph), intent(in) :: k_ranks, q_proj
  real(kind=wp), dimension(n_so), intent(in) :: spins, projs
  complex(kind=wp), intent(in) :: RED(len_sph,n_sf,n_sf)
  complex(kind=wp), intent(out) :: rho_back(n_so,n_so)
  integer(kind=iwp) :: i, j, ii, jj, k, q, l
  real(kind=wp) :: T(n_so,n_so)

  rho_back = cZero

  do l=1,len_sph
    k = k_ranks(l)
    q = q_proj(l)
    ! calculate matrix of ITO
    call ITO(n_so,k,q,spins,projs,T)
    do i=1,n_so
      do j=1,n_so
        ! determine sf indices
        ii = so_sf(i)
        jj = so_sf(j)
        rho_back(i,j) = rho_back(i,j) + RED(l,ii,jj)*T(i,j)
      enddo
    enddo
  enddo
end subroutine WERDM_back

subroutine WERSO(vso,n_so,n_sf,so_sf,spins,projs,redvso)
  ! calculates matrix elements of Wigner-Eckart reduced spin-orbit Hamiltonian
  ! <iSM|Vso|jS'M'> = sqrt(3) sum_{m=0,+-1} (-1)^{S-M+m} 3j{S1S'-MmM'} <iS||Vso||jS'>
  integer(kind=iwp), intent(in) :: n_so, n_sf
  integer(kind=iwp), dimension(n_so), intent(in) :: so_sf
  real(kind=wp), dimension(n_so), intent(in) :: spins, projs
  complex(kind=wp), intent(in) :: vso(n_so,n_so)
  complex(kind=wp), intent(out) :: redvso(n_sf,n_sf,3)
  real(kind=wp) :: s1, s2, m1, m2, threejsymb
  integer(kind=iwp) :: i, j, ii, jj, m

  redvso = cZero

  do i=1,n_so
    do j=1,n_so
      s1 = spins(i)
      s2 = spins(j)
      m1 = projs(i)
      m2 = projs(j)
      ! determine sf indices
      ii = so_sf(i)
      jj = so_sf(j)
      do m=0,2,1
        threejsymb = W3J(s1,One,s2,-m1,dble(m-1),m2)
        if (threejsymb/=Zero) then ! this condition should be checked carefully
          redvso(ii,jj,m+1) = (-1)**(nint(s1-m1+m-1))*vso(i,j)/sqrt(3.0_wp)/threejsymb
        endif
      enddo
    enddo
  enddo
  return
end subroutine WERSO

subroutine WERSO_back(redvso,n_so,n_sf,so_sf,spins,projs,vso)
  ! calculates matrix elements of Wigner-Eckart reduced spin-orbit Hamiltonian
  ! <iSM|Vso|jS'M'> = sqrt(3) sum_{m=0,+-1} (-1)^{S-M+m} 3j{S1S'-MmM'} <iS||Vso||jS'>
  integer(kind=iwp), intent(in) :: n_so, n_sf
  integer(kind=iwp), dimension(n_so), intent(in) :: so_sf
  real(kind=wp), dimension(n_so), intent(in) :: spins, projs
  complex(kind=wp), intent(in) :: redvso(n_sf,n_sf,3)
  complex(kind=wp), intent(out) :: vso(n_so,n_so)
  real(kind=wp) :: s1, s2, m1, m2
  integer(kind=iwp) :: i, j, ii, jj, m

  vso = cZero
  do i=1,n_so
    do j=1,n_so
      s1 = spins(i)
      s2 = spins(j)
      m1 = projs(i)
      m2 = projs(j)
      ! determine sf indices
      ii = so_sf(i)
      jj = so_sf(j)
      do m=0,2,1
        vso(i,j) = vso(i,j) + (-1)**(nint(s1-m1+m-1)) * sqrt(3.0_wp) * W3J(s1,One,s2,-m1,dble(m-1),m2)* redvso(ii,jj,m+1)
      end do
    end do
  end do
end subroutine WERSO_back

subroutine print_c_matrix(A, n, header)
  integer(kind=iwp), intent(in) :: n
  complex(kind=wp), intent(in) :: A(n,n)
  character(len=*), intent(in) :: header
  integer(kind=iwp) :: i,j
  call dashes()
  write(u6,*) header
  do i = 1,n
      write(u6,*) (A(i,j),j=1,n)
  end do
end subroutine print_c_matrix

subroutine check_hermicity(A,n,A_name,thrs)
! Check whether matrix A is hermitian
  integer(kind=iwp), intent(in) :: n
  complex(kind=wp), intent(in) :: A(n,n)
  character(len=*), intent(in) :: A_name
  real(kind=wp), intent(in) :: thrs
  real(kind=wp) :: diff, abserror
  integer(kind=iwp) :: i, j

  abserror = Zero
  do i=1,n
    do j=1,i
      diff = abs(real(A(i,j))-real(A(j,i)))
      if ((diff >= thrs) .and. (diff >= abserror)) then
        abserror = diff
      end if
      diff = abs(aimag(A(i,j))+aimag(A(j,i)))
      if ((diff >= thrs) .and. (diff >= abserror)) then
        abserror = diff
      end if
    end do
  end do
  if (abserror >= thrs) then
    call WarningMessage(1,'Non-hermitian matrix obtained!')
    write(u6,'(a,1x,a,1x,a,1x,g28.16)') 'Matrix', A_name, 'Abs Error =', abserror
  end if
end subroutine check_hermicity

integer(kind=iwp) function get_kq_order(k_prime,q_prime)
  integer(kind=iwp), intent(in) :: k_prime, q_prime
  integer(kind=iwp) :: k, q
  k = 0
  get_kq_order = 0
  do while (k<=k_prime)
    do q=-k,k,1
      get_kq_order = get_kq_order + 1
      if (q == q_prime .and. k == k_prime) return
    end do
    k = k + 1
  end do
end function get_kq_order

!--------------------------------------------------------------------------------------------------------------------------------
Real*8 Function fct(n)
  integer, intent(in) :: n
  integer             :: i
  real(kind=8)        :: xct
  ! this function provides correct answer till n=169 only
  xct = One
  fct = One
  if (n < 0) then
    write(u6,'(A,i0)') 'FCT:  N<0 !'
    write(u6,'(A,i0)') 'N = ', N
    write(u6,'(A   )') 'It is an impossible case.'
    fct=-9.d99
    return

  else if (n == 0) then
    return

  else if (n <= 169) then
    do i=1,n
      xct=xct*DBLE(i)
    end do

  else
    write(u6,'(A,i0)') 'FCT:   N = ',N
    write(u6,'(A)') 'Factorial of N>169 overflows on x86_64'
    write(u6,'(A)') 'Use higher numerical precision, or rethink your algorithm.'
  end if

  fct=xct

  return
end function fct

real*8 function W6J(a,b,c,d,e,f)
! c Calculates a Wigner 6-j symbol. Argument a-f are positive Integer
! c and are twice the true value of the 6-j's arguments, in the form
! c { a b c }
! c { d e f }
  integer, intent(in) :: a,b,c,d,e,f
  integer             :: n,nlow,nhig
  real(kind=8)        :: sum,isum
  W6J=0.0d0
  if (MOD(a+b,2) /= MOD(c,2)) Return
  if (MOD(c+d,2) /= MOD(e,2)) Return
  if (MOD(a+e,2) /= MOD(f,2)) Return
  if (MOD(b+d,2) /= MOD(f,2)) Return
  if ((ABS(a-b) > c) .or. (a+b < c)) Return
  if ((ABS(c-d) > e) .or. (c+d < e)) Return
  if ((ABS(a-e) > f) .or. (a+e < f)) Return
  if ((ABS(b-d) > f) .or. (b+d < f)) Return
  if (check_triangle(a,b,c).eqv. .false.) Return
  if (check_triangle(c,d,e).eqv. .false.) Return
  if (check_triangle(a,e,f).eqv. .false.) Return
  if (check_triangle(b,d,f).eqv. .false.) Return
  nlow=0
  nhig=0
  nlow = MAX( (a+b+c)/2, (c+d+e)/2, (b+d+f)/2, (a+e+f)/2 )
  nhig = MIN( (a+b+d+e)/2, (b+c+e+f)/2, (a+c+d+f)/2)
  sum =0.0d0
  do n=nlow,nhig
    isum = DBLE((-1)**n)*fct(n+1)&
           /fct(  ( a+c+d+f)/2-n)&
           /fct(  ( b+c+e+f)/2-n)&
           /fct(n-( a+b+c  )/2  )&
           /fct(n-( c+d+e  )/2  )&
           /fct(n-( a+e+f  )/2  )&
           /fct(n-( b+d+f  )/2  )&
           /fct(  ( a+b+d+e)/2-n)
    sum=sum+isum
  end do
  W6J = dlt(a,b,c)*dlt(c,d,e)*dlt(a,e,f)*dlt(b,d,f)*sum
  return
end function W6J

Real*8 function dlt(a,b,c)
! c  calculates the delta(a,b,c) function using the formula 8.2.1. from:
! c    D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
! c    "Quantum Theory of Angular Momentum", World ScientIfic, 1988.
! c
! c  a,b,c are positive Integer numbers,
! c  their values are DoUBLE than their original value
  integer, intent(in) :: a,b,c
  dlt = Zero
  if ((ABS(a-b)>c) .or. (a+b<c)) return
  if ((ABS(b-c)>a) .or. (b+c<a)) return
  if ((ABS(c-a)>b) .or. (c+a<b)) return
  if (MOD((a+b-c),2) == 1) return
  if (MOD((a-b+c),2) == 1) return
  if (MOD((-a+b+c),2) == 1) return
  if (MOD((a+b+c),2) == 1) return
  if (check_triangle(a,b,c) .eqv. .false.) return
  ! special cases:
  if (a == 0) dlt = One/SQRT(DBLE(b+1))
  if (b == 0) dlt = One/SQRT(DBLE(a+1))
  if (c == 0) dlt = One/SQRT(DBLE(a+1))
  dlt=SQRT(fct((a+b-c)/2)*fct((a-b+c)/2)*fct((-a+b+c)/2)/fct((a+b+c)/2+1))
  return
end function dlt

logical function check_triangle(a,b,c)
  integer, intent(in) :: a, b, c
  check_triangle = .false.
  if (((a+b)>=c) .and. ((b+c)>=a) .and. ((c+a)>=b)) check_triangle = .true.
  return
end function check_triangle

end module rhodyn_utils
