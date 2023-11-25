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

use linalg_mod, only: mult
use stdalloc, only: mma_allocate, mma_deallocate
#ifdef _ADDITIONAL_RUNTIME_CHECK_
use linalg_mod, only: abort_
#endif
use Constants, only: Zero, One, Two, Three, cZero
use Definitions, only: wp, iwp, u6

implicit none
private

public :: check_hermicity, compare_matrices, dashes, get_kq_order, print_c_matrix, removeColumn, removeLineAndColumn, sortci, &
          transform, W3J, W6J, WERDM, WERDM_back, WERSO, WERSO_back

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
  ! cut complex matrix a(n,m), i.e. delete rows and columns not given by array remLCarray(k)
  ! thus the result is matrix a of size (k,k)

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
  ! cut real matrix a(n,m), i.e. delete rows and columns not given by array remLCarray(k)
  ! thus the result is matrix a of size (k,k)

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
  ! cut complex matrix a(n,m), i.e. delete columns not given by array remCarray(k)
  ! thus the result is matrix a of size (n,k)

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
  ! cut real matrix a(n,m), i.e. delete columns not given by array remCarray(k)
  ! thus the result is matrix a of size (n,k)

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

subroutine compare_matrices(A,B,n,header,thrs)
  ! check if complex matrices A and B are equal within given threshold thrs

  integer(kind=iwp), intent(in) :: n
  complex(kind=wp), intent(in) :: A(n,n), B(n,n)
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
      write(u6,*) "error in",i
      exit
    end if
  end do
  if (AB_equal) write(u6,*) "matrices are equal"
  call dashes()

end subroutine compare_matrices

function DCLEBS(XJ1,XJ2,XJ3,XM1,XM2,XM3)
  ! DCLEBS: real Clebsch-Gordan coefficients
  ! From a modification of Racah''s formula. Coded: Malmqvist 1998
  ! Note carefully: The input values XJ1..XM3 are REAL, not integers. Half-integer spins are allowed
  ! Half-integers are assumed exactly represented

  real(kind=wp) :: DCLEBS
  real(kind=wp), intent(in) :: XJ1, XJ2, XJ3, XM1, XM2, XM3
  integer, parameter :: MAXJ = 10, MAXF = 3*MAXJ+1
  integer(kind=iwp), save :: icall = 0
  real(kind=wp), save :: DFACT(0:MAXF)
  real(kind=wp) :: DF, den, PRE, PRE2, SUMMA, TERM, XJSUM
  integer(kind=iwp) :: i, IA1, IA2, IA3, IB1, IB2, IB3, IX, IX1, IX2, IY, IY0, JSUM

  if (icall == 0) then
    icall = icall+1
    DF = One
    DFACT(0) = DF
    do i=1,MAXF
      DF = real(i,kind=wp)*DF
      DFACT(i) = DF
    end do
  end if
  DCLEBS = Zero
  XJSUM = XJ1+XJ2+XJ3
  JSUM = nint(XJSUM)
  if (XJSUM /= real(JSUM,kind=wp)) return
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
  IX = max(0,IX1,IX2)
  IY = min(IY0,IB1,IA2)
  SUMMA = Zero
  do I=IX,IY
    DEN = DFACT(I)*DFACT(I-IX1)*DFACT(I-IX2)*DFACT(IY0-I)*DFACT(IB1-I)*DFACT(IA2-I)
    TERM = One/DEN
    SUMMA = SUMMA+real((-1)**I,kind=wp)*TERM
  end do
  DCLEBS = PRE*SUMMA

end function DCLEBS

function W3J(j1,j2,j3,m1,m2,m3)
  ! Calculates a Wigner 3-j symbol in the form
  ! { j1 j2 j3 }
  ! { m1 m2 m3 }

  real(kind=wp) :: W3j
  real(kind=wp), intent(in) :: j1, j2, j3, m1, m2, m3
  integer(kind=iwp) :: i

  W3J = DCLEBS(j1,j2,j3,m1,m2,-m3)
  if (W3J == Zero) return
  i = nint(j1-j2-m3)
  if (i /= (i/2)*2) W3J = -W3J
  W3J = W3J/sqrt(Two*j3+One)

end function W3J

subroutine ITO(n,k,q,spins,projs,T)
  ! calculates the matrix <SM|T^K_Q|S'M'> of irreducible tensor operator

  integer(kind=iwp), intent(in) :: n, k, q
  real(kind=wp), intent(in) :: spins(n), projs(n)
  real(kind=wp), intent(out) :: T(n,n)
  integer(kind=iwp) :: i, j
  real(kind=wp) :: s1, s2, m1, m2, fact

  do i=1,n
    do j=1,n
      s1 = spins(i)
      s2 = spins(j)
      m1 = projs(i)
      m2 = projs(j)
      fact = sqrt(real((2*k+1),kind=wp))
      if (mod(int(s1-m1),2) == 1) fact = -fact
      T(i,j) = fact*W3J(s1,real(k,kind=wp),s2,-m1,real(q,kind=wp),m2)
    end do
  end do

end subroutine ITO

subroutine WERDM(rho,n_so,n_sf,k,q,spins,projs,so_sf,RED)
  ! calculates the elements of Wigner-Eckart reduced density matrix for fixed values of K, Q

  integer(kind=iwp), intent(in) :: n_so, n_sf, k, q, so_sf(n_so)
  complex(kind=wp), intent(in) :: rho(n_so,n_so)
  real(kind=wp), intent(in) :: spins(n_so), projs(n_so)
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
      RED(ii,jj) = RED(ii,jj)+rho(i,j)*T(i,j)
    end do
  end do

end subroutine WERDM

subroutine WERDM_back(RED,n_so,n_sf,len_sph,k_ranks,q_proj,spins,projs,so_sf,rho_back)
  ! calculates the density matrix in SO basis from the elements of
  ! Wigner-Eckart reduced density matrix, which are stored in 3d-matrix
  ! (len_sph,n_sf,n_sf)

  integer(kind=iwp), intent(in) :: n_so, n_sf, len_sph, k_ranks(len_sph), q_proj(len_sph), so_sf(n_so)
  complex(kind=wp), intent(in) :: RED(len_sph,n_sf,n_sf)
  real(kind=wp), intent(in) :: spins(n_so), projs(n_so)
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
        rho_back(i,j) = rho_back(i,j)+RED(l,ii,jj)*T(i,j)
      end do
    end do
  end do

end subroutine WERDM_back

subroutine WERSO(vso,n_so,n_sf,so_sf,spins,projs,redvso)
  ! calculates matrix elements of Wigner-Eckart reduced spin-orbit Hamiltonian
  ! <iSM|Vso|jS'M'> = sqrt(3) sum_{m=0,+-1} (-1)^{S-M+m} 3j{S1S'-MmM'} <iS||Vso||jS'>

  integer(kind=iwp), intent(in) :: n_so, n_sf, so_sf(n_so)
  complex(kind=wp), intent(in) :: vso(n_so,n_so)
  real(kind=wp), intent(in) :: spins(n_so), projs(n_so)
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
        threejsymb = W3J(s1,One,s2,-m1,real((m-1),kind=wp),m2)
        if (threejsymb /= Zero) then ! this condition should be checked carefully
          !!!!! test
          !if (redvso(ii,jj,m+1) /= cZero) then
          !  write(u6,*) 'Matrix red VSOC element filled twice', redvso(ii,jj,m+1)
          !  write(u6,*) 'ii,jj = ', ii,jj
          !  write(u6,*) 'i,j   = ', i,j
          !  write(u6,*) 'm     = ', m+1
          !  write(u6,*) 's1,s2 = ', s1,s2
          !  write(u6,*) 'now = ', (-1)**(nint(s1-m1+m-1))*vso(i,j)/sqrt(Three)/threejsymb
          !end if
          redvso(ii,jj,m+1) = (-1)**(nint(s1-m1+m-1))*vso(i,j)/sqrt(Three)/threejsymb
        end if
      end do
    end do
  end do

end subroutine WERSO

subroutine WERSO_back(redvso,n_so,n_sf,so_sf,spins,projs,vso)
  ! calculates matrix elements of Wigner-Eckart reduced spin-orbit Hamiltonian
  ! <iSM|Vso|jS'M'> = sqrt(3) sum_{m=0,+-1} (-1)^{S-M+m} 3j{S1S'-MmM'} <iS||Vso||jS'>

  integer(kind=iwp), intent(in) :: n_so, n_sf, so_sf(n_so)
  complex(kind=wp), intent(in) :: redvso(n_sf,n_sf,3)
  real(kind=wp), intent(in) :: spins(n_so), projs(n_so)
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
        vso(i,j) = vso(i,j)+(-1)**(nint(s1-m1+m-1))*sqrt(Three)*W3J(s1,One,s2,-m1,real((m-1),kind=wp),m2)*redvso(ii,jj,m+1)
      end do
    end do
  end do

end subroutine WERSO_back

subroutine print_c_matrix(A,n,header)
  ! prints any complex matrix

  integer(kind=iwp), intent(in) :: n
  complex(kind=wp), intent(in) :: A(n,n)
  character(len=*), intent(in) :: header
  integer(kind=iwp) :: i, j

  call dashes()
  write(u6,*) header
  do i=1,n
    write(u6,*) (A(i,j),j=1,n)
  end do

end subroutine print_c_matrix

subroutine check_hermicity(A,n,A_name,thrs)
  ! check whether matrix A is hermitian
  ! pops a warning only if hermicity is violated

  integer(kind=iwp), intent(in) :: n
  complex(kind=wp), intent(in) :: A(n,n)
  character(len=*), intent(in) :: A_name
  real(kind=wp), intent(in) :: thrs
  real(kind=wp) :: abserror, diff
  integer(kind=iwp) :: i, j

  abserror = Zero
  do i=1,n
    do j=1,i
      diff = abs(real(A(i,j))-real(A(j,i)))
      if ((diff >= thrs) .and. (diff >= abserror)) abserror = diff
      diff = abs(aimag(A(i,j))+aimag(A(j,i)))
      if ((diff >= thrs) .and. (diff >= abserror)) abserror = diff
    end do
  end do
  if (abserror >= thrs) then
    call WarningMessage(1,'Non-hermitian matrix obtained!')
    write(u6,'(a,1x,a,1x,a,1x,g28.16)') 'Matrix',A_name,'Abs Error =',abserror
  end if

end subroutine check_hermicity

function get_kq_order(k_prime,q_prime)
  integer(kind=iwp) :: get_kq_order
  integer(kind=iwp), intent(in) :: k_prime, q_prime
  integer(kind=iwp) :: k, q

  k = 0
  get_kq_order = 0
  do while (k <= k_prime)
    do q=-k,k,1
      get_kq_order = get_kq_order+1
      if (q == q_prime .and. k == k_prime) return
    end do
    k = k+1
  end do

end function get_kq_order

function fct(n)
  ! this function provides correct answer till n=169 only

  real(kind=wp) :: fct
  integer(kind=iwp), intent(in) :: n
  integer(kind=iwp) :: i
  real(kind=wp) :: xct

  xct = One
  fct = One
  if (n < 0) then
    write(u6,'(A,i0)') 'FCT:  N<0 !'
    write(u6,'(A,i0)') 'N = ',N
    write(u6,'(A   )') 'It is an impossible case.'
    fct = -9.0e99_wp
    return
  else if (n == 0) then
    return
  else if (n <= 169) then
    do i=1,n
      xct = xct*real(i,kind=wp)
    end do
  else
    write(u6,'(A,i0)') 'FCT:   N = ',N
    write(u6,'(A)') 'Factorial of N>169 overflows on x86_64'
    write(u6,'(A)') 'Use higher numerical precision, or rethink your algorithm.'
  end if
  fct = xct

end function fct

function W6J(a,b,c,d,e,f)
  ! Calculates Wigner 6j symbol. Arguments a-f are positive integers
  ! and are twice the true value of 6j's arguments, in the form
  ! { a b c }
  ! { d e f }

  real(kind=wp) :: W6j
  integer(kind=iwp), intent(in) :: a, b, c, d, e, f
  integer(kind=iwp) :: n, nlow, nhig
  real(kind=wp) :: isum, rsum

  W6J = Zero
  if (mod(a+b,2) /= mod(c,2)) return
  if (mod(c+d,2) /= mod(e,2)) return
  if (mod(a+e,2) /= mod(f,2)) return
  if (mod(b+d,2) /= mod(f,2)) return
  if ((abs(a-b) > c) .or. (a+b < c)) return
  if ((abs(c-d) > e) .or. (c+d < e)) return
  if ((abs(a-e) > f) .or. (a+e < f)) return
  if ((abs(b-d) > f) .or. (b+d < f)) return
  if (.not. check_triangle(a,b,c)) return
  if (.not. check_triangle(c,d,e)) return
  if (.not. check_triangle(a,e,f)) return
  if (.not. check_triangle(b,d,f)) return
  nlow = 0
  nhig = 0
  nlow = max((a+b+c)/2,(c+d+e)/2,(b+d+f)/2,(a+e+f)/2)
  nhig = min((a+b+d+e)/2,(b+c+e+f)/2,(a+c+d+f)/2)
  rsum = Zero
  do n=nlow,nhig
    isum = real(((-1)**n),kind=wp)*fct(n+1) &
           /fct((a+c+d+f)/2-n) &
           /fct((b+c+e+f)/2-n) &
           /fct(n-(a+b+c)/2) &
           /fct(n-(c+d+e)/2) &
           /fct(n-(a+e+f)/2) &
           /fct(n-(b+d+f)/2) &
           /fct((a+b+d+e)/2-n)
    rsum = rsum+isum
  end do
  W6J = dlt(a,b,c)*dlt(c,d,e)*dlt(a,e,f)*dlt(b,d,f)*rsum

end function W6J

function dlt(a,b,c)
  ! calculates the delta(a,b,c) function using the formula 8.2.1. from:
  ! D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
  ! "Quantum Theory of Angular Momentum", World ScientIfic, 1988.
  ! a,b,c are positive integers, their values are double than their original value

  real(kind=wp) :: dlt
  integer(kind=iwp), intent(in) :: a, b, c

  dlt = Zero
  if ((abs(a-b) > c) .or. (a+b < c)) return
  if ((abs(b-c) > a) .or. (b+c < a)) return
  if ((abs(c-a) > b) .or. (c+a < b)) return
  if (mod((a+b-c),2) == 1) return
  if (mod((a-b+c),2) == 1) return
  if (mod((-a+b+c),2) == 1) return
  if (mod((a+b+c),2) == 1) return
  if (.not. check_triangle(a,b,c)) return
  ! special cases:
  if (a == 0) dlt = One/sqrt(real((b+1),kind=wp))
  if (b == 0) dlt = One/sqrt(real((a+1),kind=wp))
  if (c == 0) dlt = One/sqrt(real((a+1),kind=wp))
  dlt = sqrt(fct((a+b-c)/2)*fct((a-b+c)/2)*fct((-a+b+c)/2)/fct((a+b+c)/2+1))

end function dlt

function check_triangle(a,b,c)
  ! boolean function, checks if arguments fulfill triangular rule

  logical :: check_triangle
  integer(kind=iwp), intent(in) :: a, b, c

  check_triangle = .false.
  if (((a+b) >= c) .and. ((b+c) >= a) .and. ((c+a) >= b)) check_triangle = .true.

end function check_triangle

end module rhodyn_utils
