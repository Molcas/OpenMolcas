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

subroutine mean_field_exch(N,H,X,Y,Z,zJ,T,W,thrs,dM,sM,ST)
! this Subroutine computes the mean field of neighboring spins ST(3)
! for zJ /= 0.0
! using ONLY Zeeman basis (N)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: H, X, Y, Z, zJ, T, W(N), thrs
complex(kind=wp), intent(in) :: DM(3,N,N), SM(3,N,N)
real(kind=wp), intent(out) :: ST(3)
integer(kind=iwp) :: i, iter, l
real(kind=wp) :: S(3), SCHK, SL(3), ZB
logical(kind=iwp) :: Conv
real(kind=wp), allocatable :: RWORK(:), WM(:)
complex(kind=wp), allocatable :: HZEE(:), SZ(:,:,:), TMP(:,:), W_c(:), WORK(:), ZM(:,:)
integer(kind=iwp), parameter :: mxIter = 100
real(kind=wp), parameter :: THRS2 = 1.0e-12_wp ! FIXME: overriding input thrs
#ifdef _DEBUGPRINT_
# define _DBG_ .true.
#else
# define _DBG_ .false.
#endif

#include "macros.fh"
unused_var(thrs)

SCHK = Zero
S(:) = Zero
ST(:) = Zero

call mma_allocate(WM,N,label='WM')
call mma_allocate(SZ,3,N,N,label='SZ')
call mma_allocate(ZM,N,N,label='ZM')
call mma_allocate(TMP,N,N,label='TMP')

! temporary arrays used in ZEEM_SA:
call mma_allocate(RWORK,3*N-2,'ZEEM_RWORK')
call mma_allocate(HZEE,nTri_Elem(N),'ZEEM_HZEE')
call mma_allocate(WORK,2*N-1,'ZEEM_WORK')
call mma_allocate(W_c,N,'ZEEM_W_c')

! zero everything:
RWORK(:) = Zero
HZEE(:) = cZero
WORK(:) = cZero
W_c(:) = cZero
! determine first the average spin of neighboring
! molecules for each temperature point
Conv = .false.
do iter=1,mxIter
  WM(:) = Zero
  ZM(:,:) = cZero
  ! build and diagonalize the Zeeman Hamiltonian (size N x N)
  ! for the field direction (X,Y,Z) and strength (H)
  call ZEEM_SA(N,H,X,Y,Z,W,dM,sM,ST,zJ,WM,ZM,_DBG_,RWORK,HZEE,WORK,W_c)
  !write(u6,'(A,3ES16.8)') 'WM:',(WM(l),l=1,N)
  !write(u6,'(A,9(2ES16.8,2x))') 'ZM:',((ZM(l,i),l=1,N),i=1,N)
  ! transform the spin momenta to the Zeeman eigenstate basis
  SZ(:,:,:) = cZero
  call UTMU(N,N,ZM,SM,SZ)
  ! compute the spin magnetization vector at this temperature (T):
  if (iter == mxIter) SL(:) = S(:)
  S(:) = Zero
  do l=1,3
    ZB = Zero
    TMP(:,:) = SZ(l,:,:)
    call calcmagn1(N,WM,TMP,T,S(l),ZB)
  end do

  ! check if average spin is converged
  SCHK = sum(abs(S(:)-ST(:)))

# ifdef _DEBUGPRINT_
  write(u6,'(A,i4,1x,A,3ES20.10,2x,A,3ES20.10)') 'ST:   End of iteration',iter,':',(ST(l),l=1,3),'DIFF:',(S(l)-ST(l),l=1,3)
  write(u6,'(A,i4,1x,A,3ES20.10,2x,A,3ES20.10)') 'SCHK: End of iteration',iter,':',SCHK
# endif

  ! decide to continue the iterative process or exit
  if (SCHK < THRS2) then
    Conv = .true.
    exit
  end if
  ST(:) = S(:)
end do ! iter

if (.not. Conv) then
  write(u6,'(A, ES24.14)') 'This message shows that the average spin did NOT converge after 100 iterations. Temp.(in K)=',T
  write(u6,'(A,4ES24.14)') 'Field: (X, Y, Z), and (H):',X,Y,Z,H
  write(u6,'(A,4ES24.14)') 'Last values of the average spin: (Sx,Sy,Sz):',(ST(i),i=1,3)
  write(u6,'(A,4ES24.14)') 'Last values of the deviation:              :',(SL(i)-ST(i),i=1,3)
  write(u6,'(A,4ES24.14)') 'Absolute value of the deviation:           :',SCHK
  write(u6,'(A,4ES24.14)') 'Convergence threshold:    THRS =           :',THRS2
  write(u6,'(A         )') 'The program will continue, using the last value of the average spin'
end if

! deallocate temporary data:
call mma_deallocate(WM)
call mma_deallocate(SZ)
call mma_deallocate(ZM)
call mma_deallocate(TMP)
call mma_deallocate(RWORK)
call mma_deallocate(HZEE)
call mma_deallocate(WORK)
call mma_deallocate(W_c)

return

end subroutine mean_field_exch
