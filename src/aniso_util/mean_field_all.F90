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

subroutine mean_field_all(EXCH,N,H,X,Y,Z,zJ,T,W,thrs,dM,SM,ST)
! this Subroutine computes the mean field of neighboring spins ST(3)
! for zJ /= 0.0
! using ONLY Zeeman basis (N)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: EXCH, N
real(kind=wp), intent(in) :: H, X, Y, Z, zJ, T, W(EXCH), thrs
real(kind=wp), intent(out) :: ST(3)
complex(kind=wp), intent(in) :: DM(3,EXCH,EXCH), SM(3,EXCH,EXCH)
integer(kind=iwp) :: i, iter, l
real(kind=wp) :: S(3), SCHK, SL(3), ZB
logical(kind=iwp) :: Conv
real(kind=wp), allocatable :: RWORK(:), WM(:)
complex(kind=wp), allocatable :: DM_TMP(:,:,:), HZEE(:), SM_TMP(:,:,:), SZ(:,:,:), TMP(:,:), W_c(:), WORK(:), ZM(:,:)
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

call mma_allocate(WM,EXCH,label='WM')
call mma_allocate(SZ,3,EXCH,EXCH,label='SZ')
call mma_allocate(ZM,N,N,label='ZM')

! temporary arrays used in ZEEM_SA:
call mma_allocate(RWORK,(3*N-2),'ZEEM_RWORK')
call mma_allocate(HZEE,nTri_Elem(N),'ZEEM_HZEE')
call mma_allocate(WORK,2*N-1,'ZEEM_WORK')
call mma_allocate(W_c,N,'ZEEM_W_c')
if (N == EXCH) then
  call mma_allocate(TMP,EXCH,EXCH,label='TMP')
else
  call mma_allocate(DM_TMP,3,N,N,label='DM_TMP')
  call mma_allocate(SM_TMP,3,N,N,label='SM_TMP')
end if

! zero everything:
RWORK(:) = Zero
HZEE(:) = cZero
WORK(:) = cZero
W_c(:) = cZero
! determine first the average spin of neighboring
! molecules for each temperature point
Conv = .false.
do iter=1,mxIter
  WM(1:N) = Zero
  ZM(:,:) = cZero
  ! build and diagonalize the Zeeman Hamiltonian (size N x N)
  ! for the field direction (X,Y,Z) and strength (H)
  if (N == EXCH) then
    call ZEEM_SA(N,H,X,Y,Z,W,dM,sM,ST,zJ,WM,ZM,_DBG_,RWORK,HZEE,WORK,W_c)
  else
    DM_TMP(:,:,:) = dM(:,1:N,1:N)
    SM_TMP(:,:,:) = sM(:,1:N,1:N)
    call ZEEM_SA(N,H,X,Y,Z,W(1:N),DM_TMP,SM_TMP,ST,zJ,WM(1:N),ZM,_DBG_,RWORK,HZEE,WORK,W_c)
  end if
  WM(N+1:) = W(N+1:)

  ! transform the spin momenta to the Zeeman eigenstate basis
  SZ(:,:,:) = cZero
  call UTMU(EXCH,N,ZM,SM,SZ)
  ! compute the spin magnetization vector at this temperature (T):
  if (iter == mxIter) SL(:) = S(:)

  ZB = Zero
  S(:) = Zero
  if (N == EXCH) then
    do L=1,3
      TMP(:,:) = SZ(l,:,:)
      call calcmagn1(EXCH,WM,TMP,T,S(l),ZB)
    end do
  else
    do L=1,3
      call calcmagn2(EXCH,N,WM,T,H,SZ,X,Y,Z,L,S(l),ZB)
    end do
  end if

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
call mma_deallocate(RWORK)
call mma_deallocate(HZEE)
call mma_deallocate(WORK)
call mma_deallocate(W_c)
if (N == EXCH) then
  call mma_deallocate(TMP)
else
  call mma_deallocate(DM_TMP)
  call mma_deallocate(SM_TMP)
end if

return

end subroutine mean_field_all
