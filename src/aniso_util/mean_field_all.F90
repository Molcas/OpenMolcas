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

subroutine mean_field_all(EXCH,N,H,X,Y,Z,zJ,T,thrs,W,dM,SM,ST)
! this Subroutine computes the mean field of neighboring spins ST(3)
! for zJ /= 0.0_wp
! using ONLY Zeeman basis (N)

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: EXCH, N
real(kind=8), intent(in) :: H, X, Y, Z, zJ, T, W(EXCH)
complex(kind=8), intent(in) :: DM(3,EXCH,EXCH), SM(3,EXCH,EXCH)
! output
real(kind=8), intent(out) :: ST(3)
#include "stdalloc.fh"
! local variables:
integer :: i, l, mxIter, iter
logical :: DBG
real(kind=8) :: WM(EXCH), S(3), ZB, SCHK, THRS, SL(3)
complex(kind=8) :: ZM(N,N), SZ(3,EXCH,EXCH)
real(kind=8), allocatable :: RWORK(:)
complex(kind=8), allocatable :: HZEE(:), WORK(:), W_c(:)

DBG = .false.

MxIter = 100
THRS = 1.0D-12
SCHK = 0.0_wp
S = 0.0_wp
ST = 0.0_wp

! temporary arrays used in ZEEM_SA:
call mma_allocate(RWORK,(3*N-2),'ZEEM_RWORK')
call mma_allocate(HZEE,(N*(N+1)/2),'ZEEM_HZEE')
call mma_allocate(WORK,(2*N-1),'ZEEM_WORK')
call mma_allocate(W_c,N,'ZEEM_W_c')

! zero everything:
call dcopy_(3*N-2,[0.0_wp],0,RWORK,1)
call zcopy_(N*(N+1)/2,[(0.0_wp,0.0_wp)],0,HZEE,1)
call zcopy_(2*N-1,[(0.0_wp,0.0_wp)],0,WORK,1)
call zcopy_(N,[(0.0_wp,0.0_wp)],0,W_c,1)
! determine first the average spin of neighboring
! molecules for each temperature point
do iter=1,mxIter
  WM(1:N) = 0.0_wp
  ZM(1:N,1:N) = (0.0_wp,0.0_wp)
  ! build and diagonalize the Zeeman Hamiltonian (size N x N)
  ! for the field direction (X,Y,Z) and strength (H)
  call ZEEM_SA(N,H,X,Y,Z,W(1:N),dM(1:3,1:N,1:N),sM(1:3,1:N,1:N),ST(1:3),zJ,WM(1:N),ZM(1:N,1:N),DBG,RWORK,HZEE,WORK,W_c)
  if (N /= EXCH) then
    do i=N+1,EXCH
      WM(i) = W(i)
    end do
  end if

  ! transform the spin momenta to the Zeeman eigenstate basis
  call zcopy_(3*EXCH*EXCH,[(0.0_wp,0.0_wp)],0,SZ,1)
  call UTMU(EXCH,N,ZM(1:N,1:N),SM,SZ)
  ! compute the spin magnetization vector at this temperature (T):
  if (iter == mxIter) then
    SL = 0.0_wp
    SL(1) = S(1)
    SL(2) = S(2)
    SL(3) = S(3)
  end if

  ZB = 0.0_wp
  S(1:3) = 0.0_wp
  if (N == EXCH) then
    do L=1,3
      call calcmagn1(EXCH,WM,SZ(l,1:EXCH,1:EXCH),T,S(l),ZB)
    end do
  else
    do L=1,3
      call calcmagn2(EXCH,N,WM,T,H,SZ,X,Y,Z,L,S(l),ZB)
    end do
  end if

  ! check if average spin is converged
  SCHK = 0.0_wp
  do L=1,3
    SCHK = SCHK+sqrt((S(L)-ST(L))*(S(L)-ST(L)))
  end do

  if (DBG) then
    write(6,'(A,i4,1x,A,3ES20.10,2x,A,3ES20.10)') 'ST:   End of iteration',iter,':',(ST(l),l=1,3),'DIFF:',(S(l)-ST(l),l=1,3)
    write(6,'(A,i4,1x,A,3ES20.10,2x,A,3ES20.10)') 'SCHK: End of iteration',iter,':',SCHK
  end if

  ! decide to continue the iterative process or exit
  if (SCHK < THRS) then
    go to 1001
  else
    do l=1,3
      ST(l) = 0.0_wp
      ST(l) = S(l)
    end do
  end if
end do ! iter

write(6,'(A, ES24.14)') 'This message shows that the average spin did NOT converge after 100 iterations. Temp.(in K)=',T
write(6,'(A,4ES24.14)') 'Field: (X, Y, Z), and (H):',X,Y,Z,H
write(6,'(A,4ES24.14)') 'Last values of the average spin: (Sx,Sy,Sz):',(ST(i),i=1,3)
write(6,'(A,4ES24.14)') 'Last values of the deviation:              :',(SL(i)-ST(i),i=1,3)
write(6,'(A,4ES24.14)') 'Absolute value of the deviation:           :',SCHK
write(6,'(A,4ES24.14)') 'Convergence threshold:    THRS =           :',THRS
write(6,'(A         )') 'The program will continue, using the last value of the average spin'

return

1001 continue

! deallocate temporary data:
call mma_deallocate(RWORK)
call mma_deallocate(HZEE)
call mma_deallocate(WORK)
call mma_deallocate(W_c)

return

end subroutine mean_field_all
