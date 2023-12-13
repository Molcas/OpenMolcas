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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!***********************************************************************

subroutine HCSCE(N,H,S,C,E,M)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Solve the secular equations HC=SCE. This routine is part of the  *
!     Davidson diagonalization procedure used to find the lowest roots *
!     of the CI-Hamiltonian. Because of that it is crucial to use      *
!     Schmidt orthogonalization such that the latest vector (in chro-  *
!     nological order) remains unchanged and the previous are ortho-   *
!     gonalized relative to it. It is also important that the input    *
!     data remain unchanged.                                           *
!                                                                      *
!     calling arguments:                                               *
!     N       : Type integer, input.                                   *
!               Dimensions of the secular equations.                   *
!     H       : Type double precision real, input.                     *
!               Hamiltonian in matrix representation.                  *
!     S       : Type double precision real, input.                     *
!               Overlap matrix.                                        *
!     C       : Type double precision real, output.                    *
!               Matrix containing the eigenvectors.                    *
!     E       : Type double precision real, output.                    *
!               Vector containing the eigenvalues.                     *
!     M       : Type integer, output.                                  *
!               This is the number of linearly independent basis       *
!               vectors that span space the metric given by the        *
!               overlap matrix.                                        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher, University of Lund, Sweden, 1993                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: H(N*(N+1)/2), S(N*(N+1)/2)
real(kind=wp), intent(out) :: C(N,N), E(N)
integer(kind=iwp), intent(inout) :: M
integer(kind=iwp) :: INFO, MMAX, NSCRATCH
real(kind=wp) :: dum1, dum2, dum3, WGronk(2)
real(kind=wp), allocatable :: Scratch(:), Temp1(:,:), Temp2(:,:), Temp3(:,:), Temp4(:)
!character(len=12) :: method
#include "WrkSpc.fh"
#include "timers.fh"

call Timing(Longines_1,dum1,dum2,dum3)

! PAM 2009: On input, M=max possible orthonormal solutions to HC=SCE
! Save it.
MMAX = M

! allocate temporary work space
call mma_allocate(Temp1,N,N,label='Temp1')
call mma_allocate(Temp2,N,N,label='Temp2')
call mma_allocate(Temp3,N,N,label='Temp3')
call mma_allocate(Temp4,N,label='Temp4')

! make local copies of H and S
call Square(S,Temp1,1,N,N)
call Square(H,Temp2,1,N,N)

! Schmidt orthogonalization
call unitmat(C,N)

!write(u6,*) ' HCSCE calling Schmidt.'
!call Schmidt(N,Temp1,C,Temp4,M)
!write(u6,*) ' HCSCE back from Schmidt. M=',M
! PAM 2009: It seems that no provision is made for the case that
! the returned M, = nr of ON vectors produced, is smaller than N?!
! Also, the whole thing looks very inefficient. But I just make
! some provisional changes now (090216).
!write(u6,*) ' HCSCE check eigenvalues. N=',N
!call eigv(N,Temp1)
!write(u6,*) ' HCSCE calling NewGS.'
call NewGS(N,Temp1,C,Temp4,M)
call mma_deallocate(Temp1)
!write(u6,*) ' HCSCE back from NewGS. M=',M
! Possibly in very difficult cases, NewGS produced too many vectors:
M = min(M,MMAX)

! transform H to an orthogonal basis
! PAM 2009: Rewritten, use only M orthogonal vectors
call DGEMM_('N','N',N,M,N,One,Temp2,N,C,N,Zero,Temp3,N)
call DGEMM_('T','N',M,M,N,One,C,N,Temp3,N,Zero,Temp2,M)

! PAM 2009: Replace by DSYEV call.
!method = 'Householder'
!method = 'Jacobi'

! diagonalize and extract eigenvalues
!if (method == 'Jacobi') then
!  call mma_allocate(Temp1,N*(N+1)/2,label='Temp1')
!  do i=1,N
!    do j=1,i
!      Temp1(j+i*(i-1)/2) = Temp2(j,i)
!    end do
!  end do
!  call Jacob(Temp1,C,N,N)
!  call JacOrd(Temp1,C,N,N)
!  do i=1,N
!    E(i) = Temp1(i*(i+1)/2)
!  end do
!  call mma_deallocate(Temp1)
!else if (method == 'Householder') then
!  call Eigen_Molcas(N,Temp2,E,Temp4)
!  call DGEMM_('N','N',N,N,N,One,C,N,Temp2,N,Zero,Temp3,N)
!  call C(:,:) = Temp3(:,:)
!end if
! PAM 2009 Equivalent, DSYEV, note now use just M, not all N:
INFO = 0
call dsyev_('V','L',M,Temp2,M,E,WGRONK,-1,INFO)
NSCRATCH = int(WGRONK(1))
call mma_allocate(Scratch,NSCRATCH,label='SCRATCH')
call dsyev_('V','L',M,Temp2,M,E,Scratch,NSCRATCH,INFO)
call mma_deallocate(Scratch)
call DGEMM_('N','N',N,M,M,One,C,N,Temp2,M,Zero,Temp3,N)
C(:,1:M) = Temp3(:,1:M)

! deallocate temporary work space
call mma_deallocate(Temp2)
call mma_deallocate(Temp3)
call mma_deallocate(Temp4)

call Timing(Longines_2,dum1,dum2,dum3)
Longines_2 = Longines_2-Longines_1
Longines_3 = Longines_3+Longines_2

return

end subroutine HCSCE
