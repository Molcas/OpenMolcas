************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine DIAG_R2(MATRIX,N,INFO,W,Z)
C
C   THIS ROUTINE PERFORMS THE DIAGONALIZATION OF A Real SQUARE
C   MATRIX WITH THE DIMENSION NBTOT. THE EIGENVALUES OF THE DIAGONALIZATION
C   ARE DIRECTED INTO W1 AND THE Real EIGENVECTORS ARE WRITTEN TO Z1.
C

      Implicit None
#include "stdalloc.fh"
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)       :: N
      Integer, intent(out)      :: INFO
      Real(kind=wp), intent(in) :: MATRIX(N,N)
      Real(kind=wp), intent(out):: W(N), Z(N,N)
      ! local variables:
      Integer       :: I,J
      Real(kind=wp), allocatable :: AP(:)   !(N*(N+1)/2)
      Real(kind=wp), allocatable :: WORK(:) !(3*N)
      Real(kind=wp), allocatable :: W1(:)   !(N)
      Real(kind=wp), allocatable :: Z1(:,:) !(N,N)

      Call qEnter('diag_r2')
C initializations
      INFO=0
      If(N<1) Return

      Call dcopy_(        N,[0.0_wp],0,  W ,1)
      Call dcopy_(      N*N,[0.0_wp],0,  Z ,1)

      Call mma_allocate(AP,(N*(N+1)/2), 'AP')
      Call mma_allocate(WORK,3*N,'WORK')
      Call mma_allocate(W1,N,'W1')
      Call mma_allocate(Z1,N,N,'Z1')
      Call dcopy_(N*(N+1)/2, [0.0_wp],0,  AP,1)
      Call dcopy_(      3*N, [0.0_wp],0,WORK,1)
      Call dcopy_(        N, [0.0_wp],0,  W1,1)
      Call dcopy_(      N*N, [0.0_wp],0,  Z1,1)

      Do j=1,N
        Do i=1,j
          AP(i+(j-1)*j/2)=MATRIX(i,j)
        End Do
      End Do
      ! diagonalize:
      Call DSPEV_('V','U',N,AP,W1,Z1,N,WORK,INFO)

      Call dcopy_(  N,W1,1,W,1)
      Call dcopy_(N*N,Z1,1,Z,1)

      Call mma_deallocate(AP)
      Call mma_deallocate(WORK)
      Call mma_deallocate(W1)
      Call mma_deallocate(Z1)
      Call qExit('diag_r2')
      Return
      End Subroutine DIAG_R2
