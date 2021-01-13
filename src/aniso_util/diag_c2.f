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
      Subroutine diag_c2(matrix,n,info,w,z)
c
c   this routine performs the diagonalization of a Complex square
c   matrix with the dimension nbtot. the eigenvalues of the diagonalization
c   are directed into w1 and the Complex eigenvectors are written to z1.
c

      Implicit None
#include "stdalloc.fh"
      Integer, parameter          ::  wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer        :: info,i,j,n
      Complex(kind=8),intent(in)  :: matrix(n,n)
      Complex(kind=8),intent(out) :: z(n,n)
      Real(kind=8), intent(out)   :: w(n)
      ! local variables:
      Real(kind=8), allocatable    :: rwork(:) !rwork(3*n-2)
      Real(kind=8), allocatable    :: w1(:)    !w1(n)
      Complex(kind=8), allocatable :: ap(:)    !ap(n*(n+1)/2)
      Complex(kind=8), allocatable :: work(:)  !work(2*n-1)
      Complex(kind=8), allocatable :: z1(:,:)  !z1(n,n)
      Real(kind=8), external       :: dznrm2_
      Real(kind=8)                 :: RM

      info=0
      Call zcopy_(       N*N,[(0.0_wp,0.0_wp)],0,   Z ,1)
      Call dcopy_(         N,        [0.0_wp] ,0,   W ,1)

      RM=0.0_wp
      RM=dznrm2_(n*n,matrix,1)

      If(RM>0.0_wp) Then
         Call mma_allocate(ap   ,(n*(n+1)/2),'ap'   )
         Call mma_allocate(work ,    (2*n-1),'work' )
         Call mma_allocate(z1   ,        n,n,'work' )
         Call zcopy_( N*(N+1)/2,[(0.0_wp,0.0_wp)],0,   AP,1)
         Call zcopy_(     2*N-1,[(0.0_wp,0.0_wp)],0, work,1)
         Call zcopy_(       N*N,[(0.0_wp,0.0_wp)],0,   Z1,1)

         Call mma_allocate(rwork,(3*n-2)    ,'rwork')
         Call mma_allocate(w1   ,          n,'w1'   )
         Call dcopy_(     3*N-2,        [0.0_wp] ,0,rwork,1)
         Call dcopy_(         N,        [0.0_wp] ,0,   W1,1)

         Do j=1,n
            Do i=1,j
              ap(i+(j-1)*j/2)=matrix(i,j)
            End Do
         End Do
         ! diagonalize:
         Call zhpev_('v','u',n,ap,w1,z1,n,work,rwork,info)
         !save results:
         Call dcopy_(   N, W1,1, W,1)
         Call zcopy_( N*N, Z1,1, Z,1)

         Call mma_deallocate(rwork)
         Call mma_deallocate(w1)
         Call mma_deallocate(ap)
         Call mma_deallocate(work)
         Call mma_deallocate(z1)
      Else
         ! return dummy results:
         Do i=1,n
            w(i)  = 0.0_wp
            z(i,i)=(1.0_wp,0.0_wp)
         End Do
      End If
      Return
      End
