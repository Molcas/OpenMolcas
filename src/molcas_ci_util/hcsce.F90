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

implicit real*8(A-H,O-Z)
dimension H(N*(N+1)/2), S(N*(N+1)/2), C(N,N), E(N)
dimension WGronk(2)
#include "WrkSpc.fh"
#include "timers.fh"
!character*12 method

call Timing(Longines_1,Swatch,Swatch,Swatch)

! PAM 2009: On input, M=max possible orthonormal solutions to HC=SCE
! Save it.
MMAX = M

! allocate temporary work space
call GetMem('Temp1','Allo','Real',lw1,N*N)
call GetMem('Temp2','Allo','Real',lw2,N*N)
call GetMem('Temp3','Allo','Real',lw3,N*N)
call GetMem('Temp4','Allo','Real',lw4,N)

! make local copies of H and S
call Square(S,Work(lw1),1,N,N)
call Square(H,Work(lw2),1,N,N)

! Schmidt orthogonalization
call dcopy_(N*N,[0.0d0],0,C,1)
call dcopy_(N,[1.0d0],0,C,N+1)

!write(6,*) ' HCSCE calling Schmidt.'
!call Schmidt(N,Work(lw1),C,Work(lw4),M)
!write(6,*) ' HCSCE back from Schmidt. M=',M
! PAM 2009: It seems that no provision is made for the case that
! the returned M, = nr of ON vectors produced, is smaller than N?!
! Also, the whole thing looks very inefficient. But I just make
! some provisional changes now (090216).
!write(6,*) ' HCSCE check eigenvalues. N=',N
!call eigv(N,Work(lw1))
!write(6,*) ' HCSCE calling NewGS.'
call NewGS(N,Work(lw1),C,Work(lw4),M)
!write(6,*) ' HCSCE back from NewGS. M=',M
! Possibly in very difficult cases, NewGS produced too many vectors:
M = min(M,MMAX)

! transform H to an orthogonal basis
! PAM 2009: Rewritten, use only M orthogonal vectors
call DGEMM_('N','N',N,M,N,1.0d0,Work(lw2),N,C,N,0.0d0,Work(lw3),N)
call DGEMM_('T','N',M,M,N,1.0d0,C,N,Work(lw3),N,0.0d0,Work(lw2),M)

! PAM 2009: Replace by DSYEV call.
!method = 'Householder'
!method = 'Jacobi'

! diagonalize and extract eigenvalues
!if (method == 'Jacobi') then
!  do i=1,N
!    do j=1,i
!      Work(j+i*(i-1)/2+lw2-1) = Work(j+(i-1)*N+lw2-1)
!    end do
!  end do
!  call Jacob(Work(lw2),C,N,N)
!  call JacOrd(Work(lw2),C,N,N)
!  do i=1,N
!    E(i) = Work(i*(i+1)/2+lw2-1)
!  end do
!else if (method == 'Householder') then
!  call Eigen_Molcas(N,Work(lw2),E,Work(lw4))
!  call DGEMM_('N','N',N,N,N,1.0d0,C,N,Work(lw2),N,0.0d0,Work(lw3),N)
!  call dcopy_(N*N,Work(lw3),1,C,1)
!end if
! PAM 2009 Equivalent, DSYEV, note now use just M, not all N:
INFO = 0
call dsyev_('V','L',M,Work(lw2),M,E,WGRONK,-1,INFO)
NSCRATCH = int(WGRONK(1))
call GETMEM('SCRATCH','ALLO','REAL',LSCRATCH,NSCRATCH)
call dsyev_('V','L',M,WORK(lw2),M,E,WORK(LSCRATCH),NSCRATCH,INFO)
call GETMEM('SCRATCH','FREE','REAL',LSCRATCH,NSCRATCH)
call DGEMM_('N','N',N,M,M,1.0d0,C,N,Work(lw2),M,0.0d0,Work(lw3),N)
call dcopy_(N*M,Work(lw3),1,C,1)

! deallocate temporary work space
call GetMem('Temp4','Free','Real',lw4,N)
call GetMem('Temp3','Free','Real',lw3,N*N)
call GetMem('Temp2','Free','Real',lw2,N*N)
call GetMem('Temp1','Free','Real',lw1,N*N)

call Timing(Longines_2,Swatch,Swatch,Swatch)
Longines_2 = Longines_2-Longines_1
Longines_3 = Longines_3+Longines_2

return

end subroutine HCSCE
