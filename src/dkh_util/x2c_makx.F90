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
!
!----------------------------------------------------------------------|
!
      subroutine x2c_makx(m,n,f,s,x)
!
! Make X matrix from m-dimensional (m=2n) Fork(f) and Overlap(s) matrix
!
!   X is the relation(transfer) matrix of Large--Small component coefficients
!   of electron solutions (positive energy solutions)
!
      implicit none
#include "WrkSpc.fh"
! Input
      integer m,n
      Real*8 f(m,m),s(m,m)
! Output
      Real*8 x(n,n)
! Temp
      integer i,j,k
      integer lwork,info,itmp,iw,itF,itS
!
      lwork = 8*m
      call getmem('TmpF ','ALLOC','REAL', itF, m*m+4 )
      call getmem('TmpS ','ALLOC','REAL', itS, m*m+4 )
      call getmem('Eig  ','ALLOC','REAL',  iw, m+4 )
      call getmem('Work ','ALLOC','REAL', itmp, lwork+4 )
!
! Copy Fock and Overlap matrix to temp arrays
!
      k = 0
      do i=1,m
        do j=1,m
          Work(itF+k) = f(j,i)
          Work(itS+k) = s(j,i)
          k = k + 1
        end do
      end do
!
! Diagonalization of Fock matrix with given overlap matrix
!
      call dsygv_(1,'V','L',m,Work(itF),m,Work(itS),m,                  &
     &           Work(iw),Work(itmp),lwork,info)
!
! Calculate the X matrix from electron solutions
!
      k = 0
      do i=1,n
        do j=1,n
!          ! store large component coefficients of electron solutions (matrix A)
          Work(itF+k) = Work(itF-1+j+(i+n-1)*m)
!          ! store small component coefficients of electron solutions (matrix B)
          Work(itS+k) = Work(itF-1+j+n+(i+n-1)*m)
          k = k + 1
        end do
      end do
!     ! compute X=BA^{-1}
      call XDR_dmatinv(Work(itF),n)
      call dmxma(n,'N','N',Work(itS),Work(itF),x,1.d0)
!
! Free temp memories
!
      call getmem('TmpF ','FREE','REAL', itF, m*m+4 )
      call getmem('TmpS ','FREE','REAL', itS, m*m+4 )
      call getmem('Eig  ','FREE','REAL',  iw, m+4 )
      call getmem('Work ','FREE','REAL', itmp, lwork+4 )
      return
      end
