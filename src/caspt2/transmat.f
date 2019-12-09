************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2019, Stefano Battaglia                                *
************************************************************************
      subroutine transmat(A,U,N)
      implicit none
#include "stdalloc.fh"
***
* This subroutine carries out the following transformation:
*       U^T * A * U
* where A and U are NxN matrices
***
      integer N
      real(8) A(N,N)
      real(8) U(N,N)
      real(8),allocatable :: B(:,:)

* Allocate temporary array B
      call mma_allocate(B,N,N)
      B=0.0d0

* B = U^T * A
      call dgemm_('T', 'N', N, N, N, 1.0d0, U, N, A, N, 0.0d0, B, N)
* A = B * U
      call dgemm_('N', 'N', N, N, N, 1.0d0, B, N, U, N, 0.0d0, A, N)

      call mma_deallocate(B)

      return
      end
