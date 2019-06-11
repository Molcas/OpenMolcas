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
C
C     Inverts a square matrix
C
      Subroutine MatInvert(A,n)
      Implicit None
#include "stdalloc.fh"
      Real*8 :: A(*)
      Integer :: n,err,nw
      Integer, Dimension(:), Allocatable :: ipiv
      Real*8, Dimension(:), Allocatable :: wrk
      Real*8 :: dum(1)
      Call mma_allocate(ipiv,n)
      Call dGeTRF_(n,n,A,n,ipiv,       err)
      Call dGeTRI_(n,  A,n,ipiv,dum,-1,err)
      nw=Int(dum(1))
      Call mma_allocate(wrk,nw)
      Call dGeTRI_(n,  A,n,ipiv,wrk,nw,err)
      Call mma_deallocate(ipiv)
      Call mma_deallocate(wrk)
      Return
      End
