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
* This subroutine diagonalizes a real symmetric matrix A using
* the Jacobi algorithm
**
      subroutine eigen(A,U,N)
      implicit real(8) (A-H,O-Z)
#include "WrkSpc.fh"

      real(8) A(N,N)
      real(8) U(N,N)

      NSCR=(N*(N+1))/2
      call getmem('SCR','ALLO','REAL',LSCR,NSCR)
      ! call getmem('EVEC','ALLO','REAL',LEVEC,Nstate**2)

      IJ=0
      do I=1,N
        do J=1,I
          IJ=IJ+1
          WORK(LSCR+IJ-1)=A(I,J)
        end do
      end do

* Initialize U as the identity matrix
      U=0.0d0
      call dcopy_(N,[1.0D0],0,U,N+1)

* Call Jacobi algorithm
      call JACOB(WORK(LSCR),U,N,N)

      call getmem('SCR','FREE','REAL',LSCR,N)

      return
      end
