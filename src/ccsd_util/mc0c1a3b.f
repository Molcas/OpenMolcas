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
       subroutine mc0c1a3b
     & (rowa,cola,rowb,colb,rowc,colc,
     & row,sum,col,a,b,c)
C
C     C = C + A*B
C
#include "ccsd1.fh"
       integer rowa,cola,rowb,colb,rowc,colc
       integer row,sum,col
       real*8  A(1:rowa,1:cola)
       real*8  B(1:rowb,1:colb)
       real*8  C(1:rowc,1:colc)
C
C      help variables
C
       integer i,j,k
c
       if (mhkey.eq.1) then
c      ESSL
       call DGEMM_('N','N',row,col,sum,1.0d0,a,rowa,b,rowb,
     & 1.0d0,c,rowc)
C
       else
c      Fortran matrix handling
C
       do 50 j=1, col
       do 40 k=1, sum
       do 30 i=1, row
       c(i,j) = c(i,j) + a(i,k)*b(k,j)
 30    continue
 40    continue
 50    continue
C
       end if
c
       return
       end
