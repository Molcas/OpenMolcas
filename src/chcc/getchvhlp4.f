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
        subroutine GetChVHlp4 (H,T1,dima,adda)
c
c       this routine do:
c       Extract H(i,a') <- T1(a,i) for given aGrp
c
c       parameter description:
c       H       - Output file (O)
c       T1      - T1 amplitudes (I)
c       dima    - dimension of given Group (I)
c       adda    - shift of a' in full a set (I)
c
c       N.B. Kvajt odflaknute, je to ExtT1 len transponovane
c
        implicit none
#include "chcc1.fh"
        integer dima,adda
        real*8 T1(1:nv,1:no)
        real*8 H(1:no,1:dima)
c
c       help variables
        integer a,aa,i
c
        do a=1,dima
        aa=adda+a
          do i=1,no
            H(i,a)=T1(aa,i)
          end do
        end do
c
        return
        end
