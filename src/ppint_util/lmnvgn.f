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
      subroutine lmnvgn_molcas(lmn1u,lmnv)
c     # generate cartesian gaussian exponent array.
c     # lmnv(*,*) = exponents of the cartesian gaussian basis functions.
c     #             s  p  d   f   g   h   i
c     #       lmn = 0  1  2   3   4   5   6
c     #    numxyz = 1, 3, 6, 10, 15  21  28      =((lmn+1)*(lmn+2))/2
      implicit logical (a-z)
      integer lmn, lmn1u, lmnv(3,*), ndx
      Integer ix,iy,iz
c
      ndx=0
      Do lmn = 0,lmn1u-1
         Do ix = lmn, 0 , -1
            Do iy = lmn-ix, 0, -1
               iz = lmn-ix-iy
               ndx=ndx+1
               lmnv(1,ndx)=ix
               lmnv(2,ndx)=iy
               lmnv(3,ndx)=iz
            End Do
         End Do
      End Do
*
      return
      end
