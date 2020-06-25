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
      subroutine genprexyz11(preY)
      implicit real*8(a-h,o-z)
#include "para.fh"
#include "Molcas.fh"
      Dimension preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
cbs #####################################################################
cbs   additional (-) signs from the (-i) factors  in the
cbs   (-) linear combinations   (see tosigX(Y,Z).f)
cbs #####################################################################
cbs   + +  - -  >   -
      do M4=-Lmax,-1
      do M3=-Lmax,-1
      do M2=0,Lmax
c     do M1=0,Lmax
            call dscal_(Lmax+1,-1d0,preY(0,m2,m3,m4),1)
c      preY(m1,m2,m3,m4)=-preY(m1,m2,m3,m4)
c      enddo
      enddo
      enddo
      enddo
      return
      end
