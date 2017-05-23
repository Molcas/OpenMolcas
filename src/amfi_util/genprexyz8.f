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
      subroutine genprexyz8(preXZ)
      implicit real*8(a-h,o-z)
#include "para.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
      Dimension preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
cbs #####################################################################
cbs   additional (-) signs from the (-i) factors  in the
cbs   (-) linear combinations   (see tosigX(Y,Z).f)
cbs #####################################################################
cbs   - + - -   =>   minus
      do M4=-Lmax,-1
      do M3=-Lmax,-1
         do M2= 0,Lmax
c      do M1=-Lmax,-1
c           preXZ(m1,m2,m3,m4)=-preXZ(m1,m2,m3,m4)
            call dscal_(Lmax,-1d0,preXZ(-Lmax,m2,m3,m4),1)
c      enddo
         enddo
      enddo
      enddo
      return
      end
