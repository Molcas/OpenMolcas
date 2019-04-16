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
      subroutine genprexyz14(icheckz,interxyz)
      implicit real*8(a-h,o-z)
#include "para.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
      integer mcheckz
      dimension icheckz(0:Lmax,0:Lmax,0:Lmax,0:Lmax),
     *interxyz(16,0:Lmax,0:Lmax,0:Lmax,0:Lmax)
cbs #####################################################################
cbs   some quick decision for interaction
cbs #####################################################################
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=0,Lmax
      do M1=0,Lmax
              icheckz(m1,m2,m3,m4)=mcheckz(m1,m2,m3,m4)
      enddo
      enddo
      enddo
      enddo
cbs #####################################################################
cbs   there are at most 16 possible combinations of signs ( 2**4)
cbs #####################################################################
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=0,Lmax
      do M1=0,Lmax
      do irun=1,16
      interxyz(irun,m1,m2,m3,m4)=0
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
