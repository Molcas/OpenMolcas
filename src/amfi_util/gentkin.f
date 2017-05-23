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
      subroutine gentkin(L,TKIN,nprims,exponents,rootOVLPinv)
      implicit real*8 (a-h,o-z)
#include "para.fh"
cbs   subroutine to generate the kinetic energy
      dimension TKIN(nprims,nprims),exponents(*),
     *dummy(MxprimL,MxprimL),dummy2(MxprimL,MxprimL),
     *rootOVLPinv(MxprimL,MxprimL)
cbs   one triangular part of the matrix
      do irun2=1,nprims
      do irun1=1,irun2
              dummy(irun1,irun2)=
     *        Tkinet(l,exponents(irun1),
     *        exponents(irun2))
      enddo
      enddo
cbs   copy to the other triangular part....
      do irun2=1,nprims-1
      do irun1=irun2+1,nprims
              dummy(irun1,irun2)=dummy(irun2,irun1)
      enddo
      enddo
cbs   now transform by rootovlp*dummy*rootovlp
      do jrun=1,nprims
      do irun=1,nprims
              TKIN(irun,jrun)=0d0
              dummy2(irun,jrun)=0d0
      enddo
      enddo
      do irun=1,nprims
      do jrun=1,nprims
      do krun=1,nprims
              dummy2(irun,jrun)=dummy2(irun,jrun)+
     *        dummy(irun,krun)*rootovlpinv(krun,jrun)
      enddo
      enddo
      enddo
      do irun=1,nprims
      do jrun=1,nprims
      do krun=1,nprims
              Tkin(irun,jrun)=Tkin(irun,jrun)+
     *        dummy2(krun,jrun)*rootovlpinv(irun,krun)
      enddo
      enddo
      enddo
      return
      end
