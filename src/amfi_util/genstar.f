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
      subroutine genstar(Lhigh)
      implicit real*8 (a-h,o-z)
cbs   purpose: generate start adresses of contraction coeffs on
cbs   contrarray for the different L-Blocks
#include "para.fh"
#include "param.fh"
      istart=1
      do L=0,Lhigh
      inc=nprimit(L)*ncontrac(L)
      iaddori(L)=istart
      istart=istart+inc
      iaddtyp1(L)=istart
      istart=istart+inc
      iaddtyp2(L)=istart
      istart=istart+inc
      iaddtyp3(L)=istart
      istart=istart+inc
      iaddtyp4(L)=istart
      istart=istart+inc
      enddo
      return
      end
