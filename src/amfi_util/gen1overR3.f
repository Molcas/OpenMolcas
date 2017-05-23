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
      subroutine gen1overR3(Lhigh,oneoverR3)
      implicit real*8 (a-h,o-z)
cbs   generates the radial integrals  for the one electron spin orbit integrals
cbs   taken the 1/r**3 formula from the documentation and included additional
cbs   factors for normalization
#include "para.fh"
#include "param.fh"
#include "dofuc.fh"
#include "real.fh"
      dimension oneoverR3((MxprimL*MxprimL+MxprimL)/2,Lmax)
      do L=1,Lhigh
      icount=0
      do iprim2=1,nprimit(L)
      alpha2=exponents(iprim2,L)
      do iprim1=1,iprim2
      alpha1=exponents(iprim1,L)
      icount=icount+1
      oneoverR3(icount,L)=sqrt(2d0/pi)*
     *(df(L+L-2)*DBLE(2**(L+3))*
     *(alpha1*alpha2)**(0.25d0*
     *DBLE(L+L+3)))/((alpha1+alpha2)**L*df(L+L+1))
      enddo
      enddo
      enddo
      return
      end
