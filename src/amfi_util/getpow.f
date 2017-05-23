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
      subroutine getpow(max,quot,quotpow,
     *nprim1,nprim2,nprim3,nprim4)
cbs   generates some powers of for the prefactors of cfunct(X)
cbs   look out for details there and in initfrac
      implicit real*8 (a-h,o-z)
#include "para.fh"
      dimension quotpow(nprim1,nprim2,
     *nprim3,nprim4),
     *quot(nprim1,nprim2,nprim3,nprim4)
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      quotpow(irun1,irun2,irun3,irun4)=
     *sqrt(quot(irun1,irun2,irun3,irun4))
      enddo
      enddo
      enddo
      enddo
      if (max.eq.1) return
cbs
      do irun=2,max
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      quotpow(irun1,irun2,irun3,irun4)=
     *quotpow(irun1,irun2,irun3,irun4)*
     *quot(irun1,irun2,irun3,irun4)
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
