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
      subroutine initfrac(nprimit1,nprimit2,
     *nprimit3,nprimit4,
     *quot1,quot2,expo1,expo2,
     *expo3,expo4)
cbs   initialize some arrays with factors  needed for cfunct(x)
      implicit real*8(a-h,o-z)
      dimension expo1(*),expo2(*),expo3(*),expo4(*),
     *quot1(nprimit1,nprimit2,nprimit3,nprimit4),
     *quot2(nprimit1,nprimit2,nprimit3,nprimit4)
      do irun4=1,nprimit4
      do irun3=1,nprimit3
      do irun2=1,nprimit2
              sum24=expo2(irun2)+expo4(irun4)
                      do irun1=1,nprimit1
                      quot1(irun1,irun2,irun3,irun4)=
     *                1d0/(1d0+(expo1(irun1)+expo3(irun3))/
     *                sum24)
                      enddo
      enddo
      enddo
      enddo
      do irun4=1,nprimit4
      do irun3=1,nprimit3
      do irun2=1,nprimit2
              sum24=expo2(irun2)+expo4(irun4)
                      do irun1=1,nprimit1
                      quot2(irun1,irun2,irun3,irun4)=
     *                1d0/(1d0+sum24/
     *                (expo1(irun1)+expo3(irun3)))
                      enddo
      enddo
      enddo
      enddo
      return
      end
