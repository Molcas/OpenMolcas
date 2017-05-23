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
      subroutine contcasb2OO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4OO)
cbs   contraction for powers (0)  with alpha3
cbs   this is one of the cases b in the documentation
cbs   use averaged integrals by interchanging kinematic factors
      implicit real*8 (a-h,o-z)
#include "para.fh"
#include "param.fh"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*)
     *,cont4OO(*)
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo
      call contract(
     *contrarray(iaddtyp1(l1)),
     *contrarray(iaddtyp3(l2)),
     *contrarray(iaddtyp4(l3)),
     *contrarray(iaddtyp1(l4)),
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo
      call contract(
     *contrarray(iaddtyp3(l1)),
     *contrarray(iaddtyp3(l2)),
     *contrarray(iaddtyp2(l3)),
     *contrarray(iaddtyp1(l4)),
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo
      call contract(
     *contrarray(iaddtyp1(l1)),
     *contrarray(iaddtyp1(l2)),
     *contrarray(iaddtyp4(l3)),
     *contrarray(iaddtyp3(l4)),
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo
      call contract(
     *contrarray(iaddtyp3(l1)),
     *contrarray(iaddtyp1(l2)),
     *contrarray(iaddtyp2(l3)),
     *contrarray(iaddtyp3(l4)),
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
      return
      end
