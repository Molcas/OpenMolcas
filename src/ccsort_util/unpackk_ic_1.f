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
       subroutine unpackk_ic_1 (i,vint,ndimv1,ndimv2,ndimv3,
     c                       Vic,ndimvi)
c
c     this routine vint(j,k,l) = <i,j|k,l>
c     for given i from incore (nonreduced) expanded block Vic
c
c     i      - value of pivot index (I)
c     vint   - array of integrals (O)
c     ndimv1 - first dimension of vint (norb(symj)) (I)
c     ndimv2 - second dimension of vint (norb(symk)) (I)
c     ndimv3 - third dimension of vint (norb(syml)) (I)
c     Vic    - incore expanded block of integrals (I)
c     ndimvi - first dimension of Vic norb(symi) (I)
c
#include "reorg.fh"

#include "SysDef.fh"
       integer i,ndimv1,ndimv2,ndimv3,ndimvi
       real*8 vint(1:ndimv1*ndimv2*ndimv3)
       real*8 Vic(1:ndimvi,1:ndimv1*ndimv2*ndimv3)
c
c     help variables
c
      integer jkl
c
c
        do jkl=1,ndimv1*ndimv2*ndimv3
        vint(jkl)=Vic(i,jkl)
        end do
c
       return
       end
