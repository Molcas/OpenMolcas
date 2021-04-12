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
       subroutine unpackk (i,vint,ndimv1,ndimv2,ndimv3,key)
c
c      unpackk process control routine
c
c     i      - value of pivot index (I)
c     vint   - array of integrals (O)
c     ndimv1 - first dimension of vint (norb(symj)) (I)
c     ndimv2 - second dimension of vint (norb(symk)) (I)
c     ndimv3 - third dimension of vint (norb(syml)) (I)
c     key    - reduced storing key (I)
c     = 0 if symj is not syml
c     = 1 if symj = syml
c
#include "reorg.fh"
       integer i,ndimv1,ndimv2,ndimv3,key
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
c
       if (zrkey.eq.1) then
         call unpackk_zr (i,vint,ndimv1,ndimv2,ndimv3,key)
       else
         call unpackk_pck (i,vint,ndimv1,ndimv2,ndimv3,key)
       end if
c
       return
       end
