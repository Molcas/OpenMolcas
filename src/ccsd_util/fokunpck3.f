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
       subroutine fokunpck3 (fok,fai,dimfok,dimfa,dimfi)
c
c     this routine distribute (Fok - dp) to Fai
c     fok    - Fok matrix (I)
c     fai    - Fai matrix (O)
c     dimfok - dimension for Fok matrix - norb (I)
c     dimfa  - dimension of virtuals - nv (I)
c     dimfi  - dimension of occupied - no (I)
c
       integer dimfok,dimfa,dimfi
c
       real*8 fok(1:dimfok,1:dimfok)
       real*8 fai(1:dimfa,1:dimfi)
c
c     help variables
c
       integer a,i
c
c1    distribute Fok to Fai
       do 300 i=1,dimfi
       do 301 a=1,dimfa
       fai(a,i)=fok(dimfi+a,i)
 301    continue
 300    continue
c
       return
       end
