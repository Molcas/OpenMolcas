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
       subroutine fokunpck2 (fok,faa,dimfok,dimfa,shift)
c
c     this routine distribute (Fok - dp) to Faa
c     fok    - Fok matrix (I)
c     faa    - Faa matrix (O)
c     dimfok - dimension for Fok matrix - norb (I)
c     dimfa  - dimension of virtuals - nv (I)
c
       integer dimfok,dimfa,shift
c
       real*8 fok(1:dimfok,1:dimfok)
       real*8 faa(1:dimfa,1:dimfa)
c
c     help variables
c
       integer a,b
c
c1    distribute Fok to Faa
       do 200 b=1,dimfa
       do 201 a=1,dimfa
       faa(a,b)=fok(shift+a,shift+b)
 201    continue
 200    continue
c
       return
       end
