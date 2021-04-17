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
       subroutine percentzero (wrk,wrksize,
     & mapd,pz)
c
c     this routine test % of small elements in meditate, decribed by mpd
c
c     mapd - direct map of required mediate (I)
c
#include "wrk.fh"
       integer mapd(0:512,1:6)
       real*8  pz
c
c     help variables
c
       integer poss,length
       integer nhelp,nzero
       real*8 zerolim
c
c     def length, poss, zerolim
c
       poss=mapd(1,1)
       nhelp=mapd(0,5)
       length=mapd(nhelp,1)+mapd(nhelp,2)-mapd(1,1)
       zerolim=1.0d-6
c
       if (length.gt.0) then
       nzero=0
       do 100 nhelp=poss,poss+length-1
       if (abs(wrk(nhelp)).lt.zerolim) then
       nzero=nzero+1
       end if
 100    continue
       pz = dble(100*nzero)/dble(length)
       else
       pz=1.0d0
       end if
c
       return
       end
