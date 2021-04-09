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
       subroutine multy0 (wrk,wrksize,
     & mvec,ix,mapdy,key)
c
c     This routine realize multiplying according mvec
c     for Y=A*B
c     N.B. if key=0, Y file is not vanished (ie can be used for
c     adding to some existing file)
c
#include "ccsd1.fh"
#include "wrk.fh"
       integer mvec(1:4096,1:7)
       integer ix,key
       integer mapdy(0:512,1:6)

c     help variables
c
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5
       integer iix,iy
c
c1    set C=0
c
       if (key.eq.1) then
c
c     Y vector must be vanished
c
       do 100 iy=1,mapdy(0,5)
       nhelp1=mapdy(iy,1)
       nhelp2=mapdy(iy,2)
       call mv0zero (nhelp2,nhelp2,wrk(nhelp1))
 100    continue
c
       end if
c
c2    Y=Y+A*B
c
       if (ix.eq.0) then
       return
       end if
c
       do 200 iix=1,ix
c
c     skip this sumation if yes/no=0
       if (mvec(iix,1).eq.0) goto 200
c
c     realize individial sumation
c
c     def possitions of A,B,Y
       nhelp1=mvec(iix,2)
       nhelp2=mvec(iix,3)
       nhelp3=mvec(iix,4)
c
c     def rowA(rowY), colA(sum)
       nhelp4=mvec(iix,5)
       nhelp5=mvec(iix,6)
c
       call mv0v1a3u (nhelp4,nhelp5,nhelp5,nhelp4,
     & nhelp4,nhelp5,1,1,
     & wrk(nhelp1),wrk(nhelp2),wrk(nhelp3))
c
 200    continue
c
       return
       end
