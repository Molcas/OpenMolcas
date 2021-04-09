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
       subroutine stz (wrk,wrksize,
     & mapda)
c
c     this routine vanish A
c     A = 0
c
c     mapda  - direct map of A m(I/O)
c
c     N.B. this routine should be done using matrix operations

#include "wrk.fh"
       integer mapda(0:512,1:6)
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp3
c
c
c1    def the length of the mediate
       nhelp1=mapda(0,5)
       nhelp3=mapda(nhelp1,1)+mapda(nhelp1,2)-mapda(1,1)
c
c2    def initial possition
       nhelp2=mapda(1,1)
c
c3    refactoring
       do 100 nhelp1=nhelp2,nhelp2+nhelp3-1
       wrk(nhelp1)=0.0d0
 100    continue
c
       return
       end
