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
       subroutine max5h2 (wrk,wrksize,
     & nind,mapd,mapi,rmax,imax,text)
c
c     this routine write:
c     a) note
c     b) 5 maximal elements with their indexes in given vector V
c     c) euclidian norm
c
c     nind  - number of indexes in V (I)
c     mapd  - direct map of V (I)
c     mapi  - inverese map of V (I)
c     rmax  - store of maximal values (I)
c     imax  - store of corr. indexes (I)
c     text  - notice (I)
c
#include "ccsd1.fh"
#include "wrk.fh"
c
       integer nind
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
       integer imax(1:8,1:5)
       real*8 rmax (1:5)
       character*8 text
c
c     help variables
c
       integer nhelp1,nhelp2,rc
       real*8 scalar
c
c1    write notice
c
       write(6,101) text
 101    format (' Five largest amplitudes of :',a8)
c
c2    write 5 maximal amplitudes
c
       write(6,102)
 102    format ('  SYMA   SYMB   SYMI   SYMJ     A      B',
     &          '      I      J     VALUE')
       do 10 nhelp1=1,5
       write(6,103) (imax(nhelp2,nhelp1),nhelp2=1,8),rmax(nhelp1)
 103    format (8(2x,i3,2x),f15.10)
 10     continue
c
c3    write euclidian norm
c
c3.1  calc euclidian norm
       call multdot (wrk,wrksize,
     & nind,mapd,mapi,1,mapd,mapi,1,scalar,rc)
       scalar=sqrt(scalar)
c
       write(6,104) scalar
 104    format (' Euclidian norm is :',f17.10)
c
       write(6,*)
c
       return
       end
