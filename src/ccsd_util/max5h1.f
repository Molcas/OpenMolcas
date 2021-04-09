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
       subroutine max5h1 (imax,rmax,symp,symq,symr,syms,p,q,r,s,value)
c
c     this routine add new max amplitude and skip smallest
c
c     imax - store of indexes of 5 max (I/O)
c     rmax - store of values of 5 max (I/O)
c     symp - symmetry of p index (I)
c     symq - symmetry of q index (I)
c     symr - symmetry of r index (I)
c     syms - symmetry of s index (I)
c     p    - value of p index (I)
c     q    - value of q index (I)
c     r    - value of r index (I)
c     s    - value of s index (I)
c     value- value of amplitude (I)
c
       integer imax(1:8,1:5)
       real*8 rmax (1:5)
       integer symp,symq,symr,syms,p,q,r,s
       real*8 value
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp3
c
c1    find possition of this value
c
       do 10 nhelp1=1,5
       if (abs(value).ge.abs(rmax(nhelp1))) then
       goto 20
       end if
 10     continue
c
c2    1-(nhelp1-1) stay untought
c
c3    push other records if necc.
c
 20     if (nhelp1.lt.5) then
       do 30 nhelp2=4,nhelp1,-1
       rmax(nhelp2+1)=rmax(nhelp2)
       do 31 nhelp3=1,8
       imax(nhelp3,nhelp2+1)=imax(nhelp3,nhelp2)
 31     continue
 30     continue
       end if
c
c4    add new one
c
       rmax(nhelp1)=value
       imax(1,nhelp1)=symp
       imax(2,nhelp1)=symq
       imax(3,nhelp1)=symr
       imax(4,nhelp1)=syms
       imax(5,nhelp1)=p
       imax(6,nhelp1)=q
       imax(7,nhelp1)=r
       imax(8,nhelp1)=s
c
       return
       end
