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
       subroutine max5 (wrk,wrksize,
     & nind,mapd,mapi,text)
c
c     this routine find and type:
c     a) note
c     b) 5 maximal elements with their indexes in given vector V
c     c) euclidian norm
c
c     nind  - number of indexes in V (I)
c     mapd  - direct map of V (I)
c     mapi  - inverese map of V (I)
c     text  - notice (I)
c
#include "ccsd1.fh"
#include "wrk.fh"
c
       integer nind
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
       character*8 text
c
c     help variables
c
       integer nhelp1,nhelp2,i,j,a,b,it
       real*8 value
       integer imax(1:8,1:5)
       real*8 rmax (1:5)
c
c0    set rmax,imax=0
c
       do 10 nhelp1=1,5
       rmax(nhelp1)=0.0d0
       do 11 nhelp2=1,8
       imax(nhelp2,nhelp1)=0
 11     continue
 10     continue
c
c
       if (nind.eq.2) then
c
c1    T1aa or T1bb amplitudes
c
c1.1  find 5 max
c
       nhelp1=mapd(1,1)
       do 100 it=1,mapd(0,5)
       do 50 i=1,dimm(mapd(0,2),mapd(it,4))
       do 51 a=1,dimm(mapd(0,1),mapd(it,3))
       value=wrk(nhelp1)
       if (abs(value).ge.abs(rmax(5))) then
c     write this amplitude
       call max5h1 (imax,rmax,mapd(it,3),0,mapd(it,4),0,a,0,i,0,
     & value)
       end if
       nhelp1=nhelp1+1
 51     continue
 50     continue
 100    continue
c
c1.2  write output
       If (fullprint.ge.0) call max5h2 (wrk,wrksize,
     & nind,mapd,mapi,rmax,imax,text)
c
       else if (mapd(0,6).eq.0) then
c
c2    T2abab amplitudes
c
c2.1  find 5 max
c
       nhelp1=mapd(1,1)
       do 200 it=1,mapd(0,5)
       do 150 j=1,dimm(mapd(0,4),mapd(it,6))
       do 151 i=1,dimm(mapd(0,3),mapd(it,5))
       do 152 b=1,dimm(mapd(0,2),mapd(it,4))
       do 153 a=1,dimm(mapd(0,1),mapd(it,3))
       value=wrk(nhelp1)
       if (abs(value).ge.abs(rmax(5))) then
c     write this amplitude
       call max5h1 (imax,rmax,mapd(it,3),mapd(it,4),
     & mapd(it,5),mapd(it,6),a,b,i,j,value)
       end if
       nhelp1=nhelp1+1
 153    continue
 152    continue
 151    continue
 150    continue
 200    continue
c
c2.2  write output
       If (fullprint.ge.0) call max5h2 (wrk,wrksize,
     & nind,mapd,mapi,rmax,imax,text)
c
       else
c
c3    T2aaaa or T2bbbb amplitudes
c
c3.1  find 5 max
c
       nhelp1=mapd(1,1)
       do 300 it=1,mapd(0,5)
c
       if (mapd(it,3).eq.mapd(it,4)) then
c     case syma=symb, symi=symj
c
       do 230 i=2,dimm(mapd(0,3),mapd(it,5))
       do 231 j=1,i-1
       do 232 a=2,dimm(mapd(0,1),mapd(it,3))
       do 233 b=1,a-1
       value=wrk(nhelp1)
       if (abs(value).ge.abs(rmax(5))) then
c     write this amplitude
       call max5h1 (imax,rmax,mapd(it,3),mapd(it,4),
     & mapd(it,5),mapd(it,6),a,b,i,j,value)
       end if
       nhelp1=nhelp1+1
 233    continue
 232    continue
 231    continue
 230    continue
c
       else
c     case syma>symb, symi> symj
       do 250 j=1,dimm(mapd(0,4),mapd(it,6))
       do 251 i=1,dimm(mapd(0,3),mapd(it,5))
       do 252 b=1,dimm(mapd(0,2),mapd(it,4))
       do 253 a=1,dimm(mapd(0,1),mapd(it,3))
       value=wrk(nhelp1)
       if (abs(value).ge.abs(rmax(5))) then
c     write this amplitude
       call max5h1 (imax,rmax,mapd(it,3),mapd(it,4),
     & mapd(it,5),mapd(it,6),a,b,i,j,value)
       end if
       nhelp1=nhelp1+1
 253    continue
 252    continue
 251    continue
 250    continue
       end if
c
 300    continue
c
c3.2  write output
       If (fullprint.ge.0) call max5h2 (wrk,wrksize,
     & nind,mapd,mapi,rmax,imax,text)
c
       end if
c
       return
       end
