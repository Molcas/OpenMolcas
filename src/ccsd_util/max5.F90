!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine max5 (wrk,wrksize,                                    &
     & nind,mapd,mapi,text)
!
!     this routine find and type:
!     a) note
!     b) 5 maximal elements with their indexes in given vector V
!     c) euclidian norm
!
!     nind  - number of indexes in V (I)
!     mapd  - direct map of V (I)
!     mapi  - inverese map of V (I)
!     text  - notice (I)
!
#include "ccsd1.fh"
#include "wrk.fh"
!
       integer nind
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
       character*8 text
!
!     help variables
!
       integer nhelp1,nhelp2,i,j,a,b,it
       real*8 value
       integer imax(1:8,1:5)
       real*8 rmax (1:5)
!
!0    set rmax,imax=0
!
       do 10 nhelp1=1,5
       rmax(nhelp1)=0.0d0
       do 11 nhelp2=1,8
       imax(nhelp2,nhelp1)=0
 11     continue
 10     continue
!
!
       if (nind.eq.2) then
!
!1    T1aa or T1bb amplitudes
!
!1.1  find 5 max
!
       nhelp1=mapd(1,1)
       do 100 it=1,mapd(0,5)
       do 50 i=1,dimm(mapd(0,2),mapd(it,4))
       do 51 a=1,dimm(mapd(0,1),mapd(it,3))
       value=wrk(nhelp1)
       if (abs(value).ge.abs(rmax(5))) then
!     write this amplitude
       call max5h1 (imax,rmax,mapd(it,3),0,mapd(it,4),0,a,0,i,0,        &
     & value)
       end if
       nhelp1=nhelp1+1
 51     continue
 50     continue
 100    continue
!
!1.2  write output
       If (fullprint.ge.0) call max5h2 (wrk,wrksize,                    &
     & nind,mapd,mapi,rmax,imax,text)
!
       else if (mapd(0,6).eq.0) then
!
!2    T2abab amplitudes
!
!2.1  find 5 max
!
       nhelp1=mapd(1,1)
       do 200 it=1,mapd(0,5)
       do 150 j=1,dimm(mapd(0,4),mapd(it,6))
       do 151 i=1,dimm(mapd(0,3),mapd(it,5))
       do 152 b=1,dimm(mapd(0,2),mapd(it,4))
       do 153 a=1,dimm(mapd(0,1),mapd(it,3))
       value=wrk(nhelp1)
       if (abs(value).ge.abs(rmax(5))) then
!     write this amplitude
       call max5h1 (imax,rmax,mapd(it,3),mapd(it,4),                    &
     & mapd(it,5),mapd(it,6),a,b,i,j,value)
       end if
       nhelp1=nhelp1+1
 153    continue
 152    continue
 151    continue
 150    continue
 200    continue
!
!2.2  write output
       If (fullprint.ge.0) call max5h2 (wrk,wrksize,                    &
     & nind,mapd,mapi,rmax,imax,text)
!
       else
!
!3    T2aaaa or T2bbbb amplitudes
!
!3.1  find 5 max
!
       nhelp1=mapd(1,1)
       do 300 it=1,mapd(0,5)
!
       if (mapd(it,3).eq.mapd(it,4)) then
!     case syma=symb, symi=symj
!
       do 230 i=2,dimm(mapd(0,3),mapd(it,5))
       do 231 j=1,i-1
       do 232 a=2,dimm(mapd(0,1),mapd(it,3))
       do 233 b=1,a-1
       value=wrk(nhelp1)
       if (abs(value).ge.abs(rmax(5))) then
!     write this amplitude
       call max5h1 (imax,rmax,mapd(it,3),mapd(it,4),                    &
     & mapd(it,5),mapd(it,6),a,b,i,j,value)
       end if
       nhelp1=nhelp1+1
 233    continue
 232    continue
 231    continue
 230    continue
!
       else
!     case syma>symb, symi> symj
       do 250 j=1,dimm(mapd(0,4),mapd(it,6))
       do 251 i=1,dimm(mapd(0,3),mapd(it,5))
       do 252 b=1,dimm(mapd(0,2),mapd(it,4))
       do 253 a=1,dimm(mapd(0,1),mapd(it,3))
       value=wrk(nhelp1)
       if (abs(value).ge.abs(rmax(5))) then
!     write this amplitude
       call max5h1 (imax,rmax,mapd(it,3),mapd(it,4),                    &
     & mapd(it,5),mapd(it,6),a,b,i,j,value)
       end if
       nhelp1=nhelp1+1
 253    continue
 252    continue
 251    continue
 250    continue
       end if
!
 300    continue
!
!3.2  write output
       If (fullprint.ge.0) call max5h2 (wrk,wrksize,                    &
     & nind,mapd,mapi,rmax,imax,text)
!
       end if
!
       return
       end
