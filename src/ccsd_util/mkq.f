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
       subroutine mkq (wrk,wrksize,
     & mapdt2,mapit2,mapdt11,mapit11,mapdt12,mapit12,
     &                 fact,rc)
c
c     this routine do:
c     t2(a,b,i,j) = fact . t2(a,b,i,j) +t11(ai).t12(bj)
c     for T2aaaa, T2bbbb and T2abab but they must be expanded (typ=0)
c
c     mapdt2  - direct map of T2 (I)
c     mapit2        - inverse map of T2 (I)
c     mapdt11        - direct map of T1 (I)
c     mapit11        - inverse map of T1 (I)
c     mapdt12        - direct map of T11 (I)
c     mapit12        - inverse map of T12 (I)
c     fact    - numerical factor (I)
c     rc        - return (error) code
c
c
#include "ccsd1.fh"
#include "wrk.fh"
       integer rc
       real*8 fact
c
       integer mapdt2(0:512,1:6)
       integer mapit2(1:8,1:8,1:8)
c
       integer mapdt11(0:512,1:6)
       integer mapit11(1:8,1:8,1:8)
c
       integer mapdt12(0:512,1:6)
       integer mapit12(1:8,1:8,1:8)
c
c     help variables
c
       integer posst2,posst11,posst12
       integer dimi,dimj,dima,dimb,syma,symb,symi,symj
       integer iit2,iit11,iit12
c
       rc=0
c
       if (mapdt2(0,6).eq.0) then
cI.1  typ of t2 is 0 (T2 is expanded)

       do 100 iit2=1,mapdt2(0,5)
c
       posst2=mapdt2(iit2,1)
       syma=mapdt2(iit2,3)
       symb=mapdt2(iit2,4)
       symi=mapdt2(iit2,5)
       symj=mapdt2(iit2,6)
       dima=dimm(mapdt2(0,1),syma)
       dimb=dimm(mapdt2(0,2),symb)
       dimi=dimm(mapdt2(0,3),symi)
       dimj=dimm(mapdt2(0,4),symj)
       iit11=mapit11(syma,1,1)
       iit12=mapit12(symb,1,1)
       posst11=mapdt11(iit11,1)
       posst12=mapdt12(iit12,1)
c
       if ((syma.eq.symi).and.(symb.eq.symj).and.(mapdt2(iit2,2).gt.0))
     & then
       call mkqhelp1 (wrk(posst2),wrk(posst11),wrk(posst12),
     & dima,dimb,dimi,dimj,noa(symi),nob(symj),fact)
       else if (mapdt2(iit2,2).gt.0) then
       call mkqhelp2 (wrk(posst2),mapdt2(iit2,2),mapdt2(iit2,2),fact)
       end if
c
 100    continue
c
c
       else
cI.4  RC=1 : typ of T2 is not 0
       rc=1
       return
       end if
c
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(mapit2)
       end
