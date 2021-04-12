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
       subroutine addinta (wrk,wrksize,
     & syma,ammap)
c
c     this routine do for all a in syma
c     1- reconstruct #2 <_a,m,p,q> from TEMPDA2 file
c     2- prepair corresponding <_am p q> (like <amef>aaaa) to #3
c     and write it to opened INTA1-4
c     N.B.  this routine use followuing foreign routines:
c     wrtmap
c     wri
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer syma
       integer ammap(1:mbas,1:8,1:8)
c
c     help variables
c
       integer lenefaaaa,lenefbaab,lenefbbbb,lenefabab
       integer lenejaaaa,lenejbaab,lenejbaba,lenejbbbb,lenejabab,
     & lenejabba
       integer posst,rc,a
c
c*    mapd2 and mapi2 of #2 <_a,m|p,q> are prepaired
c
c*    make required mapd3 and mapi3 and write them to INTA1-4
c     define lengths of this mediates
c
c*1   to INTA1 <m,_a||ef>aaaa, <m,_a||ef>baab
       call ccsort_grc0(3,2,1,3,3,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenefaaaa)
       call dawrtmap (luna1,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,2,3,4,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenefbaab)
       call dawrtmap (luna1,mapd3,mapi3,rc)
c
c*2   to INTA2 <m,_a||ef>bbbb, <m,_a||ef>abab
       call ccsort_grc0(3,2,2,4,4,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenefbbbb)
       call dawrtmap (luna2,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,1,3,4,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenefabab)
       call dawrtmap (luna2,mapd3,mapi3,rc)
c
c*3   to INTA3 <m,_a||ej>aaaa, <m,_a||ej>baab, <m,_a||ej>baba
       call ccsort_grc0(3,0,1,3,1,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenejaaaa)
       call dawrtmap (luna3,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,2,3,2,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenejbaab)
       call dawrtmap (luna3,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,2,4,1,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenejbaba)
       call dawrtmap (luna3,mapd3,mapi3,rc)
c
c*4   to INTA4 <m,_a||ej>bbbb, <m,_a||ej>abba, <m,_a||ej>abab
       call ccsort_grc0(3,0,2,4,2,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenejbbbb)
       call dawrtmap (luna4,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,1,4,1,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenejabba)
       call dawrtmap (luna4,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,1,3,2,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenejabab)
       call dawrtmap (luna4,mapd3,mapi3,rc)
c
c
c*    cycle over a
c
       do 1000 a=1,nvb(syma)
c
c*    reconstruct #2 <_a,m,p,q> for given _a
       call mkampq (wrk,wrksize,
     & a,ammap)
c
c*    get contributions to INTA2 <m,_a||ef>bbbb, <m,_a||ef>abab
c     and wtite it there
c
       if (lenefbbbb.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,2,2,4,4,1,1)
       call dawri (luna2,lenefbbbb,wrk(mapd3(1,1)))
       end if
c
       if (lenefabab.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,1,3,4,1,0)
       call dawri (luna2,lenefabab,wrk(mapd3(1,1)))
       end if
c
c*    get contributions to INTA4 <m,_a||ej>bbbb, <m,_a||ej>abba, <m,_a||ej>abab
c     and wtite it there
c
       if (lenejbbbb.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,2,4,2,1,1)
       call dawri (luna4,lenejbbbb,wrk(mapd3(1,1)))
       end if
c
       if (lenejabba.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,1,4,1,0,1)
       call dawri (luna4,lenejabba,wrk(mapd3(1,1)))
       end if
c
       if (lenejabab.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,1,3,2,1,0)
       call dawri (luna4,lenejabab,wrk(mapd3(1,1)))
       end if
c
       if (a.gt.(nvb(syma)-nva(syma))) then
c     contributions to INTA1 and INTA3 only for a-alfa
c
c*    get contributions to INTA1 <m,_a||ef>aaaa, <m,_a||ef>baab if any
c     and wtite it there
c
       if (lenefaaaa.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,2,1,3,3,1,1)
       call dawri (luna1,lenefaaaa,wrk(mapd3(1,1)))
       end if
c
       if (lenefbaab.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,2,3,4,0,1)
       call dawri (luna1,lenefbaab,wrk(mapd3(1,1)))
       end if
c
c*    get contributions to INTA3 <m,_a||ej>aaaa, <m,_a||ej>baab, <m,_a||ej>baba
c     and wtite it there
c
       if (lenejaaaa.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,1,3,1,1,1)
       call dawri (luna3,lenejaaaa,wrk(mapd3(1,1)))
       end if
c
       if (lenejbaab.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,2,3,2,0,1)
       call dawri (luna3,lenejbaab,wrk(mapd3(1,1)))
       end if
c
       if (lenejbaba.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,2,4,1,1,0)
       call dawri (luna3,lenejbaba,wrk(mapd3(1,1)))
       end if
c
       end if
c
 1000   continue
c
       return
       end
