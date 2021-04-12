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
       subroutine contt12 (wrk,wrksize)
c
c     this routine do T12 contribution:
c     t1n(a,i) <- sum(e) [FI(a,e) . T1o(e,i)
c
c     N.B. use and destroy : M1
c     N.B. Parallel : in the case where idaaaa.ne.idbaab and
c          idaaaa.ne.idbaab this routine runs contributions to
c          T1n also on idaaaa and idbbbb nodes, since on these
c          nodes there is a specific
c          part of contributions F13 (see notes in sumoverb routine)
c           which are not presented on pilot nodes
c
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
c
c     help variables
c
       integer rc,ssc
c
c1    t1n(a,i)aa <- sum(e-a) [F1(a,e)aa . T1o(e,i)aa ]
c
cpar
      if ((myRank.eq.idbaab).or.(myRank.eq.idaaaa)) then
c1.1  mult M1(a,i) <= F1(a,e)aa . T1o(e,i)aa
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdf11,mapif11,1,mapdt11,mapit11,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c1.2  add t1n(a,i)aa <- M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdt13,mapit13,1,rc)
      end if
c
c
c2    t1n(a,i)bb <- sum(e-a) [F1(a,e)bb . T1o(e,i)bb ]
c
cpar
      if ((myRank.eq.idaabb).or.(myRank.eq.idbbbb)) then
c2.1  mult M1(a,i) <= F1(a,e)bb . T1o(e,i)bb
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdf12,mapif12,1,mapdt12,mapit12,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c2.2  add t1n(a,i)bb <- M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdt14,mapit14,1,rc)
      end if
c
       return
       end
