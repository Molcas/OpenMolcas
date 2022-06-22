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
       subroutine contt13 (wrk,wrksize)
c
c     this routine do T13 contribution:
c     t1n(a,i) <- sum(m) [T1o(a,m) . FII(m,i)]
c
c     N.B. use and destroy : M1
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
c1    t1n(a,i)aa <- - sum(m-a) [T1o(a,m)aa  . FII(m,i)aa]
c
cpar
      if (myRank.eq.idbaab) then
c1.1  mult M1(a,i) <= T1o(a,m)aa  . FII(m,i)aa
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdt11,mapit11,1,mapdf21,mapif21,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c1.2  add t1n(a,i)aa <- M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdt13,mapit13,1,rc)
      end if
c
c
c
c2    t1n(a,i)bb <- - sum(m-b) [T1o(a,m)bb  . FII(m,i)bb]
c
cpar
      if (myRank.eq.idaabb) then
c2.1  mult M1(a,i) <= T1o(a,m)bb  . FII(m,i)bb
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdt12,mapit12,1,mapdf22,mapif22,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c2.2  add t1n(a,i)bb <- M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdt14,mapit14,1,rc)
      end if
c
       return
       end
