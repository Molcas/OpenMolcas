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
       subroutine contf12 (wrk,wrksize)
c
c     this routine do
c     FI2   f1(a,e) <- -0.5 sum(m) [t1o(a,m) . fok(e,m)]
c
c     N.B. use and destroy: M1,M2
c
        use Para_Info, only: MyRank
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
c
c     help variables
c
       integer rc,posst,ssc
c
c1    f1(a,e)aa <- sum(m-a) [T1o(a,m)aa . fok(e,m)aa]
c
cpar
      if (myRank.eq.idbaab) then
c
c1.1  map M1(m,e) <- fok(e,m)aa
       call map (wrk,wrksize,
     & 2,2,1,0,0,mapdfk3,mapifk3,1,mapdm1,mapim1,possm10,
     &           posst,rc)
c1.2  mult M2(a,e) = t1o(a,m)aa . M1(m,e)
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdt11,mapit11,1,mapdm1,mapim1,1,mapdm2,
     &            mapim2,ssc,possm20,rc)
c1.3  add f1(a,e)aa <- -0.5 M2(a,e)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-0.5d0,mapdm2,1,mapdf11,mapif11,1,rc)
c
       end if
c
c
c
c2    f1(a,e)bb <- sum(m-b) [T1o(a,m)bb . fok(e,m)bb]
c
cpar
      if (myRank.eq.idaabb) then
c
c2.1  map M1(m,e) <- fok(e,m)bb
       call map (wrk,wrksize,
     & 2,2,1,0,0,mapdfk4,mapifk4,1,mapdm1,mapim1,possm10,
     &           posst,rc)
c2.2  mult M2(a,e) = t1o(a,m)bb . M1(m,e)
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdt12,mapit12,1,mapdm1,mapim1,1,mapdm2,
     &            mapim2,ssc,possm20,rc)
c2.3  add f1(a,e)bb <- -0.5 M2(a,e)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-0.5d0,mapdm2,1,mapdf12,mapif12,1,rc)
c
        end if
c
       return
       end
