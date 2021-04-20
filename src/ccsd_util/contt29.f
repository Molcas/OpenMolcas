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
       subroutine contt29 (wrk,wrksize)
c
c     this routine do:
c
c     T29
c1    Q(c,d,ij)aaaa   <= - sum(m-a)  [ T1o(c,m)aa . <ij||md>aaaa ]
c     T2n(ab,ij)aaaa   <- Q(a,b,ij)aaaa - Q(b,a,ij)aaaa
c2    Q(c,d,ij)bbbb   <= - sum(m-b)  [ T1o(c,m)bb . <ij||md>bbbb ]
c     T2n(ab,ij)bbbb   <- Q(a,b,ij)bbbb - Q(b,a,ij)bbbb
c3    T2n(a,b,i,j)abab <- - sum(m-a)  [ T1o(a,m)aa . <ij||mb)abab ]
c4    T2n(a,b,i,j)abab <- + sum(m-b)  [ T1o(b,m)bb . <ij||ma>abba ]
c
c     N.B. use and destroy : V1,V2
c
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
c
c     help variables
c
       integer posst,rc,ssc
c
cpar
        if (myRank.eq.idfin) then
c
c1    T2n(ab,ij)aaaa <- -P(a,b) sum(m-a) [T1o(a,m)aa . <md||ij>aaaa]
c
c1.1  mult V1(c,d,ij) <= T1o(c,m)aa . <md||ij>aaaa
       call mult (wrk,wrksize,
     & 2,4,4,1,mapdt11,mapit11,1,mapdw11,mapiw11,1,mapdv1,
     &            mapiv1,ssc,possv10,rc)
c
c1.2  pack V2(cd,ij) <= V1(c,d,ij) - V1(d,c,ij)
       call fack (wrk,wrksize,
     & 4,4,mapdv1,1,mapiv1,mapdv2,mapiv2,possv20,rc)
c
c1.3  add t2n(ab,ij)aaaa <- - V2(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv2,1,mapdt21,mapit21,1,rc)
c
c
c
c2    T2n(ab,ij)bbbb <- -P(a,b) sum(m-b) [T1o(a,m)bb . <md||ij>bbbb]
c
c2.1  mult V1(c,d,ij) <= T1o(c,m)bb . <md||ij>bbbb
       call mult (wrk,wrksize,
     & 2,4,4,1,mapdt12,mapit12,1,mapdw12,mapiw12,1,mapdv1,
     &            mapiv1,ssc,possv10,rc)
c
c2.2  pack V2(cd,ij) <= V1(c,d,ij) - V1(d,c,ij)
       call fack (wrk,wrksize,
     & 4,4,mapdv1,1,mapiv1,mapdv2,mapiv2,possv20,rc)
c
c2.3  add t2n(ab,ij)bbbb <- - V2(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv2,1,mapdt22,mapit22,1,rc)
c
c
c
c3    T2n(a,b,i,j)abab <- - sum(m-a)  [ T1o(a,m)aa . <ij||mb)abab ]
c
c3.1  mult V1(a,b,i,j) <= T1o(a,m)aa . <mb||ij>abab
       call mult (wrk,wrksize,
     & 2,4,4,1,mapdt11,mapit11,1,mapdw13,mapiw13,1,mapdv1,
     &            mapiv1,ssc,possv10,rc)
c
c3.2  add t2n(a,b,i,j)abab <- - V1(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv1,1,mapdt23,mapit23,1,rc)
c
c
c
c4    T2n(a,b,i,j)abab <- + sum(m-b)  [ T1o(b,m)bb . <ij||ma>abba ]
c
c4.1  mult V1(b,a,i,j) <= T1o(b,m)bb . <ma||ij>baab
       call mult (wrk,wrksize,
     & 2,4,4,1,mapdt12,mapit12,1,mapdw14,mapiw14,1,mapdv1,
     &            mapiv1,ssc,possv10,rc)
c
c4.2  map V2(a,b,i,j) <= V1(b,a,i,j)
       call map (wrk,wrksize,
     & 4,2,1,3,4,mapdv1,mapiv1,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
c
c4.3  add t2n(a,b,i,j)abab <-  V2(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv2,1,mapdt23,mapit23,1,rc)
c
cpar
        end if
c
       return
       end
