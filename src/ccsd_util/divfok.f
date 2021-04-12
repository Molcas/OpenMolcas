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
       subroutine divfok (wrk,wrksize,
     & mapdfa,mapifa,possfa0,mapdfb,mapifb,possfb0,
     & mapdfk1,mapifk1,possfk10,mapdfk2,mapifk2,possfk20,
     & mapdfk3,mapifk3,possfk30,mapdfk4,mapifk4,possfk40,
     & mapdfk5,mapifk5,possfk50,mapdfk6,mapifk6,possfk60,
     & mapddp1,mapidp1,possdp10,mapddp2,mapidp2,possdp20,rc)
c
c     this routine divide fok(p,q) -> fk(a,b) + fk(a,i) + f(i,j) + dp(p)
c     to diagonal part and rest
c
c     mapd and mapi for:
c     fa,fb - fok(p,q)aa,bb
c     fk1-6 - f(ab)aa,f(ab)bb,f(ai)aa,f(ai)bb,f(ij)aa,f(ij)bb
c     dp1,2 - diagonal part dp(p)a,b
c     rc    - return (error) code
c
c
#include "ccsd1.fh"
#include "wrk.fh"
       integer rc
c
c1    maps for FOKA,FOKB
c
       integer mapdfa(0:512,1:6)
       integer mapifa(1:8,1:8,1:8)
       integer possfa0
c
       integer mapdfb(0:512,1:6)
       integer mapifb(1:8,1:8,1:8)
       integer possfb0
c
c2    maps for FK
c     FK1 - f(a,b)aa
c     FK2 - f(a,b)bb
c     FK3 - f(a,i)aa
c     FK4 - f(a,i)bb
c     FK5 - f(i,j)aa
c     FK6 - f(i,j)bb
c
       integer mapdfk1(0:512,1:6)
       integer mapifk1(1:8,1:8,1:8)
       integer possfk10
c
       integer mapdfk2(0:512,1:6)
       integer mapifk2(1:8,1:8,1:8)
       integer possfk20
c
       integer mapdfk3(0:512,1:6)
       integer mapifk3(1:8,1:8,1:8)
       integer possfk30
c
       integer mapdfk4(0:512,1:6)
       integer mapifk4(1:8,1:8,1:8)
       integer possfk40
c
       integer mapdfk5(0:512,1:6)
       integer mapifk5(1:8,1:8,1:8)
       integer possfk50
c
       integer mapdfk6(0:512,1:6)
       integer mapifk6(1:8,1:8,1:8)
       integer possfk60
c
c
c3    maps for DP - diagonal part
c     DP1 - dp(p)a
c     DP2 - dp(p)b
c
       integer mapddp1(0:512,1:6)
       integer mapidp1(1:8,1:8,1:8)
       integer possdp10
c
       integer mapddp2(0:512,1:6)
       integer mapidp2(1:8,1:8,1:8)
       integer possdp20
c
c     help variables
c
       integer symp,rc1
       integer iifoka,iifokb,iifok,iifaa,iifai,iifii,iidpa,iidpb,iidp
       integer possfoka,possfokb,possfok,possfaa,possfai,possfii
       integer possdpa,possdpb,possdp
c
       rc=0
c
c1    define dp
c
       do 500 symp=1,nsym
c
       iidpa=mapidp1(symp,1,1)
       possdpa=mapddp1(iidpa,1)
       iidpb=mapidp2(symp,1,1)
       possdpb=mapddp2(iidpb,1)
       iifoka=mapifa(symp,1,1)
       possfoka=mapdfa(iifoka,1)
       iifokb=mapifb(symp,1,1)
       possfokb=mapdfb(iifokb,1)
c
       if (norb(symp).gt.0) then
       call fokunpck5 (symp,wrk(possfoka),wrk(possfokb),
     & wrk(possdpa),wrk(possdpb),norb(symp),rc1)
       end if
c
 500    continue
c
c2    define faa,fai,fii
c
       do 1000 symp=1,nsym
       if (norb(symp).eq.0) goto 1000
c
c2.1  alpha case
c
       iifok=mapifa(symp,1,1)
       iifaa=mapifk1(symp,1,1)
       iifai=mapifk3(symp,1,1)
       iifii=mapifk5(symp,1,1)
       iidp=mapidp1(symp,1,1)
c
       possfok=mapdfa(iifok,1)
       possfaa=mapdfk1(iifaa,1)
       possfai=mapdfk3(iifai,1)
       possfii=mapdfk5(iifii,1)
       possdp=mapddp1(iidp,1)
c
       call fokunpck1 (wrk(possfok),wrk(possdp),norb(symp))
       if (nva(symp).gt.0) then
       call fokunpck2 (wrk(possfok),wrk(possfaa),norb(symp),nva(symp),
     &                 noa(symp))
       end if
       if ((noa(symp)*nva(symp)).gt.0) then
       call fokunpck3 (wrk(possfok),wrk(possfai),norb(symp),nva(symp),
     &                 noa(symp))
       end if
       if (noa(symp).gt.0) then
       call fokunpck4 (wrk(possfok),wrk(possfii),norb(symp),noa(symp))
       end if
c
c2.2  alpha case
c
       iifok=mapifb(symp,1,1)
       iifaa=mapifk2(symp,1,1)
       iifai=mapifk4(symp,1,1)
       iifii=mapifk6(symp,1,1)
       iidp=mapidp2(symp,1,1)
c
       possfok=mapdfb(iifok,1)
       possfaa=mapdfk2(iifaa,1)
       possfai=mapdfk4(iifai,1)
       possfii=mapdfk6(iifii,1)
       possdp=mapddp2(iidp,1)
c
       call fokunpck1 (wrk(possfok),wrk(possdp),norb(symp))
       if (nvb(symp).gt.0) then
       call fokunpck2 (wrk(possfok),wrk(possfaa),norb(symp),nvb(symp),
     &                 nob(symp))
       end if
       if ((nob(symp)*nvb(symp)).gt.0) then
       call fokunpck3 (wrk(possfok),wrk(possfai),norb(symp),nvb(symp),
     &                 nob(symp))
       end if
       if (nob(symp).gt.0) then
       call fokunpck4 (wrk(possfok),wrk(possfii),norb(symp),nob(symp))
       end if
c
 1000   continue
c
       return
c Avoid unused argument warnings
       if (.false.) then
         call Unused_integer(possfa0)
         call Unused_integer(possfb0)
         call Unused_integer(possfk10)
         call Unused_integer(possfk20)
         call Unused_integer(possfk30)
         call Unused_integer(possfk40)
         call Unused_integer(possfk50)
         call Unused_integer(possfk60)
         call Unused_integer(possdp10)
         call Unused_integer(possdp20)
       end if
       end
