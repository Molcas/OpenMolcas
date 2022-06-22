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
       subroutine getw3 (wrk,wrksize,
     & lunw3xxxx,nxxxx)
c
c     This routine reconstruct W3(m,e,a,j)xxxx from lunw3xxxx file
c     to V1(m,e,a,j) and define corresponding mapdv1 and mapiv1
c     This routine also close lunw3xxxx file
c
c     lunw3xxxx - lun of opened w3xxxx file
c     nxxxx     - xxxx identifier
c     1 - aaaa
c     2 - bbbb
c     3 - aabb
c     4 - abba
c     5 - baab
c     6 - bbaa
c
c     the structure of lunw3xxxx is:
c
c     do syma=1,nsym
c     map of H _a(m,e,j)bbbb (W3(m,e,a,j))
c     skip cycle over a if length of all files is zero
c     do a=1,nvx(syma) [ x is a or b ]
c     if (h1length.gt.0) then
c     write H(m,e,j)
c     end if
c     end do
c     end do
c
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
c
       integer lunw3xxxx,nxxxx
c
c     help variables
c
       integer rc,syma,a,h1length,posst,aup,aalfayes,iiv1,v1length
c
c0    def aalfayes
c
       if ((nxxxx.eq.1).or.(nxxxx.eq.5).or.(nxxxx.eq.6)) then
       aalfayes=1
       else
       aalfayes=0
       end if
c
c1.1  def maps of V1(m,e,a,j)
c
       if (nxxxx.eq.1) then
       call grc0 (4,0,1,3,3,1,1,possv10,posst,mapdv1,mapiv1)
       else if (nxxxx.eq.2) then
       call grc0 (4,0,2,4,4,2,1,possv10,posst,mapdv1,mapiv1)
       else if (nxxxx.eq.3) then
       call grc0 (4,0,1,3,4,2,1,possv10,posst,mapdv1,mapiv1)
       else if (nxxxx.eq.4) then
       call grc0 (4,0,1,4,4,1,1,possv10,posst,mapdv1,mapiv1)
       else if (nxxxx.eq.5) then
       call grc0 (4,0,2,3,3,2,1,possv10,posst,mapdv1,mapiv1)
       else if (nxxxx.eq.6) then
       call grc0 (4,0,2,4,3,1,1,possv10,posst,mapdv1,mapiv1)
       end if
c
c1.2  vanish V1
       iiv1=mapdv1(0,5)
       v1length=mapdv1(iiv1,1)+mapdv1(iiv1,2)-possv10
       call mv0zero (v1length,v1length,wrk(possv10))
c
c2    rewind tape lunw3xxxx
       call filemanager (2,lunw3xxxx,rc)
c
c3    loop over symA
c
       do 3000 syma=1,nsym
c
c3.1  get map of H _a(m,e,j) to mapd,i H1
       call getmap (lunw3xxxx,possh10,h1length,mapdh1,mapih1,rc)
c
c3.2  skip cycle over a if length of H1 is 0
       if (h1length.eq.0) goto 3000
c
c3.3  loop over all a in this symmetry
c
       if (aalfayes.eq.1) then
       aup=nva(syma)
       else
       aup=nvb(syma)
       end if
c
       do 2500 a=1,aup
c
       if (h1length.gt.0) then
c
c3.3.1read H1 if any
       call rea (lunw3xxxx,h1length,wrk(possh10))
c
c3.3.2insert H1 into V1 for given a and syma
       call add (wrk,wrksize,
     & 3,4,1,3,a,0,syma,1,1.0d0,mapdh1,syma,mapdv1,mapiv1,1,
     &           rc)
c
       end if
c
 2500   continue
c
 3000   continue
c
c4    close lunw3xxxx file
       call filemanager (3,lunw3xxxx,rc)
c
       return
       end
