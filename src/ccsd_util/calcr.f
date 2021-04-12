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
       subroutine calcr (wrk,wrksize,
     & lune)
c
c     this routine calc difference vector Tn = Tn-E
c     Tn=(T12,T22,T23,T13,T14)
c
c     lune - lun of file, where E is stored (I)
c
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
       integer lune
c
c     help variables
c
       integer rc
c
c1    rewind lune
       call filemanager (2,lune,rc)
c
c     T2aaaa
       call getmediate (wrk,wrksize,
     & lune,possv10,mapdv1,mapiv1,rc)
       call calcrh1 (wrk,wrksize,
     & mapdt21,mapdv1)
c     T2bbbb
       call getmediate (wrk,wrksize,
     & lune,possv10,mapdv1,mapiv1,rc)
       call calcrh1 (wrk,wrksize,
     & mapdt22,mapdv1)
c     T2abab
       call getmediate (wrk,wrksize,
     & lune,possv10,mapdv1,mapiv1,rc)
       call calcrh1 (wrk,wrksize,
     & mapdt23,mapdv1)
c     T1aa
       call getmediate (wrk,wrksize,
     & lune,possv10,mapdv1,mapiv1,rc)
       call calcrh1 (wrk,wrksize,
     & mapdt13,mapdv1)
c     T1bb
       call getmediate (wrk,wrksize,
     & lune,possv10,mapdv1,mapiv1,rc)
       call calcrh1 (wrk,wrksize,
     & mapdt14,mapdv1)
c
       return
       end
