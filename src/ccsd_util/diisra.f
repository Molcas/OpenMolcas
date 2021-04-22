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
       subroutine diisra (wrk,wrksize,
     & diispoint,num,
     & mapd1,mapi1,poss10,mapd2,mapi2,poss20,
     & mapd3,mapi3,poss30,mapd4,mapi4,poss40)
c
c     this routine read num vectors of amplitudes
c     form prepaired diispoint(1) - diispoint(num) files
c
c     diispoint - stack of lun numbers (I)
c     num       - number of vertors to be read (1-4) (I)
c     mapd1     - direct matrix of vector 1 (I)
c     mapi1     - inverse matrix of vector 1 (I)
c     poss10    - initial possition of vector 1 (I)
c     mapd2     - direct matrix of vector 2 (I)
c     mapi2     - inverse matrix of vector 2 (I)
c     poss20    - initial possition of vector 2 (I)
c     mapd3     - direct matrix of vector 3 (I)
c     mapi3     - inverse matrix of vector 3 (I)
c     poss30    - initial possition of vector 3 (I)
c     mapd4     - direct matrix of vector 4 (I)
c     mapi4     - inverse matrix of vector 4 (I)
c     poss40    - initial possition of vector 4 (I)
c     if there is less tahn 4 vecrors required
c     use any map's and poss's
c
#include "wrk.fh"
c
       integer diispoint(1:4)
       integer num
c
       integer mapd1(0:512,1:6)
       integer mapd2(0:512,1:6)
       integer mapd3(0:512,1:6)
       integer mapd4(0:512,1:6)
       integer mapi1(1:8,1:8,1:8)
       integer mapi2(1:8,1:8,1:8)
       integer mapi3(1:8,1:8,1:8)
       integer mapi4(1:8,1:8,1:8)
       integer poss10,poss20,poss30,poss40
c
c     help variables
c
       integer lun,rc
c
c
       if (num.eq.1) then
       lun=diispoint(1)
       call getmediate (wrk,wrksize,
     & lun,poss10,mapd1,mapi1,rc)
       else if (num.eq.2) then
       lun=diispoint(1)
       call getmediate (wrk,wrksize,
     & lun,poss10,mapd1,mapi1,rc)
       lun=diispoint(2)
       call getmediate (wrk,wrksize,
     & lun,poss20,mapd2,mapi2,rc)
       else if (num.eq.3) then
       lun=diispoint(1)
       call getmediate (wrk,wrksize,
     & lun,poss10,mapd1,mapi1,rc)
       lun=diispoint(2)
       call getmediate (wrk,wrksize,
     & lun,poss20,mapd2,mapi2,rc)
       lun=diispoint(3)
       call getmediate (wrk,wrksize,
     & lun,poss30,mapd3,mapi3,rc)
       else if (num.eq.4) then
       lun=diispoint(1)
       call getmediate (wrk,wrksize,
     & lun,poss10,mapd1,mapi1,rc)
       lun=diispoint(2)
       call getmediate (wrk,wrksize,
     & lun,poss20,mapd2,mapi2,rc)
       lun=diispoint(3)
       call getmediate (wrk,wrksize,
     & lun,poss30,mapd3,mapi3,rc)
       lun=diispoint(4)
       call getmediate (wrk,wrksize,
     & lun,poss40,mapd4,mapi4,rc)
       end if
c
       return
       end
