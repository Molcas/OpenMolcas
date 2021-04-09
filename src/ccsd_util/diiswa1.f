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
       subroutine diiswa1 (wrk,wrksize,
     & diispoint)
c
c     this routine
c     a) upgrade diispoint
c     b) write write new amplitudes T21,T22,T23,T13,T14 into
c     proper possition
c
c     diispoint - array of lun's where N-1, N-2 .. amplitudes
c     are stored (I/O)
c
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
c
       integer diispoint(1:4)
c
c     help variables
c
       integer lun1,rc,p
c
c
c1    upgrade diispoint
c
       lun1=diispoint(cycext)
c
       do 10 p=cycext-1,1,-1
       diispoint(p+1)=diispoint(p)
 10     continue
c
       diispoint(1)=lun1
c
c
c2    save aplitudes
c
c2.1  rewind lun1 file
       lun1=diispoint(1)
       call filemanager (2,lun1,rc)
c
c2.2  write T21
       call wrtmediate (wrk,wrksize,
     & lun1,mapdt21,mapit21,rc)
c
c2.3  write T22
       call wrtmediate (wrk,wrksize,
     & lun1,mapdt22,mapit22,rc)
c
c2.4  write T23
       call wrtmediate (wrk,wrksize,
     & lun1,mapdt23,mapit23,rc)
c
c2.5  write T13
       call wrtmediate (wrk,wrksize,
     & lun1,mapdt13,mapit13,rc)
c
c2.6  write T14
       call wrtmediate (wrk,wrksize,
     & lun1,mapdt14,mapit14,rc)
c
c2.1  rewind lun1 file
       call filemanager (2,lun1,rc)
c
       return
       end
