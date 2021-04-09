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
       subroutine diis (wrk,wrksize,
     & diispointt,diispointr,key)
c
c     1) increment key
c     2) if key >=firstext do:
c     Tn = DIIS (previous cycext)
c     Tn=(T21,T22,T23,T13,T14)
c     if key < firstext
c     Tn=Tn(prev)
c
c     diispointt - pointer of T stack (I)
c     diispointr - pointer of R stack (I)
c     key        - manipulation key (I/O)
c
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
c
       integer diispointt(1:4)
       integer diispointr(1:4)
       integer key
c
c     help variables
c
       real*8 rdiis1(1:4,1:4)
       real*8 cdiis(1:4)
       integer rc,lun1,nhelp
culf
      do i=1,4
         cdiis(i)=0.0
         do j=1,4
            rdiis1(i,j)=0.0
         enddo
      enddo
c1    increment key
       key=key+1
c
       if (key.lt.firstext) then
c
c     get Tn from last possition
c
c     get lun number
       lun1=diispointt(1)
c     rewind lun1
       call filemanager (2,lun1,rc)
c     T2aaaa
       call getmediate (wrk,wrksize,
     & lun1,posst210,mapdt21,mapit21,rc)
c     T2bbbb
       call getmediate (wrk,wrksize,
     & lun1,posst220,mapdt22,mapit22,rc)
c     T2abab
       call getmediate (wrk,wrksize,
     & lun1,posst230,mapdt23,mapit23,rc)
c     T1aa
       call getmediate (wrk,wrksize,
     & lun1,posst130,mapdt13,mapit13,rc)
c     T1bb
       call getmediate (wrk,wrksize,
     & lun1,posst140,mapdt14,mapit14,rc)
c     rewind lun1
       call filemanager (2,lun1,rc)
c
       return
       end if
c
c2.1  make overlap matrix
c
c2.1.1rewind R-files
       call diisrf (diispointr,cycext)
c
c2.1.2.1  read R-T21
       call diisra (wrk,wrksize,
     & diispointr,cycext,
     & mapdv1,mapiv1,possv10,mapdv2,mapiv2,possv20,
     & mapdv3,mapiv3,possv30,mapdv4,mapiv4,possv40)
c
c2.1.2.2  add overlap mtx
       call diish1 (wrk,wrksize,
     & 4,rdiis1,mapdv1,mapdv2,mapdv3,mapdv4,
     & mapiv1,mapiv2,mapiv3,mapiv4,cycext,1)
c
c
c2.1.3.1  read R-T22
       call diisra (wrk,wrksize,
     & diispointr,cycext,
     & mapdv1,mapiv1,possv10,mapdv2,mapiv2,possv20,
     & mapdv3,mapiv3,possv30,mapdv4,mapiv4,possv40)
c
c2.1.3.2  add overlap mtx
       call diish1 (wrk,wrksize,
     & 4,rdiis1,mapdv1,mapdv2,mapdv3,mapdv4,
     & mapiv1,mapiv2,mapiv3,mapiv4,cycext,0)
c
c
c2.1.4.1  read R-T23
       call diisra (wrk,wrksize,
     & diispointr,cycext,
     & mapdv1,mapiv1,possv10,mapdv2,mapiv2,possv20,
     & mapdv3,mapiv3,possv30,mapdv4,mapiv4,possv40)
c
c2.1.4.2  add overlap mtx
       call diish1 (wrk,wrksize,
     & 4,rdiis1,mapdv1,mapdv2,mapdv3,mapdv4,
     & mapiv1,mapiv2,mapiv3,mapiv4,cycext,0)
c
c
c2.1.5.1  read R-T13
       call diisra (wrk,wrksize,
     & diispointr,cycext,
     & mapdv1,mapiv1,possv10,mapdv2,mapiv2,possv20,
     & mapdv3,mapiv3,possv30,mapdv4,mapiv4,possv40)
c
c2.1.5.2  add overlap mtx
       call diish1 (wrk,wrksize,
     & 2,rdiis1,mapdv1,mapdv2,mapdv3,mapdv4,
     & mapiv1,mapiv2,mapiv3,mapiv4,cycext,0)
c
c
c2.1.6.1  read R-T14
       call diisra (wrk,wrksize,
     & diispointr,cycext,
     & mapdv1,mapiv1,possv10,mapdv2,mapiv2,possv20,
     & mapdv3,mapiv3,possv30,mapdv4,mapiv4,possv40)
c
c2.1.6.2  add overlap mtx
       call diish1 (wrk,wrksize,
     & 2,rdiis1,mapdv1,mapdv2,mapdv3,mapdv4,
     & mapiv1,mapiv2,mapiv3,mapiv4,cycext,0)
c
c
c2.2.1calc DIIS coeficients
c
       call diish2 (rdiis1,cycext,cdiis,rc)
c2.2.2write DIIS coeficients
       if (fullprint.gt.1) then
       write(6,'(6X,A,4(F9.5,2X))') 'DIIS coeficients    :',
     & (cdiis(nhelp),nhelp=1,cycext)
       end if
c
c
c2.3  make new vector
c
c2.3.1rewind T-files
       call diisrf (diispointt,cycext)
c
c2.3.2.1  read T21
       call diisra (wrk,wrksize,
     & diispointt,cycext,
     & mapdv1,mapiv1,possv10,mapdv2,mapiv2,possv20,
     & mapdv3,mapiv3,possv30,mapdv4,mapiv4,possv40)
c
c2.3.2.2  make new T21
       call diish3 (wrk,wrksize,
     & mapdt21,mapdv1,mapdv2,mapdv3,mapdv4,cdiis,cycext)
c
c2.3.3.1  read T22
       call diisra (wrk,wrksize,
     & diispointt,cycext,
     & mapdv1,mapiv1,possv10,mapdv2,mapiv2,possv20,
     & mapdv3,mapiv3,possv30,mapdv4,mapiv4,possv40)
c
c2.3.3.2  make new T22
       call diish3 (wrk,wrksize,
     & mapdt22,mapdv1,mapdv2,mapdv3,mapdv4,cdiis,cycext)
c
c2.3.4.1  read T23
       call diisra (wrk,wrksize,
     & diispointt,cycext,
     & mapdv1,mapiv1,possv10,mapdv2,mapiv2,possv20,
     & mapdv3,mapiv3,possv30,mapdv4,mapiv4,possv40)
c
c2.3.4.2  make new T23
       call diish3 (wrk,wrksize,
     & mapdt23,mapdv1,mapdv2,mapdv3,mapdv4,cdiis,cycext)
c
c2.3.5.1  read T13
       call diisra (wrk,wrksize,
     & diispointt,cycext,
     & mapdv1,mapiv1,possv10,mapdv2,mapiv2,possv20,
     & mapdv3,mapiv3,possv30,mapdv4,mapiv4,possv40)
c
c2.3.5.2  make new T13
       call diish3 (wrk,wrksize,
     & mapdt13,mapdv1,mapdv2,mapdv3,mapdv4,cdiis,cycext)
c
c2.3.6.1  read T14
       call diisra (wrk,wrksize,
     & diispointt,cycext,
     & mapdv1,mapiv1,possv10,mapdv2,mapiv2,possv20,
     & mapdv3,mapiv3,possv30,mapdv4,mapiv4,possv40)
c
c2.3.6.2  make new T14
       call diish3 (wrk,wrksize,
     & mapdt14,mapdv1,mapdv2,mapdv3,mapdv4,cdiis,cycext)
c
       return
       end
