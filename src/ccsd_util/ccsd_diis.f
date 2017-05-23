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
c
c     this file contains following routines:
c     diis
c     calcr
c     calcrh1
c     diiswa1
c     diisra
c     diisof
c     diiscf
c     diisrf
c     diish1
c     diish2
c     diish3
c
c     gauss
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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
c
c     ------------------------------
c
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
c
c     --------------
c
       subroutine calcrh1 (wrk,wrksize,
     & mapd1,mapd2)
c
c     this routine calc V1 = V1-V2
c
c     mapd1 - direct map of vecotr 1 (I)
c     mapd2 - direct map of vecotr 2 (I)
c
c     N.B. it is assumed, that V1 and V2 are of the same type
c
#include "wrk.fh"
c
       integer mapd1(0:512,1:6)
       integer mapd2(0:512,1:6)
c
c     help variables
c
       integer poss1,poss2,lenght,ii
c
c1    calc lenght
       ii=mapd1(0,5)
       lenght=mapd1(ii,1)+mapd1(ii,2)-mapd1(1,1)
c
c2    realize substract
       if (lenght.gt.0) then
       poss1=mapd1(1,1)
       poss2=mapd2(1,1)
       do 10 ii=0,lenght-1
       wrk(poss1+ii)=wrk(poss1+ii)-wrk(poss2+ii)
 10     continue
       end if
c
       return
       end
c
c     ------------------------------
c
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
c
c     ------------------------------
c
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
c
c     ------------------------------
c
       subroutine diisof (diispoint,ndiis)
c
c     this routine open 1-ndiis Temp files for DIIS procedure
c     and store lun's in stack (diispoint)
c
c     diispoint - stack of lun numbers (I)
c     ndiis     - size of diis procedure (I)
c
       integer ndiis
       integer diispoint(1:4)
c
c     help variables
c
       integer lun,rc
c
c
       if (ndiis.gt.0) then
       call filemanager (1,lun,rc)
       diispoint(1)=lun
       end if
c
       if (ndiis.gt.1) then
       call filemanager (1,lun,rc)
       diispoint(2)=lun
       end if
c
       if (ndiis.gt.2) then
       call filemanager (1,lun,rc)
       diispoint(3)=lun
       end if
c
       if (ndiis.gt.3) then
       call filemanager (1,lun,rc)
       diispoint(4)=lun
       end if
c
       return
       end
c
c     ------------------------------
c
       subroutine diisrf (diispoint,ndiis)
c
c     this routine rewind 1-ndiis Temp files for DIIS procedure
c     lun's are stored in stack (diispoint)
c
c     diispoint - stack of lun numbers (I)
c     ndiis     - size of diis procedure (I)
c
       integer ndiis
       integer diispoint(1:4)
c
c     help variables
c
       integer lun,rc
c
c
       if (ndiis.gt.0) then
       lun=diispoint(1)
       call filemanager (2,lun,rc)
       end if
c
       if (ndiis.gt.1) then
       lun=diispoint(2)
       call filemanager (2,lun,rc)
       end if
c
       if (ndiis.gt.2) then
       lun=diispoint(3)
       call filemanager (2,lun,rc)
       end if
c
       if (ndiis.gt.3) then
       lun=diispoint(4)
       call filemanager (2,lun,rc)
       end if
c
       return
       end
c
c     ------------------------------
c
       subroutine diiscf (diispoint,ndiis)
c
c     this routine close 1-ndiis Temp files for DIIS procedure
c     lun's are stored in stack (diispoint)
c
c     diispoint - stack of lun numbers (I)
c     ndiis     - size of diis procedure (I)
c
       integer ndiis
       integer diispoint(1:4)
c
c     help variables
c
       integer lun,rc
c
c
       if (ndiis.gt.0) then
       lun=diispoint(1)
       call filemanager (3,lun,rc)
       end if
c
       if (ndiis.gt.1) then
       lun=diispoint(2)
       call filemanager (3,lun,rc)
       end if
c
       if (ndiis.gt.2) then
       lun=diispoint(3)
       call filemanager (3,lun,rc)
       end if
c
       if (ndiis.gt.3) then
       lun=diispoint(4)
       call filemanager (3,lun,rc)
       end if
c
c
       return
       end
c
c     ------------------------------
c
       subroutine diish1 (wrk,wrksize,
     & nind,rdiis1,mapd1,mapd2,mapd3,mapd4,
     & mapi1,mapi2,mapi3,mapi4,ndiis,szkey)
c
c     this routine upgrade rdiis1(p,q) maptix
c     rdiis1(p,q) = szkey*rdiis1(p,q) + (xp|xq)
c
c     nind   - number of indexes in vectors (i)
c     rdiis1 - overlap matrix of 1-5 vectors (O)
c     mapd1  - direct map of vector 1 (i)
c     mapd2  - direct map of vector 2 (i)
c     mapd3  - direct map of vector 3 (i)
c     mapd4  - direct map of vector 4 (i)
c     mapi1  - inverse map of vector 1 (i)
c     mapi2  - inverse map of vector 2 (i)
c     mapi3  - inverse map of vector 3 (i)
c     mapi4  - inverse map of vector 4 (i)
c     if tere is less than 5 vectors, use any map
c     ndiis  - dimension of DIIS (2-4) (I)
c     szkey  - 0 - no vanishing rdiis1 at the beggining
c     1 - vanishing rdiis1 at the beggining
c
#include "wrk.fh"
#include "ccsd1.fh"
c
       real*8 rdiis1(1:4,1:4)
       integer mapd1(0:512,1:6)
       integer mapd2(0:512,1:6)
       integer mapd3(0:512,1:6)
       integer mapd4(0:512,1:6)
       integer mapi1(1:8,1:8,1:8)
       integer mapi2(1:8,1:8,1:8)
       integer mapi3(1:8,1:8,1:8)
       integer mapi4(1:8,1:8,1:8)
c
       integer nind,ndiis,szkey
c
c     help variables
c
       integer nhelp
       real*8 scalar
       integer rc,num
c
c
       num=ndiis+1
c
       if (szkey.eq.1) then
       nhelp=4*4
       call mv0zero(nhelp,nhelp,rdiis1)
       end if
c
       if (num.gt.0) then
c     calc X11
       call multdot (wrk,wrksize,
     & nind,mapd1,mapi1,1,mapd1,mapi1,1,scalar,rc)
       rdiis1(1,1)=rdiis1(1,1)+scalar
       end if
c
       if (num.gt.1) then
c     X21
       call multdot (wrk,wrksize,
     & nind,mapd2,mapi2,1,mapd1,mapi1,1,scalar,rc)
       rdiis1(2,1)=rdiis1(2,1)+scalar
       rdiis1(1,2)=rdiis1(1,2)+scalar
c     X22
       call multdot (wrk,wrksize,
     & nind,mapd2,mapi2,1,mapd2,mapi2,1,scalar,rc)
       rdiis1(2,2)=rdiis1(2,2)+scalar
       end if
c
       if (num.gt.2) then
c     X31
       call multdot (wrk,wrksize,
     & nind,mapd3,mapi3,1,mapd1,mapi1,1,scalar,rc)
       rdiis1(3,1)=rdiis1(3,1)+scalar
       rdiis1(1,3)=rdiis1(1,3)+scalar
c     X32
       call multdot (wrk,wrksize,
     & nind,mapd3,mapi3,1,mapd2,mapi2,1,scalar,rc)
       rdiis1(3,2)=rdiis1(3,2)+scalar
       rdiis1(2,3)=rdiis1(2,3)+scalar
c     X33
       call multdot (wrk,wrksize,
     & nind,mapd3,mapi3,1,mapd3,mapi3,1,scalar,rc)
       rdiis1(3,3)=rdiis1(3,3)+scalar
       end if
c
       if (num.gt.3) then
c     X41
       call multdot (wrk,wrksize,
     & nind,mapd4,mapi4,1,mapd1,mapi1,1,scalar,rc)
       rdiis1(4,1)=rdiis1(4,1)+scalar
       rdiis1(1,4)=rdiis1(1,4)+scalar
c     X42
       call multdot (wrk,wrksize,
     & nind,mapd4,mapi4,1,mapd2,mapi2,1,scalar,rc)
       rdiis1(4,2)=rdiis1(4,2)+scalar
       rdiis1(2,4)=rdiis1(2,4)+scalar
c     X43
       call multdot (wrk,wrksize,
     & nind,mapd4,mapi4,1,mapd3,mapi3,1,scalar,rc)
       rdiis1(4,3)=rdiis1(4,3)+scalar
       rdiis1(3,4)=rdiis1(3,4)+scalar
c     X44
       call multdot (wrk,wrksize,
     & nind,mapd4,mapi4,1,mapd4,mapi4,1,scalar,rc)
       rdiis1(4,4)=rdiis1(4,4)+scalar
       end if
c
       return
       end
c
c     ----------------------
c
       subroutine diish2 (rdiis1,ndiis,cdiis,rc)
c
c     this rouine calculate diis coeficients by solving
c
c     B-1   c  =  0
c     -1 0   l    -1
c
c     and consequent normalization of coeficients
c
c
c     r1diis  - matrix of amp. overlap of ndiid+1 iterations (I)
c     ndiis   - size of diis (2-4) (I)
c     cdiis   - final diis coeficients (O)
c     rc      - return (error) code (O)
c     0 - OK
c     1 - singular DIIS matrix
c
c
       real*8 rdiis1(1:4,1:4)
       real*8 cdiis (1:4)
       integer ndiis,rc
c
c     help variables
c
       integer p,q
       real*8 scalar
       real*8 rdiis2(1:5,1:5)
       real*8 bb(1:5)
       real*8 ci(1:5)
c
c1    vanish rdiis2 file
       call mv0zero (25,25,rdiis2)
c
c2.1  make rdiis2 matrix
c
       do 10 p=1,ndiis
       do 10 q=1,ndiis
       rdiis2(p,q)=rdiis1(p,q)
 10     continue
c
       do 11 p=1,ndiis
       rdiis2(p,ndiis+1)=-1.0d0
       rdiis2(ndiis+1,p)=-1.0d0
       bb(p)=0.0d0
 11     continue
c
       bb(ndiis+1)=-1.0d0
c
c2.2  modify rdiis2 matrix
c
c     scale matrix
       scalar=sqrt(rdiis2(1,1)*rdiis2(ndiis,ndiis))
       do 12 q=1,ndiis
       do 12 p=1,ndiis
       rdiis2(p,q)=rdiis2(p,q)/scalar
 12     continue
c
c     add penalty function
       scalar=0.01d0*rdiis2(ndiis,ndiis)
c     rdiis2(ndiis,ndiis)=rdiis2(ndiis,ndiis)+scalar
c     bb(ndiis)=scalar
c
c
c3    solve SLE
c
c3.1  vanish ci
       do 20 p=1,ndiis+1
       ci(p)=0.0d0
 20     continue
c
c3.2  get cdiis
       call gauss (ndiis+1,5,rdiis2,ci,bb)
c
c4    final modification of cdiis coeficients
c
CFUE   if (rc.eq.1) then
c     matrix R2 was singular, no extrapolation
CFUE   write(6,*) ' SINGULAR DIIS MATRIX, NO EXTRAPOLATION'
CFUE   cdiis(1)=1.0d0
CFUE   do p=2,ndiis
CFUE   cdiis(p)=0.0d0
CFUE   end do
c
CFUE   else
c     DIIS procedure was succesfull, renormalize coef.
c
       scalar=0.0d0
       do 40 p=1,ndiis
       scalar=scalar+ci(p)
 40     continue
c
       do 50 p=1,ndiis
       cdiis(p)=ci(p)/scalar
 50     continue
c
       scalar=0.0d0
       do 60 p=1,ndiis
       scalar=scalar+cdiis(p)
 60     continue
c
CFUE   end if
c
CFUE   write(6,*) cdiis(1),cdiis(2),cdiis(3),cdiis(4),scalar
c51     format (5(i2,d12.7))
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(rc)
       end
c
c     -------------------------------------
c
       subroutine diish3 (wrk,wrksize,
     & mapd0,mapd1,mapd2,mapd3,mapd4,cdiis,ndiis)
c
c     this routine produce new vector
c     v(0) = sum(i) [cdiis(i) . V(i)]
c
c     mapd0  - direct map of V0 vector (I)
c     mapd1  - direct map of V1 vector (I)
c     mapd2  - direct map of V2 vector (I)
c     mapd3  - direct map of V3 vector (I)
c     mapd4  - direct map of V4 vector (I)
c     cdiis  - vector of diis coeficients (I)
c     ndiis  - size of diis (2-4) (I)
c
#include "wrk.fh"
c
       integer mapd0(0:512,1:6)
       integer mapd1(0:512,1:6)
       integer mapd2(0:512,1:6)
       integer mapd3(0:512,1:6)
       integer mapd4(0:512,1:6)
       real*8 cdiis(1:4)
       integer ndiis
c
c     help variables
c
       integer poss0,poss1,poss2,poss3,poss4,nhelp,lenght
c
c
       if (ndiis.eq.2) then
c     2 dimensional DIIS
c
       poss0=mapd0(1,1)
       poss1=mapd1(1,1)
       poss2=mapd2(1,1)
c
       nhelp=mapd1(0,5)
       lenght=mapd1(nhelp,1)+mapd1(nhelp,2)-mapd1(1,1)
c
       if (lenght.gt.0) then
       do 20 nhelp=0,lenght-1
       wrk(poss0+nhelp)=cdiis(1)*wrk(poss1+nhelp)
     & +cdiis(2)*wrk(poss2+nhelp)
 20     continue
       end if
c
       else if (ndiis.eq.3) then
c     3 dimensional DIIS
c
       poss0=mapd0(1,1)
       poss1=mapd1(1,1)
       poss2=mapd2(1,1)
       poss3=mapd3(1,1)
c
       nhelp=mapd1(0,5)
       lenght=mapd1(nhelp,1)+mapd1(nhelp,2)-mapd1(1,1)
c
       if (lenght.gt.0) then
       do 30 nhelp=0,lenght-1
       wrk(poss0+nhelp)=cdiis(1)*wrk(poss1+nhelp)
     & +cdiis(2)*wrk(poss2+nhelp)
     & +cdiis(3)*wrk(poss3+nhelp)
 30     continue
       end if
c
       else if (ndiis.eq.4) then
c     4 dimensional DIIS
c
       poss0=mapd0(1,1)
       poss1=mapd1(1,1)
       poss2=mapd2(1,1)
       poss3=mapd3(1,1)
       poss4=mapd4(1,1)
c
       nhelp=mapd1(0,5)
       lenght=mapd1(nhelp,1)+mapd1(nhelp,2)-mapd1(1,1)
c
       if (lenght.gt.0) then
       do 40 nhelp=0,lenght-1
       wrk(poss0+nhelp)=cdiis(1)*wrk(poss1+nhelp)
     & +cdiis(2)*wrk(poss2+nhelp)
     & +cdiis(3)*wrk(poss3+nhelp)
     & +cdiis(4)*wrk(poss4+nhelp)
 40     continue
       end if
c
       end if
c
       return
       end
c
