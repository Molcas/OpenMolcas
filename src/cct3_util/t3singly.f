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
c
c     t3sgl
c     t3sgl11
c     t3sgl121
c     t3sgl122
c     t3sgl131
c     t3sgl132
c     t3sgl141
c     t3sgl142
c     t3sgl143
c     t3sgl211
c     t3sgl212
c     t3sgl221
c     t3sgl222
c     t3sgl223
c     t3sgl311
c     t3sgl312
c     t3sgl321
c     t3sgl322
c     t3sgl323
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine t3sgl (wrk,wrksize,
     & mapdw,ssw,mapds1,mapis1,mapds2,mapis2,
     & mapdd1,mapid1,mapdd2,mapid2,
     & typdiv,i,j,k,symi,symj,symk,rc1,
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)

c
c     mapdw  - direct map matrix of W (Input)
c     ssw    - overall symmetry state of matrix W (V) (Input)
c     mapds1 - direct map matrix of S1 (Input)
c     mapis1 - inverse map matrix of S1 (Input)
c     mapds2 - direct map matrix of S2 (Input)
c     mapis2 - inverse map matrix of S2 (Input)
c     mapdd1 - direct map matrix of D1 (Input)
c     mapid1 - inverse map matrix of D1 (Input)
c     mapdd2 - direct map matrix of D2 (Input)
c     mapid2 - inverse map matrix of D2 (Input)
c     (order is aa>ab>bb,
c     if there is only one spin, use any map's for 2)
c     typsgl - typ of operation (see Table) (Input)
c     i      - value of occupied index i (Inlut)
c     j      - value of occupied index j (Inlut)
c     k      - value of occupied index k (Inlut)
c     symi   - symmetry of index i (Input)
c     symj   - symmetry of index j (Input)
c     symk   - symmetry of index k (Input)
c     rc     - return (error) code (Output)
c     mapd,mapi,poss - parameterrs for M1-3,H1-3 files (I)
c
c     this routine add contributions from diconnected
c     singles, namely:
c
c     W_ijk(abc) <- P [ U1(_i,a) . U2(_jk,bc) ]
c
c     for following types of W
c
c     typsgl         Operation                   Implemented
c     1     W(abc) <- +s1(_i,a) . d1(_jk,bc)
c     -s1(_i,b) . d1(_jk,ac)
c     +s1(_i,c) . d1(_jk,ab)
c
c     -s1(_j,a) . d1(_ik,bc)
c     +s1(_j,b) . d1(_ik,ac)
c     -s1(_j,c) . d1(_ik,ab)
c
c     +s1(_k,a) . d1(_ij,bc)
c     -s1(_k,b) . d1(_ij,ac)
c     +s1(_k,c) . d1(_ij,ab)        Yes
c
c
c     2     W(ab,c)<- +s1(_i,a) . d2(_j,_k,b,c)
c     -s1(_i,b) . d2(_j,_k,a,c)
c
c     -s1(_j,a) . d2(_i,_k,b,c)
c     +s1(_j,b) . d2(_i,_k,a,c)
c
c     +s2(_k,c) . d1(_ij,ab)        Yes
c
c
c     3     W(a,bc)<- +s1(_i,a) . d2(_jk,b,c)
c
c     +s1(_j,b) . d2(_i,_k,a,c)
c     -s1(_j,c) . d2(_i,_k,a,b)
c
c     -s1(_k,b) . d2(_i,_j,a,c)
c     +s1(_k,c) . d2(_i,_j,a,b)     Yes
c
c
c
c     N.B. spin combinations aaa,bbb for 1; aab for 2; and abb for 3
c     are authomatically assumed
c
#include "t31.fh"
#include "wrk.fh"
c
       integer ssw,typdiv,i,j,k,symi,symj,symk,rc
       integer mapdw(0:512,1:6)
       integer mapds1(0:512,1:6)
       integer mapds2(0:512,1:6)
       integer mapdd1(0:512,1:6)
       integer mapdd2(0:512,1:6)
       integer mapis1(1:8,1:8,1:8)
       integer mapis2(1:8,1:8,1:8)
       integer mapid1(1:8,1:8,1:8)
       integer mapid2(1:8,1:8,1:8)
c
       integer possm10
       integer mapdm1(0:512,1:6)
       integer mapim1(1:8,1:8,1:8)
       integer possm20
       integer mapdm2(0:512,1:6)
       integer mapim2(1:8,1:8,1:8)
       integer possm30
       integer mapdm3(0:512,1:6)
       integer mapim3(1:8,1:8,1:8)
       integer possh10
       integer mapdh1(0:512,1:6)
       integer mapih1(1:8,1:8,1:8)
       integer possh20
       integer mapdh2(0:512,1:6)
       integer mapih2(1:8,1:8,1:8)
       integer possh30
       integer mapdh3(0:512,1:6)
       integer mapih3(1:8,1:8,1:8)
c
c     help variables
c
       integer iw,possw
       integer id1,id2,id3,possd1,possd2,possd3
       integer is1,is2,is3,posss1,posss2,posss3
       integer syma,symb,symc,dima,dimb,dimc
       integer nhelp1,nhelp2
       integer rc1,ssh1,ssh2,ssh3,ssm1,ssm2,ssm3
       integer symij,symik,symjk,symab,symac,symbc
c
c
c0.*  some tests
c
c
c
       if (typdiv.eq.1) then
c
c1    case W(pqr)
c
c1.*  ext H1(a) <= S1(a,i) for given i
       call ext(wrk,wrksize,
     & 2,2,i,0,0,symi,0,0,mapds1,mapis1,1,
     & possh10,mapdh1,mapih1,ssh1,rc1)
c
c1.*  ext H2(a) <= S1(a,j) for given j
       call ext(wrk,wrksize,
     & 2,2,j,0,0,symj,0,0,mapds1,mapis1,1,
     & possh20,mapdh2,mapih2,ssh2,rc1)
c
c1.*  ext H3(a) <= S1(a,k) for given k
       call ext(wrk,wrksize,
     & 2,2,k,0,0,symk,0,0,mapds1,mapis1,1,
     & possh30,mapdh3,mapih3,ssh3,rc1)
c
c1.*  ext M1(bc) <= D1(bc,jk) for given jk
       call ext(wrk,wrksize,
     & 4,7,j,k,0,symj,symk,0,mapdd1,mapid1,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c1.*  ext M2(bc) <= D1(bc,ik) for given ik
       call ext(wrk,wrksize,
     & 4,7,i,k,0,symi,symk,0,mapdd1,mapid1,1,
     & possm20,mapdm2,mapim2,ssm2,rc1)
c
c1.*  ext M3(bc) <= D1(bc,ij) for given ij
       call ext(wrk,wrksize,
     & 4,7,i,j,0,symi,symj,0,mapdd1,mapid1,1,
     & possm30,mapdm3,mapim3,ssm3,rc1)
c
c
       do 100 iw=1,mapdw(0,5)
c
c1.*  def possition of W
       possw=mapdw(iw,1)
c
c1.*  def symmetry status
       syma=mapdw(iw,3)
       symb=mapdw(iw,4)
       symc=mapdw(iw,5)
c
c1.*  def dimensions
       dima=dimm(mapdw(0,1),syma)
       dimb=dimm(mapdw(0,2),symb)
       dimc=dimm(mapdw(0,3),symc)
c1.*
       symij=mmul(symi,symj)
       symik=mmul(symi,symk)
       symjk=mmul(symj,symk)
       symab=mmul(syma,symb)
       symac=mmul(syma,symc)
       symbc=mmul(symb,symc)
c
c1.*  realize packing
c
       if (syma.eq.symc) then
c1.a  case syma=symb=symc
c
c
c1.a.1----  first triade   W(abc) <- +s1(_i,a) . d1(_jk,bc)
c     -s1(_i,b) . d1(_jk,ac)
c     +s1(_i,c) . d1(_jk,ab)
c
       if ((symi.eq.syma).and.(symjk.eq.symbc)) then
c1.a.1.*      find address for s1,d1
       is1=mapih1(1,1,1)
       id1=mapim1(symb,1,1)
c
c1.a.1.*      def possition of s1,d1
       posss1=mapdh1(is1,1)
       possd1=mapdm1(id1,1)
c
c1.a.1.*      def additional dimensions
       nhelp1=dima*(dima-1)/2
       nhelp2=dima*(dima-1)*(dima-2)/6
c
c1.a.1.*      add singly
       call t3sglh11 (wrk(possw),dima,nhelp1,nhelp2,
     & wrk(posss1),wrk(possd1),1)
       end if
c
c
c1.a.2----  2-nd. triade   W(abc) <- -s1(_j,a) . d1(_ik,bc)
c     +s1(_j,b) . d1(_ik,ac)
c     -s1(_j,c) . d1(_ik,ab)
c
       if ((symj.eq.syma).and.(symik.eq.symbc)) then
c1.a.2.*      find address for s1,d1
       is1=mapih2(1,1,1)
       id1=mapim2(symb,1,1)
c
c1.a.2.*      def possition of s1,d1
       posss1=mapdh2(is1,1)
       possd1=mapdm2(id1,1)
c
c1.a.2.*      def additional dimensions
       nhelp1=dima*(dima-1)/2
       nhelp2=dima*(dima-1)*(dima-2)/6
c
c1.a.2.*      add singly
       call t3sglh11 (wrk(possw),dima,nhelp1,nhelp2,
     & wrk(posss1),wrk(possd1),-1)
       end if
c
c
c1.a.3----  3-rd. triade   W(abc) <- +s1(_k,a) . d1(_ij,bc)
c     -s1(_k,b) . d1(_ij,ac)
c     +s1(_k,c) . d1(_ij,ab)
c
       if ((symk.eq.syma).and.(symij.eq.symbc)) then
c
c1.a.3.*      find address for s1,d1
       is1=mapih3(1,1,1)
       id1=mapim3(symb,1,1)
c
c1.a.3.*      def possition of s1,d1
       posss1=mapdh3(is1,1)
       possd1=mapdm3(id1,1)
c
c1.a.3.*      def additional dimensions
       nhelp1=dima*(dima-1)/2
       nhelp2=dima*(dima-1)*(dima-2)/6
c
c1.a.3.*      add singly
       call t3sglh11 (wrk(possw),dima,nhelp1,nhelp2,
     & wrk(posss1),wrk(possd1),1)
       end if
c
c
       else if (syma.eq.symb) then
c1.b  case syma=symb.ne.symc
c
c
c1.b.1----  first triade   W(abc) <- +s1(_i,a) . d1(_jk,bc)
c     -s1(_i,b) . d1(_jk,ac)
c     +s1(_i,c) . d1(_jk,ab)
c
       if ((symi.eq.syma).and.(symjk.eq.symbc)) then
c     if syma=symi, then obviously symc.ne.symi
c1.b.1.*      find address for s1,d1
       is1=mapih1(1,1,1)
       id1=mapim1(syma,1,1)
c
c1.b.1.*      def possition of s1,d1
       posss1=mapdh1(is1,1)
       possd1=mapdm1(id1,1)
c
c1.b.1.*      def additional dimensions
       nhelp1=dima*(dima-1)/2
c
c1.b.1.*      add singly
       call t3sglh121 (wrk(possw),dima,nhelp1,dimc,
     & wrk(posss1),wrk(possd1),1)
       else if ((symc.eq.symi).and.(symjk.eq.symab)) then
c1.b.1.*      find address for s3,d3
       is3=mapih1(1,1,1)
       id3=mapim1(syma,1,1)
c
c1.b.1.*      def possition of s1,d3
       posss3=mapdh1(is3,1)
       possd3=mapdm1(id3,1)
c
c1.b.1.*      def additional dimensions
       nhelp1=dima*(dima-1)/2
c
c1.b.1.*      add singly
       call t3sglh122 (wrk(possw),dima,nhelp1,dimc,
     & wrk(posss3),wrk(possd3),1)
       end if
c
c
c1.b.2----  2-nd. triade   W(abc) <- -s1(_j,a) . d1(_ik,bc)
c     +s1(_j,b) . d1(_ik,ac)
c     -s1(_j,c) . d1(_ik,ab)
c
       if ((symj.eq.syma).and.(symik.eq.symbc)) then
c     if syma=symj, then obviously symc.ne.symj
c1.b.2.*      find address for s1,d1
       is1=mapih2(1,1,1)
       id1=mapim2(syma,1,1)
c
c1.b.2.*      def possition of s1,d1
       posss1=mapdh2(is1,1)
       possd1=mapdm2(id1,1)
c
c1.b.2.*      def additional dimensions
       nhelp1=dima*(dima-1)/2
c
c1.b.2.*      add singly
       call t3sglh121 (wrk(possw),dima,nhelp1,dimc,
     & wrk(posss1),wrk(possd1),-1)
       else if ((symc.eq.symj).and.(symik.eq.symab)) then
c1.b.2.*      find address for s3,d3
       is3=mapih2(1,1,1)
       id3=mapim2(syma,1,1)
c
c1.b.2.*      def possition of s1,d3
       posss3=mapdh2(is3,1)
       possd3=mapdm2(id3,1)
c
c1.b.2.*      def additional dimensions
       nhelp1=dima*(dima-1)/2
c
c1.b.2.*      add singly
       call t3sglh122 (wrk(possw),dima,nhelp1,dimc,
     & wrk(posss3),wrk(possd3),-1)
       end if
c
c
c1.b.3----  3-rd. triade   W(abc) <- +s1(_k,a) . d1(_ij,bc)
c     -s1(_k,b) . d1(_ij,ac)
c     +s1(_k,c) . d1(_ij,ab)
c
       if ((symk.eq.syma).and.(symij.eq.symbc)) then
c     if syma=symk, then obviously symc.ne.symk
c1.b.3.*      find address for s1,d1
       is1=mapih3(1,1,1)
       id1=mapim3(syma,1,1)
c
c1.b.3.*      def possition of s1,d1
       posss1=mapdh3(is1,1)
       possd1=mapdm3(id1,1)
c
c1.b.3.*      def additional dimensions
       nhelp1=dima*(dima-1)/2
c
c1.b.3.*      add singly
       call t3sglh121 (wrk(possw),dima,nhelp1,dimc,
     & wrk(posss1),wrk(possd1),1)
       else if ((symc.eq.symk).and.(symij.eq.symab)) then
c1.b.3.*      find address for s3,d3
       is3=mapih3(1,1,1)
       id3=mapim3(syma,1,1)
c
c1.b.3.*      def possition of s1,d3
       posss3=mapdh3(is3,1)
       possd3=mapdm3(id3,1)
c
c1.b.3.*      def additional dimensions
       nhelp1=dima*(dima-1)/2
c
c1.b.3.*      add singly
       call t3sglh122 (wrk(possw),dima,nhelp1,dimc,
     & wrk(posss3),wrk(possd3),1)
       end if
c
c
       else if (symb.eq.symc) then
c1.c  case syma.ne.symb=symc
c
c1.c.1----  first triade   W(abc) <- +s1(_i,a) . d1(_jk,bc)
c     -s1(_i,b) . d1(_jk,ac)
c     +s1(_i,c) . d1(_jk,ab)
c
       if ((symi.eq.syma).and.(symjk.eq.symbc)) then
c     if syma=symi, then obviously symb(symc).ne.symi
c1.c.1.*      find address for s1,d1
       is1=mapih1(1,1,1)
       id1=mapim1(symb,1,1)
c
c1.c.1.*      def possition of s1,d1
       posss1=mapdh1(is1,1)
       possd1=mapdm1(id1,1)
c
c1.c.1.*      def additional dimensions
       nhelp1=dimb*(dimb-1)/2
c
c1.c.1.*      add singly
       call t3sglh131 (wrk(possw),dima,dimb,nhelp1,
     & wrk(posss1),wrk(possd1),1)
       else if ((symb.eq.symi).and.(symjk.eq.symac)) then
c1.c.1.*      find address for s3,d3
       is2=mapih1(1,1,1)
       id2=mapim1(syma,1,1)
c
c1.c.1.*      def possition of s1,d3
       posss2=mapdh1(is2,1)
       possd2=mapdm1(id2,1)
c
c1.c.1.*      def additional dimensions
       nhelp1=dimb*(dimb-1)/2
c
c1.c.1.*      add singly
       call t3sglh132 (wrk(possw),dima,dimb,nhelp1,
     & wrk(posss2),wrk(possd2),1)
       end if
c
c
c1.c.2----  2-nd. triade   W(abc) <- -s1(_j,a) . d1(_ik,bc)
c     +s1(_j,b) . d1(_ik,ac)
c     -s1(_j,c) . d1(_ik,ab)
c
       if ((symj.eq.syma).and.(symik.eq.symbc)) then
c     if syma=symj, then obviously symb(symc).ne.symi
c1.c.2.*      find address for s1,d1
       is1=mapih2(1,1,1)
       id1=mapim2(symb,1,1)
c
c1.c.2.*      def possition of s1,d1
       posss1=mapdh2(is1,1)
       possd1=mapdm2(id1,1)
c
c1.c.2.*      def additional dimensions
       nhelp1=dimb*(dimb-1)/2
c
c1.c.2.*      add singly
       call t3sglh131 (wrk(possw),dima,dimb,nhelp1,
     & wrk(posss1),wrk(possd1),-1)
       else if ((symb.eq.symj).and.(symik.eq.symac)) then
c1.c.2.*      find address for s3,d3
       is2=mapih2(1,1,1)
       id2=mapim2(syma,1,1)
c
c1.c.2.*      def possition of s1,d3
       posss2=mapdh2(is2,1)
       possd2=mapdm2(id2,1)
c
c1.c.2.*      def additional dimensions
       nhelp1=dimb*(dimb-1)/2
c
c1.c.2.*      add singly
       call t3sglh132 (wrk(possw),dima,dimb,nhelp1,
     & wrk(posss2),wrk(possd2),-1)
       end if
c
c
c1.c.3----  3-rd. triade   W(abc) <- +s1(_k,a) . d1(_ij,bc)
c     -s1(_k,b) . d1(_ij,ac)
c     +s1(_k,c) . d1(_ij,ab)
c
       if ((symk.eq.syma).and.(symij.eq.symbc)) then
c     if syma=symk, then obviously symb(symc).ne.symk
c1.c.3.*      find address for s1,d1
       is1=mapih3(1,1,1)
       id1=mapim3(symb,1,1)
c
c1.c.3.*      def possition of s1,d1
       posss1=mapdh3(is1,1)
       possd1=mapdm3(id1,1)
c
c1.c.3.*      def additional dimensions
       nhelp1=dimb*(dimb-1)/2
c
c1.c.3.*      add singly
       call t3sglh131 (wrk(possw),dima,dimb,nhelp1,
     & wrk(posss1),wrk(possd1),1)
       else if ((symb.eq.symk).and.(symij.eq.symac)) then
c1.c.3.*      find address for s3,d3
       is2=mapih3(1,1,1)
       id2=mapim3(syma,1,1)
c
c1.c.3.*      def possition of s1,d3
       posss2=mapdh3(is2,1)
       possd2=mapdm3(id2,1)
c
c1.c.3.*      def additional dimensions
       nhelp1=dimb*(dimb-1)/2
c
c1.c.3.*      add singly
       call t3sglh132 (wrk(possw),dima,dimb,nhelp1,
     & wrk(posss2),wrk(possd2),1)
       end if
c
c
       else
c1.d  case syma.ne.symb.ne.symc
c
c1.d.1----  first triade   W(abc) <- +s1(_i,a) . d1(_jk,bc)
c     -s1(_i,b) . d1(_jk,ac)
c     +s1(_i,c) . d1(_jk,ab)
c
       if ((symi.eq.syma).and.(symjk.eq.symbc)) then
c1.d.1.*      find address for s1,d1
       is1=mapih1(1,1,1)
       id1=mapim1(symb,1,1)
c
c1.d.1.*      def possition of s1,d1
       posss1=mapdh1(is1,1)
       possd1=mapdm1(id1,1)
c
c1.d.1.*      add singly
       call t3sglh141 (wrk(possw),dima,dimb,dimc,
     & wrk(posss1),wrk(possd1),1)
       else if ((symi.eq.symb).and.(symjk.eq.symac)) then
c1.d.1.*      find address for s2,d2
       is2=mapih1(1,1,1)
       id2=mapim1(syma,1,1)
c
c1.d.1.*      def possition of s1,d1
       posss2=mapdh1(is2,1)
       possd2=mapdm1(id2,1)
c
c1.d.1.*      add singly
       call t3sglh142 (wrk(possw),dima,dimb,dimc,
     & wrk(posss2),wrk(possd2),1)
       else if ((symi.eq.symc).and.(symjk.eq.symab)) then
c1.d.1.*      find address for s1,d1
       is3=mapih1(1,1,1)
       id3=mapim1(syma,1,1)
c
c1.d.1.*      def possition of s1,d1
       posss3=mapdh1(is3,1)
       possd3=mapdm1(id3,1)
c
c1.d.1.*      add singly
       call t3sglh143 (wrk(possw),dima,dimb,dimc,
     & wrk(posss3),wrk(possd3),1)
       end if
c
c
c1.d.2----  2-nd. triade   W(abc) <- -s1(_j,a) . d1(_ik,bc)
c     +s1(_j,b) . d1(_ik,ac)
c     -s1(_j,c) . d1(_ik,ab)
c
       if ((symj.eq.syma).and.(symik.eq.symbc)) then
c1.d.2.*      find address for s1,d1
       is1=mapih2(1,1,1)
       id1=mapim2(symb,1,1)
c
c1.d.2.*      def possition of s1,d1
       posss1=mapdh2(is1,1)
       possd1=mapdm2(id1,1)
c
c1.d.2.*      add singly
       call t3sglh141 (wrk(possw),dima,dimb,dimc,
     & wrk(posss1),wrk(possd1),-1)
       else if ((symj.eq.symb).and.(symik.eq.symac)) then
c1.d.2.*      find address for s2,d2
       is2=mapih2(1,1,1)
       id2=mapim2(syma,1,1)
c
c1.d.2.*      def possition of s1,d1
       posss2=mapdh2(is2,1)
       possd2=mapdm2(id2,1)
c
c1.d.2.*      add singly
       call t3sglh142 (wrk(possw),dima,dimb,dimc,
     & wrk(posss2),wrk(possd2),-1)
       else if ((symj.eq.symc).and.(symik.eq.symab)) then
c1.d.2.*      find address for s1,d1
       is3=mapih2(1,1,1)
       id3=mapim2(syma,1,1)
c
c1.d.2.*      def possition of s1,d1
       posss3=mapdh2(is3,1)
       possd3=mapdm2(id3,1)
c
c1.d.2.*      add singly
       call t3sglh143 (wrk(possw),dima,dimb,dimc,
     & wrk(posss3),wrk(possd3),-1)
       end if
c
c
c1.d.3----  3-rd. triade   W(abc) <- +s1(_k,a) . d1(_ij,bc)
c     -s1(_k,b) . d1(_ij,ac)
c     +s1(_k,c) . d1(_ij,ab)
c
       if ((symk.eq.syma).and.(symij.eq.symbc)) then
c1.d.3.*      find address for s1,d1
       is1=mapih3(1,1,1)
       id1=mapim3(symb,1,1)
c
c1.d.3.*      def possition of s1,d1
       posss1=mapdh3(is1,1)
       possd1=mapdm3(id1,1)
c
c1.d.3.*      add singly
       call t3sglh141 (wrk(possw),dima,dimb,dimc,
     & wrk(posss1),wrk(possd1),1)
       else if ((symk.eq.symb).and.(symij.eq.symac)) then
c1.d.3.*      find address for s2,d2
       is2=mapih3(1,1,1)
       id2=mapim3(syma,1,1)
c
c1.d.3.*      def possition of s1,d1
       posss2=mapdh3(is2,1)
       possd2=mapdm3(id2,1)
c
c1.d.3.*      add singly
       call t3sglh142 (wrk(possw),dima,dimb,dimc,
     & wrk(posss2),wrk(possd2),1)
       else if ((symk.eq.symc).and.(symij.eq.symab)) then
c1.d.3.*      find address for s1,d1
       is3=mapih3(1,1,1)
       id3=mapim3(syma,1,1)
c
c1.d.3.*      def possition of s1,d1
       posss3=mapdh3(is3,1)
       possd3=mapdm3(id3,1)
c
c1.d.3.*      add singly
       call t3sglh143 (wrk(possw),dima,dimb,dimc,
     & wrk(posss3),wrk(possd3),1)
       end if
c
c
       end if
c
 100    continue
c
c
       else if (typdiv.eq.2) then
c2    case W(pq,r)
c
c2.*  ext H1(a) <= S1(a,i) for given i
       call ext(wrk,wrksize,
     & 2,2,i,0,0,symi,0,0,mapds1,mapis1,1,
     & possh10,mapdh1,mapih1,ssh1,rc1)
c
c2.*  ext H2(a) <= S1(a,j) for given j
       call ext(wrk,wrksize,
     & 2,2,j,0,0,symj,0,0,mapds1,mapis1,1,
     & possh20,mapdh2,mapih2,ssh2,rc1)
c
c2.*  ext H3(a) <= S2(a,k) for given k
       call ext(wrk,wrksize,
     & 2,2,k,0,0,symk,0,0,mapds2,mapis2,1,
     & possh30,mapdh3,mapih3,ssh3,rc1)
c
c2.*  ext M1(bc) <= D2(bc,jk) for given jk
       call ext(wrk,wrksize,
     & 4,7,j,k,0,symj,symk,0,mapdd2,mapid2,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c2.*  ext M2(bc) <= D2(bc,ik) for given ik
       call ext(wrk,wrksize,
     & 4,7,i,k,0,symi,symk,0,mapdd2,mapid2,1,
     & possm20,mapdm2,mapim2,ssm2,rc1)
c
c2.*  ext M3(bc) <= D1(bc,ij) for given ij
       call ext(wrk,wrksize,
     & 4,7,i,j,0,symi,symj,0,mapdd1,mapid1,1,
     & possm30,mapdm3,mapim3,ssm3,rc1)
c
c
       do 200 iw=1,mapdw(0,5)
c
c2.*  def possition of W
       possw=mapdw(iw,1)
c
c2.*  def symmetry status
       syma=mapdw(iw,3)
       symb=mapdw(iw,4)
       symc=mapdw(iw,5)
c
c2.*  def dimensions
       dima=dimm(mapdw(0,1),syma)
       dimb=dimm(mapdw(0,2),symb)
       dimc=dimm(mapdw(0,3),symc)
c2.*
       symij=mmul(symi,symj)
       symik=mmul(symi,symk)
       symjk=mmul(symj,symk)
       symab=mmul(syma,symb)
       symac=mmul(syma,symc)
       symbc=mmul(symb,symc)
c
c2.*  add singles
       if (syma.eq.symb) then
c2.a  case syma=symb,symc
c
c2.a.11.st.diade  W(ab,c)<- +s1(_i,a) . d2(_j,_k,b,c)
c     -s1(_i,b) . d2(_j,_k,a,c)
c
       if ((syma.eq.symi).and.(symjk.eq.symbc)) then
c2.a.1.*      find address for s1,d1
       is1=mapih1(1,1,1)
       id1=mapim1(syma,1,1)
c
c2.a.1.*      def possition of s1,d1
       posss1=mapdh1(is1,1)
       possd1=mapdm1(id1,1)
c
c2.a.1.*      def additional dimensions
       nhelp1=dima*(dima-1)/2
c
c2.a.1.*      add singly
       call t3sglh211 (wrk(possw),dima,nhelp1,dimc,
     & wrk(posss1),wrk(possd1),1)
       end if
c
c2.a.22.nd.diade        W(ab,c)<- -s1(_j,a) . d2(_i,_k,b,c)
c     +s1(_j,b) . d2(_i,_k,a,c)
c
       if ((syma.eq.symj).and.(symik.eq.symbc)) then
c2.a.2.*      find address for s1,d1
       is1=mapih2(1,1,1)
       id1=mapim2(syma,1,1)
c
c2.a.2.*      def possition of s1,d1
       posss1=mapdh2(is1,1)
       possd1=mapdm2(id1,1)
c
c2.a.1.*      def additional dimensions
       nhelp1=dima*(dima-1)/2
c
c2.a.2.*      add singly
       call t3sglh211 (wrk(possw),dima,nhelp1,dimc,
     & wrk(posss1),wrk(possd1),-1)
       end if
c
c2.a.33.rd. part  W(ab,c)<- +s2(_k,c) . d1(_ij,ab)
c
       if ((symc.eq.symk).and.(symij.eq.symab)) then
c2.a.3.*      find address for s1,d1
       is1=mapih3(1,1,1)
       id1=mapim3(syma,1,1)
c
c2.a.3.*      def possition of s1,d1
       posss1=mapdh3(is1,1)
       possd1=mapdm3(id1,1)
c
c2.a.3.*      def additional dimensions
       nhelp1=dima*(dima-1)/2
c
c2.a.3.*      add singly
       call t3sglh212 (wrk(possw),dima,nhelp1,dimc,
     & wrk(posss1),wrk(possd1),1)
       end if
c
       else
c2.b  case syma>symb,symc
c
c2.b.11.st.diade  W(a>b,c)<- +s1(_i,a) . d2(_j,_k,b,c)
c     -s1(_i,b) . d2(_j,_k,a,c)
c
       if ((syma.eq.symi).and.(symjk.eq.symbc)) then
c2.b.1.*      find address for s1,d1
       is1=mapih1(1,1,1)
       id1=mapim1(symb,1,1)
c
c2.b.1.*      def possition of s1,d1
       posss1=mapdh1(is1,1)
       possd1=mapdm1(id1,1)
c
c2.b.1.*      add singly
       call t3sglh221 (wrk(possw),dima,dimb,dimc,
     & wrk(posss1),wrk(possd1),1)
       else if ((symb.eq.symi).and.(symjk.eq.symac)) then
c2.b.1.*      find address for s2,d2
       is2=mapih1(1,1,1)
       id2=mapim1(syma,1,1)
c
c2.b.1.*      def possition of s1,d1
       posss2=mapdh1(is2,1)
       possd2=mapdm1(id2,1)
c
c2.b.1.*      add singly
       call t3sglh222 (wrk(possw),dima,dimb,dimc,
     & wrk(posss2),wrk(possd2),1)
       end if
c
c2.a.22.nd.diade        W(ab,c)<- -s1(_j,a) . d2(_i,_k,b,c)
c     +s1(_j,b) . d2(_i,_k,a,c)
c
       if ((syma.eq.symj).and.(symik.eq.symbc)) then
c2.b.2.*      find address for s1,d1
       is1=mapih2(1,1,1)
       id1=mapim2(symb,1,1)
c
c2.b.2.*      def possition of s1,d1
       posss1=mapdh2(is1,1)
       possd1=mapdm2(id1,1)
c
c2.b.2.*      add singly
       call t3sglh221 (wrk(possw),dima,dimb,dimc,
     & wrk(posss1),wrk(possd1),-1)
       else if ((symb.eq.symj).and.(symik.eq.symac)) then
c2.b.2.*      find address for s2,d2
       is2=mapih2(1,1,1)
       id2=mapim2(syma,1,1)
c
c2.b.2.*      def possition of s1,d1
       posss2=mapdh2(is2,1)
       possd2=mapdm2(id2,1)
c
c2.b.2.*      add singly
       call t3sglh222 (wrk(possw),dima,dimb,dimc,
     & wrk(posss2),wrk(possd2),-1)
       end if
c
c2.a.33.rd. part  W(ab,c)<- +s2(_k,c) . d1(_ij,ab)
c
       if ((symc.eq.symk).and.(symij.eq.symab)) then
c2.b.3.*      find address for s1,d1
       is3=mapih3(1,1,1)
       id3=mapim3(syma,1,1)
c
c2.b.3.*      def possition of s1,d1
       posss3=mapdh3(is3,1)
       possd3=mapdm3(id3,1)
c
c2.b.3.*      add singly
       call t3sglh223 (wrk(possw),dima,dimb,dimc,
     & wrk(posss3),wrk(possd3),1)
       end if
c
       end if
c
 200    continue
c
c
       else if (typdiv.eq.3) then
c3    case B(p,qr)
c
c3.*  ext H1(a) <= S1(a,i) for given i
       call ext(wrk,wrksize,
     & 2,2,i,0,0,symi,0,0,mapds1,mapis1,1,
     & possh10,mapdh1,mapih1,ssh1,rc1)
c
c3.*  ext H2(a) <= S2(a,j) for given j
       call ext(wrk,wrksize,
     & 2,2,j,0,0,symj,0,0,mapds2,mapis2,1,
     & possh20,mapdh2,mapih2,ssh2,rc1)
c
c3.*  ext H3(a) <= S2(a,k) for given k
       call ext(wrk,wrksize,
     & 2,2,k,0,0,symk,0,0,mapds2,mapis2,1,
     & possh30,mapdh3,mapih3,ssh3,rc1)
c
c3.*  ext M1(bc) <= D2(bc,jk) for given jk
       call ext(wrk,wrksize,
     & 4,7,j,k,0,symj,symk,0,mapdd2,mapid2,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c3.*  ext M2(bc) <= D1(bc,ik) for given ik
       call ext(wrk,wrksize,
     & 4,7,i,k,0,symi,symk,0,mapdd1,mapid1,1,
     & possm20,mapdm2,mapim2,ssm2,rc1)
c
c3.*  ext M3(bc) <= D1(bc,ij) for given ij
       call ext(wrk,wrksize,
     & 4,7,i,j,0,symi,symj,0,mapdd1,mapid1,1,
     & possm30,mapdm3,mapim3,ssm3,rc1)
c
c
       do 300 iw=1,mapdw(0,5)
c
c3.*  def possition of W,V
       possw=mapdw(iw,1)
c
c3.*  def symmetry status
       syma=mapdw(iw,3)
       symb=mapdw(iw,4)
       symc=mapdw(iw,5)
c
c3.*  def dimensions
       dima=dimm(mapdw(0,1),syma)
       dimb=dimm(mapdw(0,2),symb)
       dimc=dimm(mapdw(0,3),symc)
c
c3.*  realize packing
c
       if (symb.eq.symc) then
c3.a  case syma,symb=symc
c
c3.a.11.st  part  W(a,bc)<- +s1(_i,a) . d2(_jk,b,c)
c
       if (symi.eq.syma) then
c3.a.1.*      find address for s1,d1
       is1=mapih1(1,1,1)
       id1=mapim1(symb,1,1)
c
c3.a.1.*      def possition of s1,d1
       posss1=mapdh1(is1,1)
       possd1=mapdm1(id1,1)
c
c3.a.1.*      def additional dimensions
       nhelp1=dimb*(dimb-1)/2
c
c3.a.1.*      add singly
       call t3sglh312 (wrk(possw),dima,dimb,nhelp1,
     & wrk(posss1),wrk(possd1),1)
       end if
c
c3.a.22.nd.diade        W(a,bc)<- +s1(_j,b) . d2(_i,_k,a,c)
c     -s1(_j,c) . d2(_i,_k,a,b)
c
       if (symj.eq.symb) then
c3.a.2.*      find address for s1,d1
       is1=mapih2(1,1,1)
       id1=mapim2(syma,1,1)
c
c3.a.2.*      def possition of s1,d1
       posss1=mapdh2(is1,1)
       possd1=mapdm2(id1,1)
c
c3.a.2.*      def additional dimensions
       nhelp1=dimb*(dimb-1)/2
c
c3.a.2.*      add singly
       call t3sglh311 (wrk(possw),dima,dimb,nhelp1,
     & wrk(posss1),wrk(possd1),1)
       end if
c
c
c3.a.33.rd.diade  W(a,bc)<- -s1(_k,b) . d2(_i,_j,a,c)
c     +s1(_k,c) . d2(_i,_j,a,b)
c
       if (symk.eq.symb) then
c3.a.3.*      find address for s1,d1
       is1=mapih3(1,1,1)
       id1=mapim3(syma,1,1)
c
c3.a.3.*      def possition of s1,d1
       posss1=mapdh3(is1,1)
       possd1=mapdm3(id1,1)
c
c3.a.3.*      def additional dimensions
       nhelp1=dimb*(dimb-1)/2
c
c3.a.3.*      add singly
       call t3sglh311 (wrk(possw),dima,dimb,nhelp1,
     & wrk(posss1),wrk(possd1),-1)
       end if
c
       else
c3.b  case syma,symb.ne.symc
c
c3.b.11.st  part  W(a,bc)<- +s1(_i,a) . d2(_jk,b,c)
c
       if (symi.eq.syma) then
c3.b.1.*      find address for s1,d1
       is1=mapih1(1,1,1)
       id1=mapim1(symb,1,1)
c
c3.b.1.*      def possition of s1,d1
       posss1=mapdh1(is1,1)
       possd1=mapdm1(id1,1)
c
c3.b.1.*      add singly
       call t3sglh323 (wrk(possw),dima,dimb,dimc,
     & wrk(posss1),wrk(possd1),1)
       end if
c
c3.b.22.nd.diade        W(a,bc)<- +s1(_j,b) . d2(_i,_k,a,c)
c     -s1(_j,c) . d2(_i,_k,a,b)
c
       if (symj.eq.symb) then
c3.b.2.*      find address for s1,d1
       is1=mapih2(1,1,1)
       id1=mapim2(syma,1,1)
c
c3.b.2.*      def possition of s1,d1
       posss1=mapdh2(is1,1)
       possd1=mapdm2(id1,1)
c
c3.b.2.*      add singly
       call t3sglh321 (wrk(possw),dima,dimb,dimc,
     & wrk(posss1),wrk(possd1),1)
       else if (symj.eq.symc) then
c3.b.2.*      find address for s1,d1
       is1=mapih2(1,1,1)
       id1=mapim2(syma,1,1)
c
c3.b.2.*      def possition of s1,d1
       posss1=mapdh2(is1,1)
       possd1=mapdm2(id1,1)
c
c3.b.2.*      add singly
       call t3sglh322 (wrk(possw),dima,dimb,dimc,
     & wrk(posss1),wrk(possd1),1)
       end if
c
c
c3.b.33.rd.diade  W(a,bc)<- -s1(_k,b) . d2(_i,_j,a,c)
c     +s1(_k,c) . d2(_i,_j,a,b)
c
       if (symk.eq.symb) then
c3.b.3.*      find address for s1,d1
       is1=mapih3(1,1,1)
       id1=mapim3(syma,1,1)
c
c3.b.3.*      def possition of s1,d1
       posss1=mapdh3(is1,1)
       possd1=mapdm3(id1,1)
c
c3.b.3.*      add singly
       call t3sglh321 (wrk(possw),dima,dimb,dimc,
     & wrk(posss1),wrk(possd1),-1)
       else if (symk.eq.symc) then
c3.b.3.*      find address for s1,d1
       is1=mapih3(1,1,1)
       id1=mapim3(syma,1,1)
c
c3.b.3.*      def possition of s1,d1
       posss1=mapdh3(is1,1)
       possd1=mapdm3(id1,1)
c
c3.b.3.*      add singly
       call t3sglh322 (wrk(possw),dima,dimb,dimc,
     & wrk(posss1),wrk(possd1),-1)
       end if
c
c
       end if
c
 300    continue
c
c
       else
c     RC=1 , typdiv is not 1,2,3 (NCI)
       rc=1
       return
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(ssw)
       end
c
c     -----------------------------
c
       subroutine t3sglh11 (w,dima,dimab,dimabc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma=symb=symc
c
c     W(abc)  <-  + S1 _i(a) . D1 _jk(bc)
c     - S1 _i(b) . D1 _jk(ac)
c     + S1 _i(c) . D1 _jk(ab)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a (b,c) index (I)
c     dimab  - dimension of ab (ac,bc) index (I)
c     dimabc - dimension of abc index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
#include "t31.fh"
       integer dima,dimab,dimabc,ns
       real*8 w(1:dimabc)
       real*8 s1(1:dima)
       real*8 d1(1:dimab)
c
c     help variables
c
       integer a,b,c,ab0,ac0,bc0,abc
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       abc=0
       do 100 a=3,dima
       s=s1(a)
       do 101 b=2,a-1
       bc0=nshf(b)
       do 102 c=1,b-1
       abc=abc+1
       w(abc)=w(abc)+d1(bc0+c)*s
 102    continue
 101    continue
 100    continue
c
       abc=0
       do 110 a=3,dima
       ac0=nshf(a)
       do 111 b=2,a-1
       s=s1(b)
       do 112 c=1,b-1
       abc=abc+1
       w(abc)=w(abc)-d1(ac0+c)*s
 112    continue
 111    continue
 110    continue
c
       abc=0
       do 120 a=3,dima
       ab0=nshf(a)
       do 121 b=2,a-1
       s=d1(ab0+b)
       do 122 c=1,b-1
       abc=abc+1
       w(abc)=w(abc)+s1(c)*s
 122    continue
 121    continue
 120    continue
c
       else
c     phase -1
c
       abc=0
       do 200 a=3,dima
       s=s1(a)
       do 201 b=2,a-1
       bc0=nshf(b)
       do 202 c=1,b-1
       abc=abc+1
       w(abc)=w(abc)-d1(bc0+c)*s
 202    continue
 201    continue
 200    continue
c
       abc=0
       do 210 a=3,dima
       ac0=nshf(a)
       do 211 b=2,a-1
       s=s1(b)
       do 212 c=1,b-1
       abc=abc+1
       w(abc)=w(abc)+d1(ac0+c)*s
 212    continue
 211    continue
 210    continue
c
       abc=0
       do 220 a=3,dima
       ab0=nshf(a)
       do 221 b=2,a-1
       s=d1(ab0+b)
       do 222 c=1,b-1
       abc=abc+1
       w(abc)=w(abc)-s1(c)*s
 222    continue
 221    continue
 220    continue
c
       end if
c
       return
       end
c
c     ---------------
c
       subroutine t3sglh121 (w,dima,dimab,dimc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma=symb > symc
c
c     W(ab,c) <-  + S1 _i(a) . D1 _jk(b,c)
c     - S1 _i(b) . D1 _jk(a,c)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a (b) index (I)
c     dimab   - dimension of ab index (I)
c     dimc   - dimension of c index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimab,dimc,ns
       real*8 w(1:dimab,1:dimc)
       real*8 s1(1:dima)
       real*8 d1(1:dima,1:dimc)
c
c     help variables
c
       integer a,b,c,ab
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 100 c=1,dimc
       ab=0
       do 101 a=2,dima
       s=s1(a)
       do 102 b=1,a-1
       ab=ab+1
       w(ab,c)=w(ab,c)+d1(b,c)*s
 102    continue
 101    continue
 100    continue
c
       do 110 c=1,dimc
       ab=0
       do 111 a=2,dima
       s=d1(a,c)
       do 112 b=1,a-1
       ab=ab+1
       w(ab,c)=w(ab,c)-s1(b)*s
 112    continue
 111    continue
 110    continue
c
       else
c     phase -1
c
       do 200 c=1,dimc
       ab=0
       do 201 a=2,dima
       s=s1(a)
       do 202 b=1,a-1
       ab=ab+1
       w(ab,c)=w(ab,c)-d1(b,c)*s
 202    continue
 201    continue
 200    continue
c
       do 210 c=1,dimc
       ab=0
       do 211 a=2,dima
       s=d1(a,c)
       do 212 b=1,a-1
       ab=ab+1
       w(ab,c)=w(ab,c)+s1(b)*s
 212    continue
 211    continue
 210    continue
c
       end if
c
       return
       end
c
c     ---------------
c
       subroutine t3sglh122 (w,dima,dimab,dimc,s3,d3,ns)
c
c     this routine add following contribution to W
c     for syma=symb > symc
c
c     W(ab,c) <-  + S3 _i(c) . D3 _jk(ab)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a (b) index (I)
c     dimab   - dimension of ab index (I)
c     dimc   - dimension of c index (I)
c     s3     - S3 matrix (I)
c     d3     - D3 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimab,dimc,ns
       real*8 w(1:dimab,1:dimc)
       real*8 s3(1:dimc)
       real*8 d3(1:dimab)
c
c     help variables
c
       integer c,ab
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 120 c=1,dimc
       s=s3(c)
       do 121 ab=1,dimab
       w(ab,c)=w(ab,c)+d3(ab)*s
 121    continue
 120    continue
c
       else
c     phase -1
c
       do 220 c=1,dimc
       s=s3(c)
       do 221 ab=1,dimab
       w(ab,c)=w(ab,c)-d3(ab)*s
 221    continue
 220    continue
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dima)
       end
c
c     ---------------
c
       subroutine t3sglh131 (w,dima,dimb,dimbc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma=symb > symc
c
c     W(a,bc) <-  + S1 _i(a) . D1 _jk(bc)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a (b) index (I)
c     dimab   - dimension of ab index (I)
c     dimc   - dimension of c index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimbc,ns
       real*8 w(1:dima,1:dimbc)
       real*8 s1(1:dima)
       real*8 d1(1:dimbc)
c
c     help variables
c
       integer a,bc
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 100 bc=1,dimbc
       s=d1(bc)
       do 101 a=1,dima
       w(a,bc)=w(a,bc)+s1(a)*s
 101    continue
 100    continue
c
       else
c     phase -1
c
       do 200 bc=1,dimbc
       s=d1(bc)
       do 201 a=1,dima
       w(a,bc)=w(a,bc)-s1(a)*s
 201    continue
 200    continue
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dimb)
       end
c
c     ---------------
c
       subroutine t3sglh132 (w,dima,dimb,dimbc,s2,d2,ns)
c
c     this routine add following contribution to W
c     for syma=symb > symc
c
c     W(a,bc) <-  - S2 _i(b) . D2 _jk(a,c)
c     + S2 _i(c) . D2 _jk(a,b)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a (b) index (I)
c     dimab   - dimension of ab index (I)
c     dimc   - dimension of c index (I)
c     s2     - S2 matrix (I)
c     d2     - D2 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimbc,ns
       real*8 w(1:dima,1:dimbc)
       real*8 s2(1:dimb)
       real*8 d2(1:dima,1:dimb)
c
c     help variables
c
       integer a,b,c,bc
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       bc=0
       do 110 b=2,dimb
       s=s2(b)
       do 111 c=1,b-1
       bc=bc+1
       do 112 a=1,dima
       w(a,bc)=w(a,bc)-d2(a,c)*s
 112    continue
 111    continue
 110    continue
c
       bc=0
       do 120 b=2,dimb
       do 121 c=1,b-1
       s=s2(c)
       bc=bc+1
       do 122 a=1,dima
       w(a,bc)=w(a,bc)+d2(a,b)*s
 122    continue
 121    continue
 120    continue
c
       else
c     phase -1
c
       bc=0
       do 210 b=2,dimb
       s=s2(b)
       do 211 c=1,b-1
       bc=bc+1
       do 212 a=1,dima
       w(a,bc)=w(a,bc)+d2(a,c)*s
 212    continue
 211    continue
 210    continue
c
       bc=0
       do 220 b=2,dimb
       do 221 c=1,b-1
       s=s2(c)
       bc=bc+1
       do 222 a=1,dima
       w(a,bc)=w(a,bc)-d2(a,b)*s
 222    continue
 221    continue
 220    continue
c
       end if
c
       return
       end
c
c     ---------------
c
       subroutine t3sglh141 (w,dima,dimb,dimc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma>symb>symc
c
c     W(a,b,c) <- + S1 _i(a) . D1 _jk(b,c)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b index (I)
c     dimc   - dimension of c index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimc,ns
       real*8 w(1:dima,1:dimb,1:dimc)
       real*8 s1(1:dima)
       real*8 d1(1:dimb,1:dimc)
c
c     help variables
c
       integer a,b,c
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 100 c=1,dimc
       do 101 b=1,dimb
       s=d1(b,c)
       do 102 a=1,dima
       w(a,b,c)=w(a,b,c)+s1(a)*s
 102    continue
 101    continue
 100    continue
c
       else
c     phase = -1
c
       do 200 c=1,dimc
       do 201 b=1,dimb
       s=d1(b,c)
       do 202 a=1,dima
       w(a,b,c)=w(a,b,c)-s1(a)*s
 202    continue
 201    continue
 200    continue
c
       end if
c
       return
       end
c
c     ---------------
c
       subroutine t3sglh142 (w,dima,dimb,dimc,s2,d2,ns)
c
c     this routine add following contribution to W
c     for syma>symb>symc
c
c     W(a,b,c) <- - S2 _i(b) . D2 _jk(a,c)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b index (I)
c     dimc   - dimension of c index (I)
c     s2     - S2 matrix (I)
c     d2     - D2 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimc,ns
       real*8 w(1:dima,1:dimb,1:dimc)
       real*8 s2(1:dimb)
       real*8 d2(1:dima,1:dimc)
c
c     help variables
c
       integer a,b,c
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 110 c=1,dimc
       do 111 b=1,dimb
       s=s2(b)
       do 112 a=1,dima
       w(a,b,c)=w(a,b,c)-d2(a,c)*s
 112    continue
 111    continue
 110    continue
c
       else
c     phase = -1
c
       do 210 c=1,dimc
       do 211 b=1,dimb
       s=s2(b)
       do 212 a=1,dima
       w(a,b,c)=w(a,b,c)+d2(a,c)*s
 212    continue
 211    continue
 210    continue
c
       end if
c
       return
       end
c
c     ---------------
c
       subroutine t3sglh143 (w,dima,dimb,dimc,s3,d3,ns)
c
c     this routine add following contribution to W
c     for syma>symb>symc
c
c     W(a,b,c) <- + S3 _i(c) . D3 _jk(a,b)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b index (I)
c     dimc   - dimension of c index (I)
c     s3     - S3 matrix (I)
c     d3     - D3 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimc,ns
       real*8 w(1:dima,1:dimb,1:dimc)
       real*8 s3(1:dimc)
       real*8 d3(1:dima,1:dimb)
c
c     help variables
c
       integer a,b,c
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 120 c=1,dimc
       s=s3(c)
       do 121 b=1,dimb
       do 122 a=1,dima
       w(a,b,c)=w(a,b,c)+d3(a,b)*s
 122    continue
 121    continue
 120    continue
c
       else
c     phase = -1
c
       do 220 c=1,dimc
       s=s3(c)
       do 221 b=1,dimb
       do 222 a=1,dima
       w(a,b,c)=w(a,b,c)-d3(a,b)*s
 222    continue
 221    continue
 220    continue
c
       end if
c
       return
       end
c
c     -----------------------------
c
       subroutine t3sglh211 (w,dima,dimab,dimc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma=symb;symc
c
c     W(ab,c) <-  + S1 _i(a) . D1 _jk(b,c)
c     - S1 _i(b) . D1 _jk(a,c)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a (b) index (I)
c     dimab  - dimension of ab  index (I)
c     dimc   - dimension of c index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimab,dimc,ns
       real*8 w(1:dimab,1:dimc)
       real*8 s1(1:dima)
       real*8 d1(1:dima,1:dimc)
c
c     help variables
c
       integer a,b,c,ab
       real*8 s
c
       if (ns.eq.1) then
c     phase + 1
c
       do 100 c=1,dimc
       ab=0
       do 101 a=2,dima
       s=s1(a)
       do 102 b=1,a-1
       ab=ab+1
       w(ab,c)=w(ab,c)+d1(b,c)*s
 102    continue
 101    continue
 100    continue
c
       do 110 c=1,dimc
       ab=0
       do 111 a=2,dima
       s=d1(a,c)
       do 112 b=1,a-1
       ab=ab+1
       w(ab,c)=w(ab,c)-s1(b)*s
 112    continue
 111    continue
 110    continue
c
       else
c     phase - 1
c
       do 200 c=1,dimc
       ab=0
       do 201 a=2,dima
       s=s1(a)
       do 202 b=1,a-1
       ab=ab+1
       w(ab,c)=w(ab,c)-d1(b,c)*s
 202    continue
 201    continue
 200    continue
c
       do 210 c=1,dimc
       ab=0
       do 211 a=2,dima
       s=d1(a,c)
       do 212 b=1,a-1
       ab=ab+1
       w(ab,c)=w(ab,c)+s1(b)*s
 212    continue
 211    continue
 210    continue
c
       end if
c
       return
       end
c
c     ---------------
c
       subroutine t3sglh212 (w,dima,dimab,dimc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma=symb;symc
c
c     W(ab,c) <-  + S1 _i(c) . D1 _jk(ab)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a (b) index (I)
c     dimab  - dimension of ab (ac,bc) index (I)
c     dimc   - dimension of c index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimab,dimc,ns
       real*8 w(1:dimab,1:dimc)
       real*8 s1(1:dimc)
       real*8 d1(1:dimab)
c
c     help variables
c
       integer c,ab
       real*8 s
c
       if (ns.eq.1) then
c     phase + 1
c
       do 100 c=1,dimc
       s=s1(c)
       do 101 ab=1,dimab
       w(ab,c)=w(ab,c)+d1(ab)*s
 101    continue
 100    continue
c
       else
c     phase - 1
c
       do 200 c=1,dimc
       s=s1(c)
       do 201 ab=1,dimab
       w(ab,c)=w(ab,c)-d1(ab)*s
 201    continue
 200    continue
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dima)
       end
c
c     ---------------
c
       subroutine t3sglh221 (w,dima,dimb,dimc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma>symb;symc
c
c     W(a,b;c) <- + S1 _i(a) . D1 _jk(b,c)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b index (I)
c     dimc   - dimension of c index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimc,ns
       real*8 w(1:dima,1:dimb,1:dimc)
       real*8 s1(1:dima)
       real*8 d1(1:dimb,1:dimc)
c
c     help variables
c
       integer a,b,c
       real*8 s
c
       if (ns.eq.1) then
c     phase + 1
c
       do 100 c=1,dimc
       do 101 b=1,dimb
       s=d1(b,c)
       do 102 a=1,dima
       w(a,b,c)=w(a,b,c)+s1(a)*s
 102    continue
 101    continue
 100    continue
c
       else
c     phase - 1
c
       do 200 c=1,dimc
       do 201 b=1,dimb
       s=d1(b,c)
       do 202 a=1,dima
       w(a,b,c)=w(a,b,c)-s1(a)*s
 202    continue
 201    continue
 200    continue
c
       end if
c
       return
       end
c
c     ---------------
c
       subroutine t3sglh222 (w,dima,dimb,dimc,s2,d2,ns)
c
c     this routine add following contribution to W
c     for syma>symb;symc
c
c     W(a,b;c) <- - S2 _i(b) . D2 _jk(a,c)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b index (I)
c     dimc   - dimension of c index (I)
c     s2     - S2 matrix (I)
c     d2     - D2 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimc,ns
       real*8 w(1:dima,1:dimb,1:dimc)
       real*8 s2(1:dimb)
       real*8 d2(1:dima,1:dimc)
c
c     help variables
c
       integer a,b,c
       real*8 s
c
       if (ns.eq.1) then
c     phase + 1
c
       do 110 c=1,dimc
       do 111 b=1,dimb
       s=s2(b)
       do 112 a=1,dima
       w(a,b,c)=w(a,b,c)-d2(a,c)*s
 112    continue
 111    continue
 110    continue
c
       else
c     phase - 1
c
       do 210 c=1,dimc
       do 211 b=1,dimb
       s=s2(b)
       do 212 a=1,dima
       w(a,b,c)=w(a,b,c)+d2(a,c)*s
 212    continue
 211    continue
 210    continue
c
       end if
c
       return
       end
c
c     ---------------
c
       subroutine t3sglh223 (w,dima,dimb,dimc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma>symb;symc
c
c     W(a,b;c) <- + S1 _i(c) . D1 _jk(a,b)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b index (I)
c     dimc   - dimension of c index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimc,ns
       real*8 w(1:dima,1:dimb,1:dimc)
       real*8 s1(1:dimc)
       real*8 d1(1:dima,1:dimb)
c
c     help variables
c
       integer a,b,c
       real*8 s
c
       if (ns.eq.1) then
c     phase + 1
c
       do 100 c=1,dimc
       s=s1(c)
       do 101 b=1,dimb
       do 102 a=1,dima
       w(a,b,c)=w(a,b,c)+d1(a,b)*s
 102    continue
 101    continue
 100    continue
c
       else
c     phase - 1
c
       do 200 c=1,dimc
       s=s1(c)
       do 201 b=1,dimb
       do 202 a=1,dima
       w(a,b,c)=w(a,b,c)-d1(a,b)*s
 202    continue
 201    continue
 200    continue
c
       end if
c
       return
       end
c
c     ---------------
c
       subroutine t3sglh311 (w,dima,dimb,dimbc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma;symb=symc
c
c     W(a;bc)  <- + S1 _i(b) . D1 _jk(a,c)
c     - S1 _i(c) . D2 _jk(a,b)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b (c) index (I)
c     dimbc  - dimension of bc index (I)
c     s1     - S1 matrix (I)
c     s2     - S2 matrix (I)
c     d1     - D1 matrix (I)
c     d2     - D2 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimbc,ns
       real*8 w(1:dima,1:dimbc)
       real*8 s1(1:dimb)
       real*8 d1(1:dima,1:dimb)
c
c     help variables
c
       integer a,b,c,bc
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       bc=0
       do 100 b=2,dimb
       s=s1(b)
       do 101 c=1,b-1
       bc=bc+1
       do 102 a=1,dima
       w(a,bc)=w(a,bc)+d1(a,c)*s
 102    continue
 101    continue
 100    continue
c
       bc=0
       do 110 b=2,dimb
       do 111 c=1,b-1
       bc=bc+1
       s=s1(c)
       do 112 a=1,dima
       w(a,bc)=w(a,bc)-d1(a,b)*s
 112    continue
 111    continue
 110    continue
c
       else
c     phase - 1
c
       bc=0
       do 200 b=2,dimb
       s=s1(b)
       do 201 c=1,b-1
       bc=bc+1
       do 202 a=1,dima
       w(a,bc)=w(a,bc)-d1(a,c)*s
 202    continue
 201    continue
 200    continue
c
       bc=0
       do 210 b=2,dimb
       do 211 c=1,b-1
       bc=bc+1
       s=s1(c)
       do 212 a=1,dima
       w(a,bc)=w(a,bc)+d1(a,b)*s
 212    continue
 211    continue
 210    continue
c
       end if
c
       return
       end
c
c     ---------------
c
       subroutine t3sglh312 (w,dima,dimb,dimbc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma;symb=symc
c
c     W(a;bc)  <- + S1 _i(a) . D1 _jk(bc)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b (c) index (I)
c     dimbc  - dimension of bc index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimbc,ns
       real*8 w(1:dima,1:dimbc)
       real*8 s1(1:dima)
       real*8 d1(1:dimbc)
c
c     help variables
c
       integer a,bc
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 100 bc=1,dimbc
       s=d1(bc)
       do 101 a=1,dima
       w(a,bc)=w(a,bc)+s1(a)*s
 101    continue
 100    continue
c
       else
c     phase - 1
c
       do 200 bc=1,dimbc
       s=d1(bc)
       do 201 a=1,dima
       w(a,bc)=w(a,bc)-s1(a)*s
 201    continue
 200    continue
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dimb)
       end
c
c     ---------------
c
       subroutine t3sglh321 (w,dima,dimb,dimc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma;symb>symc
c
c     W(a;b,c) <- + S1 _i(b) . D1 _jk(a,c)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b index (I)
c     dimc   - dimension of c index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimc,ns
       real*8 w(1:dima,1:dimb,1:dimc)
       real*8 s1(1:dimb)
       real*8 d1(1:dima,1:dimc)
c
c     help variables
c
       integer a,b,c
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 100 c=1,dimc
       do 101 b=1,dimb
       s=s1(b)
       do 102 a=1,dima
       w(a,b,c)=w(a,b,c)+d1(a,c)*s
 102    continue
 101    continue
 100    continue
c
       else
c     phase -1
c
       do 200 c=1,dimc
       do 201 b=1,dimb
       s=s1(b)
       do 202 a=1,dima
       w(a,b,c)=w(a,b,c)-d1(a,c)*s
 202    continue
 201    continue
 200    continue
c
       end if
c
       return
       end
c
c     ---------------
c
       subroutine t3sglh322 (w,dima,dimb,dimc,s2,d2,ns)
c
c     this routine add following contribution to W
c     for syma;symb>symc
c
c     W(a;b,c) <- - S2 _i(c) . D2 _jk(a,b)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b index (I)
c     dimc   - dimension of c index (I)
c     s2     - S2 matrix (I)
c     d2     - D2 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimc,ns
       real*8 w(1:dima,1:dimb,1:dimc)
       real*8 s2(1:dimc)
       real*8 d2(1:dima,1:dimb)
c
c     help variables
c
       integer a,b,c
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 110 c=1,dimc
       s=s2(c)
       do 111 b=1,dimb
       do 112 a=1,dima
       w(a,b,c)=w(a,b,c)-d2(a,b)*s
 112    continue
 111    continue
 110    continue
c
       else
c     phase -1
c
       do 210 c=1,dimc
       s=s2(c)
       do 211 b=1,dimb
       do 212 a=1,dima
       w(a,b,c)=w(a,b,c)+d2(a,b)*s
 212    continue
 211    continue
 210    continue
c
       end if
c
       return
       end
c
c     ---------------
c
       subroutine t3sglh323 (w,dima,dimb,dimc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma;symb>symc
c
c     W(a;b,c) <- + S1 _i(a) . D1 _jk(b,c)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b index (I)
c     dimc   - dimension of c index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimc,ns
       real*8 w(1:dima,1:dimb,1:dimc)
       real*8 s1(1:dima)
       real*8 d1(1:dimb,1:dimc)
c
c     help variables
c
       integer a,b,c
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 100 c=1,dimc
       do 101 b=1,dimb
       s=d1(b,c)
       do 102 a=1,dima
       w(a,b,c)=w(a,b,c)+s1(a)*s
 102    continue
 101    continue
 100    continue
c
       else
c     phase -1
c
       do 200 c=1,dimc
       do 201 b=1,dimb
       s=d1(b,c)
       do 202 a=1,dima
       w(a,b,c)=w(a,b,c)-s1(a)*s
 202    continue
 201    continue
 200    continue
c
       end if
c
       return
       end
