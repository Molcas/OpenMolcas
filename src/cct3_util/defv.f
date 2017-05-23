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
c     defv
c     defvhlp1
c     defvhlp21
c     defvhlp22
c     defvhlp3
c     defvhlp4
c     defvhlp51
c     defvhlp52
c     defvhlp53
c     defvhlp54
c     defvhlp61
c     defvhlp62
c     defvhlp7
c     defvhlp81
c     defvhlp82
c     defvhlp9
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine defv (wrk,wrksize,
     & deftyp,possv0,mapdv,mapiv,ssv,
     & mapdr,mapir,ssr,rc)
c
c     this routine define v mediate as:
c     V_i(abc) = <ab||ic>  from integrals, stored in R_i(abc)
c     R _i(a,bc)bbb is a matrix, where integrals <ab|ic> are stored
c     for b>=c
c
c     deftyp - type of definition (see table) (I)
c     possv0 - initial address of V (I)
c     mapdv  - direct map of V (O)
c     mapiv  - inverse map of V (O)
c     ssv    - overall spin of V (O)
c     mapdr  - direct map of R (I)
c     mapir  - inverse map of R (I)
c     ssr    - overall spin of R (I)
c     rc     - return (error) code (O)
c
c     Table of definitions
c
c     deftyp         Operation                 Implement.
c     1      V(ab,c)aaa  = R(abc)-R(bac)         Yes
c     2      V(ab,c)bbb  = R(abc)-R(bac)         Yes
c     3      V(a,b,c)abb = R(abc)                Yes
c     4      V(a,b,c)aba = -R(bac)                Yes
c
#include "t31.fh"
#include "wrk.fh"
c
       integer deftyp,possv0,ssv,ssr,rc
       integer mapdv(0:512,1:6),mapdr(0:512,1:6)
       integer mapiv(1:8,1:8,1:8),mapir(1:8,1:8,1:8)
c
c     help variables
c
       integer posst,possv,possr1,possr2
       integer iv,ir1,ir2,syma,symb,symc
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5
       integer nhelp6,nhelp7,nhelp8,nhelp9,nhelp10
c
c
c0.*  def mapdv,mapiv
       if (deftyp.eq.1) then
c     case V(ab,c)aaa
       call cct3_grc0 (3,1,3,3,3,0,ssr,possv0,posst,mapdv,mapiv)
       else if (deftyp.eq.2) then
c     case V(ab,c)bbb
       call cct3_grc0 (3,1,4,4,4,0,ssr,possv0,posst,mapdv,mapiv)
       else if (deftyp.eq.3) then
c     case V(a,b,c)abb
       call cct3_grc0 (3,0,3,4,4,0,ssr,possv0,posst,mapdv,mapiv)
       else if (deftyp.eq.4) then
c     case V(a,b,c)aba
       call cct3_grc0 (3,0,3,4,3,0,ssr,possv0,posst,mapdv,mapiv)
       else
c     RC=1 , deftyp out of range (1-4) (Stup)
       rc=1
       return
       end if
c
c0.*  some tests
c

c0.*  define spin
       ssv=ssr
c
       if ((deftyp.eq.1).or.(deftyp.eq.2)) then
c
c12   case V(ab,c)aaa,bbb
c
       do 12 iv=1,mapdv(0,5)
c
c12.* def symmetries
       syma=mapdv(iv,3)
       symb=mapdv(iv,4)
       symc=mapdv(iv,5)
c
       if (syma.eq.symb) then
c12.1 syma=symb
c
       if (symb.eq.symc) then
c12.1.1            case syma=symb=symc
c
c12.1.1.*     def possitions of V, R1
       possv=mapdv(iv,1)
       ir1=mapir(syma,symb,1)
       possr1=mapdr(ir1,1)
c
c12.1.1.*     def dimensions
       nhelp1=nvb(syma)
       nhelp2=nhelp1*(nhelp1+1)/2
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp3=nhelp4*(nhelp4-1)/2
       nhelp5=nvb(syma)-nhelp4
c
c12.1.1.*     realize definition of V
       call defvhlp1 (wrk(possr1),wrk(possv),
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5)
c
       else
c12.1.2            case syma=symb.ne.symc
c
c12.1.2.*     def possitions of V, R1
       possv=mapdv(iv,1)
       if (syma.ge.symc) then
       ir1=mapir(syma,symb,1)
       else
       ir1=mapir(syma,symc,1)
       end if
       possr1=mapdr(ir1,1)
c
c12.1.2.*     def dimensions
       nhelp1=nvb(syma)
       nhelp2=nvb(symc)
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp3=nhelp4*(nhelp4-1)/2
       nhelp5=dimm(mapdv(0,3),symc)
       nhelp6=nvb(syma)-nhelp4
       nhelp7=nvb(symc)-nhelp5
c
c12.1.2.*     realize definition of V
       if (syma.ge.symc) then
       call defvhlp21 (wrk(possr1),wrk(possv),nhelp1,nhelp2,
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
       else
       call defvhlp22 (wrk(possr1),wrk(possv),nhelp1,nhelp2,
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
       end if
c
       end if
c
       else
c12.2 syma > symb
c
       if (syma.eq.symc) then
c12.2.1     case syma>symb , syma=symc
c
c12.2.1.*     def possitions of V, R1, R2
       possv=mapdv(iv,1)
c     R1 is permuted, since (a=c)>b
       ir1=mapir(syma,symc,1)
       possr1=mapdr(ir1,1)
       ir2=mapir(symb,syma,1)
       possr2=mapdr(ir2,1)
c
c12.2.1.*     def dimensions
       nhelp1=nvb(syma)
       nhelp2=nvb(symb)
       nhelp3=nvb(symc)
       nhelp4=nhelp1*(nhelp1+1)/2
       nhelp5=dimm(mapdv(0,1),syma)
       nhelp6=dimm(mapdv(0,2),symb)
       nhelp7=dimm(mapdv(0,3),symc)
       nhelp8=nvb(syma)-nhelp5
       nhelp9=nvb(symb)-nhelp6
       nhelp10=nvb(symc)-nhelp7
c
c12.2.1.*     realize definition of V
       call defvhlp3 (wrk(possr1),wrk(possr2),wrk(possv),
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,
     & nhelp6,nhelp7,nhelp8,nhelp9,nhelp10)
c
       else if (symb.eq.symc) then
c12.2.2     case syma>symb , symb=symc
c
c12.2.2.*     def possitions of V, R1, R2
       possv=mapdv(iv,1)
       ir1=mapir(syma,symb,1)
       possr1=mapdr(ir1,1)
       ir2=mapir(symb,syma,1)
       possr2=mapdr(ir2,1)
c
c12.2.2.*     def dimensions
       nhelp1=nvb(syma)
       nhelp3=nvb(symb)
       nhelp4=nvb(symc)
       nhelp2=nhelp3*(nhelp3+1)/2
       nhelp5=dimm(mapdv(0,1),syma)
       nhelp6=dimm(mapdv(0,2),symb)
       nhelp7=dimm(mapdv(0,3),symc)
       nhelp8=nvb(syma)-nhelp5
       nhelp9=nvb(symb)-nhelp6
       nhelp10=nvb(symc)-nhelp7
c
c12.2.2.*     realize definition of V
       call defvhlp4 (wrk(possr1),wrk(possr2),wrk(possv),
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,
     & nhelp6,nhelp7,nhelp8,nhelp9,nhelp10)
c
       else
c12.2.3     case syma>symb , symc.ne.syma,symb
c
c
c12.2.3.*     def possitions of V, R1, R2
       possv=mapdv(iv,1)
c
       if (symb.ge.symc) then
       ir1=mapir(syma,symb,1)
       else
       ir1=mapir(syma,symc,1)
       end if
       possr1=mapdr(ir1,1)
c
       if (syma.ge.symc) then
       ir2=mapir(symb,syma,1)
       else
       ir2=mapir(symb,symc,1)
       end if
       possr2=mapdr(ir2,1)
c
c12.2.3.*     def dimensions
       nhelp1=nvb(syma)
       nhelp2=nvb(symb)
       nhelp3=nvb(symc)
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp5=dimm(mapdv(0,2),symb)
       nhelp6=dimm(mapdv(0,3),symc)
       nhelp7=nvb(syma)-nhelp4
       nhelp8=nvb(symb)-nhelp5
       nhelp9=nvb(symc)-nhelp6
c
c12.2.3.*     realize definition of V
       if ((syma.gt.symc).and.(symb.gt.symc)) then
       call defvhlp51 (wrk(possr1),wrk(possr2),wrk(possv),
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,
     & nhelp6,nhelp7,nhelp8,nhelp9)
       else if ((syma.gt.symc).and.(symb.lt.symc)) then
       call defvhlp52 (wrk(possr1),wrk(possr2),wrk(possv),
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,
     & nhelp6,nhelp7,nhelp8,nhelp9)
       else if ((syma.lt.symc).and.(symb.gt.symc)) then
       call defvhlp53 (wrk(possr1),wrk(possr2),wrk(possv),
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,
     & nhelp6,nhelp7,nhelp8,nhelp9)
       else if ((syma.lt.symc).and.(symb.lt.symc)) then
       call defvhlp54 (wrk(possr1),wrk(possr2),wrk(possv),
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,
     & nhelp6,nhelp7,nhelp8,nhelp9)
       end if
c
       end if
c
       end if
c
 12     continue
c
c
       else if (deftyp.eq.3) then
c
c3    case V(a,b,c)abb
c
       do 3 iv=1,mapdv(0,5)
c
c3.*  def symmetries
       syma=mapdv(iv,3)
       symb=mapdv(iv,4)
       symc=mapdv(iv,5)
c
       if (symb.eq.symc) then
c3.1  symb=symc
c
c3.1.*def possitions of V, R1
       possv=mapdv(iv,1)
       ir1=mapir(syma,symb,1)
       possr1=mapdr(ir1,1)
c
c3.1.*def dimensions
       nhelp1=nvb(syma)
       nhelp2=nvb(symb)
       nhelp3=nhelp2*(nhelp2+1)/2
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp5=dimm(mapdv(0,2),symb)
       nhelp6=dimm(mapdv(0,3),symc)
       nhelp7=nvb(syma)-nhelp4

c3.1.*realize definition of V
       call defvhlp7 (wrk(possr1),wrk(possv),nhelp1,nhelp2,
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
c
       else
c     case symb.ne.symc
c
c3.1.*def possitions of V, R1
       possv=mapdv(iv,1)
       if (symb.ge.symc) then
       ir1=mapir(syma,symb,1)
       else
       ir1=mapir(syma,symc,1)
       end if
       possr1=mapdr(ir1,1)
c
c3.1.*def dimensions
       nhelp1=nvb(syma)
       nhelp2=nvb(symb)
       nhelp3=nvb(symc)
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp5=dimm(mapdv(0,2),symb)
       nhelp6=dimm(mapdv(0,3),symc)
       nhelp7=nvb(syma)-nhelp4
c
c3.1.*realize definition of V
       if (symb.gt.symc) then
       call defvhlp61 (wrk(possr1),wrk(possv),nhelp1,nhelp2,
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
       else
       call defvhlp62 (wrk(possr1),wrk(possv),nhelp1,nhelp2,
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
       end if
c
       end if
c
 3      continue
c
       else if (deftyp.eq.4) then
c
c4    case V(a,b,c)aba
c
       do 4 iv=1,mapdv(0,5)
c
c4.*  def symmetries
       syma=mapdv(iv,3)
       symb=mapdv(iv,4)
       symc=mapdv(iv,5)
c
       if (syma.eq.symc) then
c4.1  syma=symc
c
c4.1.*def possitions of V, R2
       possv=mapdv(iv,1)
       ir2=mapir(symb,syma,1)
       possr2=mapdr(ir2,1)
c
c4.1.*def dimensions
       nhelp1=nvb(symb)
       nhelp2=nvb(syma)
       nhelp3=nhelp2*(nhelp2+1)/2
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp5=dimm(mapdv(0,2),symb)
       nhelp6=dimm(mapdv(0,3),symc)
       nhelp7=nvb(syma)-nhelp4
       nhelp8=nvb(symc)-nhelp6
c
c4.1.*realize definition of V
       call defvhlp9 (wrk(possr2),wrk(possv),nhelp1,nhelp2,
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8)
c
       else
c     case symb.ne.symc
c
c4.1.*def possitions of V, R2
       possv=mapdv(iv,1)
       if (syma.ge.symc) then
       ir2=mapir(symb,syma,1)
       else
       ir2=mapir(symb,symc,1)
       end if
       possr2=mapdr(ir2,1)
c
c4.1.*def dimensions
       nhelp1=nvb(symb)
       nhelp2=nvb(syma)
       nhelp3=nvb(symc)
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp5=dimm(mapdv(0,2),symb)
       nhelp6=dimm(mapdv(0,3),symc)
       nhelp7=nvb(syma)-nhelp4
       nhelp8=nvb(symc)-nhelp6
c
c4.1.*realize definition of V
       if (syma.gt.symc) then
       call defvhlp81 (wrk(possr2),wrk(possv),nhelp1,nhelp2,
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8)
       else
       call defvhlp82 (wrk(possr2),wrk(possv),nhelp1,nhelp2,
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8)
       end if
c
       end if
c
 4      continue
c
       else
c     not implemented deftyp
c
       end if
c
       return
       end
c
c     ----------------------------------------
c
       subroutine defvhlp1 (r1,v,dimr1a,dimr1bc,dimvab,dimvc,add)
c
c     this routine do
c     V(ab,c)xxx = R1(a,bc)-R1(b,ac) x=a,b
c     for syma=symb=symc
c
c     r1        - r1 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a (b,c) in R1 (I)
c     dimr1bc        - dimension of bc (ac) in R1 (I)
c     dimvab        - dimension of ab in V (I)
c     dimvc        - dimension of c (a,b) in V (I)
c     add     - additional constat (I)
c     (#of singly occ in syma for alpha, 0 for beta)
c
#include "t31.fh"
       integer dimr1a,dimr1bc,dimvab,dimvc,add
       real*8 r1(1:dimr1a,1:dimr1bc)
       real*8 v(1:dimvab,1:dimvc)
c
c     help variables
c
       integer a,b,c,ab0,ab,ar1,cr1,acr1,bk
c      integer bcr1
c
c      do 100 c=1,dimvc
c      cr1=c+add
c      do 100 a=2,dimvc
c      ar1=a+add
c      ab0=nshf(a)
c      do 100 b=1,a-1
c      bcr1=indab(b+add,cr1)
c      v(ab0+b,c)=r1(ar1,bcr1)
c100    continue
c
       do 101 c=1,dimvc
       cr1=c+add
       do 101 a=2,dimvc
       ar1=a+add
       ab0=nshf(a)
       if (c.le.(a-2)) then
c       b1 <= c1
        bk=cr1*(cr1-1)/2+add
        do b=1,c
c       bcr1=cr1*(cr1-1)/2+add+b
        v(ab0+b,c)=r1(ar1,bk+b)
        end do
c       b1 > c1
        bk=(add+c+1)*(add+c)/2+cr1
        do b=c+1,a-1
c       bcr1=(add+b)*(add+b-1)/2+cr1
        v(ab0+b,c)=r1(ar1,bk)
        bk=bk+add+b
        end do
       else
        bk=cr1*(cr1-1)/2+add
        do b=1,a-1
c       b1 <= c1
c       bcr1=cr1*(cr1-1)/2+add+b
        v(ab0+b,c)=r1(ar1,bk+b)
        end do
       end if
 101   continue
c
       do 200 c=1,dimvc
       cr1=c+add
       do 200 a=2,dimvc
       ab0=nshf(a)
c      acr1=indab(a+add,cr1)
         if ((a+add).gt.cr1) then
           acr1=(a+add)*(a+add-1)/2+cr1
         else
           acr1=cr1*(cr1-1)/2+a+add
         end if
       do 200 b=1,a-1
       ab=ab0+b
       v(ab,c)=v(ab,c)-r1(b+add,acr1)
 200    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp21 (r1,v,dimr1a,dimr1c,
     & dimvab,dimva,dimvc,adda,addc)
c
c     this routine do
c     V(ab,c)xxx = R1(a,b,c)-R1(b,a,c) x=a,b
c     for syma=symb>symc
c
c     r1        - r1 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a (b) in R1 (I)
c     dimr1c         - dimension of c in R1 (I)
c     dimvab        - dimension of ab in V (I)
c     dimva        - dimension of a (b) in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (b) (I)
c     addc    - additional constat to c (I)
c
#include "t31.fh"
       integer dimr1a,dimr1c,dimvab,dimva,dimvc,adda,addc
       real*8 r1(1:dimr1a,1:dimr1a,1:dimr1c)
       real*8 v(1:dimvab,1:dimvc)
c
c     help variables
c
       integer a,b,c,ab,ab0,ar1,cr1

c
       do 100 c=1,dimvc
       cr1=c+addc
       do 100 a=2,dimva
       ar1=a+adda
       ab0=nshf(a)
       do 100 b=1,a-1
       v(ab0+b,c)=r1(ar1,b+adda,cr1)
 100    continue
c
       do 200 c=1,dimvc
       cr1=c+addc
       do 200 a=2,dimva
       ar1=a+adda
       ab0=nshf(a)
       do 200 b=1,a-1
       ab=ab0+b
       v(ab,c)=v(ab,c)-r1(b+adda,ar1,cr1)
 200    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp22 (r1,v,dimr1a,dimr1c,
     & dimvab,dimva,dimvc,adda,addc)
c
c     this routine do
c     V(ab,c)xxx = R1(a,c,b)-R1(b,c,a) x=a,b
c     for syma=symb<symc
c
c     r1        - r1 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a (b) in R1 (I)
c     dimr1c         - dimension of c in R1 (I)
c     dimvab        - dimension of ab in V (I)
c     dimva        - dimension of a (b) in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (b) (I)
c     addc    - additional constat to c (I)
c
#include "t31.fh"
       integer dimr1a,dimr1c,dimvab,dimva,dimvc,adda,addc
       real*8 r1(1:dimr1a,1:dimr1c,1:dimr1a)
       real*8 v(1:dimvab,1:dimvc)
c
c     help variables
c
       integer a,b,c,ab,ab0,ar1,cr1

c
       do 100 c=1,dimvc
       cr1=c+addc
       do 100 a=2,dimva
       ar1=a+adda
       ab0=nshf(a)
       do 100 b=1,a-1
       v(ab0+b,c)=r1(ar1,cr1,b+adda)
 100    continue
c
       do 200 a=2,dimva
       ar1=a+adda
       ab0=nshf(a)
       do 200 c=1,dimvc
       cr1=c+addc
       do 200 b=1,a-1
       ab=ab0+b
       v(ab,c)=v(ab,c)-r1(b+adda,cr1,ar1)
 200    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp3 (r1,r2,v,dimr1a,dimr1b,dimr1c,dimr2ac,
     & dimva,dimvb,dimvc,adda,addb,addc)
c
c     this routine do
c     V(a,b,c)xxx = R1(a,c,b)-R2(b,ac) x=a,b
c     for syma.ne.symb symc.eq.syma
c
c     r1        - r1 matrix (I)
c     r2        - r2 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a in R1 (I)
c     dimr1b         - dimension of b in R1 (I)
c     dimr1c         - dimension of c in R1 (I)
c     dimr2ac - dimension of ac in R2 (I)
c     dimva         - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c     addb    - additional constat to b (I)
c     addc    - additional constat to c (I)
c
#include "t31.fh"
       integer dimr1a,dimr1b,dimr1c,dimr2ac
       integer dimva,dimvb,dimvc,adda,addb,addc
       real*8 r1(1:dimr1a,1:dimr1c,1:dimr1b)
       real*8 r2(1:dimr1b,1:dimr2ac)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c,br1,cr1,cr2,acr2
c
c
       do 100 b=1,dimvb
       br1=b+addb
       do 100 c=1,dimvc
       cr1=c+addc
       do 100 a=1,dimva
       v(a,b,c)=r1(a+adda,cr1,br1)
 100    continue
c
       do 200 c=1,dimvc
       cr2=c+addc
       do 200 a=1,dimvc
c      acr2=indab(a+adda,cr2)
         if ((a+adda).ge.cr2) then
           acr2=(a+adda)*(a+adda-1)/2+cr2
         else
           acr2=cr2*(cr2-1)/2+a+adda
         end if
       do 200 b=1,dimvb
       v(a,b,c)=v(a,b,c)-r2(b+addb,acr2)
 200    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp4 (r1,r2,v,dimr1a,dimr1bc,dimr2b,dimr2c,
     & dimva,dimvb,dimvc,adda,addb,addc)
c
c     this routine do
c     V(a,b,c)xxx = R1(a,bc)-R2(b,a,c) x=a,b
c     for syma.ne.symb symc.eq.symb
c
c     r1        - r1 matrix (I)
c     r2        - r2 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a in R1 (I)
c     dimr1bc        - dimension of bc in R1 (I)
c     dimr2b         - dimension of b in R2 (I)
c     dimr2c         - dimension of c in R2 (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c     addb    - additional constat to b (I)
c     addc    - additional constat to c (I)
c
#include "t31.fh"
       integer dimr1a,dimr1bc,dimr2b,dimr2c
       integer dimva,dimvb,dimvc,adda,addb,addc
       real*8 r1(1:dimr1a,1:dimr1bc)
       real*8 r2(1:dimr2b,1:dimr1a,dimr2c)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c,bcr1,cr1,br2,cr2
c
c
       do 100 c=1,dimvc
       cr1=c+addc
       do 100 b=1,dimvb
c      bcr1=indab(b+addb,cr1)
         if ((b+addb).gt.cr1) then
           bcr1=(b+addb)*(b+addb-1)/2+cr1
         else
           bcr1=cr1*(cr1-1)/2+b+addb
         end if
       do 100 a=1,dimva
       v(a,b,c)=r1(a+adda,bcr1)
 100    continue
c
       do 200 c=1,dimvc
       cr2=c+addc
       do 200 b=1,dimvb
       br2=b+addb
       do 200 a=1,dimva
       v(a,b,c)=v(a,b,c)-r2(br2,a+adda,cr2)
 200    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp51 (r1,r2,v,dimr1a,dimr1b,dimr1c,
     & dimva,dimvb,dimvc,adda,addb,addc)
c
c     this routine do
c     V(a,b,c)xxx = R1(a,b,c)-R2(b,a,c) x=a,b
c     for syma>symc, symb>symc, (syma>symb)
c
c     r1        - r1 matrix (I)
c     r2        - r2 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a in R (I)
c     dimr1b         - dimension of b in R (I)
c     dimr1c         - dimension of c in R (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c     addb    - additional constat to b (I)
c     addc    - additional constat to c (I)
c
#include "t31.fh"
       integer dimr1a,dimr1b,dimr1c
       integer dimva,dimvb,dimvc,adda,addb,addc
       real*8 r1(1:dimr1a,1:dimr1b,1:dimr1c)
       real*8 r2(1:dimr1b,1:dimr1a,1:dimr1c)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c,br1,cr1,br2,cr2
c
c
       do 100 c=1,dimvc
       cr1=c+addc
       do 100 b=1,dimvb
       br1=b+addb
       do 100 a=1,dimva
       v(a,b,c)=r1(a+adda,br1,cr1)
 100    continue
c
       do 200 c=1,dimvc
       cr2=c+addc
       do 200 b=1,dimvb
       br2=b+addb
       do 200 a=1,dimva
       v(a,b,c)=v(a,b,c)-r2(br2,a+adda,cr2)
 200    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp52 (r1,r2,v,dimr1a,dimr1b,dimr1c,
     & dimva,dimvb,dimvc,adda,addb,addc)
c
c     this routine do
c     V(a,b,c)xxx = R1(a,b,c)-R2(b,a,c) x=a,b
c     for syma>symc, symb<symc, (syma>symb)
c
c     r1        - r1 matrix (I)
c     r2        - r2 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a in R (I)
c     dimr1b         - dimension of b in R (I)
c     dimr1c         - dimension of c in R (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c     addb    - additional constat to b (I)
c     addc    - additional constat to c (I)
c
#include "t31.fh"
       integer dimr1a,dimr1b,dimr1c
       integer dimva,dimvb,dimvc,adda,addb,addc
       real*8 r1(1:dimr1a,1:dimr1c,1:dimr1b)
       real*8 r2(1:dimr1b,1:dimr1a,1:dimr1c)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c,br1,cr1,br2,cr2
c
c
       do 100 b=1,dimvb
       br1=b+addb
       do 100 c=1,dimvc
       cr1=c+addc
       do 100 a=1,dimva
       v(a,b,c)=r1(a+adda,cr1,br1)
 100    continue
c
       do 200 c=1,dimvc
       cr2=c+addc
       do 200 b=1,dimvb
       br2=b+addb
       do 200 a=1,dimva
       v(a,b,c)=v(a,b,c)-r2(br2,a+adda,cr2)
 200    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp53 (r1,r2,v,dimr1a,dimr1b,dimr1c,
     & dimva,dimvb,dimvc,adda,addb,addc)
c
c     this routine do
c     V(a,b,c)xxx = R1(a,b,c)-R2(b,a,c) x=a,b
c     for syma<symc, symb>symc, (syma>symb)
c
c     r1        - r1 matrix (I)
c     r2        - r2 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a in R (I)
c     dimr1b         - dimension of b in R (I)
c     dimr1c         - dimension of c in R (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c     addb    - additional constat to b (I)
c     addc    - additional constat to c (I)
c
#include "t31.fh"
       integer dimr1a,dimr1b,dimr1c
       integer dimva,dimvb,dimvc,adda,addb,addc
       real*8 r1(1:dimr1a,1:dimr1b,1:dimr1c)
       real*8 r2(1:dimr1b,1:dimr1c,1:dimr1b)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c,br1,cr1,br2,cr2
c
c
       do 100 c=1,dimvc
       cr1=c+addc
       do 100 b=1,dimvb
       br1=b+addb
       do 100 a=1,dimva
       v(a,b,c)=r1(a+adda,br1,cr1)
 100    continue
c
       do 200 c=1,dimvc
       cr2=c+addc
       do 200 b=1,dimvb
       br2=b+addb
       do 200 a=1,dimva
       v(a,b,c)=v(a,b,c)-r2(br2,cr2,a+adda)
 200    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp54 (r1,r2,v,dimr1a,dimr1b,dimr1c,
     & dimva,dimvb,dimvc,adda,addb,addc)
c
c     this routine do
c     V(a,b,c)xxx = R1(a,b,c)-R2(b,a,c) x=a,b
c     for syma>symc, symb>symc, (syma>symb)
c
c     r1        - r1 matrix (I)
c     r2        - r2 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a in R (I)
c     dimr1b         - dimension of b in R (I)
c     dimr1c         - dimension of c in R (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c     addb    - additional constat to b (I)
c     addc    - additional constat to c (I)
c
#include "t31.fh"
       integer dimr1a,dimr1b,dimr1c
       integer dimva,dimvb,dimvc,adda,addb,addc
       real*8 r1(1:dimr1a,1:dimr1c,1:dimr1b)
       real*8 r2(1:dimr1b,1:dimr1c,1:dimr1a)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c,br1,cr1,br2,cr2
c
c
       do 100 b=1,dimvb
       br1=b+addb
       do 100 c=1,dimvc
       cr1=c+addc
       do 100 a=1,dimva
       v(a,b,c)=r1(a+adda,cr1,br1)
 100    continue
c
       do 200 c=1,dimvc
       cr2=c+addc
       do 200 b=1,dimvb
       br2=b+addb
       do 200 a=1,dimva
       v(a,b,c)=v(a,b,c)-r2(br2,cr2,a+adda)
 200    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp61 (r1,v,dimr1a,dimr1b,dimr1c,
     & dimva,dimvb,dimvc,adda)
c
c     this routine do
c     V(a,b,c)abb = R1(a,b,c)
c     for symb>symc
c
c     r1        - r1 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a in R1 (I)
c     dimr1b         - dimension of b in R1 (I)
c     dimr1c         - dimension of c in R1 (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c
       integer dimr1a,dimr1b,dimr1c
       integer dimva,dimvb,dimvc,adda
       real*8 r1(1:dimr1a,1:dimr1b,1:dimr1c)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c
c
c
       do 100 c=1,dimvc
       do 100 b=1,dimvb
       do 100 a=1,dimva
       v(a,b,c)=r1(a+adda,b,c)
 100    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp62 (r1,v,dimr1a,dimr1b,dimr1c,
     & dimva,dimvb,dimvc,adda)
c
c     this routine do
c     V(a,b,c)abb = R1(a,b,c)
c     for symb<symc
c
c     r1        - r1 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a in R1 (I)
c     dimr1b         - dimension of b in R1 (I)
c     dimr1c         - dimension of c in R1 (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c
       integer dimr1a,dimr1b,dimr1c
       integer dimva,dimvb,dimvc,adda
       real*8 r1(1:dimr1a,1:dimr1c,1:dimr1b)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c
c
c
       do 100 b=1,dimvb
       do 100 c=1,dimvc
       do 100 a=1,dimva
       v(a,b,c)=r1(a+adda,c,b)
 100    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp7 (r1,v,dimr1a,dimr1b,dimr1bc,
     & dimva,dimvb,dimvc,adda)
c
c     this routine do
c     V(a,b,c)abb = R1(a,bc)
c     for symb.eq.symc
c
c     r1        - r1 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a in R1 (I)
c     dimr1b         - dimension of b (c) in R1 (I)
c     dimr1bc        - dimension of bc in R1 (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c
#include "t31.fh"
       integer dimr1a,dimr1b,dimr1bc
       integer dimva,dimvb,dimvc,adda
       real*8 r1(1:dimr1a,1:dimr1bc)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c,bcr1
c
c
       do 100 c=1,dimvc
       do 100 b=1,dimvb
c      bcr1=indab(b,c)
         if (b.ge.c) then
           bcr1=b*(b-1)/2+c
         else
           bcr1=c*(c-1)/2+b
         end if
       do 100 a=1,dimva
       v(a,b,c)=r1(a+adda,bcr1)
 100    continue
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dimr1b)
       end
c
c     ---------------
c
       subroutine defvhlp81 (r2,v,dimr2b,dimr2a,dimr2c,
     & dimva,dimvb,dimvc,adda,addc)
c
c     this routine do
c     V(a,b,c)aba = - R2(b,a,c)
c     for syma>symc
c
c     r2        - r2 matrix (I)
c     v        - v matrix (O)
c     dimr2b         - dimension of b in R2 (I)
c     dimr2a         - dimension of a in R2 (I)
c     dimr2c         - dimension of c in R2 (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c     addc    - additional constat to c (I)
c
       integer dimr2b,dimr2a,dimr2c
       integer dimva,dimvb,dimvc,adda,addc
       real*8 r2(1:dimr2b,1:dimr2a,1:dimr2c)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c,ar2,cr2
c
c
       do 100 c=1,dimvc
       cr2=c+addc
       do 100 a=1,dimva
       ar2=a+adda
       do 100 b=1,dimvb
       v(a,b,c)=-r2(b,ar2,cr2)
 100    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp82 (r2,v,dimr2b,dimr2a,dimr2c,
     & dimva,dimvb,dimvc,adda,addc)
c
c     this routine do
c     V(a,b,c)aba = - R2(b,a,c)
c     for syma<symc
c
c     r2        - r2 matrix (I)
c     v        - v matrix (O)
c     dimr2b         - dimension of b in R2 (I)
c     dimr2a         - dimension of a in R2 (I)
c     dimr2c         - dimension of c in R2 (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c     addc    - additional constat to c (I)
c
       integer dimr2b,dimr2a,dimr2c
       integer dimva,dimvb,dimvc,adda,addc
       real*8 r2(1:dimr2b,1:dimr2c,1:dimr2a)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c,ar2,cr2
c
c
       do 100 a=1,dimva
       ar2=a+adda
       do 100 c=1,dimvc
       cr2=c+addc
       do 100 b=1,dimvb
       v(a,b,c)=-r2(b,cr2,ar2)
 100    continue
c
       return
       end
c
c     ---------------
c
       subroutine defvhlp9 (r2,v,dimr2b,dimr2a,dimr2ac,
     & dimva,dimvb,dimvc,adda,addc)
c
c     this routine do
c     V(a,b,c)aba = - R2(b,ac)
c     for syma.eq.symc
c
c     r2        - r2 matrix (I)
c     v        - v matrix (O)
c     dimr2b         - dimension of b in R2 (I)
c     dimr2a         - dimension of a (c) in R2 (I)
c     dimr2ac - dimension of ac in R2 (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c     addc    - additional constat to c (I)
c
#include "t31.fh"
       integer dimr2b,dimr2a,dimr2ac
       integer dimva,dimvb,dimvc,adda,addc
       real*8 r2(1:dimr2b,1:dimr2ac)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c,cr2,acr2
c
c
       do 100 c=1,dimvc
       cr2=c+addc
       do 100 a=1,dimva
c      acr2=indab(a+adda,cr2)
         if ((a+adda).ge.cr2) then
           acr2=(a+adda)*(a+adda-1)/2+cr2
         else
           acr2=cr2*(cr2-1)/2+a+adda
         end if
       do 100 b=1,dimvb
       v(a,b,c)=-r2(b,acr2)
 100    continue
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dimr2a)
       end
c
c     ---------------
c
