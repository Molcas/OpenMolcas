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
