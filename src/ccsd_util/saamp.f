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
c     this file contains following routines
c     saamp
c     saamphlp1
c     saamphlp2
c     saamphlp3
c
c     --------------------------------------------------------
c
       subroutine saamp (wrk,wrksize,
     & key)
c
c     this routine rearrange amplitudes to be spin adapted
c     key - 0 - no adaptation
c     1 - T2 DDVV adaptation
c     2 - T2 DDVV + T1 DV adaptation
c     3 - full T1 and T2 adaptation (only for doublets)
c     4 - full T2 without SDVS (only for doublets)
c
c     amplitudes T1 are in t13 - aa, t14 - bb
c     T2 are in t21 - aaaa, t22 - bbbb, t23 - abab
c
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
c
       integer key
c
c     help variables
c
       integer symi,symj,syma,symb,syms,symij
       integer poss1,poss2,poss3,poss4,poss5,poss6,ii
c
c0    skip this routine if SA in not turn on
       if (key.eq.0) then
       return
       end if
c
c
cI    T1 adaptation
       if ((key.eq.2).or.(key.eq.3)) then
c
cI.1  def symmetry, where S orbital is situated (only for doublet states)
       syms=0
       do 10 symi=1,nsym
       if (dimm(1,symi).ne.dimm(2,symi)) then
       syms=symi
       end if
 10     continue
       if ((key.eq.2).and.(syms.eq.0)) then
       syms=1
       end if
       if (syms.eq.0) then
       write(6,*) ' Full SA is turn on and there is no S orbitals'
       Call Abend
       end if
c
cI.2  loop over symi
       do 100 symi=1,nsym
       syma=symi
c
       ii=mapit13(syma,1,1)
       poss1=mapdt13(ii,1)
       ii=mapit14(syma,1,1)
       poss2=mapdt14(ii,1)
       ii=mapit23(syma,syms,syms)
       poss3=mapdt23(ii,1)
       call saamphlp3 (wrk(poss1),wrk(poss2),wrk(poss3),
     & dimm(1,symi),dimm(2,symi),dimm(3,symi),dimm(4,symi),
     & dimm(1,syms),dimm(4,syms),key)
c
 100    continue
       end if
c
c
cII   T2 adaptation
c
       do 200 symi=1,nsym
       do 201 symj=1,symi
       symij=mmul(symi,symj)
c
       do 202 syma=1,nsym
       symb=mmul(symij,syma)
c
       if (symb.gt.syma) then
c     Meggie out
       goto 202
       end if
c
       if (symi.eq.symj) then
c     case si=sj, sa=sb
c
       ii=mapit21(syma,symb,symi)
       poss1=mapdt21(ii,1)
       ii=mapit22(syma,symb,symi)
       poss2=mapdt22(ii,1)
       ii=mapit23(syma,symb,symi)
       poss3=mapdt23(ii,1)
       call saamphlp1 (wrk(poss1),wrk(poss2),wrk(poss3),
     & dimm(1,symi),dimm(2,symi),dimm(3,syma),dimm(4,syma),key)
c
       else
c     case si>sj, sa>sb
c
       ii=mapit21(syma,symb,symi)
       poss1=mapdt21(ii,1)
       ii=mapit22(syma,symb,symi)
       poss2=mapdt22(ii,1)
       ii=mapit23(syma,symb,symi)
       poss3=mapdt23(ii,1)
       ii=mapit23(symb,syma,symj)
       poss4=mapdt23(ii,1)

       ii=mapit23(symb,syma,symi)
       poss5=mapdt23(ii,1)
       ii=mapit23(syma,symb,symj)
       poss6=mapdt23(ii,1)
       call saamphlp2 (wrk(poss1),wrk(poss2),wrk(poss3),
     & wrk(poss4),wrk(poss5),wrk(poss6),
     & dimm(1,symi),dimm(1,symj),dimm(2,symi),dimm(2,symj),
     & dimm(3,syma),dimm(3,symb),dimm(4,syma),dimm(4,symb),
     & key)
c
       end if
c
 202    continue
 201    continue
 200    continue
c
       return
       end
c
c     -----------------------------
c
       subroutine saamphlp1 (t24a,t24b,t22b,noa,nob,nva,nvb,key)
c
c     adaptation routine for T2 amplitudes for symi=symj, syma=symb
c
c     t24a  - amplitides T2(ab,ij)aaaa (I/O)
c     t24b  - amplitides T2(ab,ij)bbbb (I/O)
c     t22b  - amplitides T2(a,b,i,j)abab (I/O)
c     noa   - number of alpha occupied orbitals in symi (I)
c     nob   - number of beta occupied orbitals in symi (I)
c     nva   - number of alpha virtual orbitals in syma (I)
c     nvb   - number of beta virtual orbitals in syma (I)
c     key   - type of adaptation (I)
c     0 - no adaptation
c     1 - T2 DDVV adaptation
c     2 - T2 DDVV + T1 DV adaptation
c     3 - full T1 and T2 adaptation (only for doublets)
c     4 - full T2 without SDVS (only for doublets)
c
       integer noa,nob,nva,nvb
cT2
       real*8 t24a(1:(nva*(nva-1))/2,1:(noa*(noa-1))/2)
       real*8 t24b(1:(nvb*(nvb-1))/2,1:(nob*(nob-1))/2)
       real*8 t22b(1:nva,1:nvb,1:noa,1:nob)
c
       integer key
c
c     help variables
c
       integer nd,nv,nsi,nsa
       integer i,j,a,b,ij,ab,ab1
       real*8 taaaa,tbbbb,tabab,tabba,tbaab,tbaba,t1,t2
c
       if (key.eq.0) return
c
       nd=nob
       nv=nva
       nsi=noa-nob
       nsa=nvb-nva
c
c
c     T2 adatption
c
c     case I) DDVV i>=j,a>=b
c     turn od in any case of type of adaption
c
c     aaaa=t24a(ijab) DDVV
c     bbbb=t24b(ijab) DDVV
c     abab=t22b(ijab) DDVV
c     abba=-t22b(ijba) DDVV
c     baab=-t22b(jiab) DDVV
c     baba=t22b(ji,ba) DDVV
c
c     direct
c     t1 = (abab+baba-abba-baab)/4
c     t2 = (2aaaa+2bbbb+abab+baba+abba+baab)/12
c
c     reverse
c
c     aaaa = 2t2
c     bbbb = 2t2
c     abab = t2+t1
c     abba = t2-t1
c     baab = t2-t1
c     baba = t2+t1
c
       do 100 a=2,nv
       do 1000 b=1,a-1
       ab=(a-1)*(a-2)/2+b
       ab1=(a+nsa-1)*(a+nsa-2)/2+b+nsa
       do 1001 i=2,nd
       do 1002 j=1,i-1
       ij=(i-1)*(i-2)/2+j
c
       taaaa=t24a(ab,ij)
       tbbbb=t24b(ab1,ij)
       tabab=t22b(a,b+nsa,i,j)
       tabba=-t22b(b,a+nsa,i,j)
       tbaab=-t22b(a,b+nsa,j,i)
       tbaba=t22b(b,a+nsa,j,i)
c
       t1=(tabab+tbaba-tabba-tbaab)/4.0d0
       t2=(2.0d0*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/1.2d1
c
       t24a(ab,ij)=2.0d0*t2
       t24b(ab1,ij)=2.0d0*t2
       t22b(a,b+nsa,i,j)=t2+t1
       t22b(b,a+nsa,i,j)=t1-t2
       t22b(a,b+nsa,j,i)=t1-t2
       t22b(b,a+nsa,j,i)=t2+t1
c
 1002   continue
 1001   continue
 1000   continue
 100    continue
c
       do 101 a=2,nv
       do 1010 b=1,a-1
       ab=(a-1)*(a-2)/2+b
       ab1=(a+nsa-1)*(a+nsa-2)/2+b+nsa
       do 1011 i=1,nd
       j=i
       ij=(i-1)*(i-2)/2+j
c
c     taaaa=t24a(ab,ij)
c     tbbbb=t24b(ab1,ij)
       tabab=t22b(a,b+nsa,i,j)
       tabba=-t22b(b,a+nsa,i,j)
       tbaab=-t22b(a,b+nsa,j,i)
       tbaba=t22b(b,a+nsa,j,i)
c
       t1=(tabab+tbaba-tabba-tbaab)/4.0d0
c     t2=(2.0d0*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/1.2d1
       t2=0.0d0
c
c     t24a(ab,ij)=2.0d0*t2
c     t24b(ab1,ij)=2.0d0*t2
       t22b(a,b+nsa,i,j)=t2+t1
       t22b(b,a+nsa,i,j)=t1-t2
       t22b(a,b+nsa,j,i)=t1-t2
       t22b(b,a+nsa,j,i)=t2+t1
c
 1011   continue
 1010   continue
 101    continue
c
       do 102 a=1,nv
       b=a
       ab=(a-1)*(a-2)/2+b
       ab1=(a+nsa-1)*(a+nsa-2)/2+b+nsa
       do 1020 i=2,nd
       do 1021 j=1,i-1
       ij=(i-1)*(i-2)/2+j
c
c     taaaa=t24a(ab,ij)
c     tbbbb=t24b(ab1,ij)
       tabab=t22b(a,b+nsa,i,j)
       tabba=-t22b(b,a+nsa,i,j)
       tbaab=-t22b(a,b+nsa,j,i)
       tbaba=t22b(b,a+nsa,j,i)
c
       t1=(tabab+tbaba-tabba-tbaab)/4.0d0
c     t2=(2.0d0*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/1.2d1
       t2=0.0d0
c
c     t24a(ab,ij)=2.0d0*t2
c     t24b(ab1,ij)=2.0d0*t2
       t22b(a,b+nsa,i,j)=t2+t1
       t22b(b,a+nsa,i,j)=t1-t2
       t22b(a,b+nsa,j,i)=t1-t2
       t22b(b,a+nsa,j,i)=t2+t1
c
 1021   continue
 1020   continue
 102    continue
c
       do 103 a=1,nv
       b=a
       ab=(a-1)*(a-2)/2+b
       ab1=(a+nsa-1)*(a+nsa-2)/2+b+nsa
       do 1030 i=1,nd
       j=i
       ij=(i-1)*(i-2)/2+j
c
c     taaaa=t24a(ab,ij)
c     tbbbb=t24b(ab1,ij)
       tabab=t22b(a,b+nsa,i,j)
       tabba=-t22b(b,a+nsa,i,j)
       tbaab=-t22b(a,b+nsa,j,i)
       tbaba=t22b(b,a+nsa,j,i)
c
       t1=(tabab+tbaba-tabba-tbaab)/4.0d0
c     t2=(2.0d0*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/1.2d1
       t2=0.0d0
c
c     t24a(ab,ij)=2.0d0*t2
c     t24b(ab1,ij)=2.0d0*t2
       t22b(a,b+nsa,i,j)=t2+t1
       t22b(b,a+nsa,i,j)=t1-t2
       t22b(a,b+nsa,j,i)=t1-t2
       t22b(b,a+nsa,j,i)=t2+t1
c
 1030   continue
 103    continue
c
c
       if (((key.eq.3).or.(key.eq.4)).and.(nsa.gt.0)) then
c
c     case II) DDVS i>=j,a,b
c
c     bbbb = t24b(ij,ab)
c     abab = t22b(i,j,a,b)
c     baab =-t22b(j,i,a,b)
c
c     direct :
c
c     t1=(abab-baab)/2
c     t2=(2bbbb+abab+baab)/6
c
c     reverse:
c
c     bbbb = 2t2
c     abab = t2+t1
c     baab = t2-t1
c
       b=nsa
       do 200 a=1,nv
       ab=a*(a-1)/2+b
       do 201 i=2,nd
       do 202 j=1,i-1
       ij=(i-1)*(i-2)/2+j
c
       tbbbb=t24b(ab,ij)
       tabab=t22b(a,b,i,j)
       tbaab=-t22b(a,b,j,i)
c
       t1=(tabab-tbaab)/2.0d0
       t2=(2.0d0*tbbbb+tabab+tbaab)/6.0d0
c
       t24b(ab,ij)=2.0d0*t2
       t22b(a,b,i,j)=t2+t1
       t22b(a,b,j,i)=t1-t2
c
 202    continue
 201    continue
 200    continue
c
       b=nsa
       do 203 a=1,nv
       ab=a*(a-1)/2+b
       do 204 i=1,nd
       j=i
       ij=(i-1)*(i-2)/2+j
c
c     tbbbb=t24b(ab,ij)
       tabab=t22b(a,b,i,j)
       tbaab=-t22b(a,b,j,i)
c
       t1=(tabab-tbaab)/2.0d0
c     t2=(2.0d0*tbbbb+tabab+tbaab)/6.0d0
       t2=0.0d0
c
c     t24b(ab,ij)=2.0d0*t2
       t22b(a,b,i,j)=t2+t1
       t22b(a,b,j,i)=t1-t2
c
 204    continue
 203    continue
c
       end if
c
c
       if (((key.eq.3).or.(key.eq.4)).and.(nsi.gt.0)) then
c
c     case III) SDVV i,j,a>=b
c
c     aaaa=t24a(ij,ab)
c     abab=t22b(i,j,a,b)
c     abba=-t22b(i,j,b,a)
c
c     direct
c
c     t1=(abab-abba)/2
c     t2=(2aaaa+abab+abba)/6
c
c     reverse
c
c     aaaa=2t2
c     abab=t2+t1
c     abba=t2-t1
c
       i=nd+nsi
       do 300 a=2,nv
       do 301 b=1,a-1
       ab=(a-1)*(a-2)/2+b
       do 302 j=1,nd
       ij=(i-1)*(i-2)/2+j
c
       taaaa=t24a(ab,ij)
       tabab=t22b(a,b+nsa,i,j)
       tabba=-t22b(b,a+nsa,i,j)
c
       t1=(tabab-tabba)/2.0d0
       t2=(2.0d0*taaaa+tabab+tabba)/6.0d0
c
       t24a(ab,ij)=2.0d0*t2
       t22b(a,b+nsa,i,j)=t2+t1
       t22b(b,a+nsa,i,j)=t1-t2
c
 302    continue
 301    continue
 300    continue
c
       i=nd+nsi
       do 303 a=1,nv
       b=a
       ab=(a-1)*(a-2)/2+b
       do 304 j=1,nd
       ij=(i-1)*(i-2)/2+j
c
c     taaaa=t24a(ab,ij)
       tabab=t22b(a,b+nsa,i,j)
       tabba=-t22b(b,a+nsa,i,j)
c
       t1=(tabab-tabba)/2.0d0
c     t2=(2.0d0*taaaa+tabab+tabba)/6.0d0
       t2=0.0d0
c
c     t24a(ab,ij)=2.0d0*t2
       t22b(a,b+nsa,i,j)=t2+t1
       t22b(b,a+nsa,i,j)=t1-t2
c
 304    continue
 303    continue
c
       end if
c
       return
       end
c
c     -----------------------------
c
       subroutine saamphlp2 (t2aaaa,t2bbbb,t2abab,t2baba,t2abba,t2baab,
     & noai,noaj,nobi,nobj,nvaa,nvab,nvba,nvbb,key)
c
c     adaptation routine for T2 amplitudes for symi>symj, syma>symb
c
c     t2aaaa- amplitides T2(a,b,i,j)aaaa (I/O)
c     t2bbbb- amplitides T2(a,b,i,j)bbbb (I/O)
c     t2abab- amplitides T2(a,b,i,j)abab (I/O)
c     t2baba- amplitides T2(b,a,j,i)abab (I/O)
c     t2abba- amplitides T2(b,a,i,j)abab (I/O)
c     t2baab- amplitides T2(a,b,j,i)abab (I/O)
c     noai  - number of alpha occupied orbitals in symi (I)
c     noaj  - number of alpha occupied orbitals in symj (I)
c     nobi  - number of beta occupied orbitals in symi (I)
c     nobj  - number of beta occupied orbitals in symj (I)
c     nvaa  - number of alpha virtual orbitals in syma (I)
c     nvab  - number of alpha virtual orbitals in symb (I)
c     nvba  - number of beta virtual orbitals in syma (I)
c     nvbb  - number of beta virtual orbitals in symb (I)
c     key   - type of adaptation (I)
c     0 - no adaptation
c     1 - T2 DDVV adaptation
c     2 - T2 DDVV + T1 DV adaptation
c     3 - full T1 and T2 adaptation (only for doublets)
c     4 - full T2 without SDVS (only for doublets)
c
       integer noai,noaj,nobi,nobj,nvaa,nvab,nvba,nvbb
cT2
       real*8 t2aaaa(1:nvaa,1:nvab,1:noai,1:noaj)
       real*8 t2bbbb(1:nvba,1:nvbb,1:nobi,1:nobj)
       real*8 t2abab(1:nvaa,1:nvbb,1:noai,1:nobj)
       real*8 t2baba(1:nvab,1:nvba,1:noaj,1:nobi)
       real*8 t2abba(1:nvab,1:nvba,1:noai,1:nobj)
       real*8 t2baab(1:nvaa,1:nvbb,1:noaj,1:nobi)
c
       integer key
c
c     help variables
c
       integer ndi,ndj,nva,nvb,nsi,nsj,nsa,nsb
       integer i,j,a,b
       real*8 taaaa,tbbbb,tabab,tabba,tbaab,tbaba,t1,t2
c
       if (key.eq.0) then
       return
       end if
c
       ndi=nobi
       ndj=nobj
       nva=nvaa
       nvb=nvab
       nsi=noai-nobi
       nsj=noaj-nobj
       nsa=nvba-nvaa
       nsb=nvbb-nvab
c
c
c     T2 adatption
c
c     case I) DDVV i>=j,a>=b
c     turn od in any case of type of adaption
c
c     aaaa=t24a(ijab) DDVV
c     bbbb=t24b(ijab) DDVV
c     abab=t22b(ijab) DDVV
c     abba=-t22b(ijba) DDVV
c     baab=-t22b(jiab) DDVV
c     baba=t22b(ji,ba) DDVV
c
c     direct
c     t1 = (abab+baba-abba-baab)/4
c     t2 = (2aaaa+2bbbb+abab+baba+abba+baab)/12
c
c     reverse
c
c     aaaa = 2t2
c     bbbb = 2t2
c     abab = t2+t1
c     abba = t2-t1
c     baab = t2-t1
c     baba = t2+t1
c
       do 100 j=1,ndj
       do 101 i=1,ndi
       do 102 b=1,nvb
       do 103 a=1,nva
c
       taaaa=t2aaaa(a,b,i,j)
       tbbbb=t2bbbb(a+nsa,b+nsb,i,j)
       tabab=t2abab(a,b+nsb,i,j)
       tbaba=t2baba(b,a+nsa,j,i)
       tabba=-t2abba(b,a+nsa,i,j)
       tbaab=-t2baab(a,b+nsb,j,i)
c
       t1=(tabab+tbaba-tabba-tbaab)/4.0d0
       t2=(2.0d0*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/1.2d1
c
       taaaa=2.0d0*t2
       tbbbb=2.0d0*t2
       tabab=t1+t2
       tbaba=t1+t2
       tabba=t2-t1
       tbaab=t2-t1
c
       t2aaaa(a,b,i,j)=taaaa
       t2bbbb(a+nsa,b+nsb,i,j)=tbbbb
       t2abab(a,b+nsb,i,j)=tabab
       t2baba(b,a+nsa,j,i)=tbaba
       t2abba(b,a+nsa,i,j)=-tabba
       t2baab(a,b+nsb,j,i)=-tbaab
c
 103    continue
 102    continue
 101    continue
 100    continue
c
c
       if (((key.eq.3).or.(key.eq.4)).and.(nsb.gt.0)) then
c
c     case II) DDVS
c
c     bbbb = t24b(ij,ab)
c     abab = t22b(i,j,a,b)
c     baab =-t22b(j,i,a,b)
c
c     direct :
c
c     t1=(abab-baab)/2
c     t2=(2bbbb+abab+baab)/6
c
c     reverse:
c
c     bbbb = 2t2
c     abab = t2+t1
c     baab = t2-t1
c
       b=nsb
       do 200 j=1,ndj
       do 201 i=1,ndi
       do 202 a=1,nva
c
       tbbbb=t2bbbb(a+nsa,b,i,j)
       tabab=t2abab(a,b,i,j)
       tbaab=-t2baab(a,b,j,i)
c
       t1=(tabab-tbaab)/2.0d0
       t2=(2.0d0*tbbbb+tabab+tbaab)/6.0d0
c
       tbbbb=2.0d0*t2
       tabab=t2+t1
       tbaab=t2-t1
c
       t2bbbb(a+nsa,b,i,j)=tbbbb
       t2abab(a,b,i,j)=tabab
       t2baab(a,b,j,i)=-tbaab
c
 202    continue
 201    continue
 200    continue
c
       end if
c
c
       if (((key.eq.3).or.(key.eq.4)).and.(nsa.gt.0)) then
c
c     case III) DDSV
c
c     bbbb = - t24b(ij,ab)
c     abab = - -t22b(i,j,b,a)
c     baab = - t22b(j,i,b,a)
c
c     direct :
c
c     t1=(abab-baab)/2
c     t2=(2bbbb+abab+baab)/6
c
c     reverse:
c
c     bbbb = 2t2
c     abab = t2+t1
c     baab = t2-t1
c
       a=nsa
       do 300 j=1,ndj
       do 301 i=1,ndi
       do 302 b=1,nvb
c
       tbbbb=-t2bbbb(a,b+nsb,i,j)
       tabab=t2abba(b,a,i,j)
       tbaab=-t2baba(b,a,j,i)
c
       t1=(tabab-tbaab)/2.0d0
       t2=(2.0d0*tbbbb+tabab+tbaab)/6.0d0
c
       tbbbb=2.0d0*t2
       tabab=t2+t1
       tbaab=t2-t1
c
       t2bbbb(a,b+nsb,i,j)=-tbbbb
       t2abba(b,a,i,j)=tabab
       t2baba(b,a,j,i)=-tbaab
c
 302    continue
 301    continue
 300    continue
c
       end if
c
c
       if (((key.eq.3).or.(key.eq.4)).and.(nsi.gt.0)) then
c
c     case IV) SDVV
c
c     aaaa=t24a(ij,ab)
c     abab=t22b(i,j,a,b)
c     abba=-t22b(i,j,b,a)
c
c     direct
c
c     t1=(abab-abba)/2
c     t2=(2aaaa+abab+abba)/6
c
c     reverse
c
c     aaaa=2t2
c     abab=t2+t1
c     abba=t2-t1
c
       i=ndi+nsi
       do 400 j=1,ndj
       do 401 b=1,nvb
       do 402 a=1,nva
c
       taaaa=t2aaaa(a,b,i,j)
       tabab=t2abab(a,b+nsb,i,j)
       tabba=-t2abba(b,a+nsa,i,j)
c
       t1=(tabab-tabba)/2.0d0
       t2=(2.0d0*taaaa+tabab+tabba)/6.0d0
c
       taaaa=2.0d0*t2
       tabab=t2+t1
       tabba=t2-t1
c
       t2aaaa(a,b,i,j)=taaaa
       t2abab(a,b+nsb,i,j)=tabab
       t2abba(b,a+nsa,i,j)=-tabba
c
 402    continue
 401    continue
 400    continue
c
       end if
c
c
       if (((key.eq.3).or.(key.eq.4)).and.(nsj.gt.0)) then
c
c     case V) DSVV
c
c     aaaa= - t24a(ij,ab)
c     abab= - -t22b(j,i,a,b)
c     abba= - t22b(j,i,b,a)
c
c     direct
c
c     t1=(abab-abba)/2
c     t2=(2aaaa+abab+abba)/6
c
c     reverse
c
c     aaaa=2t2
c     abab=t2+t1
c     abba=t2-t1
c
       j=ndj+nsj
       do 500 i=1,ndi
       do 501 b=1,nvb
       do 502 a=1,nva
c
       taaaa=-t2aaaa(a,b,i,j)
       tabab=t2baab(a,b+nsb,j,i)
       tabba=-t2baba(b,a+nsa,j,i)
c
       t1=(tabab-tabba)/2.0d0
       t2=(2.0d0*taaaa+tabab+tabba)/6.0d0
c
       taaaa=2.0d0*t2
       tabab=t2+t1
       tabba=t2-t1
c
       t2aaaa(a,b,i,j)=-taaaa
       t2baab(a,b+nsb,j,i)=tabab
       t2baba(b,a+nsa,j,i)=-tabba
c
 502    continue
 501    continue
 500    continue
c
       end if
c
       return
       end
c
c     -----------------------------
c
       subroutine saamphlp3
     & (t1aa,t1bb,t2abab,noa,nob,nva,nvb,noas,nvbs,key)
c
c     this routine rearrange amplitudes to be spin adapted
c     for T1 amplitudes for given si(=sa)
c
c     t1aa - T1aa amplitudes for given si (=sa) (I/O)
c     t1bb - T1bb amplitudes for given si (=sa) (I/O)
c     t2abab - T2(arsi)abab amplitudes sa=si (sr=ss=irrep of S) (I/O)
c     noa  - number of alpha occupied in symi (I)
c     nob  - number of beta occupied in symi (I)
c     nva  - number of alpha virtuals in symi (I)
c     nvb  - number of beta virtuals in symi (I)
c     noas - number of alpha occupied in symmetry, where S orbitals is (I)
c     nvbs - number of beta virtuals in symmetry, where S orbitals is (I)
c     key  - 0 - no adaptation (I)
c     1 - T2 DDVV adaptation
c     2 - T2 DDVV + T1 DV adaptation
c     3 - full T1 and T2 adaptation (only for doublets)
c     4 - full T2 without SDVS (only for doublets)
c
       integer noa,nob,nva,nvb,nvbs,noas
cT1
       real*8 t1aa(1:nva,1:noa)
       real*8 t1bb(1:nvb,1:nob)
cT2
       real*8 t2abab(1:nva,1:nvbs,1:noas,1:nob)
c
       integer key
c
c     help variables
c
       integer nd,nv,ns
       integer i,a
       real*8 t1,t2
c
       if (key.eq.0) return
c
       nd=nob
       ns=noa-nob
       nv=nva
c
c     T1 adaption
c
c     ta=t1oaa(i,a)     DV
c     tb=t1obb(i,a)     DV
c     tc=t2o2b(p,i,a,p) SDVS
c
c     direct :
c     t1= (ta+tb)/2
c     t2= (-ta+tb+2tc)/6
c
c     reverse :
c     ta= t1-t2
c     tb= t1+t2
c     tc= 2t2
c
       if (key.eq.3) then
c
       do 10 i=1,nd
       do 11 a=1,nv
c
       t1=(t1aa(a,i)+t1bb(a+ns,i))/2.0d0
       t2=(t1bb(a+ns,i)-t1aa(a,i)+2.0d0*t2abab(a,1,noas,i))/6.0d0
c
       t1aa(a,i)=t1-t2
       t1bb(a+ns,i)=t1+t2
       t2abab(a,1,noas,i)=2.0d0*t2
c
 11     continue
 10     continue
c
       else if (key.eq.2) then
c
       do 40 i=1,nd
       do 41 a=1,nv
c
       t1=(t1aa(a,i)+t1bb(a+ns,i))/2.0d0
c
       t1aa(a,i)=t1
       t1bb(a+ns,i)=t1
c
 41     continue
 40     continue
c
c     else if (key.eq.1) then
c     no adaption in T1 turn on
       end if
c
c
       return
       end
c
c     -----------------------------
c
