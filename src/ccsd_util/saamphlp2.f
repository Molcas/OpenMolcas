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
