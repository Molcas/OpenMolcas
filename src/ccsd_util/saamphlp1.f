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
