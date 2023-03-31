!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine saamphlp1 (t24a,t24b,t22b,noa,nob,nva,nvb,key)
!
!     adaptation routine for T2 amplitudes for symi=symj, syma=symb
!
!     t24a  - amplitides T2(ab,ij)aaaa (I/O)
!     t24b  - amplitides T2(ab,ij)bbbb (I/O)
!     t22b  - amplitides T2(a,b,i,j)abab (I/O)
!     noa   - number of alpha occupied orbitals in symi (I)
!     nob   - number of beta occupied orbitals in symi (I)
!     nva   - number of alpha virtual orbitals in syma (I)
!     nvb   - number of beta virtual orbitals in syma (I)
!     key   - type of adaptation (I)
!     0 - no adaptation
!     1 - T2 DDVV adaptation
!     2 - T2 DDVV + T1 DV adaptation
!     3 - full T1 and T2 adaptation (only for doublets)
!     4 - full T2 without SDVS (only for doublets)
!
       integer noa,nob,nva,nvb
!T2
       real*8 t24a(1:(nva*(nva-1))/2,1:(noa*(noa-1))/2)
       real*8 t24b(1:(nvb*(nvb-1))/2,1:(nob*(nob-1))/2)
       real*8 t22b(1:nva,1:nvb,1:noa,1:nob)
!
       integer key
!
!     help variables
!
       integer nd,nv,nsi,nsa
       integer i,j,a,b,ij,ab,ab1
       real*8 taaaa,tbbbb,tabab,tabba,tbaab,tbaba,t1,t2
!
       if (key.eq.0) return
!
       nd=nob
       nv=nva
       nsi=noa-nob
       nsa=nvb-nva
!
!
!     T2 adatption
!
!     case I) DDVV i>=j,a>=b
!     turn od in any case of type of adaption
!
!     aaaa=t24a(ijab) DDVV
!     bbbb=t24b(ijab) DDVV
!     abab=t22b(ijab) DDVV
!     abba=-t22b(ijba) DDVV
!     baab=-t22b(jiab) DDVV
!     baba=t22b(ji,ba) DDVV
!
!     direct
!     t1 = (abab+baba-abba-baab)/4
!     t2 = (2aaaa+2bbbb+abab+baba+abba+baab)/12
!
!     reverse
!
!     aaaa = 2t2
!     bbbb = 2t2
!     abab = t2+t1
!     abba = t2-t1
!     baab = t2-t1
!     baba = t2+t1
!
       do 100 a=2,nv
       do 1000 b=1,a-1
       ab=(a-1)*(a-2)/2+b
       ab1=(a+nsa-1)*(a+nsa-2)/2+b+nsa
       do 1001 i=2,nd
       do 1002 j=1,i-1
       ij=(i-1)*(i-2)/2+j
!
       taaaa=t24a(ab,ij)
       tbbbb=t24b(ab1,ij)
       tabab=t22b(a,b+nsa,i,j)
       tabba=-t22b(b,a+nsa,i,j)
       tbaab=-t22b(a,b+nsa,j,i)
       tbaba=t22b(b,a+nsa,j,i)
!
       t1=(tabab+tbaba-tabba-tbaab)/4.0d0
       t2=(2.0d0*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/1.2d1
!
       t24a(ab,ij)=2.0d0*t2
       t24b(ab1,ij)=2.0d0*t2
       t22b(a,b+nsa,i,j)=t2+t1
       t22b(b,a+nsa,i,j)=t1-t2
       t22b(a,b+nsa,j,i)=t1-t2
       t22b(b,a+nsa,j,i)=t2+t1
!
 1002   continue
 1001   continue
 1000   continue
 100    continue
!
       do 101 a=2,nv
       do 1010 b=1,a-1
       ab=(a-1)*(a-2)/2+b
       ab1=(a+nsa-1)*(a+nsa-2)/2+b+nsa
       do 1011 i=1,nd
       j=i
       ij=(i-1)*(i-2)/2+j
!
!     taaaa=t24a(ab,ij)
!     tbbbb=t24b(ab1,ij)
       tabab=t22b(a,b+nsa,i,j)
       tabba=-t22b(b,a+nsa,i,j)
       tbaab=-t22b(a,b+nsa,j,i)
       tbaba=t22b(b,a+nsa,j,i)
!
       t1=(tabab+tbaba-tabba-tbaab)/4.0d0
!     t2=(2.0d0*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/1.2d1
       t2=0.0d0
!
!     t24a(ab,ij)=2.0d0*t2
!     t24b(ab1,ij)=2.0d0*t2
       t22b(a,b+nsa,i,j)=t2+t1
       t22b(b,a+nsa,i,j)=t1-t2
       t22b(a,b+nsa,j,i)=t1-t2
       t22b(b,a+nsa,j,i)=t2+t1
!
 1011   continue
 1010   continue
 101    continue
!
       do 102 a=1,nv
       b=a
       ab=(a-1)*(a-2)/2+b
       ab1=(a+nsa-1)*(a+nsa-2)/2+b+nsa
       do 1020 i=2,nd
       do 1021 j=1,i-1
       ij=(i-1)*(i-2)/2+j
!
!     taaaa=t24a(ab,ij)
!     tbbbb=t24b(ab1,ij)
       tabab=t22b(a,b+nsa,i,j)
       tabba=-t22b(b,a+nsa,i,j)
       tbaab=-t22b(a,b+nsa,j,i)
       tbaba=t22b(b,a+nsa,j,i)
!
       t1=(tabab+tbaba-tabba-tbaab)/4.0d0
!     t2=(2.0d0*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/1.2d1
       t2=0.0d0
!
!     t24a(ab,ij)=2.0d0*t2
!     t24b(ab1,ij)=2.0d0*t2
       t22b(a,b+nsa,i,j)=t2+t1
       t22b(b,a+nsa,i,j)=t1-t2
       t22b(a,b+nsa,j,i)=t1-t2
       t22b(b,a+nsa,j,i)=t2+t1
!
 1021   continue
 1020   continue
 102    continue
!
       do 103 a=1,nv
       b=a
       ab=(a-1)*(a-2)/2+b
       ab1=(a+nsa-1)*(a+nsa-2)/2+b+nsa
       do 1030 i=1,nd
       j=i
       ij=(i-1)*(i-2)/2+j
!
!     taaaa=t24a(ab,ij)
!     tbbbb=t24b(ab1,ij)
       tabab=t22b(a,b+nsa,i,j)
       tabba=-t22b(b,a+nsa,i,j)
       tbaab=-t22b(a,b+nsa,j,i)
       tbaba=t22b(b,a+nsa,j,i)
!
       t1=(tabab+tbaba-tabba-tbaab)/4.0d0
!     t2=(2.0d0*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/1.2d1
       t2=0.0d0
!
!     t24a(ab,ij)=2.0d0*t2
!     t24b(ab1,ij)=2.0d0*t2
       t22b(a,b+nsa,i,j)=t2+t1
       t22b(b,a+nsa,i,j)=t1-t2
       t22b(a,b+nsa,j,i)=t1-t2
       t22b(b,a+nsa,j,i)=t2+t1
!
 1030   continue
 103    continue
!
!
       if (((key.eq.3).or.(key.eq.4)).and.(nsa.gt.0)) then
!
!     case II) DDVS i>=j,a,b
!
!     bbbb = t24b(ij,ab)
!     abab = t22b(i,j,a,b)
!     baab =-t22b(j,i,a,b)
!
!     direct :
!
!     t1=(abab-baab)/2
!     t2=(2bbbb+abab+baab)/6
!
!     reverse:
!
!     bbbb = 2t2
!     abab = t2+t1
!     baab = t2-t1
!
       b=nsa
       do 200 a=1,nv
       ab=a*(a-1)/2+b
       do 201 i=2,nd
       do 202 j=1,i-1
       ij=(i-1)*(i-2)/2+j
!
       tbbbb=t24b(ab,ij)
       tabab=t22b(a,b,i,j)
       tbaab=-t22b(a,b,j,i)
!
       t1=(tabab-tbaab)/2.0d0
       t2=(2.0d0*tbbbb+tabab+tbaab)/6.0d0
!
       t24b(ab,ij)=2.0d0*t2
       t22b(a,b,i,j)=t2+t1
       t22b(a,b,j,i)=t1-t2
!
 202    continue
 201    continue
 200    continue
!
       b=nsa
       do 203 a=1,nv
       ab=a*(a-1)/2+b
       do 204 i=1,nd
       j=i
       ij=(i-1)*(i-2)/2+j
!
!     tbbbb=t24b(ab,ij)
       tabab=t22b(a,b,i,j)
       tbaab=-t22b(a,b,j,i)
!
       t1=(tabab-tbaab)/2.0d0
!     t2=(2.0d0*tbbbb+tabab+tbaab)/6.0d0
       t2=0.0d0
!
!     t24b(ab,ij)=2.0d0*t2
       t22b(a,b,i,j)=t2+t1
       t22b(a,b,j,i)=t1-t2
!
 204    continue
 203    continue
!
       end if
!
!
       if (((key.eq.3).or.(key.eq.4)).and.(nsi.gt.0)) then
!
!     case III) SDVV i,j,a>=b
!
!     aaaa=t24a(ij,ab)
!     abab=t22b(i,j,a,b)
!     abba=-t22b(i,j,b,a)
!
!     direct
!
!     t1=(abab-abba)/2
!     t2=(2aaaa+abab+abba)/6
!
!     reverse
!
!     aaaa=2t2
!     abab=t2+t1
!     abba=t2-t1
!
       i=nd+nsi
       do 300 a=2,nv
       do 301 b=1,a-1
       ab=(a-1)*(a-2)/2+b
       do 302 j=1,nd
       ij=(i-1)*(i-2)/2+j
!
       taaaa=t24a(ab,ij)
       tabab=t22b(a,b+nsa,i,j)
       tabba=-t22b(b,a+nsa,i,j)
!
       t1=(tabab-tabba)/2.0d0
       t2=(2.0d0*taaaa+tabab+tabba)/6.0d0
!
       t24a(ab,ij)=2.0d0*t2
       t22b(a,b+nsa,i,j)=t2+t1
       t22b(b,a+nsa,i,j)=t1-t2
!
 302    continue
 301    continue
 300    continue
!
       i=nd+nsi
       do 303 a=1,nv
       b=a
       ab=(a-1)*(a-2)/2+b
       do 304 j=1,nd
       ij=(i-1)*(i-2)/2+j
!
!     taaaa=t24a(ab,ij)
       tabab=t22b(a,b+nsa,i,j)
       tabba=-t22b(b,a+nsa,i,j)
!
       t1=(tabab-tabba)/2.0d0
!     t2=(2.0d0*taaaa+tabab+tabba)/6.0d0
       t2=0.0d0
!
!     t24a(ab,ij)=2.0d0*t2
       t22b(a,b+nsa,i,j)=t2+t1
       t22b(b,a+nsa,i,j)=t1-t2
!
 304    continue
 303    continue
!
       end if
!
       return
       end
