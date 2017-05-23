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
c     divt
c     divthelp1
c     divthelp2
c     divthelp3
c     mktau
c     mktauhelp1
c     mktauhelp2
c     mkq
c     mkqhelp1
c     mkqhelp2
c     max5
c     max5h1
c     max5h2
c
c     --------------------------------------------------------
c
       subroutine divt (wrk,wrksize,
     & nind,mapdt,mapit,mapddp1,mapidp1,mapddp2,
     &                  mapidp2,rc)
c
c     this routine divide T amplitudes by denominators, ie differences in dp
c
c     nint         - number of indexes in T (2 or 4) (I)
c     mapdt   - direct map of T (I)
c     mapit        - inverse map of T (I)
c     mapddp1        - direct map of dpa (I)
c     mapidp1        - inverse map of dpa (I)
c     mapddp2        - direct map of dpb (I)
c     mapidp2        - inverse map of dpb (I)
c     rc        - return (error) code
c
c
#include "ccsd1.fh"
#include "wrk.fh"
       integer nind,rc
c
       integer mapdt(0:512,1:6)
       integer mapit(1:8,1:8,1:8)
c
       integer mapddp1(0:512,1:6)
       integer mapidp1(1:8,1:8,1:8)
c
       integer mapddp2(0:512,1:6)
       integer mapidp2(1:8,1:8,1:8)
c
c     help variables
c
       integer posst,possdp,possdpa,possdpb,possdpi,possdpj
       integer dimi,dimj,dima,dimb,dimab,dimij,syma,symb,symi,symj
       integer iit,iidp,iidpa,iidpb,iidpi,iidpj
c
       rc=0
c
       if (nind.eq.2) then
cI    T1 amplitudes - 2 indexes
c
       if (mapdt(0,1).eq.3) then
cI.1  T1aa case
c
       do 100 iit=1,mapdt(0,5)
c
       posst=mapdt(iit,1)
       syma=mapdt(iit,3)
       dima=nva(syma)
       dimi=noa(syma)
       iidp=mapidp1(syma,1,1)
       possdp=mapddp1(iidp,1)
c
       if ((dima*dimi).gt.0) then
       call divthelp1 (wrk(posst),dima,dimi,wrk(possdp))
       end if
c
 100    continue
c
       else if (mapdt(0,1).eq.4) then
cI.2  T1bb case
c
       do 200 iit=1,mapdt(0,5)
c
       posst=mapdt(iit,1)
       syma=mapdt(iit,3)
       dima=nvb(syma)
       dimi=nob(syma)
       iidp=mapidp2(syma,1,1)
       possdp=mapddp2(iidp,1)
c
       if ((dima*dimi).gt.0) then
       call divthelp1 (wrk(posst),dima,dimi,wrk(possdp))
       end if
c
 200    continue
c
       else
cI.3  invalid T1 type
c     RC=1 : incorrect map for T1
       rc=1
       return
       end if
c
       else if (nind.eq.4) then
cII   T2 amplitudes - 4 indexes
c
       if (mapdt(0,6).eq.0) then
cII.1 T2abab case

       do 300 iit=1,mapdt(0,5)
c
       posst=mapdt(iit,1)
       syma=mapdt(iit,3)
       symb=mapdt(iit,4)
       symi=mapdt(iit,5)
       symj=mapdt(iit,6)
       dima=nva(syma)
       dimb=nvb(symb)
       dimi=noa(symi)
       dimj=nob(symj)
       iidpa=mapidp1(syma,1,1)
       iidpb=mapidp2(symb,1,1)
       iidpi=mapidp1(symi,1,1)
       iidpj=mapidp2(symj,1,1)
       possdpa=mapddp1(iidpa,1)
       possdpb=mapddp2(iidpb,1)
       possdpi=mapddp1(iidpi,1)
       possdpj=mapddp2(iidpj,1)
c
       if (mapdt(iit,2).gt.0) then
       call divthelp2 (wrk(posst),dima,dimb,dimi,dimj,
     & wrk(possdpa),wrk(possdpb),wrk(possdpi),wrk(possdpj),
     & noa(syma),nob(symb))
       end if
c
 300    continue
c
       else if ((mapdt(0,6).eq.4).and.(mapdt(0,1).eq.3)) then
cII.2 T2aaaa case
c
       do 400 iit=1,mapdt(0,5)
c
       posst=mapdt(iit,1)
       syma=mapdt(iit,3)
       symb=mapdt(iit,4)
       symi=mapdt(iit,5)
       symj=mapdt(iit,6)
       dima=nva(syma)
       dimb=nva(symb)
       dimi=noa(symi)
       dimj=noa(symj)
       iidpa=mapidp1(syma,1,1)
       iidpb=mapidp1(symb,1,1)
       iidpi=mapidp1(symi,1,1)
       iidpj=mapidp1(symj,1,1)
       possdpa=mapddp1(iidpa,1)
       possdpb=mapddp1(iidpb,1)
       possdpi=mapddp1(iidpi,1)
       possdpj=mapddp1(iidpj,1)
c
       if (mapdt(iit,2).eq.0) goto 400

       if (syma.ne.symb) then
c     different symmetries a,b; i,j
       call divthelp2 (wrk(posst),dima,dimb,dimi,dimj,
     & wrk(possdpa),wrk(possdpb),wrk(possdpi),wrk(possdpj),
     & noa(syma),noa(symb))
c
       else
c     same symmetries a,b; i,j
       dimab=(dima*(dima-1))/2
       dimij=(dimi*(dimi-1))/2
       call divthelp3 (wrk(posst),dimab,dimij,wrk(possdpa),wrk(possdpi),
     & dima,dimi,noa(syma))
       end if
c
 400    continue
c
       else if ((mapdt(0,6).eq.4).and.(mapdt(0,1).eq.4)) then
cII.3 T2bbbb case
c
       do 500 iit=1,mapdt(0,5)
c
       posst=mapdt(iit,1)
       syma=mapdt(iit,3)
       symb=mapdt(iit,4)
       symi=mapdt(iit,5)
       symj=mapdt(iit,6)
       dima=nvb(syma)
       dimb=nvb(symb)
       dimi=nob(symi)
       dimj=nob(symj)
       iidpa=mapidp2(syma,1,1)
       iidpb=mapidp2(symb,1,1)
       iidpi=mapidp2(symi,1,1)
       iidpj=mapidp2(symj,1,1)
       possdpa=mapddp2(iidpa,1)
       possdpb=mapddp2(iidpb,1)
       possdpi=mapddp2(iidpi,1)
       possdpj=mapddp2(iidpj,1)
c
       if (mapdt(iit,2).eq.0) goto 500

       if (syma.ne.symb) then
c     different symmetries a,b; i,j
       call divthelp2 (wrk(posst),dima,dimb,dimi,dimj,
     & wrk(possdpa),wrk(possdpb),wrk(possdpi),wrk(possdpj),
     & nob(syma),nob(symb))
c
       else
c     same symmetries a,b; i,j
       dimab=(dima*(dima-1))/2
       dimij=(dimi*(dimi-1))/2
       call divthelp3 (wrk(posst),dimab,dimij,wrk(possdpa),wrk(possdpi),
     & dima,dimi,nob(syma))
       end if
c
 500    continue
c
       else
cII.4 RC=2 : incorrect mapdt for T2
       rc=2
       return
       end if
c
       else
cIII  invalid nind
c     RC=3 : nind is not 2 or 4 (Stup)
       rc=3
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(mapit)
       end
c
c     ----------
c
       subroutine divthelp1 (t1,dima,dimi,dp)
c
c     this routine do
c     t1(a,i) = t1(a,i)/(dp(i)-dp(a))
c
c     t1        - T1 matrix (I/O)
c     dima    - v dimension of T1 (I)
c     dimi    - o domension of T1 (I)
c     dp      - diagonal part of Fok (I)
c
c     N.B. Since for T1 i and a are of the same spin, there is no reason
c     to specify spin of dp. It must be automatically the same spin ad i and a.
c
       integer dima,dimi
       real*8 t1(1:dima,1:dimi)
       real*8 dp(*)
c
c     help variables
c
       integer a,i
       real*8 dpi,den
c
       do 100 i=1,dimi
       dpi=dp(i)
       do 100 a=1,dima
c     t1(a,i)=t1(a,i)/(dpi-dp(dimi+a))
c
       den=dpi-dp(dimi+a)
       if (abs(den).lt.1.0d-7) then
       if (abs(t1(a,i)).gt.1.0d-10) then
       t1(a,i)=t1(a,i)/den
       end if
       else
       t1(a,i)=t1(a,i)/den
       end if
c
 100    continue
c
       return
       end
c
c     ----------
c
       subroutine divthelp2 (t2,dima,dimb,dimi,dimj,dpa,dpb,dpi,dpj,
     &                       shifta,shiftb)
c
c     this routine do
c     t2(a,b,i,j) = t2(a,b,i,j)/(dpi(i)+dpj(j)-dpa(a)-dpb(b))
c
c     t2        - T2 matrix (I/O)
c     dima    - 1 dimension of T2 (I)
c     dimb    - 2 dimension of T2 (I)
c     dimi    - 3 domension of T2 (I)
c     dimj    - 4 domension of T2 (I)
c     dpa     - diagonal part of Fok corresponding to irrep of a (I)
c     dpb     - diagonal part of Fok corresponding to irrep of b (I)
c     dpi     - diagonal part of Fok corresponding to irrep of i (I)
c     dpj     - diagonal part of Fok corresponding to irrep of j (I)
c     shifta        - number off occ orbitals in spin and symmetry of a (I)
c     shiftb        - number off occ orbitals in spin and symmetry of b (I)
c
c     N.B. Since for T1 i and a are of the same spin, there is no reason
c     to specify spin of dp. It must be automatically the same spin ad i and a.
c
       integer dima,dimb,dimi,dimj,shifta,shiftb
       real*8 t2(1:dima,1:dimb,1:dimi,1:dimj)
       real*8 dpa(*)
       real*8 dpb(*)
       real*8 dpi(*)
       real*8 dpj(*)
c
c     help variables
c
       integer i,j,a,b
       real*8 den,denj,denij,denijb
c
       do 100 j=1,dimj
       denj=dpj(j)
       do 100 i=1,dimi
       denij=denj+dpi(i)
       do 100 b=1,dimb
       denijb=denij-dpb(shiftb+b)
       do 100 a=1,dima
c     t2(a,b,i,j)=t2(a,b,i,j)/(denijb-dpa(shifta+a))
c
       den=denijb-dpa(shifta+a)
       if (abs(den).lt.1.0d-7) then
       if (abs(t2(a,b,i,j)).gt.1.0d-10) then
       t2(a,b,i,j)=t2(a,b,i,j)/den
       end if
       else
       t2(a,b,i,j)=t2(a,b,i,j)/den
       end if
c
 100    continue
c
       return
       end
c
c     ----------
c
       subroutine divthelp3 (t2,dimab,dimij,dpa,dpi,dima,dimi,shift)
c
c     this routine do
c     t2(ab,ij) = t2(ab,ij)/(dp(i)+dp(j)-dp(a)-dp(b))
c     for spin and symmetry of a = b
c     and spin and symmetry of i = j
c
c     t2        - T2 matrix (I/O)
c     dimab   - 1 dimension of T2 (I)
c     dimij   - 3 domension of T2 (I)
c     dpa     - diagonal part of Fok corresponding to irrep of a (I)
c     dpi     - diagonal part of Fok corresponding to irrep of i (I)
c     dima    - number of a in this spin and symm of a (I)
c     dimi    - number of a in this spin and symm of i (I)
c     shift         - number off occ orbitals in spin and symmetry of a (I)
c
c     N.B. Since for T1 i and a are of the same spin, there is no reason
c     to specify spin of dp. It must be automatically the same spin ad i and a.
c
       integer dimab,dimij,dima,dimi,shift
       real*8 t2(1:dimab,1:dimij)
       real*8 dpa(*)
       real*8 dpi(*)
c
c     help variables
c
       integer i,j,a,b,ij,ab
       real*8 den,deni,denij,denija
c
       ij=0
       do 100 i=2,dimi
       deni=dpi(i)
       do 100 j=1,i-1
       denij=deni+dpi(j)
       ij=ij+1
c
       ab=0
       do 100 a=2,dima
       denija=denij-dpa(shift+a)
       do 100 b=1,a-1
       ab=ab+1
c     t2(ab,ij)=t2(ab,ij)/(denija-dpa(shift+b))
c
       den=denija-dpa(shift+b)
       if (abs(den).lt.1.0d-7) then
       if (abs(t2(ab,ij)).gt.1.0d-10) then
       t2(ab,ij)=t2(ab,ij)/den
       end if
       else
       t2(ab,ij)=t2(ab,ij)/den
       end if
c
c
 100    continue
c
       return
       end
c
c     ---------------------------------------------------
c
       subroutine mktau (wrk,wrksize,
     & mapdt2,mapit2,mapdt1a,mapit1a,mapdt1b,mapit1b,
     &                   fact,rc)
c
c     this routine do:
c     t2(abij) = t2(abij) + fact. (t1(ai).t1(bj)-t1(bi).t1(aj))
c     N.B. T24a,4b must be of type 4, T2abab of type 0
c
c     mapdt2  - direct map of T2 (I)
c     mapit2        - inverse map of T2 (I)
c     mapdt1a        - direct map of T1aa (I)
c     mapit1a        - inverse map of T1aa (I)
c     mapdt1b        - direct map of T1bb (I)
c     mapit1b        - inverse map of T1bb (I)
c     fact    - numerical factor (I)
c     rc        - return (error) code
c
c
#include "ccsd1.fh"
#include "wrk.fh"
       integer rc
       real*8 fact
c
       integer mapdt2(0:512,1:6)
       integer mapit2(1:8,1:8,1:8)
c
       integer mapdt1a(0:512,1:6)
       integer mapit1a(1:8,1:8,1:8)
c
       integer mapdt1b(0:512,1:6)
       integer mapit1b(1:8,1:8,1:8)
c
c     help variables
c
       integer posst2,posst1a,posst1b,posst11,posst12
       integer dimi,dimj,dima,dimb,dimab,dimij,syma,symb,symi,symj
       integer iit2,iit1a,iit1b,iit11,iit12
c
       rc=0
c
       if (mapdt2(0,6).eq.0) then
cI.1  T2abab case

       do 100 iit2=1,mapdt2(0,5)
c
       posst2=mapdt2(iit2,1)
       syma=mapdt2(iit2,3)
       symb=mapdt2(iit2,4)
       symi=mapdt2(iit2,5)
       symj=mapdt2(iit2,6)
       dima=nva(syma)
       dimb=nvb(symb)
       dimi=noa(symi)
       dimj=nob(symj)
       iit1a=mapit1a(syma,1,1)
       iit1b=mapit1b(symb,1,1)
       posst1a=mapdt1a(iit1a,1)
       posst1b=mapdt1b(iit1b,1)
c
       if ((syma.eq.symi).and.(symb.eq.symj).and.(mapdt2(iit2,2).gt.0))
     & then
       call mktauhelp1 (wrk(posst2),wrk(posst1a),wrk(posst1b),
     & dima,dimb,dimi,dimj,noa(symi),nob(symj),fact)
       end if
c
 100    continue
c
       else if ((mapdt2(0,6).eq.4).and.(mapdt2(0,1).eq.3)) then
cI.2  T2aaaa case
c
       do 200 iit2=1,mapdt2(0,5)
c
       posst2=mapdt2(iit2,1)
       syma=mapdt2(iit2,3)
       symb=mapdt2(iit2,4)
       symi=mapdt2(iit2,5)
       symj=mapdt2(iit2,6)
       dima=nva(syma)
       dimb=nva(symb)
       dimi=noa(symi)
       dimj=noa(symj)
       iit11=mapit1a(syma,1,1)
       iit12=mapit1a(symb,1,1)
       posst11=mapdt1a(iit11,1)
       posst12=mapdt1a(iit12,1)
c
       if ((syma.eq.symi).and.(symb.eq.symj).
     & and.(syma.ne.symj).and.(mapdt2(iit2,2).gt.0)) then
cI.2.*case T2(sym1,sym2,sym1,sym2)
c
       call mktauhelp1 (wrk(posst2),wrk(posst11),wrk(posst12),
     & dima,dimb,dimi,dimj,noa(syma),noa(symb),fact)
c
       else if ((syma.eq.symi).and.(symb.eq.symj).
     & and.(syma.eq.symj).and.(mapdt2(iit2,2).gt.0)) then
cI.2.*case T2(sym1,sym1,sym1,sym1)
c
       dimab=(dima*(dima-1))/2
       dimij=(dimi*(dimi-1))/2
       call mktauhelp2 (wrk(posst2),wrk(posst11),
     & dimab,dimij,dima,dimi,noa(syma),fact)
c
       end if
c
 200    continue
c
       else if ((mapdt2(0,6).eq.4).and.(mapdt2(0,1).eq.4)) then
cI.3  T2bbbb case
c
       do 300 iit2=1,mapdt2(0,5)
c
       posst2=mapdt2(iit2,1)
       syma=mapdt2(iit2,3)
       symb=mapdt2(iit2,4)
       symi=mapdt2(iit2,5)
       symj=mapdt2(iit2,6)
       dima=nvb(syma)
       dimb=nvb(symb)
       dimi=nob(symi)
       dimj=nob(symj)
       iit11=mapit1b(syma,1,1)
       iit12=mapit1b(symb,1,1)
       posst11=mapdt1b(iit11,1)
       posst12=mapdt1b(iit12,1)
c
       if ((syma.eq.symi).and.(symb.eq.symj).
     & and.(syma.ne.symj).and.(mapdt2(iit2,2).gt.0)) then
cI.3.*case T2(sym1,sym2,sym1,sym2)
c
       call mktauhelp1 (wrk(posst2),wrk(posst11),wrk(posst12),
     & dima,dimb,dimi,dimj,nob(syma),nob(symb),fact)
c
       else if ((syma.eq.symi).and.(symb.eq.symj).
     & and.(syma.eq.symj).and.(mapdt2(iit2,2).gt.0)) then
cI.3.*case T2(sym1,sym1,sym1,sym1)
c
       dimab=(dima*(dima-1))/2
       dimij=(dimi*(dimi-1))/2
       call mktauhelp2 (wrk(posst2),wrk(posst11),
     & dimab,dimij,dima,dimi,nob(syma),fact)
c
       end if
c
 300    continue
c
       else
cI.4  RC=1 : incorrect mapdt for T2
       rc=1
       return
       end if
c
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(mapit2)
       end
c
c     ----------
c
       subroutine mktauhelp1 (t2,t11,t12,dima,dimb,dimi,dimj,shifta,
     &                        shiftb,fact)
c
c     this routine do
c     t2(a,b,i,j) = t2(a,b,i,j) + fact [T11(a,i).T12(b,j)]
c
c     t2        - T2 matrix (I/O)
c     t11     - T1 amplitudes corresponding to spin ia (I)
c     t12     - T1 amplitudes corresponding to spin jb (I)
c     dima    - 1 dimension of T2 (I)
c     dimb    - 2 dimension of T2 (I)
c     dimi    - 3 domension of T2 (I)
c     dimj    - 4 domension of T2 (I)
c     shifta        - number off occ orbitals in spin and symmetry of a (I)
c     shiftb        - number off occ orbitals in spin and symmetry of b (I)
c     fact    - numerical factor (I)
c
c     N.B. symi must be syma and symj must be symb
c
c
       integer dima,dimb,dimi,dimj,shifta,shiftb
       real*8 fact
       real*8 t2(1:dima,1:dimb,1:dimi,1:dimj)
       real*8 t11(1:dima,1:dimi)
       real*8 t12(1:dimb,1:dimj)
c
c     help variables
c
       integer i,j,a,b
c
       do 100 j=1,dimj
       do 100 i=1,dimi
       do 100 b=1,dimb
       do 100 a=1,dima
       t2(a,b,i,j)=t2(a,b,i,j)+fact*(t11(a,i)*t12(b,j))
 100    continue
c
       return
c Avoid unused argument warnings
       if (.false.) then
         call Unused_integer(shifta)
         call Unused_integer(shiftb)
       end if
       end
c
c     ----------
c
       subroutine mktauhelp2 (t2,t1,dimab,dimij,dima,dimi,shift,fact)
c
c     this routine do
c     t2(ab,ij) = t2(ab,ij) + fact* (T1(a,i).t1(b,j)-t1(b,i).t1(a,j))
c     for spin and symmetry of all indices equal
c
c     t2        - T2 matrix (I/O)
c     t1      - t1 matrix (I)
c     dimab   - 1 dimension of T2 (I)
c     dimij   - 3 domension of T2 (I)
c     dima    - number of a in this spin and symm of a (I)
c     dimi    - number of a in this spin and symm of i (I)
c     shift         - number off occ orbitals in spin and symmetry of a (I)
c     fact    - numerical factor (I)
c
c     N.B. Since for T1 i and a are of the same spin, there is no reason
c     to specify spin of dp. It must be automatically the same spin ad i and a.
c
       integer dimab,dimij,dima,dimi,shift
       real*8 fact
       real*8 t2(1:dimab,1:dimij)
       real*8 t1(1:dima,1:dimi)
c
c     help variables
c
       integer i,j,a,b,ij,ab
c
       ij=0
       do 100 i=2,dimi
       do 100 j=1,i-1
       ij=ij+1
c
       ab=0
       do 100 a=2,dima
       do 100 b=1,a-1
       ab=ab+1
       t2(ab,ij)=t2(ab,ij)+fact*(t1(a,i)*t1(b,j)-t1(b,i)*t1(a,j))
c
 100    continue
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(shift)
       end
c
c     ---------------------------------------------------
c
       subroutine mkq (wrk,wrksize,
     & mapdt2,mapit2,mapdt11,mapit11,mapdt12,mapit12,
     &                 fact,rc)
c
c     this routine do:
c     t2(a,b,i,j) = fact . t2(a,b,i,j) +t11(ai).t12(bj)
c     for T2aaaa, T2bbbb and T2abab but they must be expanded (typ=0)
c
c     mapdt2  - direct map of T2 (I)
c     mapit2        - inverse map of T2 (I)
c     mapdt11        - direct map of T1 (I)
c     mapit11        - inverse map of T1 (I)
c     mapdt12        - direct map of T11 (I)
c     mapit12        - inverse map of T12 (I)
c     fact    - numerical factor (I)
c     rc        - return (error) code
c
c
#include "ccsd1.fh"
#include "wrk.fh"
       integer rc
       real*8 fact
c
       integer mapdt2(0:512,1:6)
       integer mapit2(1:8,1:8,1:8)
c
       integer mapdt11(0:512,1:6)
       integer mapit11(1:8,1:8,1:8)
c
       integer mapdt12(0:512,1:6)
       integer mapit12(1:8,1:8,1:8)
c
c     help variables
c
       integer posst2,posst11,posst12
       integer dimi,dimj,dima,dimb,syma,symb,symi,symj
       integer iit2,iit11,iit12
c
       rc=0
c
       if (mapdt2(0,6).eq.0) then
cI.1  typ of t2 is 0 (T2 is expanded)

       do 100 iit2=1,mapdt2(0,5)
c
       posst2=mapdt2(iit2,1)
       syma=mapdt2(iit2,3)
       symb=mapdt2(iit2,4)
       symi=mapdt2(iit2,5)
       symj=mapdt2(iit2,6)
       dima=dimm(mapdt2(0,1),syma)
       dimb=dimm(mapdt2(0,2),symb)
       dimi=dimm(mapdt2(0,3),symi)
       dimj=dimm(mapdt2(0,4),symj)
       iit11=mapit11(syma,1,1)
       iit12=mapit12(symb,1,1)
       posst11=mapdt11(iit11,1)
       posst12=mapdt12(iit12,1)
c
       if ((syma.eq.symi).and.(symb.eq.symj).and.(mapdt2(iit2,2).gt.0))
     & then
       call mkqhelp1 (wrk(posst2),wrk(posst11),wrk(posst12),
     & dima,dimb,dimi,dimj,noa(symi),nob(symj),fact)
       else if (mapdt2(iit2,2).gt.0) then
       call mkqhelp2 (wrk(posst2),mapdt2(iit2,2),mapdt2(iit2,2),fact)
       end if
c
 100    continue
c
c
       else
cI.4  RC=1 : typ of T2 is not 0
       rc=1
       return
       end if
c
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(mapit2)
       end
c
c     ----------
c
       subroutine mkqhelp1 (t2,t11,t12,dima,dimb,dimi,dimj,shifta,
     &                      shiftb,fact)
c
c     this routine do
c     t2(a,b,i,j) = fact . t2(a,b,i,j) [T11(a,i).T12(b,j)]
c
c     t2        - T2 matrix (I/O)
c     t11     - T1 amplitudes corresponding to spin ia (I)
c     t12     - T1 amplitudes corresponding to spin jb (I)
c     dima    - 1 dimension of T2 (I)
c     dimb    - 2 dimension of T2 (I)
c     dimi    - 3 domension of T2 (I)
c     dimj    - 4 domension of T2 (I)
c     shifta        - number off occ orbitals in spin and symmetry of a (I)
c     shiftb        - number off occ orbitals in spin and symmetry of b (I)
c     fact    - numerical factor (I)
c
c     N.B. symi must be syma and symj must be symb
c
c
       integer dima,dimb,dimi,dimj,shifta,shiftb
       real*8 fact
       real*8 t2(1:dima,1:dimb,1:dimi,1:dimj)
       real*8 t11(1:dima,1:dimi)
       real*8 t12(1:dimb,1:dimj)
c
c     help variables
c
       integer i,j,a,b
c
       do 100 j=1,dimj
       do 100 i=1,dimi
       do 100 b=1,dimb
       do 100 a=1,dima
       t2(a,b,i,j)=fact*t2(a,b,i,j)+(t11(a,i)*t12(b,j))
 100    continue
c
       return
c Avoid unused argument warnings
       if (.false.) then
         call Unused_integer(shifta)
         call Unused_integer(shiftb)
       end if
       end
c
c     ----------
c
       subroutine mkqhelp2 (vector,dimv,lenght,factor)
c
c     this routine do vector = vector*factot
c     vector - multilyied vector (I/O)
c     dimv   - dimension of vecrot
c     lenght - lenght of vector to be multiplyied
c     factor - scaling factor
c
c     $N.B. this routine should be substitued by mv0s3v
c
       integer dimv,lenght
       real*8 vector(1:dimv)
       real*8 factor
c
c     help variable
c
       integer n
c
       if (lenght.gt.0) then
       do 10 n=1,lenght
       vector(n)=vector(n)*factor
 10     continue
       end if
c
       return
       end
c
c     ---------------------------------------------------
c
       subroutine max5 (wrk,wrksize,
     & nind,mapd,mapi,text)
c
c     this routine find and type:
c     a) note
c     b) 5 maximal elements with their indexes in given vector V
c     c) euclidian norm
c
c     nind  - number of indexes in V (I)
c     mapd  - direct map of V (I)
c     mapi  - inverese map of V (I)
c     text  - notice (I)
c
#include "ccsd1.fh"
#include "wrk.fh"
c
       integer nind
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
       character*8 text
c
c     help variables
c
       integer nhelp1,nhelp2,i,j,a,b,it
       real*8 value
       integer imax(1:8,1:5)
       real*8 rmax (1:5)
c
c0    set rmax,imax=0
c
       do 10 nhelp1=1,5
       rmax(nhelp1)=0.0d0
       do 10 nhelp2=1,8
       imax(nhelp2,nhelp1)=0
 10     continue
c
c
       if (nind.eq.2) then
c
c1    T1aa or T1bb amplitudes
c
c1.1  find 5 max
c
       nhelp1=mapd(1,1)
       do 100 it=1,mapd(0,5)
       do 50 i=1,dimm(mapd(0,2),mapd(it,4))
       do 50 a=1,dimm(mapd(0,1),mapd(it,3))
       value=wrk(nhelp1)
       if (abs(value).ge.abs(rmax(5))) then
c     write this amplitude
       call max5h1 (imax,rmax,mapd(it,3),0,mapd(it,4),0,a,0,i,0,
     & value)
       end if
       nhelp1=nhelp1+1
 50     continue
 100    continue
c
c1.2  write output
       If (fullprint.ge.0) call max5h2 (wrk,wrksize,
     & nind,mapd,mapi,rmax,imax,text)
c
       else if (mapd(0,6).eq.0) then
c
c2    T2abab amplitudes
c
c2.1  find 5 max
c
       nhelp1=mapd(1,1)
       do 200 it=1,mapd(0,5)
       do 150 j=1,dimm(mapd(0,4),mapd(it,6))
       do 150 i=1,dimm(mapd(0,3),mapd(it,5))
       do 150 b=1,dimm(mapd(0,2),mapd(it,4))
       do 150 a=1,dimm(mapd(0,1),mapd(it,3))
       value=wrk(nhelp1)
       if (abs(value).ge.abs(rmax(5))) then
c     write this amplitude
       call max5h1 (imax,rmax,mapd(it,3),mapd(it,4),
     & mapd(it,5),mapd(it,6),a,b,i,j,value)
       end if
       nhelp1=nhelp1+1
 150    continue
 200    continue
c
c2.2  write output
       If (fullprint.ge.0) call max5h2 (wrk,wrksize,
     & nind,mapd,mapi,rmax,imax,text)
c
       else
c
c3    T2aaaa or T2bbbb amplitudes
c
c3.1  find 5 max
c
       nhelp1=mapd(1,1)
       do 300 it=1,mapd(0,5)
c
       if (mapd(it,3).eq.mapd(it,4)) then
c     case syma=symb, symi=symj
c
       do 230 i=2,dimm(mapd(0,3),mapd(it,5))
       do 230 j=1,i-1
       do 230 a=2,dimm(mapd(0,1),mapd(it,3))
       do 230 b=1,a-1
       value=wrk(nhelp1)
       if (abs(value).ge.abs(rmax(5))) then
c     write this amplitude
       call max5h1 (imax,rmax,mapd(it,3),mapd(it,4),
     & mapd(it,5),mapd(it,6),a,b,i,j,value)
       end if
       nhelp1=nhelp1+1
 230    continue
c
       else
c     case syma>symb, symi> symj
       do 250 j=1,dimm(mapd(0,4),mapd(it,6))
       do 250 i=1,dimm(mapd(0,3),mapd(it,5))
       do 250 b=1,dimm(mapd(0,2),mapd(it,4))
       do 250 a=1,dimm(mapd(0,1),mapd(it,3))
       value=wrk(nhelp1)
       if (abs(value).ge.abs(rmax(5))) then
c     write this amplitude
       call max5h1 (imax,rmax,mapd(it,3),mapd(it,4),
     & mapd(it,5),mapd(it,6),a,b,i,j,value)
       end if
       nhelp1=nhelp1+1
 250    continue
       end if
c
 300    continue
c
c3.2  write output
       If (fullprint.ge.0) call max5h2 (wrk,wrksize,
     & nind,mapd,mapi,rmax,imax,text)
c
       end if
c
       return
       end
c
c     ---------------------
c
       subroutine max5h1 (imax,rmax,symp,symq,symr,syms,p,q,r,s,value)
c
c     this routine add new max amplitude and skip smallest
c
c     imax - store of indexes of 5 max (I/O)
c     rmax - store of values of 5 max (I/O)
c     symp - symmetry of p index (I)
c     symq - symmetry of q index (I)
c     symr - symmetry of r index (I)
c     syms - symmetry of s index (I)
c     p    - value of p index (I)
c     q    - value of q index (I)
c     r    - value of r index (I)
c     s    - value of s index (I)
c     value- value of amplitude (I)
c
       integer imax(1:8,1:5)
       real*8 rmax (1:5)
       integer symp,symq,symr,syms,p,q,r,s
       real*8 value
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp3
c
c1    find possition of this value
c
       do 10 nhelp1=1,5
       if (abs(value).ge.abs(rmax(nhelp1))) then
       goto 20
       end if
 10     continue
c
c2    1-(nhelp1-1) stay untought
c
c3    push other records if necc.
c
 20     if (nhelp1.lt.5) then
       do 30 nhelp2=4,nhelp1,-1
       rmax(nhelp2+1)=rmax(nhelp2)
       do 30 nhelp3=1,8
       imax(nhelp3,nhelp2+1)=imax(nhelp3,nhelp2)
 30     continue
       end if
c
c4    add new one
c
       rmax(nhelp1)=value
       imax(1,nhelp1)=symp
       imax(2,nhelp1)=symq
       imax(3,nhelp1)=symr
       imax(4,nhelp1)=syms
       imax(5,nhelp1)=p
       imax(6,nhelp1)=q
       imax(7,nhelp1)=r
       imax(8,nhelp1)=s
c
       return
       end
c
c     ---------------------
c
       subroutine max5h2 (wrk,wrksize,
     & nind,mapd,mapi,rmax,imax,text)
c
c     this routine write:
c     a) note
c     b) 5 maximal elements with their indexes in given vector V
c     c) euclidian norm
c
c     nind  - number of indexes in V (I)
c     mapd  - direct map of V (I)
c     mapi  - inverese map of V (I)
c     rmax  - store of maximal values (I)
c     imax  - store of corr. indexes (I)
c     text  - notice (I)
c
#include "ccsd1.fh"
#include "wrk.fh"
c
       integer nind
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
       integer imax(1:8,1:5)
       real*8 rmax (1:5)
       character*8 text
c
c     help variables
c
       integer nhelp1,nhelp2,rc
       real*8 scalar
c
c1    write notice
c
       write(6,101) text
 101    format (' Five largest amplitudes of :',a8)
c
c2    write 5 maximal amplitudes
c
       write(6,102)
 102    format ('  SYMA   SYMB   SYMI   SYMJ     A      B',
     &          '      I      J     VALUE')
       do 10 nhelp1=1,5
       write(6,103) (imax(nhelp2,nhelp1),nhelp2=1,8),rmax(nhelp1)
 103    format (8(2x,i3,2x),f15.10)
 10     continue
c
c3    write euclidian norm
c
c3.1  calc euclidian norm
       call multdot (wrk,wrksize,
     & nind,mapd,mapi,1,mapd,mapi,1,scalar,rc)
       scalar=sqrt(scalar)
c
       write(6,104) scalar
 104    format (' Euclidian norm is :',f17.10)
c
       write(6,*)
c
       return
       end
c
c     ---------------------------------------------------
c
