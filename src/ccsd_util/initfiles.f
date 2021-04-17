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
       subroutine initfiles (length,lenv,lenn)
c
c     this routine distribute work space WRK for required files
c     for fix mediates it defines also mapd and mapi, for help mediates
c     it estimates their length and distribute WRK (i.e. def poss0 parameters)
c     !N.B. This routine cannot run with +OP2 level
c
c       length  -  total length of all work space needed
c       lenv    -  length of V - type array
c       lenn    -  length of N - type array
c
#include "ccsd1.fh"
#include "ccsd2.fh"
c
       integer length,lenv,lenn
c
c     help variable
c
       integer posst,symp,symq,sympq,symr,syms
       integer lengthv,lengthm,lengthh,lengthn
       integer maxnoa,maxnvb,maxnorb
       integer maxov(1:8)
c
c1    maps and possitions for fix mediated
c
c1.0  maps for DP - diagonal part
c     N.B. DP has one degree of freedom, while other 1 index has none
c     DP1 - dp(p)a
c     DP2 - dp(p)b
c
       do symp=1,nsym
       do symq=1,nsym
       do symr=1,nsym
       mapidp1(symp,symq,symr)=0
       mapidp2(symp,symq,symr)=0
       end do
       end do
       end do
c
       posst=1
c
       possdp10=posst
       mapddp1(0,1)=5
       mapddp1(0,2)=0
       mapddp1(0,3)=0
       mapddp1(0,4)=0
       mapddp1(0,5)=nsym
       mapddp1(0,6)=0
c
       do symp=1,nsym
       mapddp1(symp,1)=posst
       mapddp1(symp,2)=norb(symp)
       mapddp1(symp,3)=symp
       mapddp1(symp,4)=1
       mapddp1(symp,5)=1
       mapddp1(symp,6)=1
       mapidp1(symp,1,1)=symp
       posst=posst+norb(symp)
       end do
c
       possdp20=posst
       mapddp2(0,1)=5
       mapddp2(0,2)=0
       mapddp2(0,3)=0
       mapddp2(0,4)=0
       mapddp2(0,5)=nsym
       mapddp2(0,6)=0
c
       do symp=1,nsym
       mapddp2(symp,1)=posst
       mapddp2(symp,2)=norb(symp)
       mapddp2(symp,3)=symp
       mapddp2(symp,4)=1
       mapddp2(symp,5)=1
       mapddp2(symp,6)=1
       mapidp2(symp,1,1)=symp
       posst=posst+norb(symp)
       end do
c
c
c1.1  maps for T1
c     T11 - t1oaa(a,i)
c     T12 - t1obb(a,i)
c     T13 - t1naa(a,i)
c     T14 - t1nbb(a,i)
c
       posst110=posst
       call grc0 (2,0,3,1,0,0,1,
     & posst110,posst,mapdt11,mapit11)
       posst120=posst
       call grc0 (2,0,4,2,0,0,1,
     & posst120,posst,mapdt12,mapit12)
       posst130=posst
       call grc0 (2,0,3,1,0,0,1,
     & posst130,posst,mapdt13,mapit13)
       posst140=posst
       call grc0 (2,0,4,2,0,0,1,
     & posst140,posst,mapdt14,mapit14)
c
c
c1.2  maps for F1
c     F11 - FI(a,e)aa
c     F12 - FI(a,e)bb
c
       possf110=posst
       call grc0 (2,0,3,3,0,0,1,
     & possf110,posst,mapdf11,mapif11)
       possf120=posst
       call grc0 (2,0,4,4,0,0,1,
     & possf120,posst,mapdf12,mapif12)
c
c
c1.3  maps for F2
c     F21 - FII(m,i)aa
c     F22 - FII(m,i)bb
c
       possf210=posst
       call grc0 (2,0,1,1,0,0,1,
     & possf210,posst,mapdf21,mapif21)
       possf220=posst
       call grc0 (2,0,2,2,0,0,1,
     & possf220,posst,mapdf22,mapif22)
c
c
c1.4  maps for F3
c     F31 - FIII(e,m)aa
c     F32 - FIII(e,m)bb
c
       possf310=posst
       call grc0 (2,0,3,1,0,0,1,
     & possf310,posst,mapdf31,mapif31)
       possf320=posst
       call grc0 (2,0,4,2,0,0,1,
     & possf320,posst,mapdf32,mapif32)
c
c
c1.5  maps for FK
c     FK1 - f(a,b)aa
c     FK2 - f(a,b)bb
c     FK3 - f(a,i)aa
c     FK4 - f(a,i)bb
c     FK5 - f(i,j)aa
c     FK6 - f(i,j)bb
c
       possfk10=posst
       call grc0 (2,0,3,3,0,0,1,
     & possfk10,posst,mapdfk1,mapifk1)
       possfk20=posst
       call grc0 (2,0,4,4,0,0,1,
     & possfk20,posst,mapdfk2,mapifk2)
       possfk30=posst
       call grc0 (2,0,3,1,0,0,1,
     & possfk30,posst,mapdfk3,mapifk3)
       possfk40=posst
       call grc0 (2,0,4,2,0,0,1,
     & possfk40,posst,mapdfk4,mapifk4)
       possfk50=posst
       call grc0 (2,0,1,1,0,0,1,
     & possfk50,posst,mapdfk5,mapifk5)
       possfk60=posst
       call grc0 (2,0,2,2,0,0,1,
     & possfk60,posst,mapdfk6,mapifk6)
c
c
c1.6  maps for T2
c     T21 - t2n(ab,ij)aaaa
c     T22 - t2n(ab,ij)bbbb
c     T33 - t2n(a,b,i,j)abab
c
       posst210=posst
       call grc0 (4,4,3,3,1,1,1,
     & posst210,posst,mapdt21,mapit21)
       posst220=posst
       call grc0 (4,4,4,4,2,2,1,
     & posst220,posst,mapdt22,mapit22)
       posst230=posst
       call grc0 (4,0,3,4,1,2,1,
     & posst230,posst,mapdt23,mapit23)
c
c
c1.7  maps for W0
c     W01 - <mn||ij>aaaa
c     W02 - <mn||ij>bbbb
c     W03 - <mn||ij>abab
c
       possw010=posst
       call grc0 (4,4,1,1,1,1,1,
     & possw010,posst,mapdw01,mapiw01)
       possw020=posst
       call grc0 (4,4,2,2,2,2,1,
     & possw020,posst,mapdw02,mapiw02)
       possw030=posst
       call grc0 (4,0,1,2,1,2,1,
     & possw030,posst,mapdw03,mapiw03)
c
c
c1.8  maps for W1
c     W11 - <ie||mn>aaaa
c     W12 - <ie||mn>bbbb
c     W13 - <ie||mn>abab
c     W14 - <ie||mn>baab
c
       possw110=posst
       call grc0 (4,3,1,3,1,1,1,
     & possw110,posst,mapdw11,mapiw11)
       possw120=posst
       call grc0 (4,3,2,4,2,2,1,
     & possw120,posst,mapdw12,mapiw12)
       possw130=posst
       call grc0 (4,0,1,4,1,2,1,
     & possw130,posst,mapdw13,mapiw13)
       possw140=posst
       call grc0 (4,0,2,3,1,2,1,
     & possw140,posst,mapdw14,mapiw14)
c
c
c2    for help files mapps are irrelevant,
c     here only estimation of maximal length is done to
c     define poss0 of help files
c     we have:
c     four V files - of vvoo type
c     four M files - of vvo  type
c     four H files - of voo  type
c     one  N file  - of nn   type
c
c2.*  def max{noa}, max{norb} ,max{nvb}, maxov(isym)=max{noa(isym),nvb(isym)}
c
       maxnoa=noa(1)
       maxnvb=nvb(1)
       maxnorb=norb(1)
       do 100 symp=1,nsym
       if (noa(symp).gt.maxnoa) then
       maxnoa=noa(symp)
       end if
       if (norb(symp).gt.maxnorb) then
       maxnorb=norb(symp)
       end if
       if (nvb(symp).gt.maxnvb) then
       maxnvb=nvb(symp)
       end if
       if (nvb(symp).gt.noa(symp)) then
       maxov(symp)=nvb(symp)
       else
       maxov(symp)=noa(symp)
       end if
 100    continue
c
c2.*  def lengths of V,M,H and N fils
c
       lengthv=0
       lengthm=0
       lengthh=0
       lengthn=0
c
       do 200 symp=1,nsym
c     symq is not known for N file
c     instead of norb(symr) maxnorb will be used so that reallength<=length
       lengthn=lengthn+norb(symp)*maxnorb
       do 201 symq=1,nsym
       sympq=mmul(symp,symq)
c     symr is not known for M and H files
c     instead of noa(symr) maxnoa will be used so that reallength<=length
       lengthm=lengthm+maxov(symp)*maxov(symq)*maxnoa
       lengthh=lengthh+maxov(symp)*noa(symq)*maxnoa
       do 202 symr=1,nsym
       syms=mmul(sympq,symr)
       lengthv=lengthv+maxov(symp)*maxov(symq)*noa(symr)*noa(syms)
 202    continue
 201    continue
 200    continue
c
c2.1  V - files
c
       possv10=posst
       posst=posst+lengthv
       possv20=posst
       posst=posst+lengthv
       possv30=posst
       posst=posst+lengthv
       possv40=posst
       posst=posst+lengthv
       lenv=lengthv
c
c2.2  M - files
c
       possm10=posst
       posst=posst+lengthm
       possm20=posst
       posst=posst+lengthm
       possm30=posst
       posst=posst+lengthm
       possm40=posst
       posst=posst+lengthm
c
c2.3  H - files
c
       possh10=posst
       posst=posst+lengthh
       possh20=posst
       posst=posst+lengthh
       possh30=posst
       posst=posst+lengthh
       possh40=posst
       posst=posst+lengthh
c
c2.4  N,P - files
c
       possn0=posst
       posst=posst+lengthn
       possp0=posst
       posst=posst+lengthn
       lenn=lengthn
c
c2.5  dedlare space for help matrix D in for matrix multiplication C=AT*B if
c     mchntyp=2
c
       if (mchntyp.eq.2) then
       possd0=posst
       if (maxnoa.le.maxnvb) then
       posst=posst+maxnoa*maxnoa*maxnvb*maxnvb
       else
       posst=posst+maxnoa*maxnoa*maxnoa*maxnoa
       end if
       end if
c
c2.6   def size of Work space
       length=posst-1
c
       return
       end
