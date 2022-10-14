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
       subroutine t3initfiles (length)
!
!     this routine distribute work space WRK for required files
!     for fix mediates it defines also mapd and mapi, for help mediates
!     it estimates their length and distribute WRK (i.e. def poss0 parameters)
!
!     length - overal requirements of work space (O)
!
!     !N.B. This routine cannot run with +OP2 level
!
       integer length
!
#include "t31.fh"
#include "t32.fh"
!
!     help variable
!
       integer posst,symp,symq,symr
       integer sizew,sizem,sizeh,sizen,sizel,sizer
       integer maxnoa,maxnvb,maxnorb
       integer nhelp1,nhelp2
!
!1    maps and possitions for fix mediated
!
!1.0  maps for DP - diagonal part
!     N.B. DP has one degree of freedom, while other 1 index has none
!     DP1 - dp(p)a
!     DP2 - dp(p)b
!
       do symp=1,nsym
       do symq=1,nsym
       do symr=1,nsym
       mapidp1(symp,symq,symr)=0
       mapidp2(symp,symq,symr)=0
       end do
       end do
       end do
!
       posst=1
!
       possdp10=posst
       mapddp1(0,1)=5
       mapddp1(0,2)=0
       mapddp1(0,3)=0
       mapddp1(0,4)=0
       mapddp1(0,5)=nsym
       mapddp1(0,6)=0
!
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
!
       possdp20=posst
       mapddp2(0,1)=5
       mapddp2(0,2)=0
       mapddp2(0,3)=0
       mapddp2(0,4)=0
       mapddp2(0,5)=nsym
       mapddp2(0,6)=0
!
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
!
!
!1.1  maps for T1
!     T11 - t1oaa(a,i)
!     T12 - t1obb(a,i)
!
       posst110=posst
       call cct3_grc0 (2,0,3,1,0,0,1,                                   &
     & posst110,posst,mapdt11,mapit11)
       posst120=posst
       call cct3_grc0 (2,0,4,2,0,0,1,                                   &
     & posst120,posst,mapdt12,mapit12)
!
!
!1.5  maps for FK
!     FK1 - f(a,b)aa
!     FK2 - f(a,b)bb
!     FK3 - f(a,i)aa
!     FK4 - f(a,i)bb
!     FK5 - f(i,j)aa
!     FK6 - f(i,j)bb
!
       possfk10=posst
       call cct3_grc0 (2,0,3,3,0,0,1,                                   &
     & possfk10,posst,mapdfk1,mapifk1)
       possfk20=posst
       call cct3_grc0 (2,0,4,4,0,0,1,                                   &
     & possfk20,posst,mapdfk2,mapifk2)
       possfk30=posst
       call cct3_grc0 (2,0,3,1,0,0,1,                                   &
     & possfk30,posst,mapdfk3,mapifk3)
       possfk40=posst
       call cct3_grc0 (2,0,4,2,0,0,1,                                   &
     & possfk40,posst,mapdfk4,mapifk4)
       possfk50=posst
       call cct3_grc0 (2,0,1,1,0,0,1,                                   &
     & possfk50,posst,mapdfk5,mapifk5)
       possfk60=posst
       call cct3_grc0 (2,0,2,2,0,0,1,                                   &
     & possfk60,posst,mapdfk6,mapifk6)
!
!
!1.6  maps for T2
!     T21 - t2o(ab,ij)aaaa
!     T22 - t2o(ab,ij)bbbb
!     T23 - t2o(a,b,i,j)abab
!
       posst210=posst
       call cct3_grc0 (4,4,3,3,1,1,1,                                   &
     & posst210,posst,mapdt21,mapit21)
       posst220=posst
       call cct3_grc0 (4,4,4,4,2,2,1,                                   &
     & posst220,posst,mapdt22,mapit22)
       posst230=posst
       call cct3_grc0 (4,0,3,4,1,2,1,                                   &
     & posst230,posst,mapdt23,mapit23)
!
!
!1.8  maps for W1
!     W11 - <ie||mn>aaaa
!     W12 - <ie||mn>bbbb
!     W13 - <ie||mn>abab
!     W14 - <ie||mn>baab
!
       possw110=posst
       call cct3_grc0 (4,3,1,3,1,1,1,                                   &
     & possw110,posst,mapdw11,mapiw11)
       possw120=posst
       call cct3_grc0 (4,3,2,4,2,2,1,                                   &
     & possw120,posst,mapdw12,mapiw12)
       possw130=posst
       call cct3_grc0 (4,0,1,4,1,2,1,                                   &
     & possw130,posst,mapdw13,mapiw13)
       possw140=posst
       call cct3_grc0 (4,0,2,3,1,2,1,                                   &
     & possw140,posst,mapdw14,mapiw14)
!
!
!1.9  maps for W2
!     W21 - <ab||ij>aaaa
!     W22 - <ab||ij>bbbb
!     W23 - <a,b|i,j>abab
!
       possw210=posst
       call cct3_grc0 (4,4,3,3,1,1,1,                                   &
     & possw210,posst,mapdw21,mapiw21)
       possw220=posst
       call cct3_grc0 (4,4,4,4,2,2,1,                                   &
     & possw220,posst,mapdw22,mapiw22)
       possw230=posst
       call cct3_grc0 (4,0,3,4,1,2,1,                                   &
     & possw230,posst,mapdw23,mapiw23)
!
!
!2    for help files mapps are irrelevant,
!     here only estimation of maximal length is done to
!     define poss0 of help files
!     we have:
!     2  W,V files - of vv2 type
!     2    L files - of vvv (vvo) type
!     3    R files - of vv2+ type
!     3    M files - of vv (vo)  type
!     3    H files - of v (o)  type
!     2  N,P files - of nn   type
!
!
       possw0=posst
!
!2.*  find maxsize of W,L,M,H
       sizew=0
       sizel=0
       sizer=0
       sizem=0
       sizeh=0
!
       do 50 nhelp1=1,nsym
!
!     W,V files
       call cct3_t3grc0 (3,2,4,4,4,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizew) then
       sizew=nhelp2
       end if
!
!     L files
       call cct3_t3grc0 (3,0,4,4,4,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizel) then
       sizel=nhelp2
       end if
       call cct3_t3grc0 (3,0,1,4,4,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizel) then
       sizel=nhelp2
       end if
!
!     R files
       call cct3_t3grc0 (3,8,4,4,4,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizer) then
       sizer=nhelp2
       end if
!
!     M files
       call cct3_t3grc0 (2,0,4,4,0,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizem) then
       sizem=nhelp2
       end if
       call cct3_t3grc0 (2,0,1,4,0,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizem) then
       sizem=nhelp2
       end if
!
!     H files
       call cct3_t3grc0 (1,0,4,0,0,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizeh) then
       sizeh=nhelp2
       end if
       call cct3_t3grc0 (1,0,1,0,0,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizeh) then
       sizeh=nhelp2
       end if
!
 50     continue
!
!
!2.*  def max{noa}, max{norb} ,max{nvb}
!
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
 100    continue
!
!
!2.*  def lengths of N fils
!
       sizen=0
!
       do 200 symp=1,nsym
!     symq is not known for N file
!     instead of norb(symr) maxnorb will be used so that reallength<=length
       sizen=sizen+norb(symp)*maxnorb
 200    continue
!
!2.1  W,V - files
!
!     possw0 is defined
       posst=possw0+sizew
       possv0=posst
       posst=posst+sizew
!
!2.2  L - files
!
       possl10=posst
       posst=posst+sizel
       possl20=posst
       posst=posst+sizel
!
!2.3  R - files
!
       possr10=posst
       posst=posst+sizer
       possr20=posst
       posst=posst+sizer
       possr30=posst
       posst=posst+sizer
!
!2.4  M - files
!
       possm10=posst
       posst=posst+sizem
       possm20=posst
       posst=posst+sizem
       possm30=posst
       posst=posst+sizem
!
!2.5  H - files
!
       possh10=posst
       posst=posst+sizeh
       possh20=posst
       posst=posst+sizeh
       possh30=posst
       posst=posst+sizeh
!
!2.6  N,P - files
!
       possn0=posst
       posst=posst+sizen
       possp0=posst
       posst=posst+sizen
!
!2.7  dedlare space for help matrix D in for matrix multiplication C=AT*B if
!     mchntyp=2
!
       if (mchntyp.eq.2) then
       possd0=posst
       if (maxnoa.le.maxnvb) then
       posst=posst+maxnoa*maxnoa*maxnvb*maxnvb
       else
       posst=posst+maxnoa*maxnoa*maxnoa*maxnoa
       end if
       end if
!
       length=posst-1
!
       return
       end
