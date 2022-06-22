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
