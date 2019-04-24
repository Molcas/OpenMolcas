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

c     this file contains routines, required for making
c     integrals for noniterative T3 calculatins:
c     t3intpck1
c     t3intpck2
c     t3reorg
c     ccsort_t3grc0 - the copy of original routine form T3
c     DefT3par
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine t3intpck1 (vint,r,dimv1,dimv2,dimv3,dima,dimbc,
     & symq,symr,syms,nob,nvb)
c
c     this routine pack integral block symi,symq,symr,syms
c     R_i(a,bc) = V_i(b,a,c)
c     for symq(b)=syms(c)
c     and write R block onto proper place of opened t3nam file - lunt3
c
c     vint  - integrals for given symetries for given i (I)
c     r     - final R_i matrix (O)
c     dimv1 - 1-st. dimension of V (I)
c     dimv2 - 2-nd. dimension of V (I)
c     dimv3 - 3-rd. dimension of V (I)
c     dima  - dimension of a in R (I)
c     dimbc - dimension of bc in R (I)
c     symq  - symmetry of q (b) (I)
c     symr  - symmetry of r (a) (I)
c     syms  - symmetry of s (c) (I)
c     nob   - number of beta occupied in each irrep (I)
c     nvb   - number of beta virtuals in each irrep (I)
c
       implicit none
#include "reorg.fh"
#include "files_ccsd.fh"
       integer symq,symr,syms
       integer dimv1,dimv2,dimv3,dima,dimbc
       integer nob(1:8)
       integer nvb(1:8)
       real*8 vint(1:dimv1,1:dimv2,1:dimv3)
       real*8 r(1:dima,1:dimbc)
c
c     help variables
c
       integer a,b,c,bc,adda,lenght,iaddr
c
c*    if there are no beta virtuals - goto write section
       if (nvb(symq)*nvb(symr)*nvb(syms).eq.0) then
       goto 200
       end if
c
c*    calc additional constant for a
       adda=nob(symr)
c
c*    do packing
c
       bc=0
       do 100 b=nob(symq)+1,nob(symq)+nvb(symq)
       do 100 c=nob(syms)+1,b
       bc=bc+1
       do 100 a=1,nvb(symr)
       r(a,bc)=vint(b,a+adda,c)
 100    continue
c
c*    write section
c
 200    lenght=dima*dimbc
       if (lenght.gt.0) then
       iaddr=daddr(lunt3)
       call ddafile (lunt3,1,r(1,1),lenght,iaddr)
       end if

       return
       end
c
c     -----------------------------
c
       subroutine t3intpck2 (vint,r,dimv1,dimv2,dimv3,dima,dimb,dimc,
     & symq,symr,syms,nob,nvb)
c
c     this routine pack integral block symi,symq,symr,syms
c     R_i(a,b,c) = V_i(b,a,c)
c     for symq(b)>syms(c)
c     and write R block onto proper place of opened t3nam file - lunt3
c
c     vint  - integrals for given symetries for given i (I)
c     r     - final R_i matrix (O)
c     dimv1 - 1-st. dimension of V (I)
c     dimv2 - 2-nd. dimension of V (I)
c     dimv3 - 3-rd. dimension of V (I)
c     dima  - dimension of a in R (I)
c     dimb  - dimension of b in R (I)
c     dimc  - dimension of c in R (I)
c     symq  - symmetry of q (b) (I)
c     symr  - symmetry of r (a) (I)
c     syms  - symmetry of s (c) (I)
c     nob   - number of beta occupied in each irrep (I)
c     nvb   - number of beta virtuals in each irrep (I)
c
       implicit none
#include "reorg.fh"
#include "files_ccsd.fh"
       integer symq,symr,syms
       integer dimv1,dimv2,dimv3,dima,dimb,dimc
       integer nob(1:8)
       integer nvb(1:8)
       real*8 vint(1:dimv1,1:dimv2,1:dimv3)
       real*8 r(1:dima,1:dimb,1:dimc)
c
c     help variables
c
       integer a,b,c,adda,addb,addc,lenght,iaddr
c
c*    if there are no beta virtuals - skip
       if (nvb(symq)*nvb(symr)*nvb(syms).eq.0) then
       return
       end if
c
c*    calc additional constants for a,b,c
       adda=nob(symr)
       addb=nob(symq)
       addc=nob(syms)

c
c*    do packing
c
       do 100 c=1,nvb(syms)
       do 100 b=1,nvb(symq)
       do 100 a=1,nvb(symr)
       r(a,b,c)=vint(b+addb,a+adda,c+addc)
 100    continue
c
c*    write section
c
        lenght=dima*dimb*dimc
       if (lenght.gt.0) then
       iaddr=daddr(lunt3)
       call ddafile (lunt3,1,r(1,1,1),lenght,iaddr)
       end if
c
       return
       end
c
c     -----------------------------
c
       subroutine t3reorg (wrk,wrksize,
     & noa,nsym)
c
c     this routine do final reorganization of t3nam file
c     and produce final form of this file
c     as it will be required in T3 and close t3nam file
c
c     noa   - array with occupation numbers
c     nsym  - actual number of irreps
c
       implicit none
#include "wrk.fh"
#include "reorg.fh"
#include "files_ccsd.fh"
       integer noa(1:8)
       integer nsym
c
c     help variables
c
       integer lenght,iri,possri
       integer posst
       integer symi,i,iaddr,iindex,iPossPack
c
c*        def iPossPack
c           iPossPack - possition of (maps+Ri) set in packed
c                        (i.e. final) of T3nam file
        iPossPack=T3IntPoss(1)
c
        iindex=0
        do symi=1,nsym
c
c0      get map's of R_i(a,bc)
        call ccsort_t3grc0
     c       (3,8,4,4,4,0,symi,possri0,posst,mapdri,mapiri)
c
        do i=1,noa(symi)
        iindex=iindex+1
c
c1        reconstruct R_i(a,bc) per blocks as in is
c         actually written in t3man file
          do iri=1,mapdri(0,5)
c
c1.1        iind address of this R_i block in t3nam file
            iaddr=T3IntPoss(iindex)+T3Off(iri,symi)
c
c1.2        def possition of of this block in R1
            possri=mapdri(iri,1)
c
c1.3        read integrals into proper possition
            lenght=mapdri(iri,2)
            if (lenght.gt.0) then
            call ddafile (lunt3,2,wrk(possri),lenght,iaddr)
            end if
c
          end do
c
c2          write into t3nam file in packed form
c            1) mapdri, mapiri
c           2) R_i
c2.1          def final (packed) address for i-th set (maps+Ri)
          T3intPoss(iindex)=iPossPack
          iaddr=T3intPoss(iindex)
c
c2.2      write maps
          call idafile (lunt3,1,mapdri,3078,iaddr)
          call idafile (lunt3,1,mapiri,512,iaddr)
c
c2.3      def actual length of Ri
          lenght=0
          do iri=1,mapdri(0,5)
          lenght=lenght+mapdri(iri,2)
          end do
c         lenght=mapdri(iri,1)+mapdri(iri,2)-mapdri(1,1)
c
c2.4          write Ri as one block
          call ddafile (lunt3,1,wrk(possri0),lenght,iaddr)
c
c2.5          save updated address as a new packed (final) possition
c         for next i
          iPossPack=iaddr
c
        end do
        end do
c
c3        store new packed (final) addreses T3IntPoss in t3nam file
c       (at the beggining)
        iaddr=0
        call idafile (lunt3,1,T3IntPoss,mbas,iaddr)
c
c4        close t3nam file
        call daclos (lunt3)
c
       return
       end
c
c     -----------------------------
c
       subroutine ccsort_t3grc0 (nind,typ,typp,typq,typr,typs,stot,
     & poss0,posst,mapd,mapi)
c
c     N.B. This routine is in principle copy of those in T3,
c     but some changes was done:
c     1) mmul is substituted by mul
c     2) dimm is added, since using of ccsd1.com is inpossible
c     3) ccsd1.com is replaced by ccsort.fh
c
c     nind   - number of indexes (I)
c     typ    - typ of mediate (I)
c     typp   - typ of index p (I)
c     typq   - typ of index q (I)
c     typr   - typ of index r (I)
c     typs   - typ of index s (I)
c     stot   - overall symetry of the mediate (I)
c     poss0  - initil possition of mediate (I)
c     posst  - final possition of the mediate (O)
c     mapd   - direct map of the mediate (O)
c     mapi   - inverse map of the mediate (O)
c
c     this routine defines mapd and mapi for given intermediat
c     it can done exactly the same maps like grc0 in CCSD
c     plus additional types of mediates are introduced:
c     type    meaning
c     5     p>q>r,s ; also p>q>r
c     6     p,q>r>s
c     7     p>=q,r,s ; also p>=q,r; p>=q
c     8          p,q>=r,s ; also p,q>=s
c     9     p,q,q>=s
c     10     p>=q,r>=s
c     11     p>=q>=r,s ; also p>=q>=r
c     12     p,q>=r>=s
c
c     currently, these new types are implemented only for nind=3
c
c     $N.B. (this routine cannot run with +OP2)
c     N.B. this routine do not test stupidities
c
c
       integer nind,typ,typp,typq,typr,typs,stot,poss0,posst
c
c@    include 'ccsd1.com'
#include "ccsort.fh"
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
       integer dimm(1:5,1:8)
c
c     help variables
c
       integer sp,sq,sr,ss,spq,spqr
       integer nsymq,nsymr
       integer poss,i,nhelp1,nhelp2,nhelp3,nhelp4
       integer rsk1,rsk2
c
c@    !!!!!!!! def dimm to je tu len terazky, lebo nemozeme pouzivat ccsd1.com !!!!
c
c      Tutok musim cosi inicializovat
       ss=0
       poss=0
       rsk1=0
       rsk2=0
       do i=1,nsym
       dimm(1,i)=noa(i)
       dimm(2,i)=nob(i)
       dimm(3,i)=nva(i)
       dimm(4,i)=nvb(i)
       dimm(5,i)=nva(i)+noa(i)
       end do
c
c@@
c     vanishing mapi files
c
       do nhelp1=1,nsym
       do nhelp2=1,nsym
       do nhelp3=1,nsym
       mapi(nhelp3,nhelp2,nhelp1)=0
       end do
       end do
       end do
c
       if (nind.eq.1) then
c
c     matrix A(p)
c
       i=1
       poss=poss0
       sp=mul(stot,1)
c
       nhelp1=dimm(typp,sp)
c
c     def mapi
       mapi(1,1,1)=i
c
c     def possition
       mapd(i,1)=poss
c
c     def lenght
       mapd(i,2)=nhelp1
c
c     def sym p,q
       mapd(i,3)=sp
       mapd(i,4)=0
       mapd(i,5)=0
       mapd(i,6)=0
c
       poss=poss+mapd(i,2)
       i=i+1
c
       else if (nind.eq.2) then
c
c     matrix A(p,q)
c
       i=1
       poss=poss0
c
       do 100 sp=1,nsym
c
       sq=mul(stot,sp)
       if ((typ.eq.1).and.(sp.lt.sq)) then
c     Meggie out
       goto 100
       end if
c
       nhelp1=dimm(typp,sp)
       nhelp2=dimm(typq,sq)
c
c     def mapi
       mapi(sp,1,1)=i
c
c     def possition
       mapd(i,1)=poss
c
c     def lenght
       if ((typ.eq.1).and.(sp.eq.sq)) then
       mapd(i,2)=nhelp1*(nhelp1-1)/2
       else
       mapd(i,2)=nhelp1*nhelp2
       end if
c
c     def sym p,q
       mapd(i,3)=sp
       mapd(i,4)=sq
       mapd(i,5)=0
       mapd(i,6)=0
c
       poss=poss+mapd(i,2)
       i=i+1
c
 100    continue
c
       else if (nind.eq.3) then
c
c     matrix A(p,q,r)
c
c     def reucion sumations keys : rsk1 for pq, rsk2 for qr
c
       if (typ.eq.0) then
       rsk1=0
       rsk2=0
       else if (typ.eq.1) then
       rsk1=1
       rsk2=0
       else if (typ.eq.2) then
       rsk1=0
       rsk2=1
       else if (typ.eq.5) then
       rsk1=1
       rsk2=1
       else if (typ.eq.7) then
       rsk1=1
       rsk2=0
       else if (typ.eq.8) then
       rsk1=0
       rsk2=1
       else if (typ.eq.11) then
       rsk1=1
       rsk2=1
       end if
c
       i=1
       poss=poss0
c
       do 200 sp=1,nsym
       if (rsk1.eq.1) then
       nsymq=sp
       else
       nsymq=nsym
       end if
c
       do 200 sq=1,nsymq
       spq=mul(sp,sq)
c
       sr=mul(stot,spq)
       if ((rsk2.eq.1).and.(sq.lt.sr)) then
c     Meggie out
       goto 200
       end if
c
       nhelp1=dimm(typp,sp)
       nhelp2=dimm(typq,sq)
       nhelp3=dimm(typr,sr)
c
c     def mapi
       mapi(sp,sq,1)=i
c
c     def possition
       mapd(i,1)=poss
c
c     def lenght
       if ((typ.eq.1).and.(sp.eq.sq)) then
       mapd(i,2)=nhelp1*(nhelp1-1)*nhelp3/2
       else if ((typ.eq.2).and.(sq.eq.sr)) then
       mapd(i,2)=nhelp1*nhelp2*(nhelp2-1)/2
       else if (typ.eq.5) then
       if (sp.eq.sr) then
       mapd(i,2)=nhelp1*(nhelp1-1)*(nhelp1-2)/6
       else if (sp.eq.sq) then
       mapd(i,2)=nhelp1*(nhelp1-1)*nhelp3/2
       else if (sq.eq.sr) then
       mapd(i,2)=nhelp1*nhelp2*(nhelp2-1)/2
       else
       mapd(i,2)=nhelp1*nhelp2*nhelp3
       end if
       else if ((typ.eq.7).and.(sp.eq.sq)) then
       mapd(i,2)=nhelp1*(nhelp1+1)*nhelp3/2
       else if ((typ.eq.8).and.(sq.eq.sr)) then
       mapd(i,2)=nhelp1*nhelp2*(nhelp2+1)/2
       else if (typ.eq.11) then
       if (sp.eq.ss) then
       mapd(i,2)=nhelp1*(nhelp1+1)*(nhelp1+2)/6
       else if (sp.eq.sq) then
       mapd(i,2)=nhelp1*(nhelp1+1)*nhelp3/2
       else if (sq.eq.sr) then
       mapd(i,2)=nhelp1*nhelp2*(nhelp2+1)/2
       else
       mapd(i,2)=nhelp1*nhelp2*nhelp3
       end if
       else
       mapd(i,2)=nhelp1*nhelp2*nhelp3
       end if
c
c     def sym p,q,r
       mapd(i,3)=sp
       mapd(i,4)=sq
       mapd(i,5)=sr
       mapd(i,6)=0
c
       poss=poss+mapd(i,2)
       i=i+1
c
 200    continue
c
       else if (nind.eq.4) then
c
c     matrix A(p,q,r,s)
c
       i=1
       poss=poss0
c
       do 300 sp=1,nsym
       if ((typ.eq.1).or.(typ.eq.4)) then
       nsymq=sp
       else
       nsymq=nsym
       end if
c
       do 300 sq=1,nsymq
       spq=mul(sp,sq)
       if (typ.eq.2) then
       nsymr=sq
       else
       nsymr=nsym
       end if
c
       do 300 sr=1,nsymr
       spqr=mul(spq,sr)
c
       ss=mul(stot,spqr)
       if (((typ.eq.3).or.(typ.eq.4)).and.(sr.lt.ss)) then
c     Meggie out
       goto 300
       end if
c
       nhelp1=dimm(typp,sp)
       nhelp2=dimm(typq,sq)
       nhelp3=dimm(typr,sr)
       nhelp4=dimm(typs,ss)
c
c     def mapi
       mapi(sp,sq,sr)=i
c
c     def possition
       mapd(i,1)=poss
c
c     def lenght
       if ((typ.eq.1).and.(sp.eq.sq)) then
       mapd(i,2)=nhelp1*(nhelp2-1)*nhelp3*nhelp4/2
       else if ((typ.eq.2).and.(sq.eq.sr)) then
       mapd(i,2)=nhelp1*nhelp2*(nhelp3-1)*nhelp4/2
       else if ((typ.eq.3).and.(sr.eq.ss)) then
       mapd(i,2)=nhelp1*nhelp2*nhelp3*(nhelp4-1)/2
       else if (typ.eq.4) then
       if ((sp.eq.sq).and.(sr.eq.ss)) then
       mapd(i,2)=nhelp1*(nhelp2-1)*nhelp3*(nhelp4-1)/4
       else if (sp.eq.sq) then
       mapd(i,2)=nhelp1*(nhelp2-1)*nhelp3*nhelp4/2
       else if (sr.eq.ss) then
       mapd(i,2)=nhelp1*nhelp2*nhelp3*(nhelp4-1)/2
       else
       mapd(i,2)=nhelp1*nhelp2*nhelp3*nhelp4
       end if
       else
       mapd(i,2)=nhelp1*nhelp2*nhelp3*nhelp4
       end if
c
c     def sym p,q,r,s
       mapd(i,3)=sp
       mapd(i,4)=sq
       mapd(i,5)=sr
       mapd(i,6)=ss
c
       poss=poss+mapd(i,2)
       i=i+1
c
 300    continue
c
       end if

c
       posst=poss
c
c     definition of other coll
c
       mapd(0,1)=typp
       mapd(0,2)=typq
       mapd(0,3)=typr
       mapd(0,4)=typs
       mapd(0,5)=i-1
       mapd(0,6)=typ
c
       return
       end
c
c     -----------------------------
c
        subroutine DefT3par (noa,nsym)
c
c        this routine do:
c        0) Open t3nam file, with lunt3 lun
c        define parameters, required for T3 integral handling, namely
c        1) def T3IndPoss(i)
c          address possitions for all occupied orbitals in t3nam file
c        2) def T3Off(ii,isym)
c           relative shifts of address for ii-th block of R_i(a,bc)
c          for each symmetry
c
c     noa   - array with occupation numbers
c     nsym  - actual number of irreps
c
        implicit none
#include "reorg.fh"
#include "files_ccsd.fh"
c
       integer noa(1:8)
       integer nsym
c
c        help variables
        integer iorb,ii,i,symi,length,posst
c
c
c0        open t3nam file
c        lunt3=1
        call daname (lunt3,t3nam)
c
c1        set address poiter to 0
        daddr(lunt3)=0
c
c2        first record in t3nam file is T3IntPoss
c       (emulate writing of T3IntPoss)
        call idafile (lunt3,0,[0],mbas,daddr(lunt3))
c
        iorb=0
c3        cycle over irreps
        do symi=1,nsym
c
c3.1      make mapd and mapi for  R_i(a,bc)
          call ccsort_t3grc0
     c         (3,8,4,4,4,0,symi,possri0,posst,mapdri,mapiri)
c
c3.2          cycle over occupied orbitals in symi
          do i=1,noa(symi)
c
c3.2.1      save initial addres for this orbital
            iorb=iorb+1
            T3IntPoss(iorb)=daddr(lunt3)
c
c3.2.2      emulate writing of mapd and mapp
            call idafile (lunt3,0,[0],513*6,daddr(lunt3))
            call idafile (lunt3,0,[0],8*8*8,daddr(lunt3))
c
c3.2.3      cycle over all blocks of R_i(a,bc), which will
c           be stored separately
            do ii=1,mapdri(0,5)
c
c3.2.3.1      def T3Off(ii,symi)
c              note, that iorb is always proper one, since only besides
c             first occ. orbital in given irrep T3Off is defined
              if (i.eq.1) then
                T3Off(ii,symi)=daddr(lunt3)-T3IntPoss(iorb)
              end if
c
c3.2.3.2      emmulate writing of each block
              length=mapdri(ii,2)
              call ddafile (lunt3,0,[0.0d0],length,daddr(lunt3))
c
            end do
          end do
        end do
c
        return
        end
