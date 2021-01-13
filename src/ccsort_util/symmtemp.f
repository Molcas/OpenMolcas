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
c     this routines are installed here temrorarry
c     for final version they will be part of SYMM routines
c
c      wrtmediate
c      wrtmap
c      rea
c      wri
c      grc0
c      ccsort_mv0zero
c      dawrtmediate
c      dawrtmap
c      dawri
c      darea
c
c     ----------------------------
c
       subroutine ccsort_wrtmediate (wrk,wrksize,
     & lun,mapd,mapi,rc)
c
c     this routine write required mediate to opened unformatted file
c     with number lun
c     it also store mapd and mapi of the given mediade
c
c     lun   - Logical unit number of file, where mediate will be stored (Input)
c     mapd  - direct map matrix corresponding to given mediate (Input)
c     mapi  - inverse map matrix corresponding to given mediate (Input)
c     rc    - return (error) code (Output)
c
c     N.B.
c     all mediates are storred as follows
c     1 - mapd, mapi
c     2 - one record with complete mediate
c
#include "wrk.fh"
       integer lun,rc
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c     help variables
c
       integer im,length,poss0
c
       rc=0
c
c1    write mapd
c
       write (lun) mapd,mapi
c
c2    calculate overall length
c
       length=0
c
       do 100 im=1,mapd(0,5)
       length=length+mapd(im,2)
 100    continue
c
c     write mediate in one block
c
       if (length.eq.0) then
c     RC=1 : there is nothing to write, length of mediate is 0
       rc=1
       return
       end if
c
       poss0=mapd(1,1)
       call ccsort_wri (lun,length,wrk(poss0))
c
       return
       end
c
c     ----------------------------
c
       subroutine ccsort_wrtmap (lun,mapd,mapi,rc)
c
c     this routine write required mapd and mapi to opened unformatted file
c     with number lun
c
c     lun   - Logical unit number of file, where mediate will be stored (Input)
c     mapd  - direct map matrix corresponding to given mediate (Input)
c     mapi  - inverse map matrix corresponding to given mediate (Input)
c     rc    - return (error) code (Output)
c
       integer lun,rc
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
       rc=0
c
c1    write mapd
c
       write (lun) mapd,mapi
c
       return
       end
c
c     ----------------------------
c
       subroutine ccsort_rea (lun,length,vector)
c
c     this routine read length-R8 numbers from opened unformatted file
c     with number lun form the given possition as one record
c
c     lun    - Logical unit number of file, where mediate is stored (Input)
c     length - # of R8 numbers to be read  (Input)
c     vector - space, where numbers are stored after reading  (Output)

c
       integer lun,length,i
       real*8 vector(1:length)
c
       read (lun) (vector(i),i=1,length)
c
       return
       end
c
c     ----------------------------
c
       subroutine ccsort_wri (lun,length,vector)
c
c     this routine write length-R8 numbers to opened unformatted file
c     with number lun at the given possition as one record
c
c     lun    - Logical unit number of file, where mediate will be stored (Input)
c     length - # of R8 numbers to be written  (Input)
c     vector - space, where numbers are stored  (Input)

c
       integer lun,length
       real*8 vector(1:length)
c
       write (lun) vector
c
       return
       end
c
c     ----------------------------
c
       subroutine ccsort_grc0 (nind,typ,typp,typq,typr,typs,stot,
     & poss0,posst,mapd,mapi)
c
c     this routine defines mapd and mapi for given intermediat
c
#include "ccsort.fh"
       integer nind,typ,typp,typq,typr,typs,stot,poss0,posst
c
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c     integer mul(1:8,1:8)
       integer dimm(1:4,1:8)
c
c     help variables
c
       integer sp,sq,sr,ss,spq,spqr
       integer nsymq,nsymr
       integer poss,i,nhelp1,nhelp2,nhelp3,nhelp4
c
c*    !!!!!!!! def dimm to je tu len terazky !!!!
c
       do i=1,nsym
       dimm(1,i)=noa(i)
       dimm(2,i)=nob(i)
       dimm(3,i)=nva(i)
       dimm(4,i)=nvb(i)
       end do
c
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
c      To get rid of tiring compiler warning
       poss=0
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
c     def length
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
c     def length
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
       i=1
       poss=poss0
c
       do 200 sp=1,nsym
       if (typ.eq.1) then
       nsymq=sp
       else
       nsymq=nsym
       end if
c
       do 201 sq=1,nsymq
       spq=mul(sp,sq)
c
       sr=mul(stot,spq)
       if ((typ.eq.2).and.(sq.lt.sr)) then
c     Meggie out
       goto 201
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
c     def length
       if ((typ.eq.1).and.(sp.eq.sq)) then
       mapd(i,2)=nhelp1*(nhelp1-1)*nhelp3/2
       else if ((typ.eq.2).and.(sq.eq.sr)) then
       mapd(i,2)=nhelp1*nhelp2*(nhelp2-1)/2
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
 201    continue
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
       do 301 sq=1,nsymq
       spq=mul(sp,sq)
       if (typ.eq.2) then
       nsymr=sq
       else
       nsymr=nsym
       end if
c
       do 302 sr=1,nsymr
       spqr=mul(spq,sr)
c
       ss=mul(stot,spqr)
       if (((typ.eq.3).or.(typ.eq.4)).and.(sr.lt.ss)) then
c     Meggie out
       goto 302
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
c     def length
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
 302    continue
 301    continue
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
       SUBROUTINE ccsort_mv0zero
     & (DD,LENGTH,MAT)
C
       INTEGER           DD
       INTEGER           LENGTH
       real*8  MAT(1:DD)
       INTEGER           INIT
       real*8  ZERO
C
       DATA              ZERO/0.0D+00/
C
C     ...loop over all elements
C
       DO 10 INIT=1,LENGTH
       MAT(INIT) = ZERO
 10     CONTINUE
C
       RETURN
       END
c
c     -----------------------------
c
       subroutine dawrtmediate (wrk,wrksize,
     & lun,mapd,mapi,rc)
c
c     this routine write required mediate to opened unformatted file
c     with number lun
c     it also store mapd and mapi of the given mediade
c
c     lun   - Logical unit number of file, where mediate will be stored (Input)
c     mapd  - direct map matrix corresponding to given mediate (Input)
c     mapi  - inverse map matrix corresponding to given mediate (Input)
c     rc    - return (error) code (Output)
c
c     N.B.
c     all mediates are storred as follows
c     1 - mapd, mapi
c     2 - one record with complete mediate
c
#include "wrk.fh"
       integer lun,rc
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c     help variables
c
       integer im,length,poss0
c
       rc=0
c
c1    write mapd
c
      call dawrtmap (lun,mapd,mapi,rc)
c
c2    calculate overall length
c
       length=0
c
       do 100 im=1,mapd(0,5)
       length=length+mapd(im,2)
 100    continue
c
c     write mediate in one block
c
       if (length.eq.0) then
c     RC=1 : there is nothing to write, length of mediate is 0
       rc=1
       return
       end if
c
       poss0=mapd(1,1)
       call dawri (lun,length,wrk(poss0))
c
       return
       end
c
c     ----------------------------
c
       subroutine dawrtmap (lun,mapd,mapi,rc)
c
c     this routine write required mapd and mapi to opened unformatted file
c     with number lun
c
c     lun   - Logical unit number of file, where mediate will be stored (Input)
c     mapd  - direct map matrix corresponding to given mediate (Input)
c     mapi  - inverse map matrix corresponding to given mediate (Input)
c     rc    - return (error) code (Output)
c
#include "files_ccsd.fh"
#include "reorg.fh"

#include "SysDef.fh"
c
       integer lun,rc
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
       rc=0
c
c1    write mapd
c
       if (iokey.eq.1) then
c      Fortran IO
       write (lun) mapd,mapi
c
       else
c      MOLCAS IO
       call idafile (lun,1,mapd,3078,daddr(lun))
       call idafile (lun,1,mapi,512,daddr(lun))
       end if
c
       return
       end
c
c     ----------------------------
c
       subroutine dawri (lun,length,vector)
c
c     this routine write length-R8 numbers to opened unformatted file
c     with number lun at the given possition as one record
c
c     lun    - Logical unit number of file, where mediate will be stored (Input)
c     length - # of R8 numbers to be written  (Input)
c     vector - space, where numbers are stored  (Input)

c
#include "files_ccsd.fh"
#include "reorg.fh"

#include "SysDef.fh"
c
       integer lun,length
       real*8 vector(1:length)
c
       if (iokey.eq.1) then
c      Fortran IO
       write (lun) vector
c
       else
c      MOLCAS IO
       call ddafile (lun,1,vector,length,daddr(lun))
       end if
c
       return
       end
c
c     ----------------------------
c
