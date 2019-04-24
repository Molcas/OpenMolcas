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
       SUBROUTINE mkadress (NOIPSB)

CFUE
CFUE   It is strictly forbidden to use any construct to fetch disk
CFUE   addresses which do not make use of subroutines provided in
CFUE   the MOLCAS utilities as otherwise transparency and portability
CFUE   is lost.
CFUE

C
c     this routine prepair information arrays
c     for integral blocks per symmetry in primitive form
c     it defines:
c     noipsb(nijkl) - number of integrals per symmetry block
c     idispsb(nijkl)- initial addreses for exach symmetry block
c     typ(pa,qa,ra) - typ of symmetry block
c     types of (ij|kl):
c     1 - si=sk, si=sj, sk=sl
c     2 - si=sk, si=sj, sk>sl
c     3 - si=sk, si>sj, sk=sl
c     4 - si=sk, si>sj, sk>sl
c     5 - si>sk, si=sj, sk=sl
c     6 - si>sk, si=sj, sk>sl
c     7 - si>sk, si>sj, sk=sl
c     8 - si>sk, si>sj, sk>sl
c
c     idis(pa,qa,ra)- initial addreses for given symmetry block
c     np(pa,qa,ra)  - possition of p index in original block (ij|kl) (on tape)
c     nq(pa,qa,ra)  - possition of q index in original block (ij|kl) (on tape)
c     nr(pa,qa,ra)  - possition of r index in original block (ij|kl) (on tape)
c     ns(pa,qa,ra)  - possition of s index in original block (ij|kl) (on tape)
c
c     N.B. typ,idis,np,nq,nr,ns are trensfered thhrough common block /edpand2/
c
c
       IMPLICIT REAL*8(A-H,O-Z)

#include "tratoc.fh"
       integer INDMAX,MXFUNC
       PARAMETER (INDMAX=nTraBuf,MXFUNC=200)


#include "SysDef.fh"
#include "ccsort.fh"
#include "reorg.fh"

       integer NOIPSB(106)
       integer idispsb(106)
c
c     help variables
c
       integer sense
       integer p,q,r,s,pa,qa,ra,sa
       integer IND,INDT,ISPQRS,NINT,NSLM,idistemp,idishelp
       integer jlow,ilow,iup,jup,kup,lup,iold,jold,kold,lold
       integer norbp,nsi,nsj,nsk,nsl,nsij,nsijk

CFUE   - pick the start addresses of each symmetry allowed integral
CFUE     block from the tranformed two electron integral file
        idistemp=0
        Call iDaFile(LUINTM,0,iTraToc,nTraToc,idistemp)
CFUE   - the option 0 in the call to dafile does not make any I/O
CFUE     but just returns the updated disk address
CFUE   - at this point idistemp points to the first record

C
       IND=0
       INDT=0
CFUE   idistemp=1
CFUE   idisadd=150
C
       do 100 NSI=1,NSYM
       do 100 NSJ=1,NSYM
       do 100 NSK=1,NSYM
       typ(NSI,NSJ,NSK)=0
 100    continue
c

       if (fullprint.gt.0) then
       Write(6,'(6X,A)') 'Transformed integral blocks:'
       Write(6,'(6X,A)') '----------------------------'
       Write(6,*)
       Write(6,'(6X,A)')
     & 'block  symmetry      no. of        no. of '
       Write(6,'(6X,A)')
     & ' no.    spec.        orbitals     integrals'
       Write(6,'(6X,A)')
     & '-------------------------------------------'
       end if
c
       ISPQRS=0
       DO 300 NSI=1,NSYM
       DO 301 NSJ=1,NSI
       NSIJ=MUL(NSI,NSJ)
       DO 302 NSK=1,NSI
       NSIJK=MUL(NSK,NSIJ)
       NSLM=NSK
       IF(NSK.EQ.NSI) NSLM=NSJ
       DO 303 NSL=1,NSLM
       IF(NSIJK.NE.NSL) GO TO 303
       NORBP=NORB(NSI)*NORB(NSJ)*NORB(NSK)*NORB(NSL)
       IF(NORBP.EQ.0)GO TO 303
       ISPQRS=ISPQRS+1
c
c     def
c
c     redefine indices from Parr to Dirac notation <p,q|r,s> <-> (i,j|k,l)
c
       p=nsi
       q=nsk
       r=nsj
       s=nsl
c
c     def. sense
c
c     type (ij|kl)
c     1 - si=sk, si=sj, sk=sl
c     2 - si=sk, si=sj, sk>sl
c     3 - si=sk, si>sj, sk=sl
c     4 - si=sk, si>sj, sk>sl
c     5 - si>sk, si=sj, sk=sl
c     6 - si>sk, si=sj, sk>sl
c     7 - si>sk, si>sj, sk=sl
c     8 - si>sk, si>sj, sk>sl
c
       if (nsi.eq.nsk) then
       if (nsi.eq.nsj) then
       if (nsk.eq.nsl) then
       sense=1
       else
       sense=2
       end if
       else
       if (nsk.eq.nsl) then
       sense=3
       else
       sense=4
       end if
       end if
       else
       if (nsi.eq.nsj) then
       if (nsk.eq.nsl) then
       sense=5
       else
       sense=6
       end if
       else
       if (nsk.eq.nsl) then
       sense=7
       else
       sense=8
       end if
       end if
       end if
c
       if (nsijk.ne.nsl) then
       sense=0
       else if (NORB(NSI)*NORB(NSJ)*NORB(NSK)*NORB(NSL).eq.0) then
       sense=0
       end if
c
c
c1:   perm <pq|rs> -> <pq|rs>
c
       pa=p
       qa=q
       ra=r
       sa=s
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=1
       nq(pa,qa,ra)=3
       nr(pa,qa,ra)=2
       ns(pa,qa,ra)=4
c
c2:   perm <pq|rs> -> <rq|ps> 1-3
c
       pa=r
       qa=q
       ra=p
       sa=s
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=2
       nq(pa,qa,ra)=3
       nr(pa,qa,ra)=1
       ns(pa,qa,ra)=4
c
c3:   perm <pq|rs> -> <ps|rq> 2-4
c
       pa=p
       qa=s
       ra=r
       sa=q
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=1
       nq(pa,qa,ra)=4
       nr(pa,qa,ra)=2
       ns(pa,qa,ra)=3
c
c4:   perm <pq|rs> -> <pq|rs> 1-3,2-4
c
       pa=r
       qa=s
       ra=p
       sa=q
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=2
       nq(pa,qa,ra)=4
       nr(pa,qa,ra)=1
       ns(pa,qa,ra)=3
c
c5:   perm <pq|rs> -> <qp|sr> 1-2,3-4
c
       pa=q
       qa=p
       ra=s
       sa=r
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=3
       nq(pa,qa,ra)=1
       nr(pa,qa,ra)=4
       ns(pa,qa,ra)=2
c
c6:   perm <pq|rs> -> <qp|sr> 1-2,3-4 -> <sp|qr> 1-3
c
       pa=s
       qa=p
       ra=q
       sa=r
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=4
       nq(pa,qa,ra)=1
       nr(pa,qa,ra)=3
       ns(pa,qa,ra)=2
c
c7:   perm <pq|rs> -> <qp|sr> 1-2,3-4 -> <qr|sp> 2-4
c
       pa=q
       qa=r
       ra=s
       sa=p
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=3
       nq(pa,qa,ra)=2
       nr(pa,qa,ra)=4
       ns(pa,qa,ra)=1
c
c8:   perm <pq|rs> -> <qp|sr> 1-2,3-4 -> <sr|qp> 1-3,2-4
c
       pa=s
       qa=r
       ra=q
       sa=p
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=4
       nq(pa,qa,ra)=2
       nr(pa,qa,ra)=3
       ns(pa,qa,ra)=1
c
c
       idispsb(ispqrs)=idistemp
       idishelp=0
C
C     ******************************************************************
C
C     LOOP OVER THE USED ORBITALS OF EACH SYMMETRY BLOCK
C     THIS LOOPING IS COPIED FROM THE MANUAL OF THE 4-INDEX
C     TRANSFORMATION PROGRAM
C
C     NINT COUNTS INTEGRAL LABELS IN THE GIVEN SYMMETRY BLOCK
C
       NINT=0
       KUP=NORB(NSK)
       DO 401 KOLD=1,KUP
C
       LUP=NORB(NSL)
       IF (NSK.EQ.NSL) LUP=KOLD
       DO 402 LOLD=1,LUP
C
       ILOW=1
       IF (NSI.EQ.NSK) ILOW=KOLD
       IUP=NORB(NSI)
       DO 403 IOLD=ILOW,IUP
C
       JLOW=1
       IF (NSI.EQ.NSK.AND.IOLD.EQ.KOLD) JLOW=LOLD
       JUP=NORB(NSJ)
       IF (NSI.EQ.NSJ) JUP=IOLD
       DO 404 JOLD=JLOW,JUP
C
       IND=IND+1
       INDT=INDT+1
       NINT=NINT+1
       idishelp=idishelp+1
       if (idishelp.gt.INDMAX) then
CFUE     - all integrals in the reord were processed, hence
CFUE       update the disk address by a dummy I/O
         Call dDaFile(LUINTM,0,[0.0d0],INDMAX,idistemp)
CFUE     idistemp=idistemp+idisadd
         idishelp=1
       end if
C
       IF (IND.LT.INDMAX) GO TO 404
       IND=0
C
 404    CONTINUE
 403    CONTINUE
 402    CONTINUE
 401    CONTINUE
       if (idishelp.gt.0) then
CFUE     - all integrals in the reord were processed, hence
CFUE       update the disk address by a dummy I/O
         Call dDaFile(LUINTM,0,[0.0d0],INDMAX,idistemp)
CFUE   idistemp=idistemp+idisadd
       end if
C
C     WRITING THE LAST RECORD OF LABELS IN THE GIVEN SYMMETRY BLOCK
C     RECORDS ON LUPACK ARE FORMATTED TO 28KB LENGTH
C
       IF(IND.NE.0) THEN
       IND=0
       ENDIF
C
       NOIPSB(ISPQRS)=NINT
C
C     ******************************************************************
C
       if (fullprint.gt.0) then
       Write(6,'(6X,I5,2X,4I2,2X,4I4,2X,I8)')
     &       ISPQRS,
     &       NSI,NSJ,NSK,NSL,
     &       IUP,JUP,KUP,LUP,
     &       NINT
       end if
C
C     ******************************************************************
C
C
 303    CONTINUE
 302    CONTINUE
 301    CONTINUE
 300    CONTINUE
       if (fullprint.gt.0) then
       Write(6,'(6X,A)')
     & '-------------------------------------------'
       end if
C
       RETURN
       END
c
