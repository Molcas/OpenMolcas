************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2004,2007, Markus Reiher                               *
*               2004, Alexander Wolf                                   *
************************************************************************
      subroutine determine_factor2(length,term,iact,nbas,posu,post,poss,
     *                       vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,scr3,
     *                       A,C,nbasp,nbaso,dkhadr,adrmem,dkh_48,
     *                       adrnext)
c
c****************************************************************************
c
c   Multiply by factor stands at position 'iact' of 'term'.
c
c   Adjust value of 'iact', such that it points at beginning of next factor.
c
c   Note: There are no brackets [] occurring in term.
c
c   input: A left factor
c          C resulting product matrix
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Markus Reiher  (ETH Zurich)
c
c   version:  2.2.0
c
c   modified: 26.02.2007  by M. Reiher (ETH Zurich)
c        * can now read in right hand side factor from disk
c
c   first version: 12.04.2004  (Theoretical Chemistry, ETH Zurich)
c
c****************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
c
      integer length,iact,nbas,posu(maxunumber),post(maxsnumber),
     *        poss(maxsnumber),snumber,tnumber,unumber,nbasp,nbaso
      character*(*) term
c
      REAL*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp)
      REAL*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),A(nbas,nbas),
     *                 C(nbas,nbas)
c
      integer idum1,dkh_char2int
      REAL*8 alpha,beta
c
      integer dkh_48,adrmem,isize,iscr,iscr2,adrnext
      integer two,dkhadr(adrmem),adr
      parameter(two=2)
      logical flag
      flag=.true.
      isize=(nbas*(nbas+1)/2)
c
c------------------------------------------------------------------------
c
      alpha=1.0d0
      beta=0.0d0
c
      if (Out_Of_Core) then
       if (term(iact:iact).eq.'S') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         idum1=1000+idum1
         iact=iact+4
         adr=dkhadr(idum1)
         isize=nbas*nbas
       else if (term(iact:iact).eq.'T') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         idum1=2000+idum1
         iact=iact+4
         adr=dkhadr(idum1)
         isize=nbas*nbas
       else if (term(iact:iact).eq.'U') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         idum1=3000+idum1
         iact=iact+4
         adr=dkhadr(idum1)
         isize=nbas*nbas
       else if (term(iact:iact).eq.'V') then
         iact=iact+1
         adr=dkhadr(1)
       else if (term(iact:iact).eq.'N') then
         iact=iact+1
         adr=dkhadr(5)
       else if (term(iact:iact).eq.'D') then
         iact=iact+1
         adr=dkhadr(2)
       else if (term(iact:iact).eq.'Y') then
         iact=iact+1
         adr=dkhadr(6)
       else if (term(iact:iact).eq.'F') then
         iact=iact+1
         adr=dkhadr(7)
       else if (term(iact:iact).eq.'G') then
         iact=iact+1
         adr=dkhadr(8)
       else if (term(iact:iact).eq.'Z') then
         call mat_times_p2b(C,A,nbas,pp)
         iact=iact+1
         flag=.false.
       else if (term(iact:iact).eq.'Q') then
         call mat_div_p2b(C,A,nbas,pp)
         iact=iact+1
         flag=.false.
       else if (term(iact:iact).eq.'X') then
         iact=iact+1
         adr=dkhadr(3)
       else if (term(iact:iact).eq.'I') then
         iact=iact+1
         adr=dkhadr(9)
       else if (term(iact:iact).eq.'J') then
         iact=iact+1
         adr=dkhadr(4)
       else if (term(iact:iact).eq.'K') then
         iact=iact+1
         adr=dkhadr(10)
       else if (term(iact:iact).eq.'L') then
         iact=iact+1
         adr=dkhadr(11)
       else if (term(iact:iact).eq.'M') then
         iact=iact+1
         adr=dkhadr(12)
       end if
       if (flag) then
CMR it would be advantageous to have a triangular and a square
CMR   scratch matrix permanently available
         call GetMem('DetFac  ','ALLO','REAL',iscr,isize+4)
         call ddafile(dkh_48,two,work(iscr),isize,adr)
         if (isize.eq.nbas*nbas) then
           call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *               nbas,work(iscr),nbas,beta,C,nbas)
         else
           call GetMem('DetFac  ','ALLO','REAL',iscr2,nbas*nbas+4)
           call mat_sq_from_t(work(iscr2),nbas,work(iscr))
           call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *               nbas,work(iscr2),nbas,beta,C,nbas)
           call GetMem('DetFac  ','FREE','REAL',iscr2,nbas*nbas+4)
         end if
         call GetMem('DetFac  ','FREE','REAL',iscr,isize+4)
       end if
      else
c
      if (term(iact:iact).eq.'S') then
        idum1=dkh_char2int(3,term(iact+1:iact+3))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,scr1(1,1,poss(idum1)),nbas,beta,
     *              C,nbas)
        iact=iact+4
      else if (term(iact:iact).eq.'T') then
        idum1=dkh_char2int(3,term(iact+1:iact+3))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,scr2(1,1,post(idum1)),nbas,beta,
     *              C,nbas)
        iact=iact+4
      else if (term(iact:iact).eq.'U') then
        idum1=dkh_char2int(3,term(iact+1:iact+3))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,scr3(1,1,posu(idum1)),nbas,beta,
     *              C,nbas)
        iact=iact+4
      else if (term(iact:iact).eq.'V') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,vv,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'N') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,nn,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'D') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,dd,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'Y') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,yy,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'F') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,ff,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'G') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,gg,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'Z') then
        call mat_times_p2b(C,A,nbas,pp)
        iact=iact+1
      else if (term(iact:iact).eq.'Q') then
        call mat_div_p2b(C,A,nbas,pp)
        iact=iact+1
      else if (term(iact:iact).eq.'X') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,xx,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'I') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,ii,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'J') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,jj,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'K') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,kk,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'L') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,ll,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'M') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,mm,nbas,beta,C,nbas)
        iact=iact+1
      else
            write (stdout,1111)
1111        format (2X,'MIRACLE in determine_factor2(): could not',
     *        ' determine 2nd factor!'//2X,'STOP.',/2X)
            write(6,*) "term(iact:iact)=",term(iact:iact)
            CALL Abend
      endif
c
      endif
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(length)
        call Unused_integer(adrnext)
      end if
      end
c
c
c
c
      subroutine determine_factor3(length,term,iact,nbas,posu,post,poss,
     *                       vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,scr3,
     *                       A,C,nbasp,nbaso,dkhadr,adrmem,dkh_48,
     *                       adrnext)
c
c****************************************************************************
c
c   Determine which factor stands at position 'iact' of 'term' on the left
c     hand side.
c
c   Adjust value of 'iact', such that it points at beginning of next factor.
c
c   Note: There are no brackets [] occurring in term.
c
c   input: A right factor
c          C resulting product matrix
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Markus Reiher  (ETH Zurich)
c
c   version:  2.2.0
c
c   modified: 26.02.2007  by M. Reiher (ETH Zurich)
c        * can now read in right hand side factor from disk
c
c   first version: 14.04.2004  (Theoretical Chemistry, ETH Zurich)
c
c****************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
c
      integer length,iact,nbas,posu(maxunumber),post(maxsnumber),
     *        poss(maxsnumber),snumber,tnumber,unumber,nbasp,nbaso
      character*(*) term
c
      REAL*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp)
      REAL*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),A(nbas,nbas),
     *                 C(nbas,nbas)
c
      integer idum1,dkh_char2int
      REAL*8 alpha,beta
c
      integer dkh_48,adrmem,isize,iscr,iscr2,adrnext
      integer two,dkhadr(adrmem),adr
      parameter(two=2)
      logical flag
      flag=.true.
      isize=(nbas*(nbas+1)/2)
c
c------------------------------------------------------------------------
c
      alpha=1.0d0
      beta=0.0d0
c
      if (Out_Of_Core) then
        if (iact.gt.3) then
         if (term(iact-3:iact-3).eq.'S') then
           idum1=dkh_char2int(3,term(iact-2:iact))
           idum1=1000+idum1
           iact=iact-4
           adr=dkhadr(idum1)
           isize=nbas*nbas
           goto 666
         else if (term(iact-3:iact-3).eq.'T') then
           idum1=dkh_char2int(3,term(iact-2:iact))
           idum1=2000+idum1
           iact=iact-4
           adr=dkhadr(idum1)
           isize=nbas*nbas
           goto 666
         else if (term(iact-3:iact-3).eq.'U') then
           idum1=dkh_char2int(3,term(iact-2:iact))
           idum1=3000+idum1
           iact=iact-4
           adr=dkhadr(idum1)
           isize=nbas*nbas
           goto 666
         end if
       end if
       if (term(iact:iact).eq.'V') then
         iact=iact-1
         adr=dkhadr(1)
       else if (term(iact:iact).eq.'N') then
         iact=iact-1
         adr=dkhadr(5)
       else if (term(iact:iact).eq.'D') then
         iact=iact-1
         adr=dkhadr(2)
       else if (term(iact:iact).eq.'Y') then
         iact=iact-1
         adr=dkhadr(6)
       else if (term(iact:iact).eq.'F') then
         iact=iact-1
         adr=dkhadr(7)
       else if (term(iact:iact).eq.'G') then
         iact=iact-1
         adr=dkhadr(8)
       else if (term(iact:iact).eq.'Z') then
         call mat_times_p2c(C,A,nbas,pp)
         iact=iact-1
         flag=.false.
       else if (term(iact:iact).eq.'Q') then
         call mat_div_p2c(C,A,nbas,pp)
         iact=iact-1
         flag=.false.
       else if (term(iact:iact).eq.'X') then
         iact=iact-1
         adr=dkhadr(3)
       else if (term(iact:iact).eq.'I') then
         iact=iact-1
         adr=dkhadr(9)
       else if (term(iact:iact).eq.'J') then
         iact=iact-1
         adr=dkhadr(4)
       else if (term(iact:iact).eq.'K') then
         iact=iact-1
         adr=dkhadr(10)
       else if (term(iact:iact).eq.'L') then
         iact=iact-1
         adr=dkhadr(11)
       else if (term(iact:iact).eq.'M') then
         iact=iact-1
         adr=dkhadr(12)
       end if
666    continue
       if (flag) then
CMR it would be advantageous to have a triangular and a square
CMR   scratch matrix permanently available
         call GetMem('DetFac  ','ALLO','REAL',iscr,isize+4)
         call ddafile(dkh_48,two,work(iscr),isize,adr)
         if (isize.eq.nbas*nbas) then
           call DGEMM_('N','N',nbas,nbas,nbas,alpha,work(iscr),
     *               nbas,A,nbas,beta,C,nbas)
         else
           call GetMem('DetFac  ','ALLO','REAL',iscr2,nbas*nbas+4)
           call mat_sq_from_t(work(iscr2),nbas,work(iscr))
           call DGEMM_('N','N',nbas,nbas,nbas,alpha,work(iscr2),
     *               nbas,A,nbas,beta,C,nbas)
           call GetMem('DetFac  ','FREE','REAL',iscr2,nbas*nbas+4)
         end if
         call GetMem('DetFac  ','FREE','REAL',iscr,isize+4)
       end if
      else
c
      if (iact.gt.3) then
       if (term(iact-3:iact-3).eq.'S') then
        idum1=dkh_char2int(3,term(iact-2:iact))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,scr1(1,1,poss(idum1)),
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-4
        return
       else if (term(iact-3:iact-3).eq.'T') then
        idum1=dkh_char2int(3,term(iact-2:iact))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,scr2(1,1,post(idum1)),
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-4
        return
       else if (term(iact-3:iact-3).eq.'U') then
        idum1=dkh_char2int(3,term(iact-2:iact))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,scr3(1,1,posu(idum1)),
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-4
        return
       end if
      end if
      if (term(iact:iact).eq.'V') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,vv,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'N') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,nn,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'D') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,dd,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'Y') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,yy,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'F') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,ff,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'G') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,gg,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'Z') then
        call mat_times_p2c(C,A,nbas,pp)
        iact=iact-1
      else if (term(iact:iact).eq.'Q') then
        call mat_div_p2c(C,A,nbas,pp)
        iact=iact-1
      else if (term(iact:iact).eq.'X') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,xx,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'I') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,ii,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'J') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,jj,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'K') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,kk,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'L') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,ll,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'M') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,mm,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else
            write (stdout,1111)
1111        format (2X,'MIRACLE in determine_factor3(): could not',
     *        ' determine 2nd factor!'//2X,'STOP.',/2X)
            CALL Abend
      endif
c
      endif
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(length)
        call Unused_integer(adrnext)
      end if
      end
c
c
c  NOTE: THE FOLLOWING ROUTINES determine_factor...b() ARE DUPLICATED
c        SO THAT THEY COULD ALSO CHECK FOR THE AUXILIARY MATRICES SET UP
c        IN evalstring()
c
c
c
      subroutine determine_factorb(length,term,iact,nbas,posu,post,poss,
     *                       vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,
     *                       scr3,factor,nbasp,nbaso,ieval,count,dkhadr,
     *                       adrmem,dkh_48,adrnext)
c
c****************************************************************************
c
c   Determine which factor stands at position 'iact' of 'term'.
c   Store its matrix representation in 'factor'.
c
c   Adjust value of 'iact', such that it points at beginning of next factor.
c
c   Note: There are no brackets [] occurring in term.
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.2.0
c
c   last modified: 26.02.2007 (M. Reiher, ETH Zurich)
c                  * can now read in right hand side factor from disk
c                  * a couple of loops were re-arranged for speed up
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c****************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
c
      integer length,iact,nbas,posu(maxunumber),post(maxsnumber),
     *        poss(maxsnumber),snumber,tnumber,unumber,nbasp,
     *        count,ieval(count),nbaso
      character*(maxlength) term
c
      REAL*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp)
      REAL*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),factor(nbas,nbas)
c
      integer idum1,dkh_char2int,i,j
c
      integer dkh_48,adrmem,isize,iscr,adrnext
      integer two,dkhadr(adrmem),adr
      parameter(two=2)
      logical flag
      flag=.true.
      isize=(nbas*(nbas+1)/2)
c
c------------------------------------------------------------------------
c
c
      if (Out_Of_Core) then
       if (term(iact:iact).eq.'S') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         idum1=1000+idum1
         iact=iact+4
         adr=dkhadr(idum1)
         isize=nbas*nbas
       else if (term(iact:iact).eq.'T') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         idum1=2000+idum1
         iact=iact+4
         adr=dkhadr(idum1)
         isize=nbas*nbas
       else if (term(iact:iact).eq.'U') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         idum1=3000+idum1
         iact=iact+4
         adr=dkhadr(idum1)
         isize=nbas*nbas
       else if (term(iact:iact).eq.'A') then
         idum1=dkh_char2int(2,term(iact+1:iact+2))
CMR         idum1=50+idum1
CMR         adr=dkhadr(idum1)
CMR         isize=nbas*nbas
         call mat_copy(factor,nbas,nbas,work(ieval(idum1)))
         iact=iact+3
         flag=.false.
       else if (term(iact:iact).eq.'V') then
         iact=iact+1
         adr=dkhadr(1)
       else if (term(iact:iact).eq.'N') then
         iact=iact+1
         adr=dkhadr(5)
       else if (term(iact:iact).eq.'D') then
         iact=iact+1
         adr=dkhadr(2)
       else if (term(iact:iact).eq.'Y') then
         iact=iact+1
         adr=dkhadr(6)
       else if (term(iact:iact).eq.'F') then
         iact=iact+1
         adr=dkhadr(7)
       else if (term(iact:iact).eq.'G') then
         iact=iact+1
         adr=dkhadr(8)
       else if (term(iact:iact).eq.'Z') then
         call mat_sq_from_d (factor,nbas,pp)
         iact=iact+1
         flag=.false.
       else if (term(iact:iact).eq.'Q') then
         call mat_sq_dev_d (factor,nbas,pp)
         iact=iact+1
         flag=.false.
       else if (term(iact:iact).eq.'X') then
         iact=iact+1
         adr=dkhadr(3)
       else if (term(iact:iact).eq.'I') then
         iact=iact+1
         adr=dkhadr(9)
       else if (term(iact:iact).eq.'J') then
         iact=iact+1
         adr=dkhadr(4)
       else if (term(iact:iact).eq.'K') then
         iact=iact+1
         adr=dkhadr(10)
       else if (term(iact:iact).eq.'L') then
         iact=iact+1
         adr=dkhadr(11)
       else if (term(iact:iact).eq.'M') then
         iact=iact+1
         adr=dkhadr(12)
       end if
       if (flag) then
        if (isize.eq.nbas*nbas) then
         call ddafile(dkh_48,two,factor,nbas*nbas,adr)
        else
CMR it would be advantageous to have a triangular scratch
CMR   matrix permanently available
         call GetMem('DetFac  ','ALLO','REAL',iscr,isize+4)
         call ddafile(dkh_48,two,work(iscr),isize,adr)
         call mat_sq_from_t(factor,nbas,work(iscr))
         call GetMem('DetFac  ','FREE','REAL',iscr,isize+4)
        end if
       end if
      else
c
      if (term(iact:iact).eq.'S') then
        idum1=dkh_char2int(3,term(iact+1:iact+3))
        do 10 i=1,nbas
          do 20 j=1,nbas
            factor(j,i)=scr1(j,i,poss(idum1))
  20      continue
  10    continue
        iact=iact+4
      else if (term(iact:iact).eq.'T') then
        idum1=dkh_char2int(3,term(iact+1:iact+3))
        do 30 i=1,nbas
          do 40 j=1,nbas
            factor(j,i)=scr2(j,i,post(idum1))
  40      continue
  30    continue
        iact=iact+4
      else if (term(iact:iact).eq.'U') then
        idum1=dkh_char2int(3,term(iact+1:iact+3))
        do 50 i=1,nbas
          do 60 j=1,nbas
            factor(j,i)=scr3(j,i,posu(idum1))
  60      continue
  50    continue
        iact=iact+4
      else if (term(iact:iact).eq.'A') then
        idum1=dkh_char2int(2,term(iact+1:iact+2))
CMR        write(*,*) "FOUND A: term=,",term(iact:iact+2),
CMR     *             " No=",idum1," work=",work(ieval(idum1))
        call mat_copy(factor,nbas,nbas,work(ieval(idum1)))
        iact=iact+3
      else if (term(iact:iact).eq.'V') then
        call mat_copy (factor,nbas,nbas,vv)
        iact=iact+1
      else if (term(iact:iact).eq.'N') then
        call mat_copy (factor,nbas,nbas,nn)
        iact=iact+1
      else if (term(iact:iact).eq.'D') then
        call mat_copy (factor,nbas,nbas,dd)
        iact=iact+1
      else if (term(iact:iact).eq.'Y') then
        call mat_copy (factor,nbas,nbas,yy)
        iact=iact+1
      else if (term(iact:iact).eq.'F') then
        call mat_copy (factor,nbas,nbas,ff)
        iact=iact+1
      else if (term(iact:iact).eq.'G') then
        call mat_copy (factor,nbas,nbas,gg)
        iact=iact+1
      else if (term(iact:iact).eq.'Z') then
        call mat_sq_from_d (factor,nbas,pp)
        iact=iact+1
      else if (term(iact:iact).eq.'Q') then
        call mat_sq_dev_d (factor,nbas,pp)
        iact=iact+1
      else if (term(iact:iact).eq.'X') then
CMR        if (nbasp.eq.1) stop "PROGRAM ERROR"
        call mat_copy (factor,nbasp,nbasp,xx)
        iact=iact+1
      else if (term(iact:iact).eq.'I') then
CMR        if (nbasp.eq.1) stop "PROGRAM ERROR"
        call mat_copy (factor,nbasp,nbasp,ii)
        iact=iact+1
      else if (term(iact:iact).eq.'J') then
CMR        if (nbasp.eq.1) stop "PROGRAM ERROR"
        call mat_copy (factor,nbasp,nbasp,jj)
        iact=iact+1
      else if (term(iact:iact).eq.'K') then
CMR        if (nbasp.eq.1) stop "PROGRAM ERROR"
        call mat_copy (factor,nbasp,nbasp,kk)
        iact=iact+1
      else if (term(iact:iact).eq.'L') then
CMR        if (nbasp.eq.1) stop "PROGRAM ERROR"
        call mat_copy (factor,nbasp,nbasp,ll)
        iact=iact+1
      else if (term(iact:iact).eq.'M') then
CMR        if (nbasp.eq.1) stop "PROGRAM ERROR"
        call mat_copy (factor,nbasp,nbasp,mm)
        iact=iact+1
      else
        write (stdout,1083)
1083        format (2X,'ERROR in determine_factorb(): could not ',
     *        'determine factor!'//2X,'STOP.',/2X)
        write(6,*) "term(iact:iact)=",term(iact:iact)
            CALL Abend
      endif
c
      endif
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(length)
        call Unused_integer(adrnext)
      end if
      end
c
c
c
c
      subroutine determine_factor2b(length,term,iact,nbas,posu,post,poss
     *                       ,vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,
     *                       scr3,A,C,nbasp,nbaso,ieval,count,dkhadr,
     *                       adrmem,dkh_48,adrnext)
c
c****************************************************************************
c
c   Determine which factor stands at position 'iact' of 'term'.
c
c   Adjust value of 'iact', such that it points at beginning of next factor.
c
c   Note: There are no brackets [] occurring in term.
c
c   input: A left factor
c          C resulting product matrix
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Markus Reiher  (ETH Zurich)
c
c   version:  2.2.0
c
c   modified: 26.02.2007  by M. Reiher (ETH Zurich)
c        * can now read in right hand side factor from disk
c
c   first version: 12.04.2004  (Theoretical Chemistry, ETH Zurich)
c
c****************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
c
      integer length,iact,nbas,posu(maxunumber),post(maxsnumber),
     *        poss(maxsnumber),snumber,tnumber,unumber,nbasp,
     *        count,ieval(count),nbaso
      character*(maxlength) term
c
      REAL*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp)
      REAL*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),A(nbas,nbas),
     *                 C(nbas,nbas)
c
      integer idum1,dkh_char2int
      REAL*8 alpha,beta
c
      integer dkh_48,adrmem,isize,iscr,iscr2,adrnext
      integer two,dkhadr(adrmem),adr
      parameter(two=2)
      logical flag
      flag=.true.
      isize=(nbas*(nbas+1)/2)
c
c------------------------------------------------------------------------
c
      alpha=1.0d0
      beta=0.0d0
c
      if (Out_Of_Core) then
       if (term(iact:iact).eq.'S') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         idum1=1000+idum1
         iact=iact+4
         adr=dkhadr(idum1)
         isize=nbas*nbas
       else if (term(iact:iact).eq.'T') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         idum1=2000+idum1
         iact=iact+4
         adr=dkhadr(idum1)
         isize=nbas*nbas
       else if (term(iact:iact).eq.'U') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         idum1=3000+idum1
         iact=iact+4
         adr=dkhadr(idum1)
         isize=nbas*nbas
       else if (term(iact:iact).eq.'A') then
         idum1=dkh_char2int(2,term(iact+1:iact+2))
         call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,work(ieval(idum1)),nbas,beta,
     *              C,nbas)
         iact=iact+3
CMR         idum1=50+idum1
CMR         adr=dkhadr(idum1)
CMR         isize=nbas*nbas
         flag=.false.
       else if (term(iact:iact).eq.'V') then
         iact=iact+1
         adr=dkhadr(1)
       else if (term(iact:iact).eq.'N') then
         iact=iact+1
         adr=dkhadr(5)
       else if (term(iact:iact).eq.'D') then
         iact=iact+1
         adr=dkhadr(2)
       else if (term(iact:iact).eq.'Y') then
         iact=iact+1
         adr=dkhadr(6)
       else if (term(iact:iact).eq.'F') then
         iact=iact+1
         adr=dkhadr(7)
       else if (term(iact:iact).eq.'G') then
         iact=iact+1
         adr=dkhadr(8)
       else if (term(iact:iact).eq.'Z') then
         call mat_times_p2b(C,A,nbas,pp)
         iact=iact+1
         flag=.false.
       else if (term(iact:iact).eq.'Q') then
         call mat_div_p2b(C,A,nbas,pp)
         iact=iact+1
         flag=.false.
       else if (term(iact:iact).eq.'X') then
         iact=iact+1
         adr=dkhadr(3)
       else if (term(iact:iact).eq.'I') then
         iact=iact+1
         adr=dkhadr(9)
       else if (term(iact:iact).eq.'J') then
         iact=iact+1
         adr=dkhadr(4)
       else if (term(iact:iact).eq.'K') then
         iact=iact+1
         adr=dkhadr(10)
       else if (term(iact:iact).eq.'L') then
         iact=iact+1
         adr=dkhadr(11)
       else if (term(iact:iact).eq.'M') then
         iact=iact+1
         adr=dkhadr(12)
       end if
       if (flag) then
CMR it would be advantageous to have a triangular and a square
CMR   scratch matrix permanently available
         call GetMem('DetFac  ','ALLO','REAL',iscr,isize+4)
         call ddafile(dkh_48,two,work(iscr),isize,adr)
         if (isize.eq.nbas*nbas) then
           call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *               nbas,work(iscr),nbas,beta,C,nbas)
         else
           call GetMem('DetFac  ','ALLO','REAL',iscr2,nbas*nbas+4)
           call mat_sq_from_t(work(iscr2),nbas,work(iscr))
           call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *               nbas,work(iscr2),nbas,beta,C,nbas)
           call GetMem('DetFac  ','FREE','REAL',iscr2,nbas*nbas+4)
         end if
         call GetMem('DetFac  ','FREE','REAL',iscr,isize+4)
       end if
      else
c
      if (term(iact:iact).eq.'S') then
        idum1=dkh_char2int(3,term(iact+1:iact+3))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,scr1(1,1,poss(idum1)),nbas,beta,
     *              C,nbas)
        iact=iact+4
      else if (term(iact:iact).eq.'T') then
        idum1=dkh_char2int(3,term(iact+1:iact+3))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,scr2(1,1,post(idum1)),nbas,beta,
     *              C,nbas)
        iact=iact+4
      else if (term(iact:iact).eq.'U') then
        idum1=dkh_char2int(3,term(iact+1:iact+3))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,scr3(1,1,posu(idum1)),nbas,beta,
     *              C,nbas)
        iact=iact+4
      else if (term(iact:iact).eq.'A') then
        idum1=dkh_char2int(2,term(iact+1:iact+2))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,work(ieval(idum1)),nbas,beta,
     *              C,nbas)
        iact=iact+3
      else if (term(iact:iact).eq.'V') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,vv,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'N') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,nn,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'D') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,dd,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'Y') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,yy,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'F') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,ff,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'G') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,gg,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'Z') then
        call mat_times_p2b(C,A,nbas,pp)
        iact=iact+1
      else if (term(iact:iact).eq.'Q') then
        call mat_div_p2b(C,A,nbas,pp)
        iact=iact+1
      else if (term(iact:iact).eq.'X') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,xx,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'I') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,ii,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'J') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,jj,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'K') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,kk,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'L') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,ll,nbas,beta,C,nbas)
        iact=iact+1
      else if (term(iact:iact).eq.'M') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,A,
     *              nbas,mm,nbas,beta,C,nbas)
        iact=iact+1
      else
            write (stdout,1111)
1111        format (2X,'MIRACLE in determine_factor2b(): could not',
     *        ' determine 2nd factor!'//2X,'STOP.',/2X)
            CALL Abend
      endif
c
      endif
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(length)
        call Unused_integer(adrnext)
      end if
      end
c
c
c
c
      subroutine determine_factor3b(length,term,iact,nbas,posu,post,poss
     *                       ,vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,scr3,
     *                       A,C,nbasp,nbaso,ieval,count,dkhadr,adrmem,
     *                       dkh_48,adrnext)
c
c****************************************************************************
c
c   Determine which factor stands at position 'iact' of 'term' on the left.
c
c   Adjust value of 'iact', such that it points at beginning of next factor.
c
c   Note: There are no brackets [] occurring in term.
c
c   input: A right factor
c          C resulting product matrix
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Markus Reiher  (ETH Zurich)
c
c   version:  2.2.0
c
c   modified: 26.02.2007  by M. Reiher (ETH Zurich)
c        * can now read in right hand side factor from disk
c
c   first version: 14.04.2004  (Theoretical Chemistry, ETH Zurich)
c
c****************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
c
      integer length,iact,nbas,posu(maxunumber),post(maxsnumber),
     *        poss(maxsnumber),snumber,tnumber,unumber,nbasp,
     *        count,ieval(count),nbaso
      character*(maxlength) term
c
      REAL*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp)
      REAL*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),A(nbas,nbas),
     *                 C(nbas,nbas)
c
      integer idum1,dkh_char2int
      REAL*8 alpha,beta
c
      integer dkh_48,adrmem,isize,iscr,iscr2,adrnext
      integer two,dkhadr(adrmem),adr
      parameter(two=2)
      logical flag
      flag=.true.
      isize=(nbas*(nbas+1)/2)
c
c------------------------------------------------------------------------
c
      alpha=1.0d0
      beta=0.0d0
c
      if (Out_Of_Core) then
       if (iact.gt.3) then
         if (term(iact-3:iact-3).eq.'S') then
           idum1=dkh_char2int(3,term(iact-2:iact))
           idum1=1000+idum1
           iact=iact-4
           adr=dkhadr(idum1)
           isize=nbas*nbas
           goto 666
         else if (term(iact-3:iact-3).eq.'T') then
           idum1=dkh_char2int(3,term(iact-2:iact))
           idum1=2000+idum1
           iact=iact-4
           adr=dkhadr(idum1)
           isize=nbas*nbas
           goto 666
         else if (term(iact-3:iact-3).eq.'U') then
           idum1=dkh_char2int(3,term(iact-2:iact))
           idum1=3000+idum1
           iact=iact-4
           adr=dkhadr(idum1)
           isize=nbas*nbas
           goto 666
         end if
       end if
       if (iact.gt.2) then
         if (term(iact-2:iact-2).eq.'A') then
           idum1=dkh_char2int(2,term(iact-1:iact))
CMR           idum1=50+idum1
CMR           adr=dkhadr(idum1)
CMR           isize=nbas*nbas
           call DGEMM_('N','N',nbas,nbas,nbas,alpha,work(ieval(idum1)),
     *              nbas,A,nbas,beta,C,nbas)
           iact=iact-3
           flag=.false.
           goto 666
         end if
       end if
       if (term(iact:iact).eq.'V') then
         iact=iact-1
         adr=dkhadr(1)
       else if (term(iact:iact).eq.'N') then
         iact=iact-1
         adr=dkhadr(5)
       else if (term(iact:iact).eq.'D') then
         iact=iact-1
         adr=dkhadr(2)
       else if (term(iact:iact).eq.'Y') then
         iact=iact-1
         adr=dkhadr(6)
       else if (term(iact:iact).eq.'F') then
         iact=iact-1
         adr=dkhadr(7)
       else if (term(iact:iact).eq.'G') then
         iact=iact-1
         adr=dkhadr(8)
       else if (term(iact:iact).eq.'Z') then
         call mat_times_p2c(C,A,nbas,pp)
         iact=iact-1
         flag=.false.
       else if (term(iact:iact).eq.'Q') then
         call mat_div_p2c(C,A,nbas,pp)
         iact=iact-1
         flag=.false.
       else if (term(iact:iact).eq.'X') then
         iact=iact-1
         adr=dkhadr(3)
       else if (term(iact:iact).eq.'I') then
         iact=iact-1
         adr=dkhadr(9)
       else if (term(iact:iact).eq.'J') then
         iact=iact-1
         adr=dkhadr(4)
       else if (term(iact:iact).eq.'K') then
         iact=iact-1
         adr=dkhadr(10)
       else if (term(iact:iact).eq.'L') then
         iact=iact-1
         adr=dkhadr(11)
       else if (term(iact:iact).eq.'M') then
         iact=iact-1
         adr=dkhadr(12)
       end if
666    continue
       if (flag) then
CMR it would be advantageous to have a triangular and a square
CMR   scratch matrix permanently available
         call GetMem('DetFac  ','ALLO','REAL',iscr,isize+4)
         call ddafile(dkh_48,two,work(iscr),isize,adr)
         if (isize.eq.nbas*nbas) then
           call DGEMM_('N','N',nbas,nbas,nbas,alpha,work(iscr),
     *               nbas,A,nbas,beta,C,nbas)
         else
           call GetMem('DetFac  ','ALLO','REAL',iscr2,nbas*nbas+4)
           call mat_sq_from_t(work(iscr2),nbas,work(iscr))
           call DGEMM_('N','N',nbas,nbas,nbas,alpha,work(iscr2),
     *               nbas,A,nbas,beta,C,nbas)
           call GetMem('DetFac  ','FREE','REAL',iscr2,nbas*nbas+4)
         end if
         call GetMem('DetFac  ','FREE','REAL',iscr,isize+4)
       end if
      else
c
      if (iact.gt.3) then
       if (term(iact-3:iact-3).eq.'S') then
        idum1=dkh_char2int(3,term(iact-2:iact))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,scr1(1,1,poss(idum1)),
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-4
        return
       else if (term(iact-3:iact-3).eq.'T') then
        idum1=dkh_char2int(3,term(iact-2:iact))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,scr2(1,1,post(idum1)),
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-4
        return
       else if (term(iact-3:iact-3).eq.'U') then
        idum1=dkh_char2int(3,term(iact-2:iact))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,scr3(1,1,posu(idum1)),
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-4
        return
       end if
      end if
      if (iact.gt.2) then
       if (term(iact-2:iact-2).eq.'A') then
        idum1=dkh_char2int(2,term(iact-1:iact))
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,work(ieval(idum1)),
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-3
        return
       end if
      end if
      if (term(iact:iact).eq.'V') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,vv,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'N') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,nn,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'D') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,dd,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'Y') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,yy,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'F') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,ff,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'G') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,gg,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'Z') then
        call mat_times_p2c(C,A,nbas,pp)
        iact=iact-1
      else if (term(iact:iact).eq.'Q') then
        call mat_div_p2c(C,A,nbas,pp)
        iact=iact-1
      else if (term(iact:iact).eq.'X') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,xx,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'I') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,ii,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'J') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,jj,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'K') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,kk,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'L') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,ll,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else if (term(iact:iact).eq.'M') then
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,mm,
     *              nbas,A,nbas,beta,C,nbas)
        iact=iact-1
      else
            write (stdout,1111)
1111        format (2X,'MIRACLE in determine_factor3b(): could not',
     *        ' determine 2nd factor!'//2X,'STOP.',/2X)
            CALL Abend
      endif
c
      endif
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(length)
        call Unused_integer(adrnext)
      end if
      end
