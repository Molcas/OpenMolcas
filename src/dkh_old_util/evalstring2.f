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
* Copyright (C) 2007, Markus Reiher                                    *
************************************************************************
      subroutine evalstring2(length,term,coeff,nbas,posu,post,poss,vv,
     *                       nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,
     *                       snumber,tnumber,unumber,scrno1,
     *                       scr1,scr2,scr3,scr5,nbasp,nbaso,darray,
     *                       dkhadr,adrmem,dkh_48,adrnext)
c
c************************************************************************
c
c   Evaluate the low-level character string 'term' of length 'length'
c     with coefficient 'coeff'. The result is stored in scr5(,,scrno1).
c
c   Note:  'term' might in general contain Sxxx, Txxx, and Uxxx expressions
c           as well as brackets [,], and V,N,D,Y,F,G,Z,Q,X,I,J,K,L,M  !!
c
c   scrno1 has to be 2 for the purposes here!
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Markus Reiher  (ETH Zurich)
c
c   version:  1.0.0
c
c   * the A matrices introduced in previous versions of this
c     routine have been eliminated
c   * this routine requires scr5(nbas,nbas,scrno1) with
c     2 (nbas X nbas)-matrices; however, scrno1=3 since
c     scr5(nbas,nbas,1) is already reserved in calc_operators()
c   * this routine now assumes that the bracket structure
c     of the term is comparatively simple
c     (either no bracket, one bracket only, contiguous brackets,
c      or only(!) nested brackets --- all other
c      cases cannot be treated here and are dealt with
c      at the level of evalstring() )
c
c   first version: 17.01.2007  (Theoretical Chemistry, ETH Zurich)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer length,nbas,posu(maxunumber),post(maxsnumber),nbasp,
     *        poss(maxsnumber),snumber,tnumber,unumber,scrno1,nbaso
      character*(maxlength) term
      REAL*8 coeff
c
      REAL*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp),e(nbas)
      REAL*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),scr5(nbas,nbas,scrno1),
     *                 darray(nbas,nbas)
c
      integer i,idum1,ileft,iright,scrleng,iposA
      character*(maxlength) scrtxt
      REAL*8 alpha,beta
      logical flag
c
      integer dkh_48,adrmem,adrnext
      integer dkhadr(adrmem)
c
      flag=.true.
      alpha=1.0d0
      beta=0.0d0
      iposA=0
c
3150  continue
c
      iright=0
      ileft=0
c
c-----------------------------------------------------------------------
c  The sequence of evaluation of 'term' is determined by the position of
c    the brackets [].
c    --> Find first bracket ']' at iright and --- if available ---
c          corresponding '[' at ileft -> yields the innermost bracket
c------------------------------------------------------------------------
c
      iright=index(term(1:length),']')
*** determine position ileft of corresponding '['
      if (iright.ne.0) then
        ileft=iright-1
1065    continue
        if (term(ileft:ileft).ne.'[') then
          ileft=ileft-1
          goto 1065
        endif
        if (ileft.eq.0 .or. ileft.ge.iright) then
          write (stdout,1010) iright,term(1:maxlength)
1010      format (2X,'Error1 in SR evalstring: unbalanced occurrence ',
     *            'of brackets:',/2X,'] at pos = ',I2,', but no [ left',
     *            ' of it in term = ',A90)
          CALL Abend
        endif
      endif
c
      if (iright.eq.0) then
c
c-----------------------------------------------------------------------
c  a) Easy case: Direct evaluation from left to right of 'term' possible
c                because NO brackets occur at all
c  -------------
*** determine left-hand-side factor if it does not exist yet
        if (flag) then
          call multiply (length,term,coeff,nbas,posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1-1),
     *                 scr5(1,1,scrno1),nbasp,nbaso,dkhadr,adrmem,
     *                 dkh_48,adrnext)
        else
*** in this case we have a left-hand factor in A00 already identified
***   and there is no additional bracket left over
*** first multiply with all matrices from the right
*** note: the coefficient has to be alpha=1, not coeff in order to avoid
***       a double multiplication
          call multiply2(length-iposA+1,term(iposA:length),coeff,nbas,
     *                 iposA,posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1-1),
     *                 scr5(1,1,scrno1),nbasp,nbaso,dkhadr,adrmem,
     *                 dkh_48,adrnext)
*** check whether there are also matrices to be multiplied from the left
          if (iposA.ne.1) then
            call multiply3(iposA-1,term(1:iposA-1),alpha,nbas,
     *                 posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1),
     *                 scr5(1,1,scrno1-1),nbasp,nbaso,dkhadr,adrmem,
     *                 dkh_48,adrnext)
            call mat_copy(scr5(1,1,scrno1),nbas,nbas,scr5(1,1,scrno1-1))
          end if
        end if
c
      else if (iright.ne.0.and.iposA.ne.0) then
c
c-----------------------------------------------------------------------
c  b) most difficult case: the previously found brackets are not embraced
c                          by the new ones but additional brackets show up
c                          directly after the first bracket on the left
c                          was closed (other cases are not permitted by evalstring() )
c                          like in ...[...]...[...] and the first
c                          brackets were already found
c
*** first check if 'old' A00 must be multiplied by matrices from the left
        if (iposA.ne.1) then
          call multiply3(iposA-1,term(1:iposA-1),alpha,nbas,
     *                 posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1-1),
     *                 scr5(1,1,scrno1),nbasp,nbaso,dkhadr,adrmem,
     *                 dkh_48,adrnext)
*** now: scr5(1,1,scrno1) contains the new A00
          call mat_copy(scr5(1,1,scrno1-1),nbas,nbas,scr5(1,1,scrno1))
*** idum1 is the length of the matrices on the left of A00
          idum1=iposA-1
          do 71 i=1,length-idum1
            term(i:i)=term(i+idum1:i+idum1)
  71      continue
*** account for new iposA
          iposA=1
          do 76 i=length-idum1+1,length
            term(i:i)=' '
  76      continue
          length=length-idum1
          ileft=ileft-idum1
          iright=iright-idum1
        end if
*** now check if there are matrices between 1st [...] and 2nd [...]
        if (iposA+3.ne.ileft) then
*** store matrix string befor '[' of length (ileft-1); note that
***   A00 is to be included because it is dealt with in multiply2()
          scrleng=ileft-1
          scrtxt(1:scrleng)=term(1:ileft-1)
*** call routine which gets A00 as input and multiplies the other
***   matrices from the left
          call multiply2(scrleng,scrtxt,alpha,nbas,
     *                 iposA,posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1-1),
     *                 scr5(1,1,scrno1),nbasp,nbaso,dkhadr,adrmem,
     *                 dkh_48,adrnext)
*** now: scr5(1,1,scrno1) contains the new A00
          call mat_copy(scr5(1,1,scrno1-1),nbas,nbas,scr5(1,1,scrno1))
*** now: scr5(1,1,scrno1-1) contains the A00 as before
*** set new length since A00 remains as it stands
          scrleng=scrleng-3
          do 72 i=4,length-scrleng
            term(i:i)=term(i+scrleng:i+scrleng)
  72      continue
          do 77 i=length-scrleng+1,length
            term(i:i)=' '
  77      continue
          length=length-scrleng
          ileft=ileft-scrleng
          iright=iright-scrleng
        end if
*** evaluate the directly following new bracket starting right after A00
*** include the energy denominators for the new bracket
*** first: determine string 'scrtxt' between these 2 brackets of 'term'
        scrleng=0
        scrleng=iright-ileft-1
c
c  note: at this stage there are no nested brackets in the second
c        bracket under consideration now as these terms would have been
c        evaluated already through evalstring()
c
*** initialize dummy string scrtxt
        scrtxt(1:scrleng)=term(ileft+1:iright-1)
        call multiply(scrleng,scrtxt,alpha,nbas,posu,post,poss,vv,nn,
     *                dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,
     *                tnumber,unumber,scr1,scr2,scr3,scr5(1,1,scrno1),
     *                darray,nbasp,nbaso,dkhadr,adrmem,dkh_48,adrnext)
        call mat_1_over_h2 (scr5(1,1,scrno1),nbas,e,darray)
*** new right-hand-side A00 is now on scr5(1,1,scrno1)
*** multiply both matrices: old A00 times new A00
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,scr5(1,1,scrno1-1),
     *              nbas,scr5(1,1,scrno1),nbas,beta,darray,nbas)
*** now: scr5(1,1,scrno1-1) contains the new A00
        call mat_copy(scr5(1,1,scrno1-1),nbas,nbas,darray)
c
*** introduce a single A00 for later multiplication purposes, which is stored
*** in scr5(1,1,scrno1-1)
*** Update 'term', and copy rest of term 'idum1' positions to the left
        iposA=1
*** determine length for string in brackets including the two brackets (iright-ileft-1)+2
        idum1=iright-ileft+1
        do 73 i=4,length-idum1
          term(i:i)=term(i+idum1:i+idum1)
  73    continue
        do 78 i=length-idum1+1,length
          term(i:i)=' '
  78    continue
        length=length-idum1
c
*** there might be another bracket term directly following the ones already obtained
        goto 3150
c
      else if (iposA.eq.0.and.iright.ne.0) then
c
c-----------------------------------------------------------------------
c  c) difficult case: Evaluate innermost brackets [..] and substitute
c  -----------------          them by a single A00 expression in 'term'
c
c
*** determine string 'scrtxt' between these 2 brackets of 'term'
        scrleng=iright-ileft-1
c
c  If there was a bracket already present, this was already evaluated and
c    coded as A00 and already stored in scr5(1,1,scrno1-1);
c    adjust its position for scrtxt
c
        if (iposA.ne.0) then
          iposA=iposA-ileft
        end if
c
*** initialize dummy string scrtxt
        scrtxt(1:scrleng)=term(ileft+1:iright-1)
c
*** Evaluate scrtxt --> result is in scr5(,,scrno1) to be transformed
***   to [...] which is then in scr5(,,scrno1-1)
        if (flag) then
          call multiply(scrleng,scrtxt,alpha,nbas,posu,post,poss,vv,nn,
     *                dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,
     *                tnumber,unumber,scr1,scr2,scr3,scr5(1,1,scrno1-1),
     *                scr5(1,1,scrno1),nbasp,nbaso,dkhadr,adrmem,dkh_48,
     *                adrnext)
*** include the energy denominators
         call mat_1_over_h2 (scr5(1,1,scrno1-1),nbas,e,scr5(1,1,scrno1))
        else
CMR ACTUALLY, THIS CASE CAN NOW NO LONGER OCCUR
*** there was already a bracket present, i.e. we already stored a A00 in scr5
*** first multiply with all matrices from the right
*** note: the coefficient has to be alpha=1, not coeff in order to avoid
***       a double multiplication
          call multiply2(scrleng-iposA+1,scrtxt(iposA:scrleng),coeff,
     *                 nbas,iposA,posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1-1),
     *                 scr5(1,1,scrno1),nbasp,nbaso,dkhadr,adrmem,
     *                 dkh_48,adrnext)
*** check whether there are also matrices to be multiplied from the left
          if (iposA.ne.1) then
            call multiply3(iposA-1,scrtxt(1:iposA-1),alpha,nbas,
     *                 posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1),
     *                 scr5(1,1,scrno1-1),nbasp,nbaso,dkhadr,adrmem,
     *                 dkh_48,adrnext)
            call mat_1_over_h2 (scr5(1,1,scrno1-1),nbas,e,
     *                          scr5(1,1,scrno1-1))
          else
*** include the energy denominators
            call mat_1_over_h2 (scr5(1,1,scrno1-1),nbas,e,
     *                          scr5(1,1,scrno1))
          end if
        end if
c
c
*** introduce a single A00 for later multiplication purposes, which is stored
*** in scr5(1,1,scrno1-1); all following multiplications start from A00 as the
*** left hand side factor and consider later eventually missing terms on the left
*** of A00.
*** Update 'term', and copy rest of term 'idum1' positions to the left
        term(ileft:ileft+2)='A00'
        iposA=ileft
        idum1=(iright-ileft+1)-3
        do 70 i=ileft+3,length-idum1
          term(i:i)=term(i+idum1:i+idum1)
  70    continue
        do 75 i=length-idum1+1,length
          term(i:i)=' '
  75    continue
        length=length-idum1
        flag=.false.
c
        goto 3150
c
      else
        write(stdout,*) "MIRACLE in evalstring2() !"
        call Abend
      endif
c
      return
      end
c
c
c
c
      subroutine evalstring2b(length,term,coeff,nbas,posu,post,poss,vv,
     *                       nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,
     *                       snumber,tnumber,unumber,scrno1,
     *                       scr1,scr2,scr3,scr5,nbasp,nbaso,darray,
     *                       ieval,count,dkhadr,adrmem,dkh_48,adrnext)
c
c************************************************************************
c
c   Evaluate the low-level character string 'term' of length 'length'
c     with coefficient 'coeff'. The result is stored in scr5(,,scrno1).
c
c   Note:  'term' might in general contain Sxxx, Txxx, and Uxxx expressions
c           as well as brackets [,], and V,N,D,Y,F,G,Z,Q,X,I,J,K,L,M  !!
c
c   This routine is essentially the same as evalstring2() with the
c     exception that it also tests for the dynamically allocated strings
c     in evalstring() !
c
c   scrno1 has to be 2 for the purposes here!
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Markus Reiher  (ETH Zurich)
c
c   version:  1.1.0
c
c   modified: 13.03.07 MR@ETH
c   first version: 18.01.2007  (Theoretical Chemistry, ETH Zurich)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
c
      integer length,nbas,posu(maxunumber),post(maxsnumber),nbasp,
     *        poss(maxsnumber),snumber,tnumber,unumber,scrno1,count,
     *        ieval(count),nbaso
      character*(maxlength) term
      REAL*8 coeff
c
      REAL*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp),e(nbas)
      REAL*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),scr5(nbas,nbas,scrno1),
     *                 darray(nbas,nbas)
c
      integer i,idum1,ileft,iright,scrleng,iposA
      character*(maxlength) scrtxt
      REAL*8 alpha,beta
      logical flag
c
      integer dkh_48,adrmem,adrnext
      integer dkhadr(adrmem)
c
      flag=.true.
      alpha=1.0d0
      beta=0.0d0
      iposA=0
c
3150  continue
c
      iright=0
      ileft=0
c
c-----------------------------------------------------------------------
c  The sequence of evaluation of 'term' is determined by the position of
c    the brackets [].
c    --> Find first bracket ']' at iright and --- if available ---
c          corresponding '[' at ileft -> yields the innermost bracket
c------------------------------------------------------------------------
c
      iright=index(term(1:length),']')
*** determine position ileft of corresponding '['
      if (iright.ne.0) then
        ileft=iright-1
1065    continue
        if (term(ileft:ileft).ne.'[') then
          ileft=ileft-1
          goto 1065
        endif
        if (ileft.eq.0 .or. ileft.ge.iright) then
          write (stdout,1010) iright,term(1:maxlength)
1010      format (2X,'Error1 in SR evalstring: unbalanced occurrence ',
     *            'of brackets:',/2X,'] at pos = ',I2,', but no [ left',
     *            ' of it in term = ',A90)
          CALL Abend
        endif
      endif
c
      if (iright.eq.0) then
c
c-----------------------------------------------------------------------
c  a) Easy case: Direct evaluation from left to right of 'term' possible
c                because NO brackets occur at all
c  -------------
*** determine left-hand-side factor if it does not exist yet
        if (flag) then
          call multiplyb(length,term,coeff,nbas,posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1-1),
     *                 scr5(1,1,scrno1),nbasp,nbaso,ieval,count,
     *                 dkhadr,adrmem,dkh_48,adrnext)
        else
*** in this case we have a left-hand factor in A00 already identified
***   and there is no additional bracket left over
*** first multiply with all matrices from the right
*** note: the coefficient has to be alpha=1, not coeff in order to avoid
***       a double multiplication
          call multiply2b(length-iposA+1,term(iposA:length),coeff,nbas,
     *                 iposA,posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1-1),
     *                 scr5(1,1,scrno1),nbasp,nbaso,ieval,count,
     *                 dkhadr,adrmem,dkh_48,adrnext)
*** check whether there are also matrices to be multiplied from the left
          if (iposA.ne.1) then
            call multiply3b(iposA-1,term(1:iposA-1),alpha,nbas,
     *                 posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1),
     *                 scr5(1,1,scrno1-1),nbasp,nbaso,ieval,count,
     *                 dkhadr,adrmem,dkh_48,adrnext)
            call mat_copy(scr5(1,1,scrno1),nbas,nbas,scr5(1,1,scrno1-1))
          end if
        end if
c
      else if (iright.ne.0.and.iposA.ne.0) then
c
c-----------------------------------------------------------------------
c  b) most difficult case: the previously found brackets are not embraced
c                          by the new ones but additional brackets show up
c                          directly after the first bracket on the left
c                          was closed (other cases are not permitted by evalstring() )
c                          like in ...[...]...[...] and the first
c                          brackets were already found
c
*** first check if 'old' A00 must be multiplied by matrices from the left
        if (iposA.ne.1) then
          call multiply3b(iposA-1,term(1:iposA-1),alpha,nbas,
     *                 posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1-1),
     *                 scr5(1,1,scrno1),nbasp,nbaso,ieval,count,
     *                 dkhadr,adrmem,dkh_48,adrnext)
*** now: scr5(1,1,scrno1) contains the new A00
          call mat_copy(scr5(1,1,scrno1-1),nbas,nbas,scr5(1,1,scrno1))
*** idum1 is the length of the matrices on the left of A00
          idum1=iposA-1
          do 71 i=1,length-idum1
            term(i:i)=term(i+idum1:i+idum1)
  71      continue
*** account for new iposA
          iposA=1
          do 76 i=length-idum1+1,length
            term(i:i)=' '
  76      continue
          length=length-idum1
          ileft=ileft-idum1
          iright=iright-idum1
        end if
*** now check if there are matrices between 1st [...] and 2nd [...]
        if (iposA+3.ne.ileft) then
*** store matrix string befor '[' of length (ileft-1); note that
***   A00 is to be included because it is dealt with in multiply2()
          scrleng=ileft-1
          scrtxt(1:scrleng)=term(iposA+3:ileft-1)
*** call routine which gets A00 as input and multiplies the other
***   matrices from the left
          call multiply2b(scrleng,scrtxt,alpha,nbas,
     *                 iposA,posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1-1),
     *                 scr5(1,1,scrno1),nbasp,nbaso,ieval,count,
     *                 dkhadr,adrmem,dkh_48,adrnext)
*** now: scr5(1,1,scrno1) contains the new A00
          call mat_copy(scr5(1,1,scrno1-1),nbas,nbas,scr5(1,1,scrno1))
*** now: scr5(1,1,scrno1-1) contains the A00 as before
*** set new length since A00 remains as it stands
          scrleng=scrleng-3
          do 72 i=4,length-scrleng
            term(i:i)=term(i+scrleng:i+scrleng)
  72      continue
          do 77 i=length-scrleng+1,length
            term(i:i)=' '
  77      continue
          length=length-scrleng
          ileft=ileft-scrleng
          iright=iright-scrleng
        end if
*** evaluate the directly following new bracket starting right after A00
*** include the energy denominators for the new bracket
*** first: determine string 'scrtxt' between these 2 brackets of 'term'
        scrleng=0
        scrleng=iright-ileft-1
c
c  note: at this stage there are no nested brackets in the second
c        bracket under consideration now as these terms would have been
c        evaluated already through evalstring()
c
*** initialize dummy string scrtxt
        scrtxt(1:scrleng)=term(ileft+1:iright-1)
        call multiplyb(scrleng,scrtxt,alpha,nbas,posu,post,poss,vv,nn,
     *                dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,
     *                tnumber,unumber,scr1,scr2,scr3,scr5(1,1,scrno1),
     *                darray,nbasp,nbaso,ieval,count,dkhadr,adrmem,
     *                dkh_48,adrnext)
        call mat_1_over_h2 (scr5(1,1,scrno1),nbas,e,darray)
*** new right-hand-side A00 is now on scr5(1,1,scrno1)
*** multiply both matrices: old A00 times new A00
        call DGEMM_('N','N',nbas,nbas,nbas,alpha,scr5(1,1,scrno1-1),
     *              nbas,scr5(1,1,scrno1),nbas,beta,darray,nbas)
*** now: scr5(1,1,scrno1-1) contains the new A00
        call mat_copy(scr5(1,1,scrno1-1),nbas,nbas,darray)
c
*** introduce a single A00 for later multiplication purposes, which is stored
*** in scr5(1,1,scrno1-1)
*** Update 'term', and copy rest of term 'idum1' positions to the left
        iposA=1
*** determine length for string in brackets including the two brackets (iright-ileft-1)+2
        idum1=iright-ileft+1
        do 73 i=4,length-idum1
          term(i:i)=term(i+idum1:i+idum1)
  73    continue
        do 78 i=length-idum1+1,length
          term(i:i)=' '
  78    continue
        length=length-idum1
c
*** there might be another bracket term directly following the ones already obtained
        goto 3150
c
      else if (iposA.eq.0.and.iright.ne.0) then
c
c-----------------------------------------------------------------------
c  c) difficult case: Evaluate innermost brackets [..] and substitute
c  -----------------          them by a single A00 expression in 'term'
c
c
*** determine string 'scrtxt' between these 2 brackets of 'term'
        scrleng=iright-ileft-1
c
c  If there was a bracket already present, this was already evaluated and
c    coded as A00 and already stored in scr5(1,1,scrno1-1);
c    adjust its position for scrtxt
c
        if (iposA.ne.0) then
          iposA=iposA-ileft
        end if
c
*** initialize dummy string scrtxt
        scrtxt(1:scrleng)=term(ileft+1:iright-1)
c
*** Evaluate scrtxt --> result is in scr5(,,scrno1) to be transformed
***   to [...] which is then in scr5(,,scrno1-1)
        if (flag) then
          call multiplyb(scrleng,scrtxt,alpha,nbas,posu,post,poss,vv,nn,
     *                dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,
     *                tnumber,unumber,scr1,scr2,scr3,scr5(1,1,scrno1-1),
     *                scr5(1,1,scrno1),nbasp,nbaso,ieval,count,dkhadr,
     *                adrmem,dkh_48,adrnext)
*** include the energy denominators
         call mat_1_over_h2 (scr5(1,1,scrno1-1),nbas,e,scr5(1,1,scrno1))
        else
CMR ACTUALLY, THIS CASE CAN NO LONGER OCCUR
*** there was already a bracket present, i.e. we already stored a A00 in scr5
*** first multiply with all matrices from the right
*** note: the coefficient has to be alpha=1, not coeff in order to avoid
***       a double multiplication
          call multiply2b(scrleng-iposA+1,scrtxt(iposA:scrleng),coeff,
     *                 nbas,iposA,posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1-1),
     *                 scr5(1,1,scrno1),nbasp,nbaso,ieval,count,
     *                 dkhadr,adrmem,dkh_48,adrnext)
*** check whether there are also matrices to be multiplied from the left
          if (iposA.ne.1) then
            call multiply3b(iposA-1,scrtxt(1:iposA-1),alpha,nbas,
     *                 posu,post,poss,vv,nn,dd,
     *                 yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,tnumber,
     *                 unumber,scr1,scr2,scr3,scr5(1,1,scrno1),
     *                 scr5(1,1,scrno1-1),nbasp,nbaso,ieval,count,
     *                 dkhadr,adrmem,dkh_48,adrnext)
            call mat_1_over_h2 (scr5(1,1,scrno1-1),nbas,e,
     *                          scr5(1,1,scrno1-1))
          else
*** include the energy denominators
            call mat_1_over_h2 (scr5(1,1,scrno1-1),nbas,e,
     *                          scr5(1,1,scrno1))
          end if
        end if
c
c
*** introduce a single A00 for later multiplication purposes, which is stored
*** in scr5(1,1,scrno1-1); all following multiplications start from A00 as the
*** left hand side factor and consider later eventually missing terms on the left
*** of A00.
*** Update 'term', and copy rest of term 'idum1' positions to the left
        term(ileft:ileft+2)='A00'
        iposA=ileft
        idum1=(iright-ileft+1)-3
        do 70 i=ileft+3,length-idum1
          term(i:i)=term(i+idum1:i+idum1)
  70    continue
        do 75 i=length-idum1+1,length
          term(i:i)=' '
  75    continue
        length=length-idum1
        flag=.false.
c
        goto 3150
c
      else
        write(stdout,*) "MIRACLE in evalstring2b() !"
        call Abend
      endif
c
      return
      end
c
c
