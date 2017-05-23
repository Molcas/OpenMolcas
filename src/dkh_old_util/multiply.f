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
* Copyright (C) 2004, Alexander Wolf                                   *
*               2004,2007, Markus Reiher                               *
************************************************************************
      subroutine multiply (length,term,coeff,nbas,posu,post,poss,vv,nn,
     *                     dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,
     *                     tnumber,unumber,scr1,scr2,scr3,A,C,nbasp,
     *                     nbaso,dkhadr,adrmem,dkh_48,adrnext)
c
c************************************************************************
c
c   Execute matrix multiplications for string 'term' of length 'length'
c     with coefficient 'coeff':
c
c      C = coeff * (A * B_1 * B_2 * B_3 ...)
c
c     (the B_i are determined in this routine according to the
c      entries in term(:))
c
c   Note:  * There are no brackets occurring in term, i.e, only Sxxx, Txxx,
c              Uxxx, Axx, V,N,D,Y,F,G,Z,Q,X,I,J,K,L,M  have to be dealt with.
c          * The result is stored in C.
c          * 'term' might contain only 1 factor.
c          * 'term' is destroyed by this procedure.
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.1.0
c
c   last modified: 13.03.2007 (M. Reiher, ETH Zurich)
c                  * latest update for IO to disk
c                  * the new version includes calls to a new version of
c                    dertermine_factor() in order to avoid unnecessary matrix
c                    copy calls
c                  * other: significant changes compared to previous versions
c                    in order to make it much more efficient
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer length,nbas,posu(maxunumber),post(maxsnumber),nbasp,
     *        poss(maxsnumber),snumber,tnumber,unumber,nbaso
      character*(*) term
      REAL*8 coeff
c
      REAL*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp)
      REAL*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),
     *                 A(nbas,nbas),C(nbas,nbas)
c
      integer iact
c
      integer dkh_48,adrmem,adrnext
      integer dkhadr(adrmem)
c
c-------------------------------------------------------------------------
c
      iact=1
      call determine_factor (length,term,iact,nbas,posu,post,poss,vv,nn,
     *                       dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,scr3,
     *                       A,nbasp,nbaso,dkhadr,adrmem,dkh_48,adrnext)
*** now: left factor is in A
c
*** check if term() contains only a single matrix and return
      if (iact.gt.length) then
        call mat_copy_c(C,nbas,nbas,A,coeff)
        goto 2000
      end if
c
1000  continue
c
*** determine right factor and multiply -> result is then in C
      call determine_factor2(length,term,iact,nbas,posu,post,poss,vv,nn,
     *                       dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,scr3,A,
     *                       C,nbasp,nbaso,dkhadr,adrmem,dkh_48,adrnext)
c
      if (iact.gt.length) then
*** last factor found: multiply and return
        call mat_copy_c(C,nbas,nbas,C,coeff)
        goto 2000
      else
*** multiply without additional copying matrices
*** now: left factor is in C
        call determine_factor2(length,term,iact,nbas,posu,post,poss,vv,
     *                       nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,scr3,C,
     *                       A,nbasp,nbaso,dkhadr,adrmem,dkh_48,adrnext)
        if (iact.gt.length) then
          call mat_copy_c(C,nbas,nbas,A,coeff)
          goto 2000
        else
          goto 1000
        end if
      end if
c
2000  return
      end
c
c
c
c
      subroutine multiply2(length,term,coeff,nbas,iposA,posu,post,poss,
     *                     vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                     snumber,tnumber,unumber,scr1,scr2,scr3,A,C,
     *                     nbasp,nbaso,dkhadr,adrmem,dkh_48,adrnext)
c
c************************************************************************
c
c   Execute matrix multiplications for string 'term' of length 'length'
c     with coefficient 'coeff':
c
c       C = coeff * (A * B_1 * B_2 * B_3 ...)
c
c     (the B_i are determined in this routine according to the
c      entries in term(:))
c
c   Note:  * There are no brackets occurring in term, i.e, only Sxxx, Txxx,
c              Uxxx, Axx, V,N,D,Y,F,G,Z,Q,X,I,J,K,L,M  have to be dealt with.
c          * The result is stored in C.
c          * 'term' might contain only 1 factor.
c          * 'term' is destroyed by this procedure.
c          * most important: this routine obtains the left-hand-side factor
c            on input left factor is in A
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Markus Reiher  (ETH Zurich)
c
c   version:  2.2.0
c
c   modified: 13.03.07 MR@ETH
c   first version: 14.01.2007  (ETH Zurich)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer length,nbas,posu(maxunumber),post(maxsnumber),nbasp,
     *        poss(maxsnumber),snumber,tnumber,unumber,iposA,nbaso
      character*(*) term
      REAL*8 coeff
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
      integer iact
c
      integer dkh_48,adrmem,adrnext
      integer dkhadr(adrmem)
c
c-------------------------------------------------------------------------
c
CMR      iact=iposA
CMR This sub assumes that "A00" starts at position 1 in term !
      iact=1
*** this routine gets A00 as a first factor so that the next factor starts at iact=4
      iact=iact+3
c
*** check if term() contains only a single matrix
      if (iact.gt.length) then
        call mat_copy_c(C,nbas,nbas,A,coeff)
        goto 2000
      end if
c
1000  continue
c
*** determine right factor and multiply -> result is then in C
      call determine_factor2(length,term,iact,nbas,posu,post,poss,vv,nn,
     *                       dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,scr3,A,
     *                       C,nbasp,nbaso,dkhadr,adrmem,dkh_48,adrnext)
c
      if (iact.gt.length) then
*** last factor found: multiply and return
        call mat_copy_c(C,nbas,nbas,C,coeff)
        goto 2000
      else
*** multiply without additional copying matrices
*** now: left factor is in C
        call determine_factor2(length,term,iact,nbas,posu,post,poss,vv,
     *                       nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,scr3,C,
     *                       A,nbasp,nbaso,dkhadr,adrmem,dkh_48,adrnext)
        if (iact.gt.length) then
          call mat_copy_c(C,nbas,nbas,A,coeff)
          goto 2000
        else
          goto 1000
        end if
      end if
c
2000  return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(iposA)
      end
c
c
c
c
      subroutine multiply3(length,term,coeff,nbas,posu,post,poss,vv,nn,
     *                     dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,snumber,
     *                     tnumber,unumber,scr1,scr2,scr3,A,C,nbasp,
     *                     nbaso,dkhadr,adrmem,dkh_48,adrnext)
c
c************************************************************************
c
c   Execute matrix multiplications for string 'term' of length 'length'
c     with coefficient 'coeff':
c
c       C = coeff * (... B_3 * B_2 * B_1 * A)
c
c     (the B_i are determined in this routine according to the
c      entries in term(:))
c
c   Note:  * There are no brackets occurring in term, i.e, only Sxxx, Txxx,
c              Uxxx, Axx, V,N,D,Y,F,G,Z,Q,X,I,J,K,L,M  have to be dealt with.
c          * The result is stored in C.
c          * 'term' might contain only 1 factor.
c          * 'term' is destroyed by this procedure.
c          * most important: this routine obtains the left-hand-side factor
c            on input left factor is in A
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Markus Reiher  (ETH Zurich)
c
c   version:  2.2.0
c
c   modified: 13.03.2007 MR@ETH
c   first version: 14.01.2007  (ETH Zurich)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer length,nbas,posu(maxunumber),post(maxsnumber),nbasp,
     *        poss(maxsnumber),snumber,tnumber,unumber,nbaso
      character*(*) term
      REAL*8 coeff
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
      integer iact
c
      integer dkh_48,adrmem,adrnext
      integer dkhadr(adrmem)
c
c-------------------------------------------------------------------------
c
*** this routine gets A00 as a first right-hand-side factor so that the next
***   factor starts at length-1
      iact=length
c
*** check if term() contains only a single matrix
      if (iact.lt.1) then
        write(6,*) "CMR: SHOULD NEVER OCCUR HERE !"
        call mat_copy_c(C,nbas,nbas,A,coeff)
        goto 2000
      end if
c
1000  continue
c
*** determine left factor and multiply -> result is then in C
      call determine_factor3(length,term,iact,nbas,posu,post,poss,vv,nn,
     *                       dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,scr3,A,
     *                       C,nbasp,nbaso,dkhadr,adrmem,dkh_48,adrnext)
c
      if (iact.lt.1) then
*** last factor found: multiply and return
        call mat_copy_c(C,nbas,nbas,C,coeff)
        goto 2000
      else
*** multiply without additional copying matrices
*** now: left factor is in C
        call determine_factor3(length,term,iact,nbas,posu,post,poss,vv,
     *                       nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,scr3,C,
     *                       A,nbasp,nbaso,dkhadr,adrmem,dkh_48,adrnext)
        if (iact.lt.1) then
          call mat_copy_c(C,nbas,nbas,A,coeff)
          goto 2000
        else
          goto 1000
        end if
      end if
c
2000  return
      end
c
c
