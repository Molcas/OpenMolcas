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
* Copyright (C) 2004,2005, Alexander Wolf                              *
*               2004,2005, Markus Reiher                               *
************************************************************************
      subroutine replace2(length,coeff,operator,wordercounter,wopsleng,
     *                    dwops,wops,sused,stimes,wstimes,s,scrchar,
     *                    scrleng,ttimes,termcounter,termleng,termleng2,
     *                    dtcoeff,dtcoeff2,
     *                    term)
c
c******************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 12.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer length
      character*(maxlength) operator
      REAL*8 coeff
c
      integer wordercounter(maxorder),wopsleng(maxuops),sused,
     *        stimes(maxsnumber),wstimes(maxuops,maxsnumber),
     *        scrleng(maxsnumber),ttimes(maxsnumber)
      REAL*8 dwops(maxuops)
      character*(maxlength) wops(maxuops)
      character*(4) s(maxsnumber)
      character*(9) scrchar(maxsnumber)
c
      integer termcounter,termleng(maxoperators),
     *        termcounter2,termleng2(maxoperators)
      REAL*8 dtcoeff(maxoperators),dtcoeff2(maxoperators)
#if defined(_MOLCAS_) || defined(MOLPRO)
#include "WrkSpc.fh"
      character*(maxlength) termstr
      integer term(*)
      integer term2,lwop,nwop,intrea
#else
      character*(maxlength) term(maxoperators),term2(maxoperators)
Cfrueher_benutzt      common/cdkhops2/ term2
#endif
c
C     integer i,k,l,m,pos,posw,poso,pose,reslengl,reslengr,dummyleng,
C    *        tmpleng,hit,istart,iposs,
C    *        idum,idum1,idum2,wstart,dkh_char2int
      integer i,k,    pos,posw,          reslengl,reslengr,dummyleng,
     *                hit,istart,
     *        idum,idum1,      wstart,dkh_char2int
      character*(maxlength) rescharl,rescharr,dummychar
c
      dummyleng=0
      do 10 k=1,maxlength
        dummychar(k:k)=' '
  10  continue
#if defined(_MOLCAS_) || defined(MOLPRO)
      lwop=8/intrea()
      nwop=(maxlength-1)/lwop+1
      call getmem('term2','Allo','Inte',term2,maxoperators*nwop)
      termstr=' '
#endif
      termcounter=1
      termcounter2=0
      do 20 i=1,maxoperators
        termleng(i)=0
        dtcoeff(i)=0.0d0
        termleng2(i)=0
        dtcoeff2(i)=0.0d0
#if defined(_MOLCAS_) || defined(MOLPRO)
        call put_dkoperators(i,termstr,iwork(term2))
        call put_dkoperators(i,termstr,term)
#else
        do 30 k=1,maxlength
          term(i)(k:k)=' '
          term2(i)(k:k)=' '
  30    continue
#endif
  20  continue
      termleng(1)=length
      dtcoeff(1)=coeff
#if defined(_MOLCAS_) || defined(MOLPRO)
      call put_dkoperators(1,operator,term)
#else
      term(1)(1:termleng(1))=operator(1:length)
#endif
2050  continue
      termcounter2=0
      hit=0
c
      do 200 i=1,termcounter
        dummyleng=2
        pos=0
        posw=0
        idum=0
        istart=1
2063    continue
#if defined(_MOLCAS_) || defined(MOLPRO)
        call get_dkoperators(i,termstr,term)
        posw=index(termstr(istart:termleng(i)),'W')
#else
        posw=index(term(i)(istart:termleng(i)),'W')
#endif
        if (posw.gt.0) then
#if defined(_MOLCAS_) || defined(MOLPRO)
          idum1=dkh_char2int(dummyleng,
     *                       termstr(istart+posw:istart+posw+1))
#else
          idum1=dkh_char2int(dummyleng,
     *                       term(i)(istart+posw:istart+posw+1))
#endif
          if (idum1.gt.idum) then
            idum=idum1
            pos=istart+posw-1
          endif
          istart=istart+pos+2
          goto 2063
        endif
        if (pos.eq.0) then
c  term(i) does not contain a 'W':  copy term --> term2 without any modification
          termcounter2=termcounter2+1
          termleng2(termcounter2)=termleng(i)
          dtcoeff2(termcounter2)=dtcoeff(i)
#if defined(_MOLCAS_) || defined(MOLPRO)
          call put_dkoperators(termcounter2,termstr,iwork(term2))
#else
          term2(termcounter2)(1:termleng2(termcounter2))
     *      =term(i)(1:termleng(i))
#endif
c
        else
c
          hit=hit+1
          wstart=1
          do 341 k=1,idum-1
            wstart=wstart+wordercounter(k)
 341      continue
          dummyleng=3
#if defined(_MOLCAS_) || defined(MOLPRO)
          call store_reschar (termleng(i),termstr,pos,dummyleng,
     *                        reslengl,reslengr,rescharl,rescharr)
#else
          call store_reschar (termleng(i),term(i),pos,dummyleng,
     *                        reslengl,reslengr,rescharl,rescharr)
#endif
          do 351 k=1,wordercounter(idum)
            termcounter2=termcounter2+1
            termleng2(termcounter2)=termleng(i)-3+wopsleng(k+wstart-1)
            if (termleng2(termcounter2).gt.maxlength) then
              write (stdout,3034) termleng2(termcounter2)
3034          format (/2X,'ERROR3 in subroutine "replace2": termleng =',
     *              I3,' is larger than maxlength.',//2X,'Increase ',
     *              'parameter maxlength in "dkhparameters.fh".',//2X,
     *              'STOP.',/)
              CALL Abend
            endif
            dtcoeff2(termcounter2)=dtcoeff(i)*dwops(k+wstart-1)
#if defined(_MOLCAS_) || defined(MOLPRO)
            call concatenate (termleng2(termcounter2),
     *                        termstr,reslengl,rescharl,
     *                        wopsleng(k+wstart-1),
     *                        wops(k+wstart-1)(1:wopsleng(k+wstart-1)),
     *                        reslengr,rescharr)
            call put_dkoperators(termcounter2,termstr,iwork(term2))
#else
            call concatenate (termleng2(termcounter2),
     *                        term2(termcounter2),reslengl,rescharl,
     *                        wopsleng(k+wstart-1),
     *                        wops(k+wstart-1)(1:wopsleng(k+wstart-1)),
     *                        reslengr,rescharr)
#endif
            call calc_stimes1 (k,wstart,sused,stimes,wstimes,ttimes,
     *                         reslengl,reslengr,rescharl,rescharr)
c
 351      continue
        endif
 200  continue
      do 400 i=1,termcounter2
        termleng(i)=termleng2(i)
        dtcoeff(i)=dtcoeff2(i)
#if defined(_MOLCAS_) || defined(MOLPRO)
        call copy_dkoperators(i,iwork(term2),i,term)
#else
        term(i)(1:termleng(i))=term2(i)(1:termleng2(i))
#endif
 400  continue
      termcounter=termcounter2
c
      if (hit.gt.0) goto 2050
2060  continue
      termcounter2=0
      hit=0
      do 1200 i=1,termcounter
        pos=0
#if defined(_MOLCAS_) || defined(MOLPRO)
        call get_dkoperators(i,termstr,term)
        pos=index(termstr(1:termleng(i)),'O01')
#else
        pos=index(term(i)(1:termleng(i)),'O01')
#endif
c
        if (pos.eq.0) then
          termcounter2=termcounter2+1
          termleng2(termcounter2)=termleng(i)
          dtcoeff2(termcounter2)=dtcoeff(i)
#if defined(_MOLCAS_) || defined(MOLPRO)
          call copy_dkoperators(i,term,termcounter2,iwork(term2))
#else
          term2(termcounter2)(1:termleng2(termcounter2))
     *      =term(i)(1:termleng(i))
#endif
        else
          hit=hit+1
          dummyleng=3
#if defined(_MOLCAS_) || defined(MOLPRO)
          call store_reschar (termleng(i),termstr,pos,dummyleng,
     *                        reslengl,reslengr,rescharl,rescharr)
#else
          call store_reschar (termleng(i),term(i),pos,dummyleng,
     *                        reslengl,reslengr,rescharl,rescharr)
#endif
          termcounter2=termcounter2+1
          termleng2(termcounter2)=termleng(i)
          dtcoeff2(termcounter2)=dtcoeff(i)
          dummyleng=3
          dummychar(1:dummyleng)='BPV'
#if defined(_MOLCAS_) || defined(MOLPRO)
          call concatenate (termleng2(termcounter2),termstr,
     *                      reslengl,rescharl,dummyleng,
     *                      dummychar(1:dummyleng),reslengr,rescharr)
          call put_dkoperators(termcounter2,termstr,iwork(term2))
#else
          call concatenate (termleng2(termcounter2),term2(termcounter2),
     *                      reslengl,rescharl,dummyleng,
     *                      dummychar(1:dummyleng),reslengr,rescharr)
#endif
          termcounter2=termcounter2+1
          termleng2(termcounter2)=termleng(i)
          dtcoeff2(termcounter2)=-1.0d0*dtcoeff(i)
          dummychar(1:dummyleng)='BVP'
#if defined(_MOLCAS_) || defined(MOLPRO)
          call concatenate (termleng2(termcounter2),termstr,
     *                      reslengl,rescharl,dummyleng,
     *                      dummychar(1:dummyleng),reslengr,rescharr)
          call put_dkoperators(termcounter2,termstr,iwork(term2))
#else
          call concatenate (termleng2(termcounter2),term2(termcounter2),
     *                      reslengl,rescharl,dummyleng,
     *                      dummychar(1:dummyleng),reslengr,rescharr)
#endif
          call calc_stimes2 (stimes,ttimes,reslengl,reslengr,rescharl,
     *                       rescharr)
c
        endif
1200  continue
      do 1400 i=1,termcounter2
        termleng(i)=termleng2(i)
        dtcoeff(i)=dtcoeff2(i)
#if defined(_MOLCAS_) || defined(MOLPRO)
        call copy_dkoperators(i,iwork(term2),i,term)
#else
        term(i)(1:termleng(i))=term2(i)(1:termleng2(i))
#endif
1400  continue
      termcounter=termcounter2
c
      if (hit.gt.0) goto 2060
c
c
3060  continue
      termcounter2=0
      hit=0
      do 1270 i=1,termcounter
        pos=0
#if defined(_MOLCAS_) || defined(MOLPRO)
        call get_dkoperators(i,termstr,term)
        pos=index(termstr(1:termleng(i)),'CO0')
#else
        pos=index(term(i)(1:termleng(i)),'CO0')
#endif
c
        if (pos.eq.0) then
          termcounter2=termcounter2+1
          termleng2(termcounter2)=termleng(i)
          dtcoeff2(termcounter2)=dtcoeff(i)
#if defined(_MOLCAS_) || defined(MOLPRO)
          call put_dkoperators(termcounter2,termstr,iwork(term2))
#else
          term2(termcounter2)(1:termleng2(termcounter2))
     *      =term(i)(1:termleng(i))
#endif
        else
          hit=hit+1
          dummyleng=3
#if defined(_MOLCAS_) || defined(MOLPRO)
          call store_reschar (termleng(i),termstr,pos,dummyleng,
     *                        reslengl,reslengr,rescharl,rescharr)
#else
          call store_reschar (termleng(i),term(i),pos,dummyleng,
     *                        reslengl,reslengr,rescharl,rescharr)
#endif
          termcounter2=termcounter2+1
          termleng2(termcounter2)=termleng(i)
          dtcoeff2(termcounter2)=dtcoeff(i)
          dummyleng=3
          dummychar(1:dummyleng)='BPX'
#if defined(_MOLCAS_) || defined(MOLPRO)
          call concatenate (termleng2(termcounter2),termstr,
     *                      reslengl,rescharl,dummyleng,
     *                      dummychar(1:dummyleng),reslengr,rescharr)
          call put_dkoperators(termcounter2,termstr,iwork(term2))
#else
          call concatenate (termleng2(termcounter2),term2(termcounter2),
     *                      reslengl,rescharl,dummyleng,
     *                      dummychar(1:dummyleng),reslengr,rescharr)
#endif
          termcounter2=termcounter2+1
          termleng2(termcounter2)=termleng(i)
          dtcoeff2(termcounter2)=-1.0d0*dtcoeff(i)
          dummychar(1:dummyleng)='BXP'
#if defined(_MOLCAS_) || defined(MOLPRO)
          call concatenate (termleng2(termcounter2),termstr,
     *                      reslengl,rescharl,dummyleng,
     *                      dummychar(1:dummyleng),reslengr,rescharr)
          call put_dkoperators(termcounter2,termstr,iwork(term2))
#else
          call concatenate (termleng2(termcounter2),term2(termcounter2),
     *                      reslengl,rescharl,dummyleng,
     *                      dummychar(1:dummyleng),reslengr,rescharr)
#endif
          call calc_stimes2 (stimes,ttimes,reslengl,reslengr,
     *                       rescharl,rescharr)
c
        endif
1270  continue
      do 1470 i=1,termcounter2
        termleng(i)=termleng2(i)
        dtcoeff(i)=dtcoeff2(i)
#if defined(_MOLCAS_) || defined(MOLPRO)
        call copy_dkoperators(i,iwork(term2),i,term)
#else
        term(i)(1:termleng(i))=term2(i)(1:termleng2(i))
#endif
1470  continue
      termcounter=termcounter2
c
      if (hit.gt.0) goto 3060
c
#if defined(_MOLCAS_) || defined(MOLPRO)
      call getmem('term2','Free','Inte',term2,maxoperators*nwop)
#endif
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_character(s)
        call Unused_character(scrchar)
        call Unused_integer_array(scrleng)
      end if
      end
