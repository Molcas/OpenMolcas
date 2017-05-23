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
      subroutine evaluate_Txxx (dkhorder,xorder,paramtype,dkhscfflg,
     *             wordercounter,wopsleng,woporder,eowops,dwops,wops,
     *             termleng,termleng2,termorder,dtcoeff,dtcoeff2,
     *             sused,snumber,scounter,stimes,wstimes,sorder,s,
     *             scrchar,scrleng,stimestot,tnumber,tcounter,ttimes,
     *             tmult,t,tscrchar,tscrleng,ttimestot,uused,unumber,
     *             utimes,uorder,umult,u,uscrchar,uscrleng,utimestot,
     *             umulttot,totmult,dkhzerothrsh)
c
c*************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 19.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c*************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
c
      integer dkhorder,xorder
      character*(3) paramtype
      logical dkhscfflg
c
C     integer wordercounter(maxorder),wopscounter,wopsleng(maxuops),
      integer wordercounter(maxorder),wopsleng(maxuops),
     *        woporder(maxuops,3),eowops(maxuops)
      real*8 dwops(maxuops)
      character*(maxlength) wops(maxuops)
c
      integer sused,snumber,scounter(maxsnumber),stimes(maxsnumber),
     *        wstimes(maxuops,maxsnumber),sorder(maxsnumber,3),
     *        scrleng(maxsnumber),stimestot
      integer tnumber,tcounter(maxsnumber),ttimes(maxsnumber),
     *        tmult(maxsnumber),tscrleng(maxsnumber),ttimestot
      integer uused,unumber,utimes(maxunumber),uorder(maxunumber,3),
     *        umult(maxunumber),uscrleng(maxunumber),utimestot,umulttot
      character*(3) uscrchar(maxunumber)
      character*(4) s(maxsnumber),t(maxsnumber),u(maxunumber)
      character*(9) scrchar(maxsnumber)
      character*(11) tscrchar(maxsnumber)
c
      integer stimes2(maxsnumber),ttimes2(maxsnumber),
     *        utimes2(maxunumber)
c
      integer idum,dummyleng
      integer tmulttot,totmult
      real*8 dkhzerothrsh(0:maxorder)
c
      integer termcounter,termleng(maxoperators),
     *        termleng2(maxoperators),
     *        termorder(maxoperators,3)
      real*8 dtcoeff(maxoperators),dtcoeff2(maxoperators),
     *                 dummycoeff
#if defined(_MOLCAS_) || defined(MOLPRO)
#include "WrkSpc.fh"
      integer term,lwop,nwop,intrea
      character*(maxlength) termstr
#else
      character*(maxlength) term(maxoperators)
#endif
      character*(maxlength) dummychar
c
      integer j,k,l
#ifdef _MOLCAS_
      Integer IsFreeUnit
#endif
CMR      write (stdout,1007)
CMR1007  format (//5X,'|',70('-'),'|',/5X,'|',6X,'STEP 5: Evaluate ',
CMR     *        'auxiliary matrices Txxx = PSxxxP',15X,'|',/5X,'|',
CMR     *        70('-'),'|',/2X)

      if (dbgflg.ge.1) then
        write (dbgunit,1009)
1009    format (/2X,'Step 5: Enter subroutine "evaluate_Txxx" ',//2X,
     *          'Evaluate even-order auxiliary matrices Txxx = PSxxxP',
     *          ', which have been used in step 3.')
      endif
#ifndef MOLPRO
#ifdef _MOLCAS_
      dkhunit4=14
      dkhunit4=IsFreeUnit(dkhunit4)
      Call Molcas_Open(dkhunit4,'dkhops.14')
#else
      open (dkhunit4, file='dkhops.14', status='UNKNOWN',
     *        form='FORMATTED')
#endif
#endif
      rewind(dkhunit4)
      write (dkhunit4,1021) dkhorder,paramtype,xorder,dkhscfflg
1021  format (50('-'),/2X,'dkhorder = ',I2,10X,A3,/2X,'xorder   = ',I2,
     *        /2X,'dkhscfflg = ',L1,/'+++',47('-'))
      write (dkhunit4,'(I3)') tnumber
#if defined(_MOLCAS_) || defined(MOLPRO)
      lwop=8/intrea()
      nwop=(maxlength-1)/lwop+1
      call getmem('term','Allo','Inte',term,maxoperators*nwop)
#endif
      idum=0
      do 10 j=1,sused
        if (ttimes(j).gt.0) then
          idum=idum+1
          do 20 k=1,maxlength
            dummychar(k:k)=' '
  20      continue
          if (tscrleng(j).eq.8) dummyleng=8
          if (tscrleng(j).eq.11) dummyleng=11
          dummycoeff=1.0d0
          dummychar(1:dummyleng)=tscrchar(j)
#if defined(_MOLCAS_) || defined(MOLPRO)
          call replace3 (dummyleng,dummycoeff,dummychar,wordercounter,
     *                   wopsleng,dwops,wops,scrchar,scrleng,ttimes,t,
     *                   tscrchar,tscrleng,termcounter,termleng,
     *                   termleng2,dtcoeff,dtcoeff2,
     *                   iwork(term))
#else
          call replace3 (dummyleng,dummycoeff,dummychar,wordercounter,
     *                   wopsleng,dwops,wops,scrchar,scrleng,ttimes,t,
     *                   tscrchar,tscrleng,termcounter,termleng,
     *                   termleng2,dtcoeff,dtcoeff2,
     *                   term)
#endif
          tcounter(j)=termcounter
          if (tcounter(j).gt.maxoperators) then
            write (stdout,1033) j,tcounter(j)
1033        format (/2X,'ERROR in SR "evaluate_Txxx": tcounter(',I3,
     *            ') = ',I9,' > maxoperators.',//2X,'Increase ',
     *            'maxoperators in "dkhparameters.fh".',//2X,'STOP.',/)
            CALL Abend
          endif
1042      continue
          do 30 k=1,termcounter
c
c  Initialize auxiliary counters stimes2, ttimes2, utimes2
            do 36 l=1,sused
              stimes2(l)=0
              ttimes2(l)=0
  36        continue
            do 37 l=1,uused
              utimes2(l)=0
  37        continue
#if defined(_MOLCAS_) || defined(MOLPRO)
            call get_dkoperators(k,termstr,iwork(term))
            call removeB1 (termleng(k),dtcoeff(k),termstr)
            call removeB2 (termleng(k),termstr)
            call insert_ri (termleng(k),termstr)
            call movebraces (termleng(k),termstr)
            call finalize (termleng(k),termstr,stimes2,ttimes2,t)
            call finalize2 (termleng(k),termstr,uused,utimes2,u,
     *                        uscrleng,uscrchar)
            termorder(k,3)=sorder(j,3)
            call calc_orders (dkhscfflg,termleng(k),termorder(k,1),
     *                        termorder(k,2),termorder(k,3),termstr,
     *                        sorder,uorder)
            call put_dkoperators(k,termstr,iwork(term))
#else
            call removeB1 (termleng(k),dtcoeff(k),term(k))
            call removeB2 (termleng(k),term(k))
            call insert_ri (termleng(k),term(k))
            call movebraces (termleng(k),term(k))
            call finalize (termleng(k),term(k),stimes2,ttimes2,t)
            call finalize2 (termleng(k),term(k),uused,utimes2,u,
     *                        uscrleng,uscrchar)
              termorder(k,3)=sorder(j,3)
              call calc_orders (dkhscfflg,termleng(k),termorder(k,1),
     *                          termorder(k,2),termorder(k,3),term(k),
     *                          sorder,uorder)
#endif
c
CMR ISN'T IT POSSIBLE TO DO THIS PRIOR TO ALL THESE CALCULATIONS ?? -> CHECK AGAIN!
            if ( (dynthrsh.eq.0.and.abs(dtcoeff(k)).lt.dkhzero) .or.
     *           (dynthrsh.eq.1.and.
     *             abs(dtcoeff(k)).lt.dkhzerothrsh(sorder(j,3))) .or.
     *           (dkhscfflg.and.termorder(k,2).gt.xorder) ) then
              do 40 l=k,termcounter-1
                termleng(l)=termleng(l+1)
                dtcoeff(l)=dtcoeff(l+1)
#if defined(_MOLCAS_) || defined(MOLPRO)
CMR s. evaluate_Sxxx.f an dieser Stelle
                call copy_dkoperators(l+1,iwork(term),l,iwork(term))
#else
                term(l)=term(l+1)
#endif
  40          continue
              termcounter=termcounter-1
              tcounter(j)=tcounter(j)-1
CMR WHY SHOULD ONE START WITH 1 AGAIN INSTEAD OF K IN LOOP 30: k=k_actual, termcounter ???
              goto 1042
            else
              do 47 l=1,sused
                stimes(l)=stimes(l)+stimes2(l)
                ttimes(l)=ttimes(l)+ttimes2(l)
  47          continue
              do 48 l=1,uused
                utimes(l)=utimes(l)+utimes2(l)
  48          continue
            endif
c
  30      continue
c
#if defined(_MOLCAS_) || defined(MOLPRO)
          call calc_smult (j,tmult(j),termcounter,termleng,iwork(term))
          call output12 (dkhunit4,tcounter(j),sorder(j,1),sorder(j,2),
     *                   sorder(j,3),t(j),
     *                   tscrchar(j),termcounter,termleng,termorder,
     &                   dtcoeff,iwork(term))
#else
          call calc_smult (j,tmult(j),termcounter,termleng,term)
          call output12 (dkhunit4,tcounter(j),sorder(j,1),sorder(j,2),
     *                   sorder(j,3),t(j),
     *                   tscrchar(j),termcounter,termleng,dtcoeff,term)
#endif
c
CMR          write (stdout,2045) idum,t(j),tcounter(j)
CMR2045      format (2X,I3,4X,A4,'  with ',I7,' low-level lines ',
CMR     *            'completed.')
        endif
  10  continue
c
#ifdef _MOLCAS_
      close (dkhunit4)
#endif
c
      tmulttot=0
      do 220 j=1,sused
        if (ttimes(j).gt.0) then
          tmulttot=tmulttot+tmult(j)
        endif
 220  continue
c
      unumber=0
      utimestot=0
      umulttot=0
      do 225 j=1,uused
        if (utimes(j).ge.1) then
          unumber=unumber+1
          utimestot=utimestot+utimes(j)
          umulttot=umulttot+umult(j)
        endif
 225  continue
c
      totmult=totmult+tmulttot
      if (dbgflg.ge.1) then
        write (dbgunit,4026) tnumber
4026    format (/2X,'All ',I3,' auxiliary matrices Txxx have been ',
     *          'reduced to their low-level form.')
        call output7 (dbgunit,dkhscfflg,sused,tcounter,ttimes,sorder,
     *                tmult,t,tscrchar,ttimestot,tmulttot)
        write (dbgunit,4027) tnumber
4027    format (2X,'Final low-level expressions for all ',I3,' Txxx',
     *          ' have been written to the file "dkhops.14".'/2X)
        call output8 (dbgunit,dkhscfflg,uused,unumber,utimes,uorder,
     *                umult,u,uscrchar,uscrleng,utimestot,umulttot)
        write (dbgunit,2769)
2769    format (/2X,'Subroutine evaluate_Txxx (step5) has now been ',
     *          'completed.',//2X,145('*'))
      endif
CMR      write (stdout,1034) tnumber,tnumber,tmulttot
CMR1034  format (/2X,'All ',I3,' high-level expressions Txxx have been ',
CMR     *        'reduced to their low-level',/2X,'form. The results ',
CMR     *         'are stored in the file "dkhops.14".',
CMR     *        //2X,'No. of low-level matrix multiplications needed',
CMR     *        /4X,'for the evaluation of all',I3,' Txxx-operators:',
CMR     *        12X,I10)
c
#if defined(_MOLCAS_) || defined(MOLPRO)
      call getmem('term','Free','Inte',term,maxoperators*nwop)
#endif
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer_array(woporder)
        call Unused_integer_array(eowops)
        call Unused_integer(snumber)
        call Unused_integer_array(scounter)
        call Unused_integer_array(wstimes)
        call Unused_character(s)
        call Unused_integer(stimestot)
      end if
      end
