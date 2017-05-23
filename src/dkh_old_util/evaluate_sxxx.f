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
      subroutine evaluate_Sxxx (dkhorder,xorder,paramtype,dkhscfflg,
     *             wordercounter,wopsleng,woporder,eowops,dwops,wops,
     *             termleng,termleng2,termorder,dtcoeff,dtcoeff2,
     *             sused,snumber,scounter,stimes,wstimes,sorder,smult,s,
     *             scrchar,scrleng,stimestot,tnumber,tcounter,ttimes,
     *             tmult,t,tscrchar,tscrleng,ttimestot,uused,unumber,
     *             utimes,uorder,umult,u,uscrchar,uscrleng,totmult,
     *             dkhzerothrsh)
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
      integer dkhorder,xorder
      character*(3) paramtype
      logical dkhscfflg
c
C     integer wordercounter(maxorder),wopscounter,wopsleng(maxuops),
      integer wordercounter(maxorder),            wopsleng(maxuops),
     *        woporder(maxuops,3),eowops(maxuops)
      real*8 dwops(maxuops)
      character*(maxlength) wops(maxuops)
c
      integer sused,snumber,scounter(maxsnumber),stimes(maxsnumber),
     *        wstimes(maxuops,maxsnumber),sorder(maxsnumber,3),
     *        smult(maxsnumber),scrleng(maxsnumber),stimestot
      integer tnumber,tcounter(maxsnumber),ttimes(maxsnumber),
     *        tmult(maxsnumber),tscrleng(maxsnumber),ttimestot
      integer uused,unumber,utimes(maxunumber),uorder(maxunumber,3),
     *        umult(maxunumber),uscrleng(maxunumber)
      character*(3) uscrchar(maxunumber)
      character*(4) s(maxsnumber),t(maxsnumber),u(maxunumber)
      character*(9) scrchar(maxsnumber)
      character*(11) tscrchar(maxsnumber)
c
      integer stimes2(maxsnumber),ttimes2(maxsnumber),
     *        utimes2(maxunumber)
c
      integer idum
      integer dummyleng,smulttot
      integer totmult
      real*8 dkhzerothrsh(0:maxorder)
      integer termcounter,termleng(maxoperators),
     *        termleng2(maxoperators),
     *        termorder(maxoperators,3)
      real*8 dtcoeff(maxoperators),dtcoeff2(maxoperators),
     *                 dummycoeff
#if defined(_MOLCAS_) || defined(MOLPRO)
#include "WrkSpc.fh"
      integer term,nwop,lwop,intrea
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
CMR      write (stdout,1017)
CMR1017  format (//5X,'|',70('-'),'|',/5X,'|',6X,'STEP 4: Evaluate ',
CMR     *        'even-order auxiliary matrices Sxxx',13X,'|',/5X,'|',
CMR     *        70('-'),'|',/)

      if (dbgflg.ge.1) then
        write (dbgunit,1019)
1019    format (/2X,'Step 4: Enter subroutine "evaluate_Sxxx" ',//2X,
     *          'Evaluate even-order auxiliary matrices Sxxx, which ',
     *          'have been used in step 2.')
      endif
#ifndef MOLPRO
#ifdef _MOLCAS_
      dkhunit3=13
      dkhunit3=IsFreeUnit(dkhunit3)
      Call Molcas_Open(dkhunit3,'dkhops.13')
#else
      open (dkhunit3, file='dkhops.13', status='UNKNOWN',
     *        form='FORMATTED')
#endif
#endif
      rewind(dkhunit3)
      write (dkhunit3,1021) dkhorder,paramtype,xorder,dkhscfflg
1021  format (50('-'),/2X,'dkhorder = ',I2,10X,A3,/2X,'xorder   = ',I2,
     *        /2X,'dkhscfflg = ',L1,/'+++',47('-'))
      write (dkhunit3,'(I3)') snumber
#if defined(_MOLCAS_) || defined(MOLPRO)
      lwop=8/intrea()
      nwop=(maxlength-1)/lwop+1
      call getmem('term_Sxxx','Allo','Inte',term,maxoperators*nwop)
#endif
      idum=0
      do 10 j=1,sused
        smult(j)=0
        if (stimes(j).gt.0) then
          idum=idum+1
          do 20 k=1,maxlength
            dummychar(k:k)=' '
  20      continue
          if (scrleng(j).eq.6) dummyleng=6
          if (scrleng(j).eq.9) dummyleng=9
          dummycoeff=1.0d0
          dummychar(1:dummyleng)=scrchar(j)
c
#if defined(_MOLCAS_) || defined(MOLPRO)
          call replace2 (dummyleng,dummycoeff,dummychar,wordercounter,
     *                   wopsleng,dwops,wops,sused,stimes,wstimes,s,
     *                   scrchar,scrleng,ttimes,termcounter,
     *                   termleng,termleng2,
     *                   dtcoeff,dtcoeff2,iwork(term))
#else
          call replace2 (dummyleng,dummycoeff,dummychar,wordercounter,
     *                   wopsleng,dwops,wops,sused,stimes,wstimes,s,
     *                   scrchar,scrleng,ttimes,termcounter,
     *                   termleng,termleng2,
     *                   dtcoeff,dtcoeff2,term)
#endif
c
          scounter(j)=termcounter
          if (scounter(j).gt.maxoperators) then
            write (stdout,1033) j,scounter(j)
1033        format (/2X,'ERROR in SR "evaluate_Sxxx": scounter(',I3,
     *            ') = ',I9,' > maxoperators.',//2X,'Increase ',
     *            'maxoperators in "dkhparameters.fh".',//2X,'STOP.',/)
            CALL Abend
          endif
1042      continue
          do 30 k=1,termcounter
c
            do 36 l=1,sused
              stimes2(l)=0
              ttimes2(l)=0
  36        continue
            do 37 l=1,uused
              utimes2(l)=0
  37        continue
c
#if defined(_MOLCAS_) || defined(MOLPRO)
            call get_dkoperators(k,termstr,iwork(term))
            call removeB1 (termleng(k),dtcoeff(k),termstr)
            call removeB2 (termleng(k),termstr)
            call insert_ri (termleng(k),termstr)
            call movebraces (termleng(k),termstr)
            call finalize (termleng(k),termstr,stimes2,ttimes2,t)
            call finalize2 (termleng(k),termstr,uused,utimes2,u,
     *                      uscrleng,uscrchar)
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
     *                      uscrleng,uscrchar)
            termorder(k,3)=sorder(j,3)
            call calc_orders (dkhscfflg,termleng(k),termorder(k,1),
     *                        termorder(k,2),termorder(k,3),term(k),
     *                        sorder,uorder)
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
CMR das funktioniert hier so vielleicht nicht, weil iwork(term) nicht da ist, sondern term
                call copy_dkoperators(l+1,iwork(term),l,iwork(term))
#else
                term(l)=term(l+1)
#endif
  40          continue
              termcounter=termcounter-1
              scounter(j)=scounter(j)-1
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
#if defined(_MOLCAS_) || defined(MOLPRO)
          call calc_smult (j,smult(j),termcounter,termleng,iwork(term))
          call output11 (dkhunit3,scounter(j),sorder(j,1),sorder(j,2),
     *                   sorder(j,3),s(j),scrchar(j),termcounter,
     *                   termleng,termorder,dtcoeff,iwork(term))
#else
          call calc_smult (j,smult(j),termcounter,termleng,term)
          call output11 (dkhunit3,scounter(j),sorder(j,1),sorder(j,2),
     *                   sorder(j,3),s(j),scrchar(j),termcounter,
     *                   termleng,termorder,dtcoeff,term)
#endif

CMR          write (stdout,2045) idum,s(j),scounter(j)
CMR2045      format (2X,I3,4X,A4,'  with ',I7,' low-level lines ',
CMR     *            'completed.')
        endif
  10  continue
c
#ifdef _MOLCAS_
      close (dkhunit3)
#endif
      tnumber=0
      do 50 j=1,sused
        if (ttimes(j).gt.0) then
          tnumber=tnumber+1
          ttimestot=ttimestot+ttimes(j)
        endif
  50  continue
      smulttot=0
      stimestot=0
      do 220 j=1,sused
        if (stimes(j).gt.0) then
          stimestot=stimestot+stimes(j)
          smulttot=smulttot+smult(j)
        endif
 220  continue
c
      totmult=totmult+smulttot
      if (dbgflg.ge.1) then
        write (dbgunit,4026) snumber
4026    format (/2X,'All ',I3,' auxiliary matrices Sxxx have been ',
     *          'reduced to their low-level form.')
        call output6 (dbgunit,dkhscfflg,sused,scounter,stimes,sorder,
     *                smult,s,scrchar,stimestot,smulttot)
        write (dbgunit,4027) snumber
4027    format (2X,'Final low-level expressions for all ',I3,' Sxxx',
     *          ' have been written to the file "dkhops.13".')
        write (dbgunit,4030) tnumber
4030    format (//2X,'Thereby,',I3,' different Txxx-expressions have',
     *          ' been used so far (steps 3 and 4).')
        if (dkhscfflg) then
          write (dbgunit,4072)
4072      format (/5X,'#',7X,'step1',4X,'-->',3X,'step3',3X,'order',3X,
     *          'no. of subs.',3X,/36X,'(tot)',5X,'(ttimes)',/)
        else
          write (dbgunit,4073)
4073      format (/5X,'#',7X,'step1',4X,'-->',3X,'step3',3X,'order',1X,
     *            'order',1X,'order',2X,'no. of subs.',3X,/37X,'(V)',3X,
     *            '(X)',2X,'(tot)',3X,'(ttimes)'/)
        endif
        call output5 (dkhscfflg,dbgunit,sused,ttimes,sorder,t,tscrchar,
     *                ttimestot)
        write (dbgunit,2769)
2769    format (/2X,'Subroutine evaluate_Sxxx (step3) has now been ',
     *          'completed.',//2X,145('*'))
      endif
CMR      write (stdout,1034) snumber
CMR1034  format (/2X,'All ',I3,' high-level expressions Sxxx have been ',
CMR     *        'reduced to their low-level',/2X,'form. The results ',
CMR     *        'are stored in the file "dkhops.13".',/2X)
CMR      write (stdout,1036) tnumber,ttimestot
CMR1036  format (2X,'No. of different Txxx expressions used so far:',19X,
CMR     *        I3,/2X,'No. of high-level subs. (Txxx) used so far:',18X,
CMR     *        I7)
CMR      write (stdout,1038) snumber,smulttot
CMR1038  format (/2X,'No. of low-level matrix multiplications needed',
CMR     *        /4X,'for the evaluation of all ',I3,' Sxxx-operators:',
CMR     *        11X,I10)
c
#if defined(_MOLCAS_) || defined(MOLPRO)
      call getmem('term_Sxxx','Free','Inte',term,maxoperators*nwop)
#endif
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer_array(woporder)
        call Unused_integer_array(eowops)
        call Unused_integer_array(tcounter)
        call Unused_integer_array(tmult)
        call Unused_integer_array(tscrleng)
        call Unused_integer(unumber)
        call Unused_integer_array(umult)
      end if
      end
