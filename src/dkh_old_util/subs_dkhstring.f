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
      subroutine subs_dkhstring (dkhorder,xorder,dkhscfflg,ordercounter,
     *          opcounter,operleng,oporder,evenodd,doperators,operators,
     *          xordercounter,xopcounter,xoperleng,xoporder,xevenodd,
     *          xdoperators,xoperators,wordercounter,wopscounter,
     *          wopsleng,woporder,eowops,dwops,wops,sused,snumber,
     *          scounter,stimes,wstimes,sorder,smult,s,scrchar,scrleng,
     *          stimestot,tnumber,tcounter,ttimes,tmult,t,
     *          tscrchar,tscrleng,ttimestot,no_s)
c
c************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.1.0
c
c   last modified: 15.01.2007 (M. Reiher, ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer dkhorder,xorder
      logical dkhscfflg,no_s
c
      integer ordercounter(0:maxorder),opcounter,operleng(maxoperators),
     *        oporder(maxoperators,3),evenodd(maxoperators),
     *        xordercounter(0:maxorder),xopcounter,
     *        xoperleng(maxoperators),xoporder(maxoperators),
     *        xevenodd(maxoperators)
      REAL*8 doperators(maxoperators),
     *                 xdoperators(maxoperators)
#if defined(_MOLCAS_) || defined(MOLPRO)
      character*(maxlength) opstring,xopstring
      integer operators(*),xoperators(*)
#else
      character*(maxlength) operators(maxoperators),
     *                      xoperators(maxoperators)
#endif
c
      integer wordercounter(maxorder),wopscounter,wopsleng(maxuops),
     *        woporder(maxuops,3),eowops(maxuops)
      REAL*8 dwops(maxuops)
      character*(maxlength) wops(maxuops)
c
      integer sused,snumber,scounter(maxsnumber),stimes(maxsnumber),
     *        wstimes(maxuops,maxsnumber),sorder(maxsnumber,3),
     *        smult(maxsnumber),scrleng(maxsnumber),stimestot
      integer tnumber,tcounter(maxsnumber),ttimes(maxsnumber),
     *        tmult(maxsnumber),tscrleng(maxsnumber),ttimestot
      character*(4) s(maxsnumber),t(maxsnumber)
      character*(9) scrchar(maxsnumber)
      character*(11) tscrchar(maxsnumber)
c
C     integer j,k,l,pos,hits,idum1,idum2,dummyleng,istart
      integer j,k,  pos,hits,idum1,idum2,dummyleng,istart
C     logical again
      character*(maxlength) dummychar
c
CMR      write (stdout,1017)
CMR1017  format (//5X,'|',70('-'),'|',/5X,'|',6X,'STEP 2: Substitution ',
CMR     *        'of even-order high-level expressions',7X,'|',/5X,'|',
CMR     *        70('-'),'|')
c
      if (dbgflg.ge.1) then
        write (dbgunit,1019)
1019    format (/2X,'Step 2: Enter subroutine "subs_dkhstring" ',//2X,
     *          'Start high-level substitution of suitable expressions',
     *          ' occurring in the',/4X,'DKH Hamiltonian (operators),',
     *          ' the property X (xoperators), and the W-operators ',
     *          '(wops).')
      endif
c
      call s_create (dkhscfflg,sused,snumber,scounter,stimes,wstimes,
     *               sorder,smult,s,scrleng,scrchar,stimestot,
     *               wopscounter,no_s)
c
      call t_create (sused,scrleng,scrchar,tnumber,tcounter,ttimes,
     *               tmult,t,tscrleng,tscrchar,ttimestot)
c
      dummyleng=0
      do 20 j=1,maxlength
        dummychar(j:j)=' '
  20  continue
      do 30 k=1,wopscounter
        dummyleng=wopsleng(k)-3
        dummychar(1:dummyleng)=wops(k)(3:wopsleng(k)-1)
        if (index(dummychar(1:dummyleng),'E00').ne.0) then
          write (stdout,1061) k
1061      format (/2X,'ERROR1 in SR "subs_dkhstring": In term no. ',I6,
     *            ' of wops occurs expression "E00".',//2X,
     *            'All "E00" should have been removed at this stage.',
     *            //2X,'STOP.',/)
          CALL Abend
        endif
        hits=0
        pos=1
  50    continue
        if (pos.lt.dummyleng) then
          idum1=index(dummychar(pos:dummyleng),'E01')
          idum2=index(dummychar(pos:dummyleng),'CE0')
          if ((idum2.lt.idum1.and.idum2.gt.0).or.idum1.eq.0) idum1=idum2
          if (idum1.ne.0) then
            hits=hits+1
            pos=pos+idum1+1
            goto 50
          endif
        endif
        if (hits.gt.1)  then
          write (stdout,1073) k
1073      format (/2X,'ERROR2 in SR "subs_dkhstring": In term no. ',I6,
     *            ' of wops occur expressions "E01" or "CE0" more than',
     *            ' once.',//2X,'This is not possible at this stage.',
     *            //2X,'STOP.',/)
          CALL Abend
        endif
        idum1=1
        call substitute_term (idum1,k,hits,dummyleng,dummychar,sused,
     *                        stimes,wstimes,s,scrleng,scrchar)
        wopsleng(k)=dummyleng+3
        wops(k)(3:wopsleng(k)-1)=dummychar(1:dummyleng)
        wops(k)(wopsleng(k):wopsleng(k))=']'
  30  continue
c
      if (dbgflg.ge.1)  then
        write (dbgunit,1233)
1233    format (/2X,'Form of stored W-operators (wops) after ',
     *          'substitution by Sxxx-expressions:',/2X,73('-'))
        call output2 (dbgunit,wopscounter,wopsleng,woporder,eowops,
     *                dwops,wops)
      endif
c
c
      if (dkhscfflg) istart=4
      if (.not.dkhscfflg) istart=3
c
      do 100 k=istart,opcounter
#if defined(_MOLCAS_) || defined(MOLPRO)
        call get_dkoperators(k,opstring,operators)
        if (index(opstring(1:operleng(k)),'E00').ne.0) then
#else
        if (index(operators(k)(1:operleng(k)),'E00').ne.0) then
#endif
          write (stdout,1621) k
1621      format (/2X,'ERROR3 in SR "subs_dkhstring": In term no. ',I6,
     *            ' occurs expression "E00".',//2X,'All "E00" should ',
     *            'have been removed at this stage.',//2X,'STOP.',/)
          CALL Abend
        endif
        hits=0
        pos=1
 122    continue
        if (pos.lt.operleng(k)) then
#if defined(_MOLCAS_) || defined(MOLPRO)
          idum1=index(opstring(pos:operleng(k)),'E01')
          idum2=index(opstring(pos:operleng(k)),'CE0')
#else
          idum1=index(operators(k)(pos:operleng(k)),'E01')
          idum2=index(operators(k)(pos:operleng(k)),'CE0')
#endif
          if ((idum2.lt.idum1.and.idum2.gt.0).or.idum1.eq.0) idum1=idum2
          if (idum1.ne.0) then
            hits=hits+1
            pos=pos+idum1+1
            goto 122
          endif
        endif
        if (hits.gt.1)  then
          write (stdout,1625) k
1625      format (/2X,'ERROR4 in SR "subs_dkhstring": In term no. ',I6,
     *            ' of operators occurs expression "E01" or "CE0" more',
     *            ' than once.',//2X,'This is not possible at this ',
     *            'stage.',//2X,'STOP.',/)
          CALL Abend
        endif
        idum1=2
#if defined(_MOLCAS_) || defined(MOLPRO)
        call substitute_term (idum1,k,hits,operleng(k),opstring,
     *                        sused,stimes,wstimes,s,scrleng,scrchar)
        call Put_dkoperators(k,opstring,operators)
#else
        call substitute_term (idum1,k,hits,operleng(k),operators(k),
     *                        sused,stimes,wstimes,s,scrleng,scrchar)
#endif

 100  continue
c
      if (dbgflg.ge.1) then
        write (dbgunit,1701)
1701    format (/2X,'Substitution of even-order high-level terms of ',
     *          'DKH Hamiltonian successful (SR "substitute_term").')
        if (dbgflg.ge.2) then
          write (dbgunit,1721) dkhorder
1721      format (/2X,'Final high-level DKH Hamiltonian up to dkhorder',
     *            ' = ',I2,':',/2X,53('-'))
          call output3 (dbgunit,opcounter,operleng,oporder,evenodd,
     *                  doperators,operators)
        endif
      endif
c
c
      if (.not.dkhscfflg) then
c
        do 200 k=2,xopcounter
          hits=0
          pos=1
 222      continue
          if (pos.lt.xoperleng(k)) then
#if defined(_MOLCAS_) || defined(MOLPRO)
          call get_dkoperators(k,xopstring,xoperators)
          idum1=index(xopstring(pos:xoperleng(k)),'CE0')
#else
          idum1=index(xoperators(k)(pos:xoperleng(k)),'CE0')
#endif
            if (idum1.ne.0) then
              hits=hits+1
              pos=pos+idum1+1
              goto 222
            endif
          endif
          if (hits.gt.1)  then
            write (stdout,1725) k
1725        format (/2X,'ERROR: In term no. ',I6,' of xoperators ',
     *              'occurs expression "CE0" more than once.',//2X,
     *              'This is not possible at this stage.',//2X,
     *              'STOP.',/)
            CALL Abend
          endif
          idum1=2
#if defined(_MOLCAS_) || defined(MOLPRO)
          call substitute_term (idum1,k,hits,xoperleng(k),xopstring,
     *                          sused,stimes,wstimes,s,scrleng,scrchar)
          call Put_dkoperators(k,xopstring,xoperators)
#else
          call substitute_term (idum1,k,hits,xoperleng(k),xoperators(k),
     *                          sused,stimes,wstimes,s,scrleng,scrchar)
#endif
c
 200    continue
c
        if (dbgflg.ge.1) then
          write (dbgunit,1801)
1801      format (/2X,'Substitution of even-order high-level terms of ',
     *            'property X successful (SR "substitute_term").')
          if (dbgflg.ge.2) then
            write (dbgunit,1821) xorder
1821        format (/2X,'Final high-level property X up to xorder',
     *              ' = ',I2,':',/2X,48('-'))
            call output3b (dbgunit,xopcounter,xoperleng,xoporder,
     *                     xevenodd,xdoperators,xoperators)
          endif
        endif
c
      endif
c
      do 120 j=1,sused
        if (stimes(j).gt.0) then
          snumber=snumber+1
          stimestot=stimestot+stimes(j)
        endif
 120  continue
c
      if (dbgflg.ge.1) then
        if (dkhscfflg) then
          write (dbgunit,1871)
1871      format (//2X,'Statistic about high-level substitution ',
     *            'process (Sxxx) in step 2:   (only for "operators"!)',
     *          /2X,65('-'),//5X,'Note:  * All values refer ',
     *          'to high-level situation of step2 and',/14X,'will ',
     *          'change (increase) dramatically in further steps.',/12X,
     *          '* All values depend strongly on dkhorder and xorder.',
     *          /12X,'* All other patterns did not occur.',/12X)
          write (dbgunit,1872)
1872      format (//5X,'#',6X,'step1',3X,'-->',3X,'step2',3X,'order',3X,
     *            'no. of subs.',3X,/34X,'(tot)',5X,'(stimes)',/)
        else
          write (dbgunit,1873)
1873      format (//2X,'Statistic about high-level substitution ',
     *            'process (Sxxx) in step 2:   (only for "operators"',
     *            ' and "xoperators"!)',
     *          /2X,65('-'),//5X,'Note:  * All values refer ',
     *          'to high-level situation of step2 and',/14X,'will ',
     *          'change (increase) dramatically in further steps.',/12X,
     *          '* All values depend strongly on dkhorder and xorder.',
     *          /12X,'* All other patterns did not occur.',/12X)
          write (dbgunit,1874)
1874      format (//5X,'#',6X,'step1',3X,'-->',3X,'step2',3X,'order',
     *            1X,'order',1X,'order',2X,'no. of subs.',3X,/35X,'(V)',
     *            3X,'(X)',2X,'(tot)',3X,'(stimes)',/)
        endif
        call output4 (dkhscfflg,dbgunit,sused,stimes,sorder,s,scrchar,
     *                stimestot)
c
        write (dbgunit,1769)
1769    format (/2X,'Subroutine subs_dkhstring (step2) has now been ',
     *          'completed.',//2X,145('*'))
      endif
CMR      if (dkhscfflg) then
CMR        write (stdout,1902) snumber,stimestot
CMR1902    format (/2X,'No. of diff. Sxxx-patterns suited for ',
CMR     *          'operator subst.:',11X,I3,/2X,'No. of high-level ',
CMR     *          'subst. (Sxxx) for operators:',14X,I8,
CMR     *           //2X,'Note: Also substitutions (Sxxx) in wops ',
CMR     *           'have been executed,',/8X,'but they have not been ',
CMR     *           'considered (stimes) for the statistic.')
CMR      else
CMR        write (stdout,1903) snumber,stimestot
CMR1903    format (/2X,'No. of diff. Sxxx-patterns suited for ',
CMR     *          '(x)operator subst.:',8X,I3,/2X,'No. of high-level ',
CMR     *          'subst. (Sxxx) for (x)operators:',11X,I8,
CMR     *           //2X,'Note: Also substitutions (Sxxx) in wops ',
CMR     *           'have been executed,',/8X,'but they have not been ',
CMR     *           'considered (stimes) for the statistic.')
CMR      endif
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer_array(ordercounter)
        call Unused_integer_array(xordercounter)
        call Unused_integer_array(wordercounter)
      end if
      end
