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
* Copyright (C) 2004,2006,2007, Markus Reiher                          *
*               2004,2006, Alexander Wolf                              *
************************************************************************
      subroutine dkhparser_symbolic (dkhorder,xorder,paramtype,
     *                               dkhscfflg,no_prop,no_s,no_u)
c
c*****************************************************************************************
c
c                  Implementation of the scalar-relativistic
c
c                       I N F I N I T E - O R D E R
c
c      G E N E R A L I Z E D   D O U G L A S - K R O L L - H E S S  (DKH)
c
c                         T R A N S F O R M A T I O N
c
c                                    F O R
c
c             H A M I L T O N I A N   A N D   P R O P E R T Y
c
c                 2006 Markus Reiher and Alexander Wolf, ETH Zurich
c                      Infinite-order DKH
c                      {markus.reiher,alexander.wolf}@phys.chem.ethz.ch
c
c
c   Reference to the infinite-order DKH method:
c
c     M. Reiher, A. Wolf, J.Chem.Phys. 121 (2004) 2037--2047   part I   [Theory]
c     M. Reiher, A. Wolf, J.Chem.Phys. 121 (2004) 10945--10956 part II  [Algorithm]
c     A. Wolf, M. Reiher, J.Chem.Phys. 124 (2006) 064102       part III (Properties-Theory
c     A. Wolf, M. Reiher, J.Chem.Phys. 124 (2006) 064103       part IV  (Properties-Implem
c
c   Reference to the generalized higher-order DKH method:
c
c     A. Wolf, M. Reiher, B.A. Hess, J. Chem. Phys. 117 (2002) 9215
c
c-----------------------------------------------------------------------------------------
c
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.1.0
c
c   last modified: 15.01.2007 (M. Reiher, ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer dkhorder,xorder
      character*(3) paramtype
      logical dkhscfflg,no_prop,no_s,no_u
c
      integer ordercounter(0:maxorder),opcounter,
     *         xordercounter(0:maxorder),xopcounter
#if defined(_MOLCAS_) || defined(MOLPRO)
#include "WrkSpc.fh"
      integer operators,lwop,nwop,intrea,doperators,
     *        operleng,oporder,evenodd,xoperleng,xoporder,
     *        xevenodd,xdoperators,xoperators,dwops,wops,
     *        wopsleng,woporder,eowops,termleng,termleng2,
     *        termorder,dtcoeff,dtcoeff2
#else
Cfrueher      common/cdkhops1/ operators,xoperators,wops
      integer operleng(maxoperators),oporder(maxoperators,3),
     *        evenodd(maxoperators),
     *        xoperleng(maxoperators),xoporder(maxoperators),
     *        xevenodd(maxoperators),
     *        termleng(maxoperators),termleng2(maxoperators),
     *        termorder(maxoperators2,3)
      REAL*8 dtcoeff(maxoperators),dtcoeff2(maxoperators)
      REAL*8 doperators(maxoperators),
     *                 xdoperators(maxoperators)
      character*(maxlength) operators(maxoperators),
     *                      xoperators(maxoperators)
cweg  W-operators
      REAL*8 dwops(maxuops)
      character*(maxlength) wops(maxuops)
      integer wopsleng(maxuops),woporder(maxuops,3),eowops(maxuops)
#endif
c
      integer wordercounter(maxorder),wopscounter
c
      integer sused,snumber,scounter(maxsnumber),stimes(maxsnumber),
     *        wstimes(maxuops,maxsnumber),sorder(maxsnumber,3),
     *        smult(maxsnumber),scrleng(maxsnumber),stimestot
      integer tnumber,tcounter(maxsnumber),ttimes(maxsnumber),
     *        tmult(maxsnumber),tscrleng(maxsnumber),ttimestot
      integer uused,unumber,utimes(maxunumber),uorder(maxunumber,3),
     *        umult(maxunumber),uscrleng(maxunumber),utimestot,umulttot
      character*(3) uscrchar(maxunumber)
      character*(4) s(maxsnumber),t(maxsnumber),u(maxunumber)
      character*(9) scrchar(maxsnumber)
      character*(11) tscrchar(maxsnumber)
c
      REAL*8 dkhzerothrsh(0:maxorder)
      integer totmult
c
c-------------------------------------------------------------------------------
c
CMR
CMR NOTE: One might hard wire all low-order expressions as an option
CMR       here - however, since the set up is fast it might not be
CMR       worth the effort.
CMR
c
c-------------------------------------------------------------------------------
c
      call output1 (dkhorder,xorder,paramtype,dkhscfflg)
c
#if defined(_MOLCAS_) || defined(MOLPRO)
      lwop=8/intrea()
      nwop=(maxlength-1)/lwop+1
      call getmem('operators','Allo','Char',operators,maxoperators*nwop)
      call getmem('doperators','Allo','Real',doperators,maxoperators)
      call getmem('operleng','Allo','Inte',operleng,maxoperators)
      call getmem('oporder','Allo','Inte',oporder,maxoperators*3)
      call getmem('evenodd','Allo','Inte',evenodd,maxoperators)
      call getmem('dwops','Allo','Real',dwops,maxuops)
      call getmem('wops','Allo','Char',wops,maxuops*maxlength) !fix this
      call getmem('wopsleng','Allo','Inte',wopsleng,maxuops)
      call getmem('woporder','Allo','Inte',woporder,maxuops*3)
      call getmem('eowops','Allo','Inte',eowops,maxuops)
      call getmem('xoperators','Allo','Char',xoperators,
     *            maxoperators*nwop)
      call getmem('xdoperators','Allo','Real',xdoperators,maxoperators)
      call getmem('xoperleng','Allo','Inte',xoperleng,maxoperators)
      call getmem('xoporder','Allo','Inte',xoporder,maxoperators)
      call getmem('xevenodd','Allo','Inte',xevenodd,maxoperators)
      call getmem('termleng','Allo','Inte',termleng,maxoperators)
      call getmem('termleng2','Allo','Inte',termleng2,maxoperators)
      call getmem('termorder','Allo','Inte',termorder,3*maxoperators)
      call getmem('dtcoeff','Allo','Real',dtcoeff,maxoperators)
      call getmem('dtcoeff2','Allo','Real',dtcoeff2,maxoperators)
*                                                                      *
************************************************************************
*                                                                      *
c
      call get_dkhstring (dkhorder,xorder,paramtype,dkhscfflg,
     *          ordercounter,opcounter,iwork(operleng),iwork(oporder),
     *          iwork(evenodd),work(doperators),cwork(operators),
     *          xordercounter,xopcounter,iwork(xoperleng),
     *          iwork(xoporder),iwork(xevenodd),work(xdoperators),
     *          cwork(xoperators),wordercounter,
     *          wopscounter,iWork(wopsleng),iWork(woporder),
     *          iWork(eowops),Work(dwops),cWork(wops),dkhzerothrsh)
c
      call subs_dkhstring (dkhorder,xorder,dkhscfflg,ordercounter,
     *          opcounter,iwork(operleng),iwork(oporder),iwork(evenodd),
     *          work(doperators),cwork(operators),
     *          xordercounter,xopcounter,iwork(xoperleng),
     *          iwork(xoporder),iwork(xevenodd),work(xdoperators),
     *          cwork(xoperators),wordercounter,wopscounter,
     *          iWork(wopsleng),iWork(woporder),iWork(eowops),
     *          Work(dwops),cWork(wops),sused,snumber,scounter,stimes,
     *          wstimes,sorder,smult,s,scrchar,scrleng,stimestot,
     *          tnumber,tcounter,ttimes,tmult,t,tscrchar,tscrleng,
     *          ttimestot,no_s)
c
      call red_dkhstring (dkhorder,xorder,paramtype,dkhscfflg,
     *             ordercounter,opcounter,iwork(operleng),
     *             iwork(oporder),iwork(evenodd),
     *             work(doperators),cwork(operators),xordercounter,
     *             xopcounter,iwork(xoperleng),iwork(xoporder),
     *             iwork(xevenodd),work(xdoperators),cwork(xoperators),
     *             wordercounter,wopscounter,iWork(wopsleng),
     *             iWork(woporder),iWork(eowops),Work(dwops),
     *             cWork(wops),
     *             iwork(termleng),iwork(termleng2),iwork(termorder),
     *             work(dtcoeff),work(dtcoeff2),
     *             sused,snumber,scounter,stimes,wstimes,
     *             sorder,smult,s,scrchar,scrleng,stimestot,tnumber,
     *             tcounter,ttimes,tmult,t,tscrchar,tscrleng,ttimestot,
     *             uused,unumber,utimes,uorder,umult,u,uscrchar,
     *             uscrleng,totmult,dkhzerothrsh,no_u)
c
      call evaluate_Sxxx (dkhorder,xorder,paramtype,dkhscfflg,
     *             wordercounter,iWork(wopsleng),iWork(woporder),
     *             iWork(eowops),Work(dwops),cWork(wops),
     *             iwork(termleng),iwork(termleng2),iwork(termorder),
     *             work(dtcoeff),work(dtcoeff2),
     *             sused,snumber,scounter,stimes,wstimes,sorder,smult,s,
     *             scrchar,scrleng,stimestot,tnumber,tcounter,ttimes,
     *             tmult,t,tscrchar,tscrleng,ttimestot,uused,unumber,
     *             utimes,uorder,umult,u,uscrchar,uscrleng,totmult,
     *             dkhzerothrsh)
c
      call evaluate_Txxx (dkhorder,xorder,paramtype,dkhscfflg,
     *             wordercounter,iWork(wopsleng),iWork(woporder),
     *             iWork(eowops),Work(dwops),cWork(wops),
     *             iwork(termleng),iwork(termleng2),iwork(termorder),
     *             work(dtcoeff),work(dtcoeff2),
     *             sused,snumber,scounter,stimes,wstimes,sorder,s,
     *             scrchar,scrleng,stimestot,tnumber,tcounter,ttimes,
     *             tmult,t,tscrchar,tscrleng,ttimestot,uused,unumber,
     *             utimes,uorder,umult,u,uscrchar,uscrleng,utimestot,
     *             umulttot,totmult,dkhzerothrsh)
c
      call evaluate_Uxxx (dkhorder,xorder,paramtype,dkhscfflg,uused,
     *                    unumber,utimes,uorder,umult,u,uscrchar,
     *                    uscrleng,utimestot,umulttot,totmult)
*                                                                      *
************************************************************************
*                                                                      *
c
      call getmem('dtcoeff2','Free','Real',dtcoeff2,maxoperators)
      call getmem('dtcoeff','Free','Real',dtcoeff,maxoperators)
      call getmem('termorder','Free','Inte',termorder,3*maxoperators)
      call getmem('termleng2','Free','Inte',termleng2,maxoperators)
      call getmem('termleng','Free','Inte',termleng,maxoperators)
      call getmem('xevenodd','Free','Inte',xevenodd,maxoperators)
      call getmem('xoporder','Free','Inte',xoporder,maxoperators)
      call getmem('xoperleng','Free','Inte',xoperleng,maxoperators)
      call getmem('xdoperators','Free','Real',xdoperators,maxoperators)
      call getmem('xoperators','Free','Char',xoperators,maxoperators*
     *            nwop)
      call getmem('eowops','Free','Inte',eowops,maxuops)
      call getmem('woporder','Free','Inte',woporder,maxuops*3)
      call getmem('wopsleng','Free','Inte',wopsleng,maxuops)
      call getmem('wops','Free','Char',wops,maxuops)
      call getmem('dwops','Free','Real',dwops,maxuops)
      call getmem('evenodd','Free','Inte',evenodd,maxoperators)
      call getmem('oporder','Free','Inte',oporder,maxoperators*3)
      call getmem('operleng','Free','Inte',operleng,maxoperators)
      call getmem('doperators','Free','Real',doperators,maxoperators)
      call getmem('operators','Free','Char',operators,maxoperators*nwop)
#else
      call get_dkhstring (dkhorder,xorder,paramtype,dkhscfflg,
     *          ordercounter,opcounter,operleng,oporder,evenodd,
     *          doperators,operators,xordercounter,xopcounter,xoperleng,
     *          xoporder,xevenodd,xdoperators,xoperators,wordercounter,
     *          wopscounter,wopsleng,woporder,eowops,dwops,wops,
     *          dkhzerothrsh)
      call subs_dkhstring (dkhorder,xorder,dkhscfflg,ordercounter,
     *          opcounter,operleng,oporder,evenodd,doperators,operators,
     *          xordercounter,xopcounter,xoperleng,xoporder,xevenodd,
     *          xdoperators,xoperators,wordercounter,wopscounter,
     *          wopsleng,woporder,eowops,dwops,wops,sused,snumber,
     *          scounter,stimes,wstimes,sorder,smult,s,scrchar,scrleng,
     *          stimestot,tnumber,tcounter,ttimes,tmult,t,tscrchar,
     *          tscrleng,ttimestot,no_s)
      call red_dkhstring (dkhorder,xorder,paramtype,dkhscfflg,
     *             ordercounter,opcounter,operleng,oporder,evenodd,
     *             doperators,operators,xordercounter,xopcounter,
     *             xoperleng,xoporder,xevenodd,xdoperators,xoperators,
     *             wordercounter,wopscounter,wopsleng,woporder,eowops,
     *             dwops,wops,
     *             termleng,termleng2,termorder,dtcoeff,dtcoeff2,
     *             sused,snumber,scounter,stimes,wstimes,
     *             sorder,smult,s,scrchar,scrleng,stimestot,tnumber,
     *             tcounter,ttimes,tmult,t,tscrchar,tscrleng,ttimestot,
     *             uused,unumber,utimes,uorder,umult,u,uscrchar,
     *             uscrleng,totmult,dkhzerothrsh,no_u)
      call evaluate_Sxxx (dkhorder,xorder,paramtype,dkhscfflg,
     *             wordercounter,wopsleng,woporder,eowops,dwops,wops,
     *             termleng,termleng2,termorder,dtcoeff,dtcoeff2,
     *             sused,snumber,scounter,stimes,wstimes,sorder,smult,s,
     *             scrchar,scrleng,stimestot,tnumber,tcounter,ttimes,
     *             tmult,t,tscrchar,tscrleng,ttimestot,uused,unumber,
     *             utimes,uorder,umult,u,uscrchar,uscrleng,totmult,
     *             dkhzerothrsh)
      call evaluate_Txxx (dkhorder,xorder,paramtype,dkhscfflg,
     *             wordercounter,wopsleng,woporder,eowops,dwops,wops,
     *             termleng,termleng2,termorder,dtcoeff,dtcoeff2,
     *             sused,snumber,scounter,stimes,wstimes,sorder,s,
     *             scrchar,scrleng,stimestot,tnumber,tcounter,ttimes,
     *             tmult,t,tscrchar,tscrleng,ttimestot,uused,unumber,
     *             utimes,uorder,umult,u,uscrchar,uscrleng,utimestot,
     *             umulttot,totmult,dkhzerothrsh)
      call evaluate_Uxxx (dkhorder,xorder,paramtype,dkhscfflg,uused,
     *                    unumber,utimes,uorder,umult,u,uscrchar,
     *                    uscrleng,utimestot,umulttot,totmult)
#endif
c
c
CMR      write (stdout,2001)
CMR2001  format (/2X,'Symbolic DKH PARSER routines "dkhparser1" ',
CMR     *        'successfully completed.',//2X)
c
#ifdef _MOLCAS_
      if (dbgflg.ge.1) close (dbgunit)
#endif
c
      If (DKH_Verbose) Then
         Write (stdout,1111)
 1111    Format (//2X,
     &        'Start evaluation of the DKH Hamiltonian by matrix ',
     &        'multiplies ...')
      End If
*
      return
c Avoid unused argument warnings
      if (.false.) call Unused_logical(no_prop)
      end
