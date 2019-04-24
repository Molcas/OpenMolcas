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
* Copyright (C) 2004,2006, Markus Reiher                               *
*               2004, Alexander Wolf                                   *
************************************************************************
      subroutine get_dkhstring (dkhorder,xorder,paramtype,dkhscfflg,
     *             ordercounter,opcounter,operleng,oporder,evenodd,
     *             doperators,operators,xordercounter,xopcounter,
     *             xoperleng,xoporder,xevenodd,xdoperators,xoperators,
     *             wordercounter,wopscounter,wopsleng,woporder,eowops,
     *             dwops,wops,dkhzerothrsh)
c
c******************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  M. Reiher and A. Wolf  (Univ. Jena)
c
c   version:  2.0.3
c
c   last modified: 19.10.2006 (MR, ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer dkhorder,xorder
      character*(3) paramtype
      logical dkhscfflg
c
      integer ordercounter(0:maxorder),opcounter,
     *        operleng(maxoperators),oporder(maxoperators,3),
     *        Evenodd(Maxoperators),
     *        xordercounter(0:maxorder),xopcounter,
     *        xoperleng(maxoperators),xoporder(maxoperators),
     *        xevenodd(maxoperators)
      REAL*8 doperators(maxoperators),
     *                 xdoperators(maxoperators)
#if defined(_MOLCAS_) || defined(MOLPRO)
#include "WrkSpc.fh"
      character operators(*),xoperators(*)
      REAL*8 dwops(maxuops)
      character*(maxlength) wops(maxuops)
      integer wordercounter(maxorder),wopscounter,wopsleng(maxuops),
     *        woporder(maxuops,3),eowops(maxuops),
     *        oddcounter
      Integer dodd, odd, oddleng, oddorder, eoodd
      Integer uops,uopsleng,uoporder,eouops,duops
      integer ucounter, puop
#else
      character*(maxlength) operators(maxoperators),
     *                        xoperators(maxoperators)
      REAL*8 dwops(maxuops),dodd(maxuops)
      character*(maxlength) wops(maxuops),odd(maxuops)
      integer wordercounter(maxorder),wopscounter,wopsleng(maxuops),
     *        woporder(maxuops,3),eowops(maxuops),
     *        oddcounter,oddleng(maxuops),oddorder(maxuops,3),
     *        eoodd(maxuops)
      integer ucounter,uopsleng(maxuops),uoporder(maxuops),
     *        eouops(maxuops),puop(maxuops)
      REAL*8 duops(maxuops)
      character*(maxlength) uops(maxuops)
#endif
c
      REAL*8 ducoeffs(maxorder)
c
      REAL*8 dkhzerothrsh(0:maxorder)
      integer i,uniorder
      intrinsic DBLE,INT
c
CMR      write (stdout,1017)
CMR1017  format (5X,'|',70('-'),'|',/5X,'|',6X,'STEP 1: Symbolic ',
CMR     *        'high-level evaluation of DKH Hamiltonians',
CMR     *        6X,'|',/5X,'|',70('-'),'|')
c
      Call GetMem('odd','Allo','Char',odd,maxuops*maxlength)
      Call GetMem('oddleng','Allo','Inte',oddleng,maxuops)
      Call GetMem('oddorder','Allo','Inte',oddorder,maxuops*3)
      Call GetMem('eoodd','Allo','Inte',eoodd,maxuops)
      Call GetMem('dodd','Allo','Real',dodd,maxuops)
      Call GetMem('uops','Allo','Char',uops,maxuops*maxlength)
      Call GetMem('uopsleng','Allo','Inte',uopsleng,maxuops)
      Call GetMem('uoporder','Allo','Inte',uoporder,maxuops)
      Call GetMem('eouops','Allo','Inte',eouops,maxuops)
      Call GetMem('duops','Allo','Real',duops,maxuops)
      Call GetMem('puop','Allo','Inte',puop,maxuops)

      call initialize1 (dkhorder,xorder,paramtype,dkhscfflg,
     *         ordercounter,opcounter,operleng,oporder,evenodd,
     *         doperators,operators,xordercounter,xopcounter,xoperleng,
     *         xoporder,xevenodd,xdoperators,xoperators,wordercounter,
     *         wopscounter,wopsleng,woporder,eowops,dwops,wops,
     *         oddcounter,iWork(oddleng),iWork(oddorder),iWork(eoodd),
     &         Work(dodd),cWork(odd),ducoeffs,
     *         dkhzerothrsh)
      if (dkhscfflg) then
        uniorder = INT(DBLE(dkhorder)/2.0d0)
      else
        uniorder = max(INT(DBLE(dkhorder)/2.0d0),xorder)
      endif
c
      do 10 i=1,uniorder
c
CMR        if (i.le.9) then
CMR          write (stdout,1021) i
CMR1021      format (/2X,'Start unitary transformation U_',I1,':'/2X,
CMR     *            33('-'))
CMR        endif
CMR        if (i.ge.10) then
CMR          write (stdout,1022) i
CMR1022      format (/2X,'Start unitary transformation U_',I2,':'/2X,
CMR     *            34('-'))
CMR        endif
        call build_U (i,uniorder,dkhorder,xorder,dkhscfflg,ucounter,
     *             iwork(uopsleng),iWork(uoporder),iWork(eouops),
     &             Work(duops),iWork(puop),cWork(uops),opcounter,
     *             operleng,oporder,evenodd,doperators,operators,
     *             wordercounter,wopscounter,wopsleng,woporder,eowops,
     *             dwops,wops,oddcounter,iWork(oddleng),iWork(oddorder),
     &             iWork(eoodd),Work(dodd),
     *             cWork(odd),ducoeffs)
        call U_operators (i,dkhorder,xorder,dkhscfflg,
     *               ordercounter,opcounter,operleng,oporder,evenodd,
     *               doperators,operators,xordercounter,xopcounter,
     *               xoperleng,xoporder,xevenodd,xdoperators,xoperators,
     *               ucounter,iwork(uopsleng),iWork(uoporder),
     &               iWork(eouops),Work(duops),cWork(uops))
        call operators_U (i,uniorder,dkhorder,xorder,dkhscfflg,
     *               ordercounter,opcounter,operleng,oporder,evenodd,
     *               doperators,operators,xordercounter,xopcounter,
     *               xoperleng,xoporder,xevenodd,xdoperators,xoperators,
     *               ucounter,iwork(uopsleng),iWork(uoporder),
     &               iWork(eouops),iWork(puop),Work(duops),cWork(uops))
        call simplify1 (i,uniorder,dkhorder,xorder,dkhscfflg,
     *               ordercounter,opcounter,operleng,oporder,evenodd,
     *               doperators,operators,xordercounter,xopcounter,
     *               xoperleng,xoporder,xevenodd,xdoperators,xoperators)
        call simplify2 (i,uniorder,dkhorder,ordercounter,opcounter,
     *               operleng,oporder,evenodd,doperators,operators,
     *               wordercounter,wopsleng,woporder,eowops,dwops,wops,
     *               oddcounter,iWork(oddleng),iWork(oddorder),
     &               iWork(eoodd),Work(dodd),cWork(odd))
        call simplify3 (i,uniorder,dkhorder,xorder,dkhscfflg,
     *               ordercounter,opcounter,operleng,oporder,evenodd,
     *               doperators,operators,xordercounter,xopcounter,
     *               xoperleng,xoporder,xevenodd,xdoperators,xoperators,
     *               wordercounter,wopsleng,woporder,eowops,dwops,wops,
     *               dkhzerothrsh)
c
  10  continue
      if (dbgflg.ge.1) then
c
        if (uniorder.le.8) then
          write (dbgunit,1029) uniorder+1,uniorder,uniorder
1029      format (/2X,'Terms of Hamiltonian H',I1,' = [U',I1,' * ',
     *          '(operators) * U',I1,'^t]:',
     *          /2X,52('-'))
        endif
        if (uniorder.eq.9) then
          write (dbgunit,1030) uniorder+1,uniorder,uniorder
1030      format (/2X,'Terms of Hamiltonian H',I2,' = [U',I1,' * ',
     *          '(operators) * U',I1,'^t]:',
     *          /2X,53('-'))
        endif
        if (uniorder.ge.10) then
          write (dbgunit,1031) uniorder+1,uniorder,uniorder
1031      format (/2X,'Terms of Hamiltonian H',I2,' = [U',I2,' * ',
     *          '(operators) * U',I2,'^t]:',
     *          /2X,55('-'))
        endif
        call output3 (dbgunit,opcounter,operleng,
     *                oporder,evenodd,doperators,operators)
c
        if (.not.dkhscfflg) then
          if (uniorder.le.8) then
            write (dbgunit,1033) uniorder+1,uniorder,uniorder
1033        format (/2X,'Terms of property X',I1,' = [U',I1,' * ',
     *              '(xoperators) * U',I1,'^t]:',/2X,50('-'))
          endif
          if (uniorder.eq.9) then
            write (dbgunit,1034) uniorder+1,uniorder,uniorder
1034        format (/2X,'Terms of property X',I2,' = [U',I1,' * ',
     *              '(xoperators) * U',I1,'^t]:',/2X,51('-'))
          endif
          if (uniorder.ge.10) then
            write (dbgunit,1035) uniorder+1,uniorder,uniorder
1035        format (/2X,'Terms of property X',I2,' = [U',I2,' * ',
     *             '(operators) * U',I2,'^t]:',/2X,53('-'))
          endif
          call output3b (dbgunit,xopcounter,xoperleng,xoporder,
     *                   xevenodd,xdoperators,xoperators)
        endif
c
        write (dbgunit,1043) paramtype
1043    format (//2X,'Number of high-level terms of DKH Hamiltonian ',
     *        'and property after step 1 (get_dkhstring):',2X,A3)
        call distribution (dbgunit,dkhorder,xorder,dkhscfflg,
     *              ordercounter,opcounter,operleng,oporder,evenodd,
     *              doperators,operators,xordercounter,xopcounter,
     *              xoperleng,xoporder,xevenodd,xdoperators,xoperators)
        write (dbgunit,1051)
1051    format (/2X,'Subroutine get_dkhstring (step1) has now been ',
     *          'completed.',//2X,145('*'))
c
      endif
c
      Call GetMem('puop','Free','Inte',puop,maxuops)
      Call GetMem('duops','Free','Real',duops,maxuops)
      Call GetMem('eouops','Free','Inte',eouops,maxuops)
      Call GetMem('uoporder','Free','Inte',uoporder,maxuops)
      Call GetMem('uopsleng','Free','Inte',uopsleng,maxuops)
      Call GetMem('uops','Free','Char',uops,maxuops*maxlength)
      call getmem('dodd','Free','Real',dodd,maxuops)
      call getmem('eoodd','Free','Inte',eoodd,maxuops)
      call getmem('oddorder','Free','Inte',oddorder,maxuops*3)
      call getmem('oddleng','Free','Inte',oddleng,maxuops)
      call getmem('odd','Free','Char',odd,maxuops)
c
      return
      end
