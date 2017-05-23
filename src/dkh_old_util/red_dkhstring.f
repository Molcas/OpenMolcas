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
      subroutine red_dkhstring (dkhorder,xorder,paramtype,dkhscfflg,
     *             ordercounter,opcounter,operleng,oporder,evenodd,
     *             doperators,operators,xordercounter,xopcounter,
     *             xoperleng,xoporder,xevenodd,xdoperators,xoperators,
     *             wordercounter,wopscounter,wopsleng,woporder,eowops,
     *             dwops,wops,termleng,termleng2,termorder,
     *             dtcoeff,dtcoeff2,
     *             sused,snumber,scounter,stimes,wstimes,
     *             sorder,smult,s,scrchar,scrleng,stimestot,
     *             tnumber,tcounter,ttimes,tmult,t,tscrchar,tscrleng,
     *             ttimestot,uused,unumber,utimes,uorder,
     *             umult,u,uscrchar,uscrleng,totmult,dkhzerothrsh,no_u)
c
c******************************************************************************
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
c******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer dkhorder,xorder
      character*(3) paramtype
      logical dkhscfflg,no_u
c
      integer ordercounter(0:maxorder),opcounter,operleng(maxoperators),
     *        oporder(maxoperators,3),evenodd(maxoperators),
     *        xordercounter(0:maxorder),xopcounter,
     *        xoperleng(maxoperators),xoporder(maxoperators),
     *        xevenodd(maxoperators),wordercounter(maxorder),
     *        wopscounter,wopsleng(maxuops),woporder(maxuops,3),
     *        eowops(maxuops)
      real*8 doperators(maxoperators),
     *                 xdoperators(maxoperators),dwops(maxuops)
#if defined(_MOLCAS_) || defined(MOLPRO)
#include "WrkSpc.fh"
      character*(maxlength) opstring,xopstring
      integer operators(*),xoperators(*)
#else
      character*(maxlength) operators(maxoperators),
     *                      xoperators(maxoperators)
#endif
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
      integer totmult
      real*8 dkhzerothrsh(0:maxorder)
c
      integer termcounter,termleng(maxoperators),
     *        termleng2(maxoperators),
     *        termorder(maxoperators,3)
      real*8 dtcoeff(maxoperators),dtcoeff2(maxoperators),
     *                 dummycoeff
#if defined(_MOLCAS_) || defined(MOLPRO)
      character*(maxlength) termstr
      integer term,nwop,lwop,intrea
#else
      character*(maxlength) term(maxoperators)
#endif
c
      integer i,j,k,l,wpos,opos,epos,cpos1,cpos2,idum1,idum2,idum3,
     *        idum4,ilow,istart,dummyleng,ctrflg(0:maxorder),
     *        matmulttot1,matmulttot2,matmultterm,opmult(0:maxorder),
     *        xopmult(0:maxorder)
      character*(maxlength) dummychar
#ifdef _MOLCAS_
      Integer IsFreeUnit
#endif
c
      If (dkhscfflg) Then
        Write (stdout,1016)
1016    format (//5X,'|',70('-'),'|',/5X,'|',6X,'STEP 3: Simplify ',
     *          'final Hamiltonian',30X,'|',/5X,'|',
     *          70('-'),'|',/2X)
CMR      Else
CMR        Write (stdout,1017)
CMR1017    format (//5X,'|',70('-'),'|',/5X,'|',6X,'STEP 3: Simplify ',
CMR     *          'final Hamiltonian and property X',15X,'|',/5X,'|',
CMR     *          70('-'),'|',/2X)
      End If
c
      If (dbgflg.ge.1) Then
        Write (dbgunit,1019)
1019    format (/2X,'Step 3: Enter subroutine "red_dkhstring" ',//2X,
     *          'Simplify (if still necessary) terms of the final ',
     *           'Hamiltonian and property X.')
      End If
      Call u_create (uused,unumber,utimes,uorder,umult,u,uscrleng,
     *               uscrchar,no_u)
      do 20 i=0,maxorder
        ctrflg(i)=0
  20  continue
c
#if defined(_MOLCAS_) || defined(MOLPRO)
      lwop=8/intrea()
      nwop=(maxlength-1)/lwop+1
      Call GetMem('term_red','Allo','Inte',term,maxoperators*nwop)
#endif
2456  continue
      do 30 i=0,dkhorder
c
        istart=1
        do 40 l=0,i-1
          istart=istart+ordercounter(l)
  40    continue
c
        do 45 j=istart,istart+ordercounter(i)-1
#if defined(_MOLCAS_) || defined(MOLPRO)
          Call get_dkoperators(j,opstring,operators)
          wpos=index(opstring(1:operleng(j)),'W')
          opos=index(opstring(1:operleng(j)),'O01')
          epos=index(opstring(1:operleng(j)),'E01')
          cpos1=index(opstring(1:operleng(j)),'CO0')
          cpos2=index(opstring(1:operleng(j)),'CE0')
#else
          wpos=index(operators(j)(1:operleng(j)),'W')
          opos=index(operators(j)(1:operleng(j)),'O01')
          epos=index(operators(j)(1:operleng(j)),'E01')
          cpos1=index(operators(j)(1:operleng(j)),'CO0')
          cpos2=index(operators(j)(1:operleng(j)),'CE0')
#endif
c
          If (wpos.gt.0 .or. opos.gt.0 .or. epos.gt.0 .or.
     *        cpos1.gt.0 .or. cpos2.gt.0) Then
            do 50 k=1,maxlength
              dummychar(k:k)=' '
  50        continue
            dummyleng=operleng(j)
            dummycoeff=doperators(j)
c
#if defined(_MOLCAS_) || defined(MOLPRO)
            dummychar=opstring(1:operleng(j))
            Call replace2 (dummyleng,dummycoeff,dummychar,wordercounter,
     *                     wopsleng,dwops,wops,sused,stimes,wstimes,s,
     *                     scrchar,scrleng,ttimes,termcounter,
     *                     termleng,termleng2,
     *                     dtcoeff,dtcoeff2,iwork(term))
#else
            dummychar(1:dummyleng)=operators(j)(1:operleng(j))
            Call replace2 (dummyleng,dummycoeff,dummychar,wordercounter,
     *                     wopsleng,dwops,wops,sused,stimes,wstimes,s,
     *                     scrchar,scrleng,ttimes,termcounter,
     *                     termleng,termleng2,dtcoeff,dtcoeff2,term)
#endif
c
            If (opcounter+termcounter.gt.maxoperators) Then
              Write (stdout,1033) j,opcounter+termcounter
1033          format (/2X,'In SR red_dkhstring: For j = ',I9,':  ',
     *                'opcounter+termcounter = ',I9,' > maxoperators.',
     *                //2X,'Increase maxoperators.',//2X,'STOP.',/)
              Call Abend
            End If
            idum2=0
            do 70 k=1,termcounter
#if defined(_MOLCAS_) || defined(MOLPRO)
              Call get_dkoperators(k,termstr,iwork(term))
              Call removeB1 (termleng(k),dtcoeff(k),termstr)
              Call removeB2 (termleng(k),termstr)
              Call insert_ri (termleng(k),termstr)
              Call movebraces (termleng(k),termstr)
              Call finalize (termleng(k),termstr,stimes,ttimes,t)
              Call finalize2 (termleng(k),termstr,uused,utimes,u,
     *                        uscrleng,uscrchar)
              termorder(k,3)=i
              Call calc_orders (dkhscfflg,termleng(k),termorder(k,1),
     *                          termorder(k,2),termorder(k,3),termstr,
     *                          sorder,uorder)
              Call put_dkoperators(k,termstr,iwork(term))
#else
              Call removeB1 (termleng(k),dtcoeff(k),term(k))
              Call removeB2 (termleng(k),term(k))
              Call insert_ri (termleng(k),term(k))
              Call movebraces (termleng(k),term(k))
              Call finalize (termleng(k),term(k),stimes,ttimes,t)
              Call finalize2 (termleng(k),term(k),uused,utimes,u,
     *                        uscrleng,uscrchar)
              termorder(k,3)=i
              Call calc_orders (dkhscfflg,termleng(k),termorder(k,1),
     *                          termorder(k,2),termorder(k,3),term(k),
     *                          sorder,uorder)
#endif
              If (dkhscfflg.and. termorder(k,2).gt.xorder) idum2=idum2+1
  70        continue
c
            idum3=termcounter-idum2
            If (idum3.ge.1) Then
              do 90 k=opcounter,j+1,-1
                operleng(k+idum3-1)=operleng(k)
                oporder(k+idum3-1,1)=oporder(k,1)
                oporder(k+idum3-1,2)=oporder(k,2)
                oporder(k+idum3-1,3)=oporder(k,3)
                evenodd(k+idum3-1)=evenodd(k)
                doperators(k+idum3-1)=doperators(k)
#if defined(_MOLCAS_) || defined(MOLPRO)
                Call copy_dkoperators(k,operators,
     *                              k+idum3-1,operators)
#else
              operators(k+idum3-1)=operators(k)
#endif
  90          continue
              idum4=0
              do 100 k=1,termcounter
                If (dynthrsh.eq.0 .and.abs(dtcoeff(k)).ge.dkhzero .and.
     *              (.not.dkhscfflg .or. termorder(k,2).le.xorder)) Then
                  idum4=idum4+1
                  operleng(j+idum4-1)=termleng(k)
                  oporder(j+idum4-1,1)=termorder(k,1)
                  oporder(j+idum4-1,2)=termorder(k,2)
                  oporder(j+idum4-1,3)=termorder(k,3)
                  evenodd(j+idum4-1)=evenodd(j)
                  doperators(j+idum4-1)=dtcoeff(k)
#if defined(_MOLCAS_) || defined(MOLPRO)
                  Call copy_dkoperators(k,iwork(term),j+idum4-1,
     *                                  operators)
#else
              operators(j+idum4-1)=term(k)
#endif
                End If
                If (dynthrsh.eq.1 .and.
     *              abs(dtcoeff(k)).ge.dkhzerothrsh(i) .and.
     *              (.not.dkhscfflg .or. termorder(k,2).le.xorder)) Then
                  idum4=idum4+1
                  operleng(j+idum4-1)=termleng(k)
                  oporder(j+idum4-1,1)=termorder(k,1)
                  oporder(j+idum4-1,2)=termorder(k,2)
                  oporder(j+idum4-1,3)=termorder(k,3)
                  evenodd(j+idum4-1)=evenodd(j)
                  doperators(j+idum4-1)=dtcoeff(k)
#if defined(_MOLCAS_) || defined(MOLPRO)
                  Call copy_dkoperators(k,iwork(term),j+idum4-1,
     *                                  operators)
#else
              operators(j+idum4-1)=term(k)
#endif
                End If
 100          continue
              opcounter=opcounter+idum3-1
              ordercounter(i)=ordercounter(i)+idum3-1
CMR: Is this necessary? Can't one start with reduced loop limits?
              goto 2456
            End If
c
            If (idum3.eq.0) Then
              do 110 k=j,opcounter-1
                operleng(k)=operleng(k+1)
                oporder(k,1)=oporder(k+1,1)
                oporder(k,2)=oporder(k+1,2)
                oporder(k,3)=oporder(k+1,3)
                evenodd(k)=evenodd(k+1)
                doperators(k)=doperators(k+1)
#if defined(_MOLCAS_) || defined(MOLPRO)
CMR Does this work?? Was adopted from above
C               Call copy_dkoperators(k+1,iwork(term),k,operators)
                Call copy_dkoperators(k+1,operators,k,operators)
#else
                operators(k)=operators(k+1)
#endif
 110          continue
              opcounter=opcounter-1
              ordercounter(i)=ordercounter(i)-1
CMR Does this work?? Was adopted from above
              goto 2456
            End If
c
          End If
c
  45    continue
c
        If (ctrflg(i).eq.0) Then
CMR          Write (stdout,1123) i,ordercounter(i)
CMR1123      format (2X,'order ',I2,':  all ',I8,' final terms ',
CMR     *            '(DKH Hamiltonian) successfully simplified.')
        End If
        ctrflg(i)=1
  30  continue
c     Write (stdout,*)
c
#ifndef MOLPRO
#ifdef _MOLCAS_
      dkhunit1=11
      dkhunit1=IsFreeUnit(dkhunit1)
      Call Molcas_Open(dkhunit1,'dkhops.11')
#else
      open (dkhunit1, file='dkhops.11', status='UNKNOWN',
     *        form='FORMATTED')
#endif
#endif
      rewind(dkhunit1)
      Write (dkhunit1,1427) dkhorder,paramtype,xorder,dkhscfflg
1427  format (50('-'),/2X,'dkhorder = ',I2,10X,A3,/2X,'xorder   = ',I2,
     *        /2X,'dkhscfflg = ',L1,/'+++',47('-'))
      ilow=ordercounter(0)+ordercounter(1)+1
      do 300 j=2,dkhorder
        Write (dkhunit1,1001) j,ordercounter(j)
1001    format ('***',/5X,I2,9X,'orders',4X,I6)
        idum1=0
        do 400 k=ilow,ilow+ordercounter(j)-1
          idum1=idum1+1
#if defined(_MOLCAS_) || defined(MOLPRO)
          Call get_dkoperators(k,opstring,operators)
          Write (dkhunit1,1467) idum1,operleng(k),oporder(k,1),
     *                          oporder(k,2),oporder(k,3),
     *                          opstring(1:operleng(k)),doperators(k)
#else
          Write (dkhunit1,1467) idum1,operleng(k),oporder(k,1),
     *                         oporder(k,2),oporder(k,3),
     *                         operators(k)(1:operleng(k)),doperators(k)
#endif
1467      format (I7,1X,I3,4X,I2,1X,I2,1X,I2,1X,A90,4X,F17.14)
 400    continue
        ilow=ilow+ordercounter(j)
 300  continue
#ifdef _MOLCAS_
      close (dkhunit1)
#endif
c
      If (dbgflg.ge.1) Then
        Write (dbgunit,1221) dkhorder,opcounter
1221    format (/2X,'Final low-level DKH Hamiltonian up to dkhorder = ',
     *        I3,':   (opcounter = ',I6,')',/2X,53('-'))
        If (dbgflg.ge.2) Call output3 (dbgunit,opcounter,operleng,
     *                             oporder,evenodd,doperators,operators)
      End If
      If (.not.dkhscfflg) Then
        do 21 i=0,maxorder
          ctrflg(i)=0
  21    continue
c
3456    continue
        do 31 i=0,xorder
c
          istart=1
          do 41 l=0,i-1
            istart=istart+xordercounter(l)
  41      continue
c
          do 245 j=istart,istart+xordercounter(i)-1
#if defined(_MOLCAS_) || defined(MOLPRO)
            Call get_dkoperators(j,xopstring,xoperators)
            wpos=index(xopstring(1:xoperleng(j)),'W')
            opos=index(xopstring(1:xoperleng(j)),'O01')
            epos=index(xopstring(1:xoperleng(j)),'E01')
            cpos1=index(xopstring(1:xoperleng(j)),'CO0')
            cpos2=index(xopstring(1:xoperleng(j)),'CE0')
#else
            wpos=index(xoperators(j)(1:xoperleng(j)),'W')
            opos=index(xoperators(j)(1:xoperleng(j)),'O01')
            epos=index(xoperators(j)(1:xoperleng(j)),'E01')
            cpos1=index(xoperators(j)(1:xoperleng(j)),'CO0')
            cpos2=index(xoperators(j)(1:xoperleng(j)),'CE0')
#endif
c
            If (wpos.gt.0.or.opos.gt.0.or.epos.gt.0.or.cpos1.gt.0
     *            .or.cpos2.gt.0) Then
              do 57 k=1,maxlength
                dummychar(k:k)=' '
  57          continue
              dummyleng=xoperleng(j)
              dummycoeff=xdoperators(j)
c
#if defined(_MOLCAS_) || defined(MOLPRO)
            dummychar=xopstring(1:xoperleng(j))
            Call replace2 (dummyleng,dummycoeff,dummychar,wordercounter,
     *                     wopsleng,dwops,wops,sused,stimes,wstimes,s,
     *                     scrchar,scrleng,ttimes,termcounter,
     *                     termleng,termleng2,
     *                     dtcoeff,dtcoeff2,iwork(term))
#else
              dummychar(1:dummyleng)=xoperators(j)(1:xoperleng(j))
              Call replace2 (dummyleng,dummycoeff,dummychar,
     *                   wordercounter,wopsleng,dwops,wops,sused,stimes,
     *                   wstimes,s,scrchar,scrleng,ttimes,termcounter,
     *                   termleng,termleng2,dtcoeff,dtcoeff2,term)
#endif
c
              If (xopcounter+termcounter.gt.maxoperators) Then
                Write (stdout,2133) j,xopcounter+termcounter
2133            format (/2X,'In SR red_xopstring: For j = ',I9,':  ',
     *                  'xopcounter+termcounter = ',I9,
     *                  ' > maxoperators.',//2X,'Increase maxoperators',
     *                  //2X,'STOP.',/)
                Call Abend
              End If
              idum2=0
              do 71 k=1,termcounter
                If (dynthrsh.eq.0 .and. abs(dtcoeff(k)).lt.dkhzero)
     *                                                     idum2=idum2+1
                If (dynthrsh.eq.1 .and.
     *               abs(dtcoeff(k)).lt.dkhzerothrsh(i))  idum2=idum2+1
c
#if defined(_MOLCAS_) || defined(MOLPRO)
              Call get_dkoperators(k,termstr,iwork(term))
              Call removeB1 (termleng(k),dtcoeff(k),termstr)
              Call removeB2 (termleng(k),termstr)
              Call insert_ri (termleng(k),termstr)
              Call movebraces (termleng(k),termstr)
              Call finalize (termleng(k),termstr,stimes,ttimes,t)
              Call finalize2 (termleng(k),termstr,uused,utimes,u,
     *                        uscrleng,uscrchar)
              Call put_dkoperators(k,termstr,iwork(term))
#else
                Call removeB1 (termleng(k),dtcoeff(k),term(k))
                Call removeB2 (termleng(k),term(k))
                Call insert_ri (termleng(k),term(k))
                Call movebraces (termleng(k),term(k))
                Call finalize (termleng(k),term(k),stimes,ttimes,t)
                Call finalize2 (termleng(k),term(k),uused,utimes,u,
     *                          uscrleng,uscrchar)
#endif
c
  71          continue
c
              idum3=termcounter-idum2
              If (idum3.ge.1) Then
                do 91 k=xopcounter,j+1,-1
                  xoperleng(k+idum3-1)=xoperleng(k)
                  xoporder(k+idum3-1)=xoporder(k)
                  xevenodd(k+idum3-1)=xevenodd(k)
                  xdoperators(k+idum3-1)=xdoperators(k)
#if defined(_MOLCAS_) || defined(MOLPRO)
                  Call copy_dkoperators(k,xoperators,
     *                                  k+idum3-1,xoperators)
#else
                xoperators(k+idum3-1)=xoperators(k)
#endif
  91            continue
                idum4=0
                do 102 k=1,termcounter
                  If (dynthrsh.eq.0 .and.
     *                                 abs(dtcoeff(k)).ge.dkhzero) Then
                    idum4=idum4+1
                    xoperleng(j+idum4-1)=termleng(k)
                    xoporder(j+idum4-1)=xoporder(j)
                    xevenodd(j+idum4-1)=xevenodd(j)
                    xdoperators(j+idum4-1)=dtcoeff(k)
#if defined(_MOLCAS_) || defined(MOLPRO)
                    Call copy_dkoperators(k,iwork(term),j+idum4-1,
     *                                    xoperators)
#else
                  xoperators(j+idum4-1)=term(k)
#endif
                  End If
                  If (dynthrsh.eq.1 .and.
     *                         abs(dtcoeff(k)).ge.dkhzerothrsh(i)) Then
                    idum4=idum4+1
                    xoperleng(j+idum4-1)=termleng(k)
                    xoporder(j+idum4-1)=xoporder(j)
                    xevenodd(j+idum4-1)=xevenodd(j)
                    xdoperators(j+idum4-1)=dtcoeff(k)
#if defined(_MOLCAS_) || defined(MOLPRO)
                    Call copy_dkoperators(k,iwork(term),j+idum4-1,
     *                                    xoperators)
#else
                  xoperators(j+idum4-1)=term(k)
#endif
                  End If
 102            continue
                xopcounter=xopcounter+idum3-1
                xordercounter(i)=xordercounter(i)+idum3-1
                goto 3456
              End If
c
              If (idum3.eq.0) Then
                do 112 k=j,xopcounter-1
                  xoperleng(k)=xoperleng(k+1)
                  xoporder(k)=xoporder(k+1)
                  xevenodd(k)=xevenodd(k+1)
                  xdoperators(k)=xdoperators(k+1)
                  xoperators(k)=xoperators(k+1)
#if defined(_MOLCAS_) || defined(MOLPRO)
CMR Does this work?? Was adopted from above
                  Call copy_dkoperators(k+1,xoperators,k,xoperators)
#else
                xoperators(k)=xoperators(k+1)
#endif
 112            continue
                xopcounter=xopcounter-1
                xordercounter(i)=xordercounter(i)-1
                goto 3456
              End If
c
            End If
 245      continue
          If (ctrflg(i).eq.0) Then
CMR            Write (stdout,1129) i,xordercounter(i)
CMR1129        format (2X,'order ',I2,':  all ',I8,' final terms ',
CMR     *              '(property X) successfully simplified.')
          End If
          ctrflg(i)=1
  31    continue
c
#ifndef MOLPRO
#ifdef _MOLCAS_
      dkhunit2=12
      dkhunit2=IsFreeUnit(dkhunit2)
      Call Molcas_Open(dkhunit2,'dkhops.12')
#else
      open (dkhunit2, file='dkhops.12', status='UNKNOWN',
#endif
#endif
        rewind(dkhunit2)
        Write (dkhunit2,3427) dkhorder,paramtype,xorder,dkhscfflg
3427    format (50('-'),/2X,'dkhorder = ',I2,10X,A3,/2X,'xorder   = ',
     *          I2,/2X,'dkhscfflg = ',L1,/'+++',47('-'))
        ilow=2
c  Loop 321 runs over all orders
        do 321 j=1,xorder
          Write (dkhunit2,1011) j,xordercounter(j)
1011      format ('***',/5X,I2,7X,'order(V)',2X,I6)
c  Loop 423 processes all terms of one specific order j
          idum1=0
          do 423 k=ilow,ilow+xordercounter(j)-1
            idum1=idum1+1
#if defined(_MOLCAS_) || defined(MOLPRO)
          Call get_dkoperators(k,xopstring,xoperators)
          Write (dkhunit2,1567) idum1,xoperleng(k),xoporder(k),
     *                          xopstring(1:xoperleng(k)),xdoperators(k)
#else
          Write (dkhunit2,1567) idum1,xoperleng(k),xoporder(k),
     *              xoperators(k)(1:xoperleng(k)),xdoperators(k)
#endif
1567        format (I7,1X,I3,5X,I2,5X,A90,4X,F17.14)
 423      continue
          ilow=ilow+xordercounter(j)
 321    continue
#ifdef _MOLCAS_
        close (dkhunit2)
#endif
c
        If (dbgflg.ge.1) Then
          Write (dbgunit,2221) xorder,xopcounter
2221      format (/2X,'Final low-level property X up to xorder = ',I3,
     *                ':   (xopcounter = ',I6,')',/2X,46('-'))
          If (dbgflg.ge.2) Call output3b (dbgunit,xopcounter,xoperleng,
     *                                    xoporder,xevenodd,xdoperators,
     *                                    xoperators)
        End If
c
      End If
c
      snumber=0
      stimestot=0
      do 120 j=1,sused
        If (stimes(j).gt.0) Then
          snumber=snumber+1
          stimestot=stimestot+stimes(j)
        End If
 120  continue
      tnumber=0
      ttimestot=0
      do 121 j=1,sused
        If (ttimes(j).gt.0) Then
          tnumber=tnumber+1
          ttimestot=ttimestot+ttimes(j)
        End If
 121  continue
      do 167 j=0,maxorder
        opmult(j)=0
        xopmult(j)=0
 167  continue
c
      ilow=1
      matmulttot1=0
      matmulttot2=0
      do 180 j=0,dkhorder
        do 200 k=ilow,ilow+ordercounter(j)-1
          matmultterm=0
#if defined(_MOLCAS_) || defined(MOLPRO)
          Call get_dkoperators(k,opstring,operators)
#endif
          do 250 l=1,operleng(k)
#if defined(_MOLCAS_) || defined(MOLPRO)
            If (opstring(l:l).eq.'V'.or.opstring(l:l).eq.'N'.or.
     *          opstring(l:l).eq.'D'.or.opstring(l:l).eq.'Y'.or.
     *          opstring(l:l).eq.'F'.or.opstring(l:l).eq.'G'.or.
     *          opstring(l:l).eq.'Z'.or.opstring(l:l).eq.'Q'.or.
     *          opstring(l:l).eq.'S'.or.opstring(l:l).eq.'T'.or.
     *          opstring(l:l).eq.'U'.or.opstring(l:l).eq.'X'.or.
     *          opstring(l:l).eq.'I'.or.opstring(l:l).eq.'J'.or.
     *          opstring(l:l).eq.'K'.or.opstring(l:l).eq.'L'.or.
     *          opstring(l:l).eq.'M') Then
#else
            If (operators(k)(l:l).eq.'V'.or.operators(k)(l:l).eq.'N'.or.
     *          operators(k)(l:l).eq.'D'.or.operators(k)(l:l).eq.'Y'.or.
     *          operators(k)(l:l).eq.'F'.or.operators(k)(l:l).eq.'G'.or.
     *          operators(k)(l:l).eq.'Z'.or.operators(k)(l:l).eq.'Q'.or.
     *          operators(k)(l:l).eq.'S'.or.operators(k)(l:l).eq.'T'.or.
     *          operators(k)(l:l).eq.'U'.or.operators(k)(l:l).eq.'X'.or.
     *          operators(k)(l:l).eq.'I'.or.operators(k)(l:l).eq.'J'.or.
     *          operators(k)(l:l).eq.'K'.or.operators(k)(l:l).eq.'L'.or.
     *          operators(k)(l:l).eq.'M') Then
#endif
              matmultterm=matmultterm+1
            End If
 250      continue
          If (matmultterm.ge.2) Then
            matmulttot1=matmulttot1+matmultterm-1
            opmult(j)=opmult(j)+matmultterm-1
          End If
 200    continue
        ilow=ilow+ordercounter(j)
 180  continue
c
c
      If (.not.dkhscfflg) Then
        ilow=1
        do 182 j=0,xorder
#if defined(_MOLCAS_) || defined(MOLPRO)
          Call get_dkoperators(k,xopstring,xoperators)
#endif
          do 202 k=ilow,ilow+xordercounter(j)-1
            matmultterm=0
            do 252 l=1,xoperleng(k)
#if defined(_MOLCAS_) || defined(MOLPRO)
            if(xopstring(l:l).eq.'V'.or.xopstring(l:l).eq.'N'.or.
     *         xopstring(l:l).eq.'D'.or.xopstring(l:l).eq.'Y'.or.
     *         xopstring(l:l).eq.'F'.or.xopstring(l:l).eq.'G'.or.
     *         xopstring(l:l).eq.'Z'.or.xopstring(l:l).eq.'Q'.or.
     *         xopstring(l:l).eq.'S'.or.xopstring(l:l).eq.'T'.or.
     *         xopstring(l:l).eq.'U'.or.xopstring(l:l).eq.'X'.or.
     *         xopstring(l:l).eq.'I'.or.xopstring(l:l).eq.'J'.or.
     *         xopstring(l:l).eq.'K'.or.xopstring(l:l).eq.'L'.or.
     *         xopstring(l:l).eq.'M') Then
#else
           if(xoperators(k)(l:l).eq.'V'.or.xoperators(k)(l:l).eq.'N'.or.
     *        xoperators(k)(l:l).eq.'D'.or.xoperators(k)(l:l).eq.'Y'.or.
     *        xoperators(k)(l:l).eq.'F'.or.xoperators(k)(l:l).eq.'G'.or.
     *        xoperators(k)(l:l).eq.'Z'.or.xoperators(k)(l:l).eq.'Q'.or.
     *        xoperators(k)(l:l).eq.'S'.or.xoperators(k)(l:l).eq.'T'.or.
     *        xoperators(k)(l:l).eq.'U'.or.xoperators(k)(l:l).eq.'X'.or.
     *        xoperators(k)(l:l).eq.'I'.or.xoperators(k)(l:l).eq.'J'.or.
     *        xoperators(k)(l:l).eq.'K'.or.xoperators(k)(l:l).eq.'L'.or.
     *        xoperators(k)(l:l).eq.'M') Then
#endif
              matmultterm=matmultterm+1
           End If
 252        continue
            If (matmultterm.ge.2) Then
              matmulttot2=matmulttot2+matmultterm-1
              xopmult(j)=xopmult(j)+matmultterm-1
            End If
 202      continue
          ilow=ilow+xordercounter(j)
 182    continue
      End If
c
      totmult=matmulttot1+matmulttot2
c
      If (dbgflg.ge.1) Then
        If (dkhscfflg) Then
          Write (dbgunit,2274)
2274      format (//2X,'Statistic about high-level substitution ',
     *            'process (Txxx) in step 3:  (only for "operators"!)',
     *            /2X,65('-'))
          Write (dbgunit,2275)
2275      format (/5X,'#',7X,'step1',4X,'-->',3X,'step3',3X,'order',3X,
     *          'no. of subs.',3X,/36X,'(tot)',5X,'(ttimes)'/)
        Else
          Write (dbgunit,2277)
2277      format (//2X,'Statistic about high-level substitution ',
     *            'process (Txxx) in step 3:  (only for "operators"',
     *            ' and "xoperators"!)',/2X,65('-'))
          Write (dbgunit,2278)
2278      format (/5X,'#',7X,'step1',4X,'-->',3X,'step3',3X,'order',1X,
     *            'order',1X,'order',2X,'no. of subs.',3X,/37X,'(V)',3X,
     *            '(X)',2X,'(tot)',3X,'(ttimes)'/)
        End If
        Call output5 (dkhscfflg,dbgunit,sused,ttimes,sorder,t,tscrchar,
     *                ttimestot)
        If (dkhscfflg) Then
          Write (dbgunit,1268)
1268      format (/2X,'No. of matrix multiplications necessary for ',
     *            'evaluation of operators:'/2X,68('-'),//9X,
     *            'order',10X,'opmult',/9X,'(tot)',/)
          do 290 j=0,dkhorder
            Write (dbgunit,1269) j,opmult(j)
1269        format (10X,I2,5X,I10)
 290      continue
          Write (dbgunit,1270) matmulttot1
1270      format (17X,10('-'),/17X,I10)
        Else
          Write (dbgunit,1278)
1278      format (/2X,'No. of matrix multiplications necessary for ',
     *            'evaluation of operators and xoperators:'/2X,83('-'),
     *            //9X,'order',10X,'opmult',13X,'xopmult',/10X,'(V)')
          do 291 j=0,max(dkhorder,xorder)
            Write (dbgunit,1279) j,opmult(j),xopmult(j)
1279        format (10X,I2,5X,I10,10X,I10)
 291      continue
          Write (dbgunit,1280) matmulttot1,matmulttot2,totmult
1280      format (17X,10('-'),10X,10('-'),/17X,I10,5X,'+',4X,I10,
     *            5X,'=',I10)
        End If
        Write (dbgunit,1043) paramtype
1043    format (//2X,'Number of terms of DKH Hamiltonian ',
     *        'and property X after step 3 (red_dkhstring):',2X,A3)
        Call distribution (dbgunit,dkhorder,xorder,dkhscfflg,
     *              ordercounter,opcounter,operleng,oporder,evenodd,
     *              doperators,operators,xordercounter,xopcounter,
     *              xoperleng,xoporder,xevenodd,xdoperators,xoperators)
c
        Write (dbgunit,1327)
1327    format (/2X,'Final operators have been written to the file ',
     *          '"dkhops.11".')
        If (.not.dkhscfflg) Write (dbgunit,1369)
1369    format (2X,'Final xoperators have been written to the file ',
     *          '"dkhops.12".')
c
        Write (dbgunit,1769)
1769    format (/2X,'Subroutine red_dkhstring (step3) has now been ',
     *          'completed.',//2X,145('*'))
      End If
c
CMR      If (dkhscfflg) Then
CMR        Write (stdout,1401) tnumber,ttimestot
CMR1401    format (/2X,'No. of different Sxxx-expressions occurring in',
CMR     *          /4X,'operators suited for substitution by Txxx = ',
CMR     *          'PSxxxP:',7X,I8,/2X,'No. of high-level subs. (Txxx) ',
CMR     *          'used in step3:',15X,I8)
CMR        Write (stdout,1402) opcounter-3,matmulttot1
CMR1402    format (/2X,'No. of low-level matrix-multipl. needed ',
CMR     *        'for the',/4X,'evaluation of all ',I8,
CMR     *        ' terms of "dkhops.11":',10X,I8)
CMR
CMR      Else
CMR        Write (stdout,1403) tnumber,ttimestot
CMR1403    format (/2X,'No. of different Sxxx-expressions occurring in',
CMR     *          /4X,'(x)operators suited for substitution by Txxx = ',
CMR     *          'PSxxxP:',4X,I8,/2X,'No. of high-level subs. (Txxx) ',
CMR     *          'used in step3:',15X,I8)
CMR        Write (stdout,1404) opcounter-2,totmult
CMR1404    format (/2X,'No. of low-level matrix-multipl. needed ',
CMR     *        'for the',/4X,'evaluation of all ',I8,
CMR     *        ' terms of "dkhops.11(12)":',6X,I8)
CMR      End If
c
#if defined(_MOLCAS_) || defined(MOLPRO)
      Call GetMem('term_red','Free','Inte',term,maxoperators*nwop)
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer(wopscounter)
        Call Unused_integer_array(woporder)
        Call Unused_integer_array(eowops)
        Call Unused_integer_array(scounter)
        Call Unused_integer_array(smult)
        Call Unused_integer_array(tcounter)
        Call Unused_integer_array(tmult)
        Call Unused_integer_array(tscrleng)
      End If
      End
