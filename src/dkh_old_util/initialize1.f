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
      subroutine initialize1 (dkhorder,xorder,paramtype,dkhscfflg,
     *         ordercounter,opcounter,operleng,oporder,evenodd,
     *         doperators,operators,xordercounter,xopcounter,xoperleng,
     *         xoporder,xevenodd,xdoperators,xoperators,wordercounter,
     *         wopscounter,wopsleng,woporder,eowops,dwops,wops,
     *         oddcounter,oddleng,oddorder,eoodd,dodd,odd,ducoeffs,
     *         dkhzerothrsh)

c
c************************************************************************
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
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer dkhorder,xorder,ordercounter(0:maxorder),opcounter,
     *        operleng(maxoperators),oporder(maxoperators,3),
     *        evenodd(maxoperators),
     *        xordercounter(0:maxorder),xopcounter,
     *        xoperleng(maxoperators),xoporder(maxoperators),
     *        xevenodd(maxoperators)
      REAL*8 doperators(maxoperators),
     *                 xdoperators(maxoperators)
#if defined(_MOLCAS_) || defined(MOLPRO)
      character operators(*),xoperators(*)
      character*(maxlength) opstring,xopstring
#else
      character*(maxlength) operators(maxoperators),
     *                        xoperators(maxoperators)
#endif
c
      integer wordercounter(maxorder),wopscounter,wopsleng(maxuops),
     *        woporder(maxuops,3),eowops(maxuops),
     *        oddcounter,oddleng(maxuops),oddorder(maxuops,3),
     *        eoodd(maxuops)
      REAL*8 dwops(maxuops),dodd(maxuops)
      character*(maxlength) wops(maxuops),odd(maxuops)
c
      character*(3) paramtype
      logical dkhscfflg
c
      integer i,j
      REAL*8 ducoeffs(maxorder),dkhzerothrsh(0:maxorder)
      intrinsic INT,DBLE
#ifdef _DEBUG_
#ifdef _MOLCAS_
      Integer IsFreeUnit,k
#endif
#endif
c
      if (DBLE(maxorder)/2.0d0.ne.DBLE(INT(DBLE(maxorder)/2.0d0))) then
        write (stdout,1002) maxorder
1002    format(5X,'Parameter  maxorder = ',I2,' has to be even.',//5X,
     *         'Adjust it in parameter file "dkhparameters.h"!',
     *         //5X,'STOP.')
        call Abend
      endif
      if (dkhorder.gt.maxorder .or. xorder.gt.maxorder) then
        write (stdout,1004) maxorder
1004    format(5X,'Parameter  maxorder = ',I2,' is too small.',//5X,
     *         'Increase it in parameter file "dkhparameters.h"!',
     *         //5X,'STOP.')
        call Abend
      endif
      if (paramtype .eq. 'EXP') then
        call exp_param (ducoeffs)
      else
        call other_param (ducoeffs,dkhorder,paramtype)
      endif
      opcounter=0
      xopcounter=0
      do 10 i=0,maxorder
        ordercounter(i)=0
        xordercounter(i)=0
  10  continue
      do 20 i=1,maxoperators
        operleng(i)=0
        do 21 j=1,3
          oporder(i,j)=0
  21    continue
#if defined(_MOLCAS_) || defined(MOLPRO)
      opstring=' '
      xopstring=' '
#endif
        evenodd(i)=0
        doperators(i)=0.0d0
        xoperleng(i)=0
        xoporder(i)=0
        xevenodd(i)=0
        xdoperators(i)=0.0d0
#if defined(_MOLCAS_) || defined(MOLPRO)
        call put_dkoperators(i,opstring,operators)
        call put_dkoperators(i,xopstring,xoperators)
#else
        do 25 k=1,maxlength
          operators(i)(k:k)=' '
          xoperators(i)(k:k)=' '
  25    continue
#endif
  20  continue
c
      wopscounter=0
      oddcounter=0
      do 40 i=1,maxorder
        wordercounter(i)=0
  40  continue
      do 50 i=1,maxuops
        wopsleng(i)=0
        eowops(i)=0
        dwops(i)=0.0d0
        oddleng(i)=0
        eoodd(i)=0
        dodd(i)=0.0d0
        do 51 j=1,3
          woporder(i,j)=0
          oddorder(i,j)=0
  51    continue
C        do 53 k=1,maxlength
CC         wops(i)(k:k)= ' '
CC         odd(i)(k:k)=  ' '
C          Write(wops(i)(k:k),'(A)') ' '
C          Write(odd(i)(k:k),'(A)')  ' '
C  53    continue
         wops(i)=' '
         odd(i)=' '
  50  continue
c
      do 60 i=0,maxorder
        dkhzerothrsh(i)=dkhzero
  60  continue
c
      dkhzerothrsh(5) =1.0d-3
      dkhzerothrsh(6) =1.0d-3
      dkhzerothrsh(7) =1.0d-2
      dkhzerothrsh(8) =1.0d-2
      dkhzerothrsh(9) =1.0d-2
      dkhzerothrsh(10)=1.0d-2
      dkhzerothrsh(11)=1.0d-2
      dkhzerothrsh(12)=1.0d-2
      dkhzerothrsh(13)=1.0d-2
      dkhzerothrsh(14)=1.0d-2
      dkhzerothrsh(15)=1.0d-2
      dkhzerothrsh(16)=1.0d-2
      dkhzerothrsh(17)=1.0d-2
      dkhzerothrsh(18)=1.0d-2
      dkhzerothrsh(19)=1.0d-2
      dkhzerothrsh(20)=1.0d-2
      dkhzerothrsh(21)=1.0d-2
      dkhzerothrsh(22)=1.0d-2
c
      opcounter=1
      ordercounter(0)=ordercounter(0)+1
      operleng(opcounter)=3
      oporder(opcounter,1)=0
      oporder(opcounter,2)=0
      oporder(opcounter,3)=oporder(opcounter,1)+oporder(opcounter,2)
      evenodd(opcounter)=+1
      doperators(opcounter)=1.0d0
#if defined(_MOLCAS_) || defined(MOLPRO)
      opstring='E00'
      call put_dkoperators(opcounter,opstring,operators)
#else
      operators(opcounter)='E00'
#endif
      opcounter=opcounter+1
      ordercounter(1)=ordercounter(1)+1
      operleng(opcounter)=3
      oporder(opcounter,1)=1
      oporder(opcounter,2)=0
      oporder(opcounter,3)=oporder(opcounter,1)+oporder(opcounter,2)
      evenodd(opcounter)=+1
      doperators(opcounter)=1.0d0
#if defined(_MOLCAS_) || defined(MOLPRO)
      opstring='E01'
      call put_dkoperators(opcounter,opstring,operators)
#else
      operators(opcounter)='E01'
#endif
      opcounter=opcounter+1
      ordercounter(1)=ordercounter(1)+1
      operleng(opcounter)=3
      oporder(opcounter,1)=1
      oporder(opcounter,2)=0
      oporder(opcounter,3)=oporder(opcounter,1)+oporder(opcounter,2)
      evenodd(opcounter)=-1
      doperators(opcounter)=1.0d0
#if defined(_MOLCAS_) || defined(MOLPRO)
      opstring='O01'
      call put_dkoperators(opcounter,opstring,operators)
#else
      operators(opcounter)='O01'
#endif
c
      if (dkhscfflg) then
        opcounter=opcounter+1
        ordercounter(1)=ordercounter(1)+1
        operleng(opcounter)=3
        oporder(opcounter,1)=0
        oporder(opcounter,2)=1
        oporder(opcounter,3)=oporder(opcounter,1)+oporder(opcounter,2)
        evenodd(opcounter)=+1
        doperators(opcounter)=1.0d0
#if defined(_MOLCAS_) || defined(MOLPRO)
      opstring='CE0'
      call put_dkoperators(opcounter,opstring,operators)
#else
      operators(opcounter)='CE0'
#endif
        opcounter=opcounter+1
        ordercounter(1)=ordercounter(1)+1
        operleng(opcounter)=3
        oporder(opcounter,1)=0
        oporder(opcounter,2)=1
        oporder(opcounter,3)=oporder(opcounter,1)+oporder(opcounter,2)
        evenodd(opcounter)=-1
        doperators(opcounter)=1.0d0
#if defined(_MOLCAS_) || defined(MOLPRO)
      opstring='CO0'
      call put_dkoperators(opcounter,opstring,operators)
#else
      operators(opcounter)='CO0'
#endif
      endif
      if (.not.dkhscfflg) then
        xopcounter=1
        xordercounter(0)=xordercounter(0)+1
        xoperleng(xopcounter)=3
        xoporder(xopcounter)=0
        xevenodd(xopcounter)=+1
        xdoperators(xopcounter)=1.0d0
#if defined(_MOLCAS_) || defined(MOLPRO)
        xopstring='CE0'
        call put_dkoperators(xopcounter,xopstring,xoperators)
#else
        xoperators(xopcounter)='CE0'
#endif
        xopcounter=xopcounter+1
        xordercounter(0)=xordercounter(0)+1
        xoperleng(xopcounter)=3
        xoporder(xopcounter)=0
        xevenodd(xopcounter)=-1
        xdoperators(xopcounter)=1.0d0
#if defined(_MOLCAS_) || defined(MOLPRO)
        xopstring='CO0'
        call put_dkoperators(xopcounter,xopstring,xoperators)
#else
        xoperators(xopcounter)='CO0'
#endif
      endif
c
#ifdef _DEBUG_
      if (dbgflg.ge.1) then
#ifdef _MOLCAS_
        dbgunit=32
        dbgunit=IsFreeUnit(dbgunit)
        Call Molcas_Open(dbgunit,'dbgout.32')
#else
        open (dbgunit, file='dbgout.32', status='UNKNOWN',
     *        form='FORMATTED')
#endif
        rewind(dbgunit)
        write (dbgunit,1012) dkhorder,xorder,dkhscfflg,paramtype,
     *                       dynthrsh
1012    format (/2X,'dkhorder  = ',I2,/2X,'xorder',4X,'= ',I2,/2X,
     *          'dkhscfflg =  ',L1,/2X,'paramtype = ',A3,/2X,'dynthrsh',
     *          '  =  ',I1,//2X,'Form of the initial fpFW ',
     *          'transformed Hamiltonian H1:',/2X,52('-'))
        call output3 (dbgunit,opcounter,operleng,oporder,evenodd,
     *                doperators,operators)
        if (.not.dkhscfflg) then
          write (dbgunit,1013)
1013      format (//2X,'Form of the initial fpFW transformed property',
     *          ' operator X1:',/2X,58('-'))
          call output3b (dbgunit,xopcounter,xoperleng,xoporder,xevenodd,
     *                  xdoperators,xoperators)
        endif
      endif

      if (dbgflg.ge.1) then
        write (dbgunit,1019)
1019    format (//2X,'Coefficients for the parametrization',10X,
     *          'zero thresholds',/4X,'of the unitary transformations:',
     *       17X,'dkhzero(k):',/2X,36('-'),8X,23('-'))
        write (dbgunit,1021)
1021    format (5X,'a_0   =   1.0000000000')
        do 80 k=1,maxorder
          if (k.le.9) then
            write (dbgunit,1025) k, ducoeffs(k), k, dkhzerothrsh(k)
1025        format (5X,'a_',I1,'   =  ',F13.10,19X,'dkhzero(',I1,
     *              ')  = ',D9.2)
          endif
          if (k.gt.9) then
            write (dbgunit,1029) k, ducoeffs(k), k, dkhzerothrsh(k)
1029        format (5X,'a_',I2,'  =  ',D13.6,19X,'dkhzero(',I2,
     *              ') = ',D9.2)
          endif
  80    continue
        write (dbgunit,1033)
1033    format(2X/)
      endif
c
      if (dbgflg.ge.1) write (dbgunit,1037)
1037  format (/2X,145('*'),//2X,'Step 1: Enter subroutine',
     *        ' "get_dkhstring"')
c
#endif
      return
      end
