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
      subroutine U_operators (i,dkhorder,xorder,dkhscfflg,
     *                ordercounter,opcounter,operleng,oporder,evenodd,
     *                doperators,operators,xordercounter,xopcounter,
     *                xoperleng,xoporder,xevenodd,xdoperators,
     *                xoperators,ucounter,uopsleng,uoporder,eouops,
     *                duops,uops)
c
c************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.3
c
c   last modified: 19.10.2006 (MR, ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer i,dkhorder,xorder
      logical dkhscfflg
c
      integer ordercounter(0:maxorder),opcounter,
     *        operleng(maxoperators),oporder(maxoperators,3),
     *        evenodd(maxoperators),
     *        xordercounter(0:maxorder),xopcounter,
     *        xoperleng(maxoperators),xoporder(maxoperators),
     *        xevenodd(maxoperators),
     *        ucounter,uopsleng(maxuops),uoporder(maxuops),
     *        eouops(maxuops)
c
      real*8 doperators(maxoperators),
     *       xdoperators(maxoperators),duops(maxuops)
      character*(maxlength) uops(maxuops)
#if defined(_MOLCAS_) || defined(MOLPRO)
      character*(maxlength) opstring,opstring1,xopstring,xopstring1
      integer operators(*),xoperators(*)
#else
      character*(maxlength) operators(maxoperators),
     *         xoperators(maxoperators)
#endif
c
      integer idum1,idum2,j,l
c
      idum1=opcounter
      idum2=xopcounter
      do 10 j=1,ucounter
        do 20 l=1,idum1
          if (uoporder(j)+oporder(l,3).gt.dkhorder) goto 9
          opcounter=opcounter+1
          if (opcounter.gt.maxoperators) then
            write (stdout,1005)
1005        format (/2X,'WARNING1: Too many operators (opcounter >',
     *              ' maxoperators) in SR "U_operators"!',/11X,
     *              'Increase maxoperators (e.g., by one order of ',
     *              'magnitude)!')
            CALL Abend
          endif
          operleng(opcounter)=uopsleng(j)+operleng(l)
          if (operleng(opcounter).gt.maxlength) then
            write (stdout,1007)
1007        format (/2X,'WARNING2:  operleng > maxlength in SR',
     *              ' "U_operators".',/2X,'STOP.',/2X)
            CALL Abend
          endif
#if defined(_MOLCAS_) || defined(MOLPRO)
          call get_dkoperators(l,opstring1,operators)
          opstring=uops(j)(1:uopsleng(j))//opstring1(1:operleng(l))
          call put_dkoperators(opcounter,opstring,operators)
#else
          operators(opcounter)(1:operleng(opcounter))
     *      =uops(j)(1:uopsleng(j))//operators(l)(1:operleng(l))
#endif
          oporder(opcounter,3)=uoporder(j)+oporder(l,3)
          if (.not.dkhscfflg) then
            oporder(opcounter,1)=oporder(opcounter,3)
            oporder(opcounter,2)=0
          endif
          ordercounter(oporder(opcounter,3))=
     *      ordercounter(oporder(opcounter,3))+1
          evenodd(opcounter)=eouops(j)*evenodd(l)
          doperators(opcounter)=duops(j)*doperators(l)
  20    continue
c
   9    continue
c
        if (.not.dkhscfflg) then
* Peng          idum2=xopcounter
          do 21 l=1,idum2
            if (uoporder(j)+xoporder(l).gt.xorder) goto 10
            xopcounter=xopcounter+1
            xoperleng(xopcounter)=uopsleng(j)+xoperleng(l)
#if defined(_MOLCAS_) || defined(MOLPRO)
          call get_dkoperators(l,xopstring1,xoperators)
          xopstring=uops(j)(1:uopsleng(j))//xopstring1(1:xoperleng(l))
          call put_dkoperators(xopcounter,xopstring,xoperators)
#else
          xoperators(xopcounter)(1:xoperleng(xopcounter))
     *      =uops(j)(1:uopsleng(j))//xoperators(l)(1:xoperleng(l))
#endif
            xoporder(xopcounter)=uoporder(j)+xoporder(l)
            xordercounter(xoporder(xopcounter))=
     *        xordercounter(xoporder(xopcounter))+1
            xevenodd(xopcounter)=eouops(j)*xevenodd(l)
            xdoperators(xopcounter)=duops(j)*xdoperators(l)
  21      continue
        endif
  10  continue
c
c
CAW------  not needed for standard runs
c      if (i.le.9) then
c        write (stdout,1011) i,opcounter
c1011    format (2X,'Total number of operators after U',I1,' * ',
c     *          '(operators):',11X,I8)
c        if (.not.dkhscfflg) write (stdout,1012) i,xopcounter
c1012    format (2X,'Total number of xoperators',
c     *          ' after U',I1,' * (xoperators):',9X,I8)
c      endif
c      if (i.ge.10) then
c        write (stdout,1013) i,opcounter
c1013    format (2X,'Total number of operators after U',I2,' * ',
c     *          '(operators):',10X,I8)
c        if (.not.dkhscfflg) write (stdout,1014) i,xopcounter
c1014    format (2X,'Total number of xoperators',
c     *          ' after U',I1,' * (xoperators):',8X,I8)
c      endif
CAW------------------------------------
c
      if (dbgflg.ge.1 .and. i.le.9) write (dbgunit,1016) i
1016  format (/2X,'Calculation of  U',I1,' * (operators) has been',
     *        ' successful.')
      if (dbgflg.ge.1 .and. i.le.9 .and. .not.dkhscfflg)
     *                                            write (dbgunit,1017) i
1017  format (2X,'Calculation of  U',I1,' * (xoperators) has been ',
     *        'successful.')
      if (dbgflg.ge.1 .and. i.ge.10) write (dbgunit,1018) i
1018  format (/2X,'Calculation of  U',I2,' * (operators) has been',
     *        ' successful.')
      if (dbgflg.ge.1 .and. i.ge.10 .and. .not.dkhscfflg)
     *                                            write (dbgunit,1019) i
1019  format (2X,'Calculation of  U',I1,' * (xoperators) has been ',
     *        'successful.')
c
      idum1=0
      idum2=0
      do 30 j=0,maxorder
        idum1=idum1+ordercounter(j)
        idum2=idum2+xordercounter(j)
  30  continue
      if (idum1.ne.opcounter) then
        write (stdout,1022) opcounter
1022    format (/2X,'ERROR in SR "U_operators": Number of operators ',
     *          '"opcounter" = ',I7,' deviates from sum over each ',
     *          'order "ordercounter(j)".',//2X,'STOP.',/)
            CALL Abend
      endif
      if (idum2.ne.xopcounter) then
        write (stdout,1023) xopcounter
1023    format (/2X,'ERROR in SR "U_operators": Number of xoperators ',
     *          '"xopcounter" = ',I7,' deviates from sum over each ',
     *          'order "xordercounter(j)".',//2X,'STOP.',/)
            CALL Abend
      endif
      if (dkhscfflg) then
        call sort_op1 (dkhorder,opcounter,operleng,oporder,
     *                 evenodd,doperators,operators)
        if (dbgflg.ge.1) write (dbgunit,1066)
1066    format (/2X,'Sorting of terms of Hamiltonian (operators) ',
     *          'in ascending order in V(X) succesful.')
      else
        call sort_op1 (dkhorder,opcounter,operleng,oporder(1,1),
     *                 evenodd,doperators,operators)
        call sort_op2 (xorder,xopcounter,xoperleng,xoporder,
     *                 xevenodd,xdoperators,xoperators)
        if (dbgflg.ge.1) write (dbgunit,1067)
1067    format (/2X,'Sorting of terms of Hamiltonian (operators) ',
     *          'in ascending order in V succesful.')
        if (dbgflg.ge.1) write (dbgunit,1068)
1068    format (2X,'Sorting of terms of property (xoperators) in',
     *          ' ascending order in V succesful.')
      endif
c
      if (dbgflg.ge.3)  then
        if (i.le.9) then
          write (dbgunit,1029) i
1029      format (/2X,'Form of intermediate Hamiltonian [U',I1,' * ',
     *          '(operators)]:',/2X,52('-'))
        endif
        if (i.ge.10) then
          write (dbgunit,1030) i
1030      format (/2X,'Form of intermediate Hamiltonian [U',I2,' * ',
     *          '(operators)]:',/2X,53('-'))
        endif
        call output3 (dbgunit,opcounter,operleng,oporder,evenodd,
     *                doperators,operators)
      endif
c
      if (.not.dkhscfflg .and. dbgflg.ge.3)  then
        if (i.le.9) then
          write (dbgunit,1039) i
1039      format (/2X,'Form of intermediate property [U',I1,' * ',
     *          '(xoperators)]:',/2X,50('-'))
        endif
        if (i.ge.10) then
          write (dbgunit,1040) i
1040      format (/2X,'Form of intermediate property [U',I2,' * ',
     *          '(xoperators)]:',/2X,51('-'))
        endif
        call output3b (dbgunit,xopcounter,xoperleng,xoporder,xevenodd,
     *                 xdoperators,xoperators)
      endif
c
      return
      end
