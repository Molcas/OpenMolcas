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
*               2004,2006, Markus Reiher                               *
************************************************************************
      subroutine operators_U (i,uniorder,dkhorder,xorder,dkhscfflg,
     *               ordercounter,opcounter,operleng,oporder,evenodd,
     *               doperators,operators,xordercounter,xopcounter,
     *               xoperleng,xoporder,xevenodd,xdoperators,xoperators,
     *               ucounter,uopsleng,uoporder,eouops,puop,duops,uops)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.3
c
c   last modified: 19.10.2006 (MR. ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer i,uniorder,dkhorder,xorder,ordercounter(0:maxorder),
     *        opcounter,operleng(maxoperators),oporder(maxoperators,3),
     *        evenodd(maxoperators),
     *        xordercounter(0:maxorder),xopcounter,
     *        xoperleng(maxoperators),xoporder(maxoperators),
     *        xevenodd(maxoperators)
      REAL*8 doperators(maxoperators),
     *                 xdoperators(maxoperators)
c
      integer ucounter,uopsleng(maxuops),uoporder(maxuops),
     *        eouops(maxuops),puop(maxuops)
      REAL*8 duops(maxuops)
      character*(maxlength) uops(maxuops)
#if defined(_MOLCAS_) || defined(MOLPRO)
      character operators(*),xoperators(*)
      character*(maxlength) opstring,opstring1,xopstring,xopstring1
#else
      character*(maxlength) operators(maxoperators)
     *                      xoperators(maxoperators)
#endif
c
      logical dkhscfflg,lasttrans
      integer idum1,idum2,j,l
      intrinsic INT,NINT,DBLE
c
      idum1=opcounter
      lasttrans=(i.eq.uniorder)
      do 10 j=1,ucounter
        do 20 l=1,idum1
          if (oporder(l,3)+uoporder(j).gt.dkhorder) goto 10
          if (lasttrans) then
            if (evenodd(l)*eouops(j).eq.-1) goto 20
          endif
          opcounter=opcounter+1
          if (opcounter.gt.maxoperators) then
            write (stdout,1001)
1001        format (/2X,'WARNING1: Too many operators (opcounter >',
     *              ' maxoperators)!',/11X,'Increase maxoperators ',
     *              '(e.g., by one order of magnitude)!',/2X,'STOP.',/)
            CALL Abend
          endif
          operleng(opcounter)=operleng(l)+uopsleng(j)
          if (operleng(opcounter).gt.maxlength) then
            write (stdout,1003)
1003        format (/2X,'WARNING2:  operleng > maxlength in SR',
     *              ' "operators_U".',/2X,'STOP.',/2X)
            CALL Abend
          endif
#if defined(_MOLCAS_) || defined(MOLPRO)
          call get_dkoperators(l,opstring1,operators)
          opstring=opstring1(1:operleng(l))//uops(j)(1:uopsleng(j))
          call put_dkoperators(opcounter,opstring,operators)
#else
          operators(opcounter)(1:operleng(opcounter))
     *        =operators(l)(1:operleng(l))//uops(j)(1:uopsleng(j))
#endif
          oporder(opcounter,3)=oporder(l,3)+uoporder(j)
          if (.not.dkhscfflg) then
            oporder(opcounter,1)=oporder(opcounter,3)
            oporder(opcounter,2)=0
          endif
          ordercounter(oporder(opcounter,3))=
     *      ordercounter(oporder(opcounter,3))+1
          evenodd(opcounter)=evenodd(l)*eouops(j)
          doperators(opcounter)=doperators(l)*duops(j)
* determine sign (because this is U_dagger !)
          if (INT(DBLE(puop(j))/2.0d0).ne.NINT(DBLE(puop(j))/2.0d0))
     *          doperators(opcounter)=-doperators(opcounter)
  20    continue
  10  continue
c
      if (.not.dkhscfflg) then
        idum2=xopcounter
        do 11 j=1,ucounter
          do 21 l=1,idum2
            if (xoporder(l)+uoporder(j).gt.xorder) goto 11
            if (lasttrans) then
              if (xevenodd(l)*eouops(j).eq.-1) goto 21
            endif
            xopcounter=xopcounter+1
            if (xopcounter.gt.maxoperators) then
              write (stdout,1005)
1005          format (/2X,'WARNING1: Too many xoperators (xopcounter >',
     *                ' maxoperators)!',/11X,'Increase maxoperators ',
     *                '(e.g., by one order of magnitude)!',/2X,'STOP.',
     *                /)
              CALL Abend
            endif
            xoperleng(xopcounter)=xoperleng(l)+uopsleng(j)
            if (xoperleng(xopcounter).gt.maxlength) then
              write (stdout,1006)
1006          format (/2X,'WARNING2:  xoperleng > maxlength in SR',
     *                ' "operators_U".',/2X,'STOP.',/2X)
              CALL Abend
            endif
#if defined(_MOLCAS_) || defined(MOLPRO)
          call get_dkoperators(l,xopstring1,xoperators)
          xopstring=xopstring1(1:xoperleng(l))//uops(j)(1:uopsleng(j))
          call put_dkoperators(xopcounter,xopstring,xoperators)
#else
            xoperators(xopcounter)(1:xoperleng(xopcounter))
     *          =xoperators(l)(1:xoperleng(l))//uops(j)(1:uopsleng(j))
#endif
            xoporder(xopcounter)=xoporder(l)+uoporder(j)
            xordercounter(xoporder(xopcounter))=
     *        xordercounter(xoporder(xopcounter))+1
            xevenodd(xopcounter)=xevenodd(l)*eouops(j)
            xdoperators(xopcounter)=xdoperators(l)*duops(j)
* determine sign (because this is U_dagger !)
            if (INT(DBLE(puop(j))/2.0d0).ne.NINT(DBLE(puop(j))/2.0d0))
     *          xdoperators(xopcounter)=-xdoperators(xopcounter)
  21      continue
  11    continue
      endif
c
      if (dbgflg.ge.1 .and. i.le.9) write (dbgunit,1011) i,i
1011  format (/2X,'Calculation of  U',I1,' * (operators) * U',I1,
     *        '^t has been successful.')
      if (dbgflg.ge.1 .and. i.le.9 .and. .not.dkhscfflg)
     *                                          write (dbgunit,1012) i,i
1012  format (2X,'Calculation of  U',I1,' * (xoperators) * U',I1,
     *        '^t has been successful.')
      if (dbgflg.ge.1 .and. i.ge.10) write (dbgunit,1013) i,i
1013  format (/2X,'Calculation of  U',I2,' * (operators) * U',I2,
     *        '^t has been successful.')
      if (dbgflg.ge.1 .and. i.ge.10 .and. .not.dkhscfflg)
     *                                          write (dbgunit,1014) i,i
1014  format (2X,'Calculation of  U',I2,' * (operators) * U',I2,
     *        '^t has been successful.')
c
      idum1=0
      idum2=0
      do 40 j=0,maxorder
        idum1=idum1+ordercounter(j)
        idum2=idum2+xordercounter(j)
  40  continue
      if (idum1.ne.opcounter) then
        write (stdout,1021) opcounter
1021    format (/2X,'WARNING: Number of operators "opcounter" = ',I7,
     *          ' deviates from sum over each order "ordercounter(j)".',
     *          //2X,'STOP in subroutine "operators_U.f".',/)
        CALL Abend
      endif
      if (idum2.ne.xopcounter) then
        write (stdout,1022) xopcounter
1022    format (/2X,'WARNING: Number of xoperators "xopcounter" = ',I7,
     *         ' deviates from sum over each order "xordercounter(j)".',
     *          //2X,'STOP in subroutine "operators_U.f".',/)
        CALL Abend
      endif
      if (dkhscfflg) then
        call sort_op1 (dkhorder,opcounter,operleng,oporder,
     &                 evenodd,doperators,operators)
        if (dbgflg.ge.1) write (dbgunit,1066)
1066    format (/2X,'Sorting of terms of Hamiltonian (operators) ',
     *          'in ascending order in V(X) succesful.')
      else
        call sort_op1 (dkhorder,opcounter,operleng,oporder(1,1),
     &                 evenodd,doperators,operators)
        call sort_op2 (xorder,xopcounter,xoperleng,xoporder,
     &                 xevenodd,xdoperators,xoperators)
        if (dbgflg.ge.1) write (dbgunit,1067)
1067    format (/2X,'Sorting of terms of Hamiltonian (operators) ',
     *          'in ascending order in V successful.')
        if (dbgflg.ge.1) write (dbgunit,1068)
1068    format (2X,'Sorting of terms of property (xoperators) in',
     *          ' ascending order in V successful.')
      endif
c
      if (dbgflg.ge.2) then
        if (i.le.9) then
          write (dbgunit,1026) i+1,i,i
1026      format (/2X,'Terms of Hamiltonian H',I1,' = [U',I1,' * ',
     *            '(operators) * U',I1,'^t]:',/2X,52('-'))
        endif
        if (i.ge.10) then
          write (dbgunit,1027) i+1,i,i
1027      format (/2X,'Terms of Hamiltonian H',I2,' = [U',I2,' * ',
     *            '(operators) * U',I2,'^t]:',/2X,55('-'))
        endif
        call output3 (dbgunit,opcounter,operleng,oporder,evenodd,
     *                doperators,operators)
      endif
c
      if (.not.dkhscfflg .and. dbgflg.ge.2) then
        if (i.le.9) then
          write (dbgunit,1028) i+1,i,i
1028      format (/2X,'Terms of property X',I1,' = [U',I1,' * ',
     *            '(xoperators) * U',I1,'^t]:',/2X,50('-'))
        endif
        if (i.ge.10) then
          write (dbgunit,1029) i+1,i,i
1029      format (/2X,'Terms of property X',I2,' = [U',I2,' * ',
     *            '(xoperators) * U',I2,'^t]:',/2X,53('-'))
        endif
        call output3b (dbgunit,xopcounter,xoperleng,xoporder,xevenodd,
     *                 xdoperators,xoperators)
      endif
c
      return
      end
