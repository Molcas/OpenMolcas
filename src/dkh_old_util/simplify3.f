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
      subroutine simplify3 (i,uniorder,dkhorder,xorder,dkhscfflg,
     *               ordercounter,opcounter,operleng,oporder,evenodd,
     *               doperators,operators,xordercounter,xopcounter,
     *               xoperleng,xoporder,xevenodd,xdoperators,xoperators,
     *               wordercounter,wopsleng,woporder,eowops,dwops,wops,
     *               dkhzerothrsh)
c
c***************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 30.06.2005
c
c   first version: 23.07.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer i,uniorder,dkhorder,xorder,ordercounter(0:maxorder),
     *        opcounter,operleng(maxoperators),oporder(maxoperators,3),
     *        evenodd(maxoperators),xordercounter(0:maxorder),
     *        xopcounter,xoperleng(maxoperators),xoporder(maxoperators),
     *        xevenodd(maxoperators),
     *        wordercounter(maxorder),wopsleng(maxuops),
     *        woporder(maxuops,3),eowops(maxuops)
      real*8 doperators(maxoperators),
     *                 xdoperators(maxoperators),dwops(maxuops)
      character*(maxlength) wops(maxuops)
#if defined(_MOLCAS_) || defined(MOLPRO)
      character*(maxlength) opstring1,opstring2,xopstring1,xopstring2
      integer operators(*),xoperators(*)
#else
      character*(maxlength) operators(maxoperators),
     *          xoperators(maxoperators)
#endif
c
      real*8 dkhzerothrsh(0:maxorder)
c
      logical dkhscfflg
      integer j,k,l,m,ilow
      intrinsic DABS
 240  continue
      ilow=1
      do 10 j=0,dkhorder
        do 20 k=ilow,ilow+ordercounter(j)-1
#if defined(_MOLCAS_) || defined(MOLPRO)
          call get_dkoperators(k,opstring1,operators)
#endif
          do 30 l=k+1,ilow+ordercounter(j)-1
#if defined(_MOLCAS_) || defined(MOLPRO)
            call get_dkoperators(l,opstring2,operators)
            if ( (operleng(k).eq.operleng(l)) .and.
     *               (oporder(k,3).eq.oporder(l,3)) .and.
     *               opstring1(1:operleng(k)).eq.
     *               opstring2(1:operleng(l)) ) then
#else
            if ( (operleng(k).eq.operleng(l)) .and.
     *               (oporder(k,3).eq.oporder(l,3)) .and.
     *               operators(k)(1:operleng(k)).eq.
     *                   operators(l)(1:operleng(l)) ) then
#endif
              doperators(k)=doperators(k)+doperators(l)
              do 40 m=l,opcounter
                operleng(m)=operleng(m+1)
                oporder(m,1)=oporder(m+1,1)
                oporder(m,2)=oporder(m+1,2)
                oporder(m,3)=oporder(m+1,3)
                evenodd(m)=evenodd(m+1)
                doperators(m)=doperators(m+1)
#if defined(_MOLCAS_) || defined(MOLPRO)
                call copy_dkoperators(m+1,operators,m,operators)
#else
                operators(m)=operators(m+1)
#endif
  40          continue
              opcounter=opcounter-1
              ordercounter(j)=ordercounter(j)-1
              if (abs(doperators(k)).lt.dkhzero) then
                do 50 m=k,opcounter
                  operleng(m)=operleng(m+1)
                  oporder(m,1)=oporder(m+1,1)
                  oporder(m,2)=oporder(m+1,2)
                  oporder(m,3)=oporder(m+1,3)
                  evenodd(m)=evenodd(m+1)
                  doperators(m)=doperators(m+1)
#if defined(_MOLCAS_) || defined(MOLPRO)
                call copy_dkoperators(m+1,operators,m,operators)
#else
                  operators(m)=operators(m+1)
#endif
  50            continue
                opcounter=opcounter-1
                ordercounter(j)=ordercounter(j)-1
              endif
              goto 240
            endif
  30      continue
  20    continue
        ilow=ilow+ordercounter(j)
  10  continue
c
      if (dynthrsh.ne.0 .and. dynthrsh.ne.1) then
        write (stdout,2045) dynthrsh
2045    format (/2X,'ERROR in SR simplify3: Parameter dynthrsh = ',I2,
     *          ' has illegal value!',//2X,'dynthrsh has to be 0 or 1.',
     *          //2X,'STOP.',/2X)
        Call Abend
      endif
      if (dynthrsh.eq.0) then
c
 440    continue
        ilow=1
        do 110 j=0,dkhorder
          do 120 k=ilow,ilow+ordercounter(j)-1
            if (abs(doperators(k)).lt.dkhzero) then
              do 130 m=k,opcounter-1
                operleng(m)=operleng(m+1)
                oporder(m,1)=oporder(m+1,1)
                oporder(m,2)=oporder(m+1,2)
                oporder(m,3)=oporder(m+1,3)
                evenodd(m)=evenodd(m+1)
                doperators(m)=doperators(m+1)
#if defined(_MOLCAS_) || defined(MOLPRO)
                call copy_dkoperators(m+1,operators,m,operators)
#else
                operators(m)=operators(m+1)
#endif
 130          continue
              opcounter=opcounter-1
              ordercounter(j)=ordercounter(j)-1
              goto 440
            endif
 120      continue
          ilow=ilow+ordercounter(j)
 110    continue
c
      endif
c
      if (dynthrsh.eq.1) then
c
 540    continue
        ilow=1
        do 210 j=0,dkhorder
          do 220 k=ilow,ilow+ordercounter(j)-1
            if (abs(doperators(k)).lt.dkhzerothrsh(j)) then
              do 230 m=k,opcounter-1
                operleng(m)=operleng(m+1)
                oporder(m,1)=oporder(m+1,1)
                oporder(m,2)=oporder(m+1,2)
                oporder(m,3)=oporder(m+1,3)
                evenodd(m)=evenodd(m+1)
                doperators(m)=doperators(m+1)
#if defined(_MOLCAS_) || defined(MOLPRO)
                call copy_dkoperators(m+1,operators,m,operators)
#else
                operators(m)=operators(m+1)
#endif
 230          continue
              opcounter=opcounter-1
              ordercounter(j)=ordercounter(j)-1
              goto 540
            endif
 220      continue
          ilow=ilow+ordercounter(j)
 210    continue
c
      endif
c
      if (.not.dkhscfflg) then
c
 640    continue
        ilow=1
        do 310 j=0,dkhorder
          do 320 k=ilow,ilow+xordercounter(j)-1
#if defined(_MOLCAS_) || defined(MOLPRO)
            call get_dkoperators(k,xopstring1,xoperators)
#endif
            do 330 l=k+1,ilow+xordercounter(j)-1
#if defined(_MOLCAS_) || defined(MOLPRO)
            call get_dkoperators(l,xopstring2,xoperators)
            if ( (xoperleng(k).eq.xoperleng(l)) .and.
     *               (xoporder(k).eq.xoporder(l)) .and.
     *               xopstring1(1:xoperleng(k)).eq.
     *               xopstring2(1:xoperleng(l)) ) then
#else
            if ( (xoperleng(k).eq.xoperleng(l)) .and.
     *               (xoporder(k).eq.xoporder(l)) .and.
     *               xoperators(k)(1:xoperleng(k)).eq.
     *                   xoperators(l)(1:xoperleng(l)) ) then
#endif
                xdoperators(k)=xdoperators(k)+xdoperators(l)
                do 340 m=l,xopcounter
                  xoperleng(m)=xoperleng(m+1)
                  xoporder(m)=xoporder(m+1)
                  xevenodd(m)=xevenodd(m+1)
                  xdoperators(m)=xdoperators(m+1)
#if defined(_MOLCAS_) || defined(MOLPRO)
                  call copy_dkoperators(m+1,xoperators,m,xoperators)
#else
                  xoperators(m)=xoperators(m+1)
#endif
 340            continue
                xopcounter=xopcounter-1
                xordercounter(j)=xordercounter(j)-1
                if (abs(xdoperators(k)).lt.dkhzero) then
                  do 350 m=k,xopcounter
                    xoperleng(m)=xoperleng(m+1)
                    xoporder(m)=xoporder(m+1)
                    xevenodd(m)=xevenodd(m+1)
                    xdoperators(m)=xdoperators(m+1)
#if defined(_MOLCAS_) || defined(MOLPRO)
                    call copy_dkoperators(m+1,xoperators,m,xoperators)
#else
                    xoperators(m)=xoperators(m+1)
#endif
 350              continue
                  xopcounter=xopcounter-1
                  xordercounter(j)=xordercounter(j)-1
                endif
                goto 640
              endif
 330        continue
 320      continue
          ilow=ilow+xordercounter(j)
 310    continue
c
        if (dynthrsh.eq.0) then
c
 740      continue
          ilow=1
          do 810 j=0,dkhorder
            do 820 k=ilow,ilow+xordercounter(j)-1
              if (abs(xdoperators(k)).lt.dkhzero) then
                do 830 m=k,xopcounter-1
                  xoperleng(m)=xoperleng(m+1)
                  xoporder(m)=xoporder(m+1)
                  xevenodd(m)=xevenodd(m+1)
                  xdoperators(m)=xdoperators(m+1)
#if defined(_MOLCAS_) || defined(MOLPRO)
                  call copy_dkoperators(m+1,xoperators,m,xoperators)
#else
                  xoperators(m)=xoperators(m+1)
#endif
 830            continue
                xopcounter=xopcounter-1
                xordercounter(j)=xordercounter(j)-1
                goto 740
              endif
 820        continue
            ilow=ilow+xordercounter(j)
 810      continue
        endif
c
        if (dynthrsh.eq.1) then
c
 940      continue
          ilow=1
          do 910 j=0,dkhorder
            do 920 k=ilow,ilow+xordercounter(j)-1
              if (abs(xdoperators(k)).lt.dkhzerothrsh(j)) then
                do 930 m=k,xopcounter-1
                  xoperleng(m)=xoperleng(m+1)
                  xoporder(m)=xoporder(m+1)
                  xevenodd(m)=xevenodd(m+1)
                  xdoperators(m)=xdoperators(m+1)
#if defined(_MOLCAS_) || defined(MOLPRO)
                  call copy_dkoperators(m+1,xoperators,m,xoperators)
#else
                  xoperators(m)=xoperators(m+1)
#endif
 930            continue
                xopcounter=xopcounter-1
                xordercounter(j)=xordercounter(j)-1
                goto 940
              endif
 920        continue
            ilow=ilow+xordercounter(j)
 910      continue
        endif
c
      endif
c
CMR      if (i.le.9) then
CMR        write (stdout,3197) i,i,opcounter
CMR3197    format (2X,'Number of operators after exploitation of ',
CMR     *          'commutator',/3X,'symmetry [W',I1,',E0] = -O',I1,
CMR     *          '  and combination of equal terms:',4X,I8)
CMR        if (.not.dkhscfflg) write (stdout,3198) xopcounter
CMR3198    format (2X,'Number of xoperators after combination of equal ',
CMR     *          'terms:',6X,I8)
CMR      endif
CMR      if (i.ge.10) then
CMR        write (stdout,3199) i,i,opcounter
CMR3199    format (2X,'Number of operators after exploitation of ',
CMR     *          'commutator',/3X,'symmetry [W',I2,',E0] = -O',I2,
CMR     *          '  and combination of equal terms:',2X,I8)
CMR        if (.not.dkhscfflg) write (stdout,3200) xopcounter
CMR3200    format (2X,'Number of xoperators after combination of equal ',
CMR     *          'terms:',6X,I8)
CMR      endif
c
      if (dbgflg.ge.1 .and. dkhscfflg) write (dbgunit,4036)
4036  format (/2X,'Terms occurring more than once have been summarized',
     *        ' for operators (SR "simplify3").')
      if (dbgflg.ge.1 .and. .not.dkhscfflg) write (dbgunit,4037)
4037  format (/2X,'Terms occurring more than once have been summarized',
     *        ' for operators and xoperators (SR "simplify3").')
c
      if (dbgflg.ge.2) then
        if (i.le.9) then
          write (dbgunit,5123) i+1,i,i
5123      format (/2X,'Terms of Hamiltonian H',I1,' = [U',I1,' * ',
     *            '(operators) * U',I1,'^t]:',/2X,52('-'))
          call output3 (dbgunit,opcounter,operleng,oporder,evenodd,
     *                  doperators,operators)
          if (.not.dkhscfflg) then
            write (dbgunit,5124) i+1,i,i
5124        format (/2X,'Terms of property X',I1,' = [U',I1,' * ',
     *              '(xoperators) * U',I1,'^t]:',/2X,50('-'))
            call output3b (dbgunit,xopcounter,xoperleng,xoporder,
     *                     xevenodd,xdoperators,xoperators)
          endif
        endif
c
        if (i.ge.10) then
          write (dbgunit,5125) i+1,i,i
5125      format (/2X,'Terms of Hamiltonian H',I2,' = [U',I2,' * ',
     *            '(operators) * U',I2,'^t]:',/2X,55('-'))
          call output3 (dbgunit,opcounter,operleng,oporder,evenodd,
     *                  doperators,operators)
          if (.not.dkhscfflg) then
            write (dbgunit,5126) i+1,i,i
5126        format (/2X,'Terms of property X',I2,' = [U',I2,' * ',
     *              '(xoperators) * U',I2,'^t]:',/2X,53('-'))
            call output3 (dbgunit,xopcounter,xoperleng,xoporder,
     *                    xevenodd,xdoperators,xoperators)
          endif
        endif
      endif
c
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(uniorder)
        call Unused_integer(xorder)
        call Unused_integer_array(wordercounter)
        call Unused_integer_array(wopsleng)
        call Unused_integer_array(woporder)
        call Unused_integer_array(eowops)
        call Unused_real_array(dwops)
        call Unused_character(wops)
      end if
      end
