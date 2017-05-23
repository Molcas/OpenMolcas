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
      subroutine simplify1 (i,uniorder,dkhorder,xorder,dkhscfflg,
     *                 ordercounter,opcounter,operleng,oporder,evenodd,
     *                 doperators,operators,xordercounter,xopcounter,
     *                 xoperleng,xoporder,xevenodd,xdoperators,
     *                 xoperators)
c
c***************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 01.07.2005
c
c   first version: 18.06.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***************************************************************************
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
#if defined(_MOLCAS_) || defined(MOLPRO)
      integer operators(*),xoperators(*)
      character*(maxlength) opstring,xopstring
#else
      character*(*) operators(maxoperators)
     *                      xoperators(maxoperators)
#endif
c
      logical dkhscfflg
C     integer l,k,idum1,idum2
      integer l,  idum1,idum2
      idum1=0
      do 10 l=1,opcounter
        if (i.eq.uniorder) then
  20      continue
          if (evenodd(l+idum1).eq.-1) then
            ordercounter(oporder(l+idum1,3))=
     *        ordercounter(oporder(l+idum1,3))-1
            idum1=idum1+1
            if (evenodd(l+idum1).eq.-1) goto 20
          endif
        else
  30      continue
          if (evenodd(l+idum1).eq.-1.and.oporder(l+idum1,3).eq.i ) then
            ordercounter(oporder(l+idum1,3))=
     *        ordercounter(oporder(l+idum1,3))-1
            idum1=idum1+1
            if (evenodd(l+idum1).eq.-1) goto 30
          endif
        endif
        if (l+idum1.gt.opcounter) goto 10
        oporder(l,1)=oporder(l+idum1,1)
        oporder(l,2)=oporder(l+idum1,2)
        oporder(l,3)=oporder(l+idum1,3)
        evenodd(l)=evenodd(l+idum1)
        doperators(l)=doperators(l+idum1)
        operleng(l)=operleng(l+idum1)
#if defined(_MOLCAS_) || defined(MOLPRO)
        call copy_dkoperators(l+idum1,operators,l,operators)
#else
        operators(l)=operators(l+idum1)
#endif
  10  continue
c
#if defined(_MOLCAS_) || defined(MOLPRO)
      opstring=' '
#endif
      do 40 l=opcounter-idum1+1,opcounter
        oporder(l,1)=0
        oporder(l,2)=0
        oporder(l,3)=0
        evenodd(l)=0
        doperators(l)=0.0d0
#if defined(_MOLCAS_) || defined(MOLPRO)
        call put_dkoperators(l,opstring,operators)
#else
        do 50 k=1,operleng(l)
          operators(l)(k:k)=' '
50      continue
#endif
        operleng(l)=0
  40  continue
c
      opcounter=opcounter-idum1
c
      if (.not.dkhscfflg) then
        idum2=0
        do 100 l=1,xopcounter
          if (i.eq.uniorder) then
 110        continue
            if (xevenodd(l+idum2).eq.-1) then
              xordercounter(xoporder(l+idum2))=
     *          xordercounter(xoporder(l+idum2))-1
              idum2=idum2+1
              if (xevenodd(l+idum2).eq.-1) goto 110
            endif
          else
 120        continue
            if ( xevenodd(l+idum2).eq.-1 .and.
     *           (xoporder(l+idum2).gt.xorder-i-1) ) then
              xordercounter(xoporder(l+idum2))=
     *          xordercounter(xoporder(l+idum2))-1
              idum2=idum2+1
              if (xevenodd(l+idum2).eq.-1) goto 120
            endif
          endif
          if (l+idum2.gt.xopcounter) goto 100
          xoporder(l)=xoporder(l+idum2)
          xevenodd(l)=xevenodd(l+idum2)
          xdoperators(l)=xdoperators(l+idum2)
          xoperleng(l)=xoperleng(l+idum2)
#if defined(_MOLCAS_) || defined(MOLPRO)
          call copy_dkoperators(l+idum2,xoperators,l,xoperators)
#else
          xoperators(l)=xoperators(l+idum2)
#endif
 100    continue
c
#if defined(_MOLCAS_) || defined(MOLPRO)
        xopstring=' '
#endif
        do 140 l=xopcounter-idum2+1,xopcounter
          xoporder(l)=0
          xevenodd(l)=0
          xdoperators(l)=0.0d0
#if defined(_MOLCAS_) || defined(MOLPRO)
        call put_dkoperators(l,xopstring,xoperators)
#else
        do 150 k=1,xoperleng(l)
          xoperators(l)(k:k)=' '
 150    continue
#endif
          xoperleng(l)=0
 140    continue
c
        xopcounter=xopcounter-idum2
      endif
c
c
CAW---------
c      write (stdout,3097) opcounter
c3097  format (2X,'Reduced no. of operators after elimination of ',
c     *        'odd terms:',4X,I8)
c      if (.not.dkhscfflg) write (stdout,3098) xopcounter
c3098  format (2X,'Reduced no. of xoperators after elimination of odd',
c     *        ' terms:',3X,I8)
CAW---------
c
c
c  DEBUGGING: write actual operators to dbgout.32
      if (dbgflg.ge.1) write (dbgunit,4025)
4025  format (/2X,'Irrelevant odd terms have been removed.',
     *        ' (SR "simplify1")')
      if (dbgflg.ge.3) then
        if (i.le.9) then
          write (dbgunit,4030) i+1,i,i
4030      format (/2X,'Terms of Hamiltonian H',I1,' = [U',I1,' * ',
     *            '(operators) * U',I1,'^t]:',/2X,52('-'))
          call output3 (dbgunit,opcounter,operleng,oporder,evenodd,
     *                  doperators,operators)
          if (.not.dkhscfflg) then
            write (dbgunit,4031) i+1,i,i
4031        format (/2X,'Terms of property X',I1,' = [U',I1,' * ',
     *              '(xoperators) * U',I1,'^t]:',/2X,50('-'))
            call output3b (dbgunit,xopcounter,xoperleng,xoporder,
     *                     xevenodd,xdoperators,xoperators)
          endif
        endif
        if (i.ge.10) then
          write (dbgunit,4032) i+1,i,i
4032      format (/2X,'Terms of Hamiltonian H',I2,' = [U',I2,' * ',
     *            '(operators) * U',I2,'^t]:',/2X,55('-'))
          call output3 (dbgunit,opcounter,operleng,oporder,evenodd,
     *                  doperators,operators)
          if (.not.dkhscfflg) then
            write (dbgunit,4033) i+1,i,i
4033        format (/2X,'Terms of property X',I2,' = [U',I2,' * ',
     *              '(xoperators) * U',I2,'^t]:',/2X,53('-'))
            call output3b (dbgunit,xopcounter,xoperleng,xoporder,
     *                     xevenodd,xdoperators,xoperators)
          endif
        endif
      endif
c
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dkhorder)
      end
