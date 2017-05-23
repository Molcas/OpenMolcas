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
      subroutine simplify2 (i,uniorder,dkhorder,ordercounter,opcounter,
     *               operleng,oporder,evenodd,doperators,operators,
     *               wordercounter,wopsleng,woporder,eowops,dwops,wops,
     *               oddcounter,oddleng,oddorder,eoodd,dodd,odd)
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
c   first version: 23.05.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer i,uniorder,dkhorder,ordercounter(0:maxorder),opcounter,
     *        operleng(maxoperators),oporder(maxoperators,3),
     *        evenodd(maxoperators),
     *        wordercounter(maxorder),wopsleng(maxuops),
     *        woporder(maxuops,3),eowops(maxuops),oddcounter,
     *        oddleng(maxuops),oddorder(maxuops,3),eoodd(maxuops)
      REAL*8 doperators(maxoperators),dwops(maxuops),
     *                 dodd(maxuops)
      character*(maxlength) wops(maxuops),odd(maxuops)
#if defined(_MOLCAS_) || defined(MOLPRO)
      character*(maxlength) opstring
      integer operators(*)
#else
      character*(maxlength) operators(maxoperators)
#endif
c
      integer scrleng1,scrleng2,reslengl1,reslengr1,reslengl2,reslengr2
      character*(maxlength) scrchar1,scrchar2,rescharl1,rescharr1,
     *                      rescharl2,rescharr2
c
      integer j,k,l,ilow,hit1,hit2
      character*(3) scrchar0, dkh_int2char
      logical match,equalstring
c
      if(i.ge.1.and.i.le.99) then
          scrchar0=dkh_int2char(i)
          scrchar0(1:1)='W'
      endif
cvv      do 10 k=1,99
cvv        if (i.eq.k) then
cvv          scrchar0=dkh_int2char(k)
cvv          scrchar0(1:1)='W'
cvv        endif
cvv 10   continue
      if (uniorder.gt.99) then
        write (stdout,1002)
1002    format (/2X,'Note: So far only 99 unitary transformations U_i',
     *          ' are possible in SR "simplify2."')
        Call Abend
      endif
      do 15 k=1,maxlength
        scrchar1(k:k)=' '
        scrchar2(k:k)=' '
  15  continue
      scrleng1=6
      scrleng2=6
      scrchar1(1:scrleng1)=scrchar0(1:3)//'E00'
      scrchar2(1:scrleng2)='E00'//scrchar0(1:3)
 412  continue
      match=.false.
      ilow=1
      do 100 j=0,dkhorder
        do 200 k=ilow,ilow+ordercounter(j)-1
#if defined(_MOLCAS_) || defined(MOLPRO)
          call get_dkoperators(k,opstring,operators)
          hit1=index(opstring(1:operleng(k)),scrchar1(1:6))
#else
          hit1=index(operators(k)(1:operleng(k)),scrchar1(1:6))
#endif
          if (hit1.ne.0) then
#if defined(_MOLCAS_) || defined(MOLPRO)
            call store_reschar (operleng(k),opstring,hit1,scrleng1,
     *                          reslengl1,reslengr1,rescharl1,rescharr1)
#else
            call store_reschar (operleng(k),operators(k),hit1,scrleng1,
     *                          reslengl1,reslengr1,rescharl1,rescharr1)
#endif
            do 300 l=ilow,ilow+ordercounter(j)-1
#if defined(_MOLCAS_) || defined(MOLPRO)
              call get_dkoperators(l,opstring,operators)
              hit2=index(opstring(1:operleng(l)),scrchar2(1:6))
#else
              hit2=index(operators(l)(1:operleng(l)),scrchar2(1:6))
#endif
              if (hit2.ne.0) then
#if defined(_MOLCAS_) || defined(MOLPRO)
                call store_reschar (operleng(l),opstring,hit2,
     *                              scrleng2,reslengl2,
     *                              reslengr2,rescharl2,rescharr2)
#else
                call store_reschar (operleng(l),operators(l),hit2,
     *                              scrleng2,reslengl2,
     *                              reslengr2,rescharl2,rescharr2)
#endif
                match=equalstring(reslengl1,reslengr1,rescharl1,
     *                            rescharr1,reslengl2,reslengr2,
     *                            rescharl2,rescharr2)
              endif
              if (match) then
                call replace1 (i,j,k,l,ordercounter,opcounter,operleng,
     *                         oporder,evenodd,doperators,operators,
     *                         reslengl1,reslengr1,rescharl1,rescharr1,
     *                         oddcounter,oddleng,oddorder,eoodd,dodd,
     *                         odd)
                goto 412
              endif
 300        continue
            do 400 l=opcounter,ilow+ordercounter(j),-1
              operleng(l+2)=operleng(l)
              oporder(l+2,1)=oporder(l,1)
              oporder(l+2,2)=oporder(l,2)
              oporder(l+2,3)=oporder(l,3)
              evenodd(l+2)=evenodd(l)
              doperators(l+2)=doperators(l)
#if defined(_MOLCAS_) || defined(MOLPRO)
              call copy_dkoperators(l,operators,l+2,operators)
#else
              operators(l+2)=operators(l)
#endif
 400        continue
c                                              and ilow+ordercounter(j)+1
            operleng(ilow+ordercounter(j))=operleng(k)
            oporder(ilow+ordercounter(j),1)=oporder(k,1)
            oporder(ilow+ordercounter(j),2)=oporder(k,2)
            oporder(ilow+ordercounter(j),3)=oporder(k,3)
            evenodd(ilow+ordercounter(j))=evenodd(k)
            doperators(ilow+ordercounter(j))=-1.d0*doperators(k)
#if defined(_MOLCAS_) || defined(MOLPRO)
            opstring=rescharl1(1:reslengl1)//
     *               scrchar2(1:scrleng2)//rescharr1(1:reslengr1)
            call put_dkoperators(ilow+ordercounter(j),opstring,
     *               operators)
#else
            operators(ilow+ordercounter(j))=rescharl1(1:reslengl1)//
     *                    scrchar2(1:scrleng2)//rescharr1(1:reslengr1)
#endif
c
            operleng(ilow+ordercounter(j)+1)=operleng(k)
            oporder(ilow+ordercounter(j)+1,1)=oporder(k,1)
            oporder(ilow+ordercounter(j)+1,2)=oporder(k,2)
            oporder(ilow+ordercounter(j)+1,3)=oporder(k,3)
            evenodd(ilow+ordercounter(j)+1)=evenodd(k)
            doperators(ilow+ordercounter(j)+1)=doperators(k)
#if defined(_MOLCAS_) || defined(MOLPRO)
            opstring=rescharl1(1:reslengl1)//
     *               scrchar2(1:scrleng2)//rescharr1(1:reslengr1)
            call put_dkoperators(ilow+ordercounter(j)+1,opstring,
     *               operators)
#else
            operators(ilow+ordercounter(j)+1)=rescharl1(1:reslengl1)//
     *                    scrchar2(1:scrleng2)//rescharr1(1:reslengr1)
#endif
            ordercounter(j)=ordercounter(j)+2
            opcounter=opcounter+2
            goto 412
          endif
 200    continue
        ilow=ilow+ordercounter(j)
 100  continue
c
      if (dbgflg.ge.1 .and. i.le.9) write (dbgunit,4026) i,i
4026  format (/2X,'Commutator symmetry [W',I1,',E0] = -O',I1,' has ',
     *        'been exploited for operators.  (SR "simplify2")')
      if (dbgflg.ge.1 .and. i.ge.10) write (dbgunit,4027) i,i
4027  format (/2X,'Commutator symmetry [W',I2,',E0] = -O',I2,' has ',
     *        'been exploited for operators.  (SR "simplify2")')
      if (dbgflg.ge.3) then
        write (dbgunit,1026) i+1,i,i
1026    format (/2X,'Terms of Hamiltonian H',I1,' = [U',I1,' * ',
     *            '(operators) * U',I1,'^t]:',/2X,52('-'))
        call output3 (dbgunit,opcounter,operleng,oporder,evenodd,
     *                doperators,operators)
      endif
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer_array(wordercounter)
        call Unused_integer_array(wopsleng)
        call Unused_integer_array(woporder)
        call Unused_integer_array(eowops)
        call Unused_real_array(dwops)
        call Unused_character(wops)
      end if
      end
