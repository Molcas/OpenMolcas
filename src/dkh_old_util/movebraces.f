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
      subroutine movebraces (length,operator)
c
c*************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 13.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c*************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
C     integer length,idum,posp1,posp2,istart,posq,pospp
      integer length,     posp1,posp2,istart,posq,pospp
      character*(maxlength) operator
CMR      logical action
cAOM<
c... AOM: logical flag is introduced
c...      there was a bug with index boundaries
      logical action,flag
cAOM>
c
      istart=1
1000  continue
      action=.false.
      if (index(operator(istart:length),'P').gt.0) then
        posp1=istart-1+index(operator(istart:length),'P')
      else
        istart=1
        goto 2000
      endif
      if (index(operator(posp1+1:length),'P').gt.0) then
        posp2=posp1+index(operator(posp1+1:length),'P')
      else
        istart=1
        goto 2000
      endif
CMR<
c      if ( operator(posp1+1:posp1+1).ne.'V' .and.
c     *       operator(posp1+1:posp1+1).ne.'E' .and.
c     *       operator(posp1+1:posp1+1).ne.'S' .and.
c     *       operator(posp1+1:posp1+2).ne.'[P' .and.
c     *       operator(posp1+1:posp1+1).ne.'X' .and.
c     *       operator(posp1+1:posp1+1).ne.'C' .and.
c     *       posp1.lt.posp2-1) then
CMR the following 3 lines may be erased once the code has been tested
CMR      if(posp1.lt.length-1)
CMR     *  write(*,*) "WARNING in movebraces: posp1.lt.length-1, posp1=",
CMR     *             posp1
      flag=.true.
      if (operator(posp1+1:posp1+1).eq.'V') flag=.false.
      if (operator(posp1+1:posp1+1).eq.'E') flag=.false.
      if (operator(posp1+1:posp1+1).eq.'S') flag=.false.
      if (operator(posp1+1:posp1+2).eq.'[P') flag=.false.
      if (operator(posp1+1:posp1+1).eq.'X') flag=.false.
      if (operator(posp1+1:posp1+1).eq.'C'.and.
     *    posp1.lt.length-1) flag=.false.
      if (posp1.ge.posp2-1) flag=.false.
      if (flag) then
        action=.true.
        operator(posp1:posp1)=operator(posp1+1:posp1+1)
        operator(posp1+1:posp1+1)='P'
      endif
CMR>

CMR<
c      if (operator(posp2-1:posp2-1).ne.'V' .and.
c     *      operator(posp2-3:posp2-1).ne.'E01' .and.
c     *      operator(posp2-2:posp2-1).ne.'P]' .and.
c     *      operator(posp2-4:posp2-4).ne.'S' .and.
c     *      operator(posp2-1:posp2-1).ne.'X' .and.
c     *      operator(posp2-3:posp2-1).ne.'CE0' .and.
c     *      posp1.lt.posp2-1) then
CMR the following line may be erased once the code has been tested
CMR      write(*,*) "WARNING in movebraces: posp2.ge.length ..., posp2=",posp2
      flag=.true.
      if (operator(posp2-1:posp2-1).eq.'V') flag=.false.
      if (operator(posp2-1:posp2-1).eq.'X') flag=.false.
      if (posp2-4.gt.0) then
        if (operator(posp2-4:posp2-4).eq.'S') flag=.false.
        if (operator(posp2-3:posp2-1).eq.'E01') flag=.false.
        if (operator(posp2-3:posp2-1).eq.'CE0') flag=.false.
        if (operator(posp2-2:posp2-1).eq.'P]') flag=.false.
      else if (posp2-3.gt.0) then
        if (operator(posp2-3:posp2-1).eq.'E01') flag=.false.
        if (operator(posp2-3:posp2-1).eq.'CE0') flag=.false.
        if (operator(posp2-2:posp2-1).eq.'P]') flag=.false.
      else if (posp2-2.gt.0) then
        if (operator(posp2-2:posp2-1).eq.'P]') flag=.false.
      end if
      if (posp1.ge.posp2-1) flag=.false.
      if (flag) then
        action=.true.
        operator(posp2:posp2)=operator(posp2-1:posp2-1)
        operator(posp2-1:posp2-1)='P'
        goto 1000
      endif
CMR>
      if (action) goto 1000
      istart=posp2+1
      if (istart.lt.length-2) goto 1000
      istart=1
2000  continue
      action=.false.
      if (index(operator(istart:length),'Q').gt.0) then
        posq=istart-1+index(operator(istart:length),'Q')
      else
        istart=1
        goto 3000
      endif
cAOM<
      flag=.true.
      if(posq.le.1) flag=.false.
CMR the following 2 lines may be erased once the code has been tested
CMR      if(posq.le.1)
CMR     *  write(*,*) "WARNING in movebraces: posq.le.1, posq=",posq
      if (flag.and.operator(posq-1:posq-1).eq.'[') then
        action=.true.
        operator(posq-1:posq-1)='Q'
        operator(posq:posq)='['
        istart=posq-1
      endif
      flag=.true.
      if(posq.ge.length) flag=.false.
CMR the following 2 lines may be erased once the code has been tested
CMR      if(posq.ge.length)
CMR     *  write(*,*) "WARNING in movebraces: posq.ge.length, posq=",posq
      if (flag.and.operator(posq+1:posq+1).eq.']') then
        action=.true.
        operator(posq+1:posq+1)='Q'
        operator(posq:posq)=']'
        istart=posq+1
      endif
cAOM>
      if (action) goto 2000
      istart=posq+2
      if (istart.lt.length-2) goto 2000
      istart=1
3000  continue
      action=.false.
      if (index(operator(istart:length),'PP').gt.0) then
        pospp=istart-1+index(operator(istart:length),'PP')
      else
        goto 4000
      endif
cAOM<
      flag=.true.
      if(pospp.le.1.or.pospp.ge.length) flag=.false.
CMR the following 3 lines may be erased once the code has been tested
CMR      if(pospp.ge.(length-1))
CMR     * write(*,*) "WARNING in movebraces: pospp.le.1.or.pospp.ge.length"
CMR     *            //", pospp=",pospp
      if (flag.and.operator(pospp-1:pospp-1).eq.'[') then
        action=.true.
        operator(pospp-1:pospp)='PP'
        operator(pospp+1:pospp+1)='['
        istart=pospp-1
      endif
      flag=.true.
      if(pospp.ge.(length-1)) flag=.false.
CMR the following 3 lines may be erased once the code has been tested
CMR      if(pospp.ge.(length-1))
CMR     *  write(*,*) "WARNING in movebraces: pospp.ge.(length-1), pospp=",
CMR     *             pospp
      if (flag.and.operator(pospp+2:pospp+2).eq.']') then
        action=.true.
        operator(pospp+1:pospp+2)='PP'
        operator(pospp:pospp)=']'
        istart=pospp+1
      endif
cAOM>
      if (action) goto 3000
      istart=pospp+2
      if (istart.lt.length-2) goto 3000
c
4000  continue
c
      return
      end
