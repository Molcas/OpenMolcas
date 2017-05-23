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
      subroutine sort_op2 (dkhorder,opcounter,operleng,oporder,
     *                     evenodd,doperators,operators)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.3
c
c   last modified: 19.10.2005  (MR, ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer dkhorder,opcounter,operleng(maxoperators),
     *        oporder(maxoperators),evenodd(maxoperators),
     *        opcounter2
      REAL*8 doperators(maxoperators)
#if defined(_MOLCAS_) || defined(MOLPRO)
#include "WrkSpc.fh"
      integer operleng2,oporder2,evenodd2,doperators2,operators2
      integer intrea,nwop,lwop
      integer operators(*)
#else
      integer operleng2(maxoperators),oporder2(maxoperators,3),
     *        evenodd2(maxoperators)
c:      character*(maxlength) operators(maxoperators),operators2(maxoperators)
      REAL*8 doperators2(maxoperators)
#endif
C     integer j,k,l
      integer j,  l
c
#if defined(_MOLCAS_) || defined(MOLPRO)
      lwop=8/intrea()
      nwop=(maxlength-1)/lwop+1
      call getmem('doperators2','Allo','Real',doperators2,maxoperators)
      call getmem('operators2','Allo','Inte',operators2,
     *             maxoperators*nwop)
      call getmem('operleng2','Allo','Inte',operleng2,maxoperators)
      call getmem('oporder2','Allo','Inte',oporder2,maxoperators)
      call getmem('evenodd2','Allo','Inte',evenodd2,maxoperators)
c
      opcounter2=0
      do j=0,dkhorder
        do l=1,opcounter
          if (oporder(l).eq.j) then
            iwork(operleng2+opcounter2)=operleng(l)
            iwork(oporder2+opcounter2)=oporder(l)
            iwork(evenodd2+opcounter2)=evenodd(l)
            work(doperators2+opcounter2)=doperators(l)
            opcounter2=opcounter2+1
            call copy_dkoperators(l,operators,
     *             opcounter2,iwork(operators2))
          endif
        end do
      end do
#else
      do 120 j=1,opcounter
        operleng2(j)=0
        oporder2(j)=0
        evenodd2(j)=0
        doperators2(j)=0.0d0
        do 130 k=1,maxlength
          operators2(j)(k:k)=' '
 130    continue
 120  continue
c
c     Write operators in ascending order in V to operators2
c
      opcounter2=0
      do 350 j=0,dkhorder
        do 360 l=1,opcounter
          if (oporder(l).eq.j) then
            opcounter2=opcounter2+1
            operleng2(opcounter2)=operleng(l)
            oporder2(opcounter2)=oporder(l)
            evenodd2(opcounter2)=evenodd(l)
            doperators2(opcounter2)=doperators(l)
            operators2(opcounter2)=operators(l)
          endif
 360    continue
 350  continue
#endif
      if (opcounter2.ne.opcounter) then
        write (stdout,1043) opcounter2,opcounter
1043    format (/2X,'ERROR in sort_op1: opcounter2 = ',I8,
     *          ' not equal to opcounter = ',I8,'.',//2X,'STOP.',/)
            CALL Abend
      endif
#if defined(_MOLCAS_) || defined(MOLPRO)
      do l=1,opcounter
        operleng(l)=iwork(operleng2-1+l)
        oporder(l)=iwork(oporder2-1+l)
        evenodd(l)=iwork(evenodd2-1+l)
        doperators(l)=work(doperators2-1+l)
        call copy_dkoperators(l,iwork(operators2),l,operators)
      end do
      call getmem('evenodd2','Free','Inte',evenodd2,maxoperators)
      call getmem('oporder2','Free','Inte',oporder2,maxoperators)
      call getmem('operleng2','Free','Inte',operleng2,maxoperators)
      call getmem('operators2','Free','Inte',operators2,
     *             maxoperators*nwop)
      call getmem('doperators2','Free','Real',doperators2,maxoperators)
#else
      do 361 l=1,opcounter
        operleng(l)=operleng2(l)
        oporder(l)=oporder2(l)
        evenodd(l)=evenodd2(l)
        doperators(l)=doperators2(l)
        operators(l)=operators2(l)
 361  continue
#endif
c
      return
      end
