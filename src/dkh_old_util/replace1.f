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
      subroutine replace1 (i,j,k,l,ordercounter,opcounter,operleng,
     *                     oporder,evenodd,doperators,operators,
     *                     reslengl,reslengr,rescharl,rescharr,
     *                     oddcounter,oddleng,oddorder,eoodd,dodd,odd)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 30.06.2005
c
c   first version: 18.06.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      real*8 one
      parameter (one=+1.0d0)
c
      integer i,j,k,l,ordercounter(0:maxorder),opcounter,
     *        operleng(maxoperators),oporder(maxoperators,3),
     *        evenodd(maxoperators),
     *        oddcounter,oddleng(maxuops),oddorder(maxuops,3),
     *        eoodd(maxuops),reslengl,reslengr
      real*8 doperators(maxoperators),dodd(maxuops),
     *                 scrcoeff1,scrcoeff2
      character*(maxlength) odd(maxuops),rescharl,rescharr
#if defined(_MOLCAS_) || defined(MOLPRO)
      character*(maxlength) opstring
      character operators(*)
#else
      character*(maxlength) operators(maxoperators)
#endif
C     integer idum1,idum2,idum3,m,n
      integer idum1,idum2,idum3,m
      logical ctrlflg
      intrinsic ABS,DBLE,SIGN
c
      idum1=j
      idum2=k
      idum3=l
      scrcoeff1=doperators(k)
      scrcoeff2=doperators(l)
      ctrlflg=.false.
c
      if (abs(doperators(idum2)+doperators(idum3)).lt.dkhzero) then
        if(idum2.eq.idum3) then
          write (stdout,3012) idum2
3012      format (/2X,'ERROR in subroutine "replace1": idum2 = idum3 =',
     *          I6,//2X,'STOP.',/)
          Call Abend
        endif
        do 20 m=idum3,opcounter
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
  20    continue
        ordercounter(idum1)=ordercounter(idum1)-1
        opcounter=opcounter-1
        if (idum2.gt.idum3) idum2=idum2-1
        do 30 m=opcounter,idum2+1,-1
          operleng(m+oddcounter-1)=operleng(m)
          oporder(m+oddcounter-1,1)=oporder(m,1)
          oporder(m+oddcounter-1,2)=oporder(m,2)
          oporder(m+oddcounter-1,3)=oporder(m,3)
          evenodd(m+oddcounter-1)=evenodd(m)
          doperators(m+oddcounter-1)=doperators(m)
#if defined(_MOLCAS_) || defined(MOLPRO)
          call copy_dkoperators(m,operators,m+oddcounter-1,operators)
#else
          operators(m+oddcounter-1)=operators(m)
#endif
  30    continue
        do 40 m=1,oddcounter
          operleng(idum2+m-1)=reslengl+oddleng(m)+reslengr
          oporder(idum2+m-1,1)=oporder(idum2,1)
          oporder(idum2+m-1,2)=oporder(idum2,2)
          oporder(idum2+m-1,3)=oporder(idum2,3)
          evenodd(idum2+m-1)=evenodd(idum2)
          doperators(idum2+m-1)=scrcoeff1*(-1.0d0)*dodd(m)
#if defined(_MOLCAS_) || defined(MOLPRO)
          call concatenate (operleng(idum2+m-1),opstring,
     *                      reslengl,rescharl,oddleng(m),odd(m),
     *                      reslengr,rescharr)
          call put_dkoperators(idum2+m-1,opstring,operators)
#else
          call concatenate (operleng(idum2+m-1),operators(idum2+m-1),
     *                      reslengl,rescharl,oddleng(m),odd(m),
     *                      reslengr,rescharr)
#endif
  40    continue
        ordercounter(idum1)=ordercounter(idum1)+oddcounter-1
        opcounter=opcounter+oddcounter-1
        ctrlflg=.true.
        goto 6734
      endif
      if ((abs(doperators(idum2))-abs(doperators(idum3)))
     *                                                 .gt.dkhzero) then
        if(idum2.eq.idum3) then
          write (stdout,3022) idum2
 3022     format (/2X,'ERROR in subroutine "replace1": idum2 = idum3 =',
     *          I6,//2X,'STOP.',/)
          Call Abend
        endif
c
        doperators(idum2)=scrcoeff1+scrcoeff2
        do 120 m=idum3,opcounter
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
 120    continue
        ordercounter(idum1)=ordercounter(idum1)-1
        opcounter=opcounter-1
        if (idum2.gt.idum3) idum2=idum2-1
        do 130 m=opcounter,idum2+1,-1
          operleng(m+oddcounter)=operleng(m)
          oporder(m+oddcounter,1)=oporder(m,1)
          oporder(m+oddcounter,2)=oporder(m,2)
          oporder(m+oddcounter,3)=oporder(m,3)
          evenodd(m+oddcounter)=evenodd(m)
          doperators(m+oddcounter)=doperators(m)
#if defined(_MOLCAS_) || defined(MOLPRO)
          call copy_dkoperators(m,operators,m+oddcounter,operators)
#else
          operators(m+oddcounter)=operators(m)
#endif
 130    continue
        do 140 m=1,oddcounter
          operleng(idum2+m)=reslengl+oddleng(m)+reslengr
          oporder(idum2+m,1)=oporder(idum2,1)
          oporder(idum2+m,2)=oporder(idum2,2)
          oporder(idum2+m,3)=oporder(idum2,3)
          evenodd(idum2+m)=evenodd(idum2)
          doperators(idum2+m)=scrcoeff2*dodd(m)
#if defined(_MOLCAS_) || defined(MOLPRO)
          call concatenate (operleng(idum2+m),opstring,
     *                      reslengl,rescharl,oddleng(m),odd(m),
     *                      reslengr,rescharr)
          call put_dkoperators(idum2+m,opstring,operators)
#else
          call concatenate (operleng(idum2+m),operators(idum2+m),
     *                      reslengl,rescharl,oddleng(m),odd(m),
     *                      reslengr,rescharr)
#endif
 140    continue
        ordercounter(idum1)=ordercounter(idum1)+oddcounter
        opcounter=opcounter+oddcounter
        ctrlflg=.true.
        goto 6734
      endif
c
      if ((abs(doperators(idum3))-abs(doperators(idum2)))
     *                                                 .gt.dkhzero) then
        if(idum2.eq.idum3) then
          write (stdout,3032) idum2
 3032     format (/2X,'ERROR in subroutine "replace1": idum2 = idum3 =',
     *          I6,//2X,'STOP.',/)
          Call Abend
        endif

        doperators(idum3)=scrcoeff1+scrcoeff2
        do 220 m=idum2,opcounter
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
 220    continue
        ordercounter(idum1)=ordercounter(idum1)-1
        opcounter=opcounter-1
        if (idum3.gt.idum2) idum3=idum3-1
        do 230 m=opcounter,idum3+1,-1
          operleng(m+oddcounter)=operleng(m)
          oporder(m+oddcounter,1)=oporder(m,1)
          oporder(m+oddcounter,2)=oporder(m,2)
          oporder(m+oddcounter,3)=oporder(m,3)
          evenodd(m+oddcounter)=evenodd(m)
          doperators(m+oddcounter)=doperators(m)
#if defined(_MOLCAS_) || defined(MOLPRO)
          call copy_dkoperators(m,operators,m+oddcounter,operators)
#else
          operators(m+oddcounter)=operators(m)
#endif
 230    continue
        do 240 m=1,oddcounter
          operleng(idum3+m)=reslengl+oddleng(m)+reslengr
          oporder(idum3+m,1)=oporder(idum3,1)
          oporder(idum3+m,2)=oporder(idum3,2)
          oporder(idum3+m,3)=oporder(idum3,3)
          evenodd(idum3+m)=evenodd(idum3)
          doperators(idum3+m)=scrcoeff1*(-1.0d0)*dodd(m)
#if defined(_MOLCAS_) || defined(MOLPRO)
          call concatenate (operleng(idum3+m),opstring,
     *                      reslengl,rescharl,oddleng(m),odd(m),
     *                      reslengr,rescharr)
          call put_dkoperators(idum3+m,opstring,operators)
#else
          call concatenate (operleng(idum3+m),operators(idum3+m),
     *                      reslengl,rescharl,oddleng(m),odd(m),
     *                      reslengr,rescharr)
#endif
 240    continue
        ordercounter(idum1)=ordercounter(idum1)+oddcounter
        opcounter=opcounter+oddcounter
        ctrlflg=.true.
        goto 6734
      endif
c
6734  continue
c
      if (.not.ctrlflg) then
        write (stdout,6327)
 6327   format (/2X,'ERROR in subroutine replace1: None of the 3 cases',
     *          ' has been detected.',//2X,'STOP.',/)
        CALL Abend
      endif
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(i)
        call Unused_integer_array(oddorder)
        call Unused_integer_array(eoodd)
      end if
      end
