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
      subroutine substitute_term (idum,knumber,hits,leng,term,sused,
     *                            stimes,wstimes,s,scrleng,scrchar)
c
c******************************************************************************
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
c******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer idum,knumber,hits,leng
c
      integer sused,stimes(maxsnumber),wstimes(maxuops,maxsnumber),
     *        scrleng(maxsnumber)
      character*(maxlength) term
      character*(4) s(maxsnumber)
      character*(9) scrchar(maxsnumber)
c
      integer i,l,m,n,pos,istart,idum2
      logical action
c
      if (idum.eq.1 .and. leng.lt.3) Then
         Write (stdout,*) 'ERROR1 in substitute_term'
         Call Abend
      End If
      if (idum.eq.2 .and. leng.lt.6) Then
         Write (stdout,*) 'ERROR2 in substitute_term'
         Call Abend
      End If
c
      if (hits.eq.0) then
        if (idum.eq.1) Then
           n=(leng-3)/6
        Else if (idum.eq.2) Then
           n=leng/6
        Else
           n=0 ! Dummy define
           Write (stdout,*) 'Substitute_term: Illegal idum value'
           Call Abend
        End If
        if ((idum.eq.1 .and. abs(DBLE(n)-DBLE(leng-3)/6.d0).gt.dkhzero)
     *          .or. (idum.eq.2 .and. abs(DBLE(n)-DBLE(leng)/6.d0)
     *                                  .gt.dkhzero) )  then
          write (stdout,1034) hits,idum,leng
1034      format(/2X,'ERROR3 in subroutine substitute_term: hits = ',I2,
     *          //2X,'idum = ',I1,//2X,'leng = ',I2,' should be ',
     *               'equal to 6*n+3 (idum=1) or 6*n (idum=2).',
     *          //'STOP.',/)
          Call Abend
        endif
        istart=1
        do 10 i=1,n
          action=.false.
          do 20 l=1,sused
            if (scrleng(l).eq.6 .and.
     *            scrchar(l).eq.term(istart:istart+5)) then
              action=.true.
              term(istart:istart+3)=s(l)
              if (idum.eq.1) then
                wstimes(knumber,l)=wstimes(knumber,l)+1
                if (stimes(l).eq.0) stimes(l)=stimes(l)+1
              endif
              if (idum.eq.2) stimes(l)=stimes(l)+1
              do 30 m=istart+4,leng-2
                term(m:m)=term(m+2:m+2)
  30          continue
              term(leng-1:leng)='  '
              leng=leng-2
              istart=istart+4
              goto 10
            endif
  20      continue
          if (.not.action) istart=istart+6
  10    continue
      endif
c
c
      if (hits.eq.1) then
        if (idum.eq.1) Then
           n=(leng-6)/6
        Else if (idum.eq.2) Then
            n=(leng-3)/6
        Else
           n=0 ! Dummy define
           Write (stdout,*) 'Substitute_term: Illegal idum value'
           Call Abend
        End If
        if ((idum.eq.1 .and. abs(DBLE(n)-DBLE(leng-6)/6.d0).gt.dkhzero)
     *          .or. (idum.eq.2 .and. abs(DBLE(n)-DBLE(leng-3)/6.d0)
     *                                  .gt.dkhzero) )  then
          write (stdout,1036) hits,idum,leng
1036      format (/2X,'ERROR in subroutine substitute_term: hits = ',I2,
     *            //2X,'idum = ',I1,//2X,'leng = ',I2,' should be ',
     *            'equal to 6*n+6 (idum=1) or 6*n+3 (idum2).',
     *            //'STOP.',/)
          Call Abend
        endif
        istart=1
        do 60 i=1,n
          action=.false.
          pos=index(term(istart:leng),'E01')
          idum2=index(term(istart:leng),'CE0')
          if ((idum2.lt.pos.and.idum2.gt.0).or.pos.eq.0) pos=idum2
c
          if (pos.ge.7 .or. pos.eq.0) then
            do 70 l=1,sused
              if (scrleng(l).eq.6 .and.
     *              scrchar(l).eq.term(istart:istart+5)) then
                action=.true.
                term(istart:istart+3)=s(l)
                if (idum.eq.1) then
                  wstimes(knumber,l)=wstimes(knumber,l)+1
                  if (stimes(l).eq.0) stimes(l)=stimes(l)+1
                endif
                if (idum.eq.2) stimes(l)=stimes(l)+1
                do 80 m=istart+4,leng-2
                  term(m:m)=term(m+2:m+2)
  80            continue
                term(leng-1:leng)='  '
                leng=leng-2
                istart=istart+4
                goto 60
              endif
  70        continue
            if (.not.action) istart=istart+6
          endif
c
          if (pos.eq.1) then
            do 90 l=1,sused
              if (scrleng(l).eq.6 .and.
     *              scrchar(l).eq.term(istart+3:istart+8)) then
                action=.true.
                term(istart+3:istart+6)=s(l)
                if (idum.eq.1) then
                  wstimes(knumber,l)=wstimes(knumber,l)+1
                  if (stimes(l).eq.0) stimes(l)=stimes(l)+1
                endif
                if (idum.eq.2) stimes(l)=stimes(l)+1
                do 100 m=istart+7,leng-2
                  term(m:m)=term(m+2:m+2)
 100            continue
                term(leng-1:leng)='  '
                leng=leng-2
                istart=istart+7
                goto 60
              endif
  90        continue
            if (.not.action) istart=istart+9
          endif
c
          if (pos.eq.4) then
            do 110 l=1,sused
              if (scrleng(l).eq.9 .and.
     *              scrchar(l).eq.term(istart:istart+8)) then
                action=.true.
                term(istart:istart+3)=s(l)
                if (idum.eq.1) then
                  wstimes(knumber,l)=wstimes(knumber,l)+1
                  if (stimes(l).eq.0) stimes(l)=stimes(l)+1
                endif
                if (idum.eq.2) stimes(l)=stimes(l)+1
                do 120 m=istart+4,leng-5
                  term(m:m)=term(m+5:m+5)
 120            continue
                term(leng-4:leng)='     '
                leng=leng-5
                istart=istart+4
                goto 60
              endif
 110        continue
            if (.not.action) istart=istart+9
          endif
c
  60    continue
      endif
c
      return
      end
