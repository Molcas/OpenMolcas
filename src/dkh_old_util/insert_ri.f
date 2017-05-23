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
      subroutine insert_ri (length,operator)
c
c***********************************************************************
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
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer length,k,istart,v1pos,v2pos,p1pos,p2pos,idum1,idum2,idum3,
     *        idum4, iMax
      character*(maxlength) operator
c
      istart=1
      iMax=0
c
 300  continue
c
      iMax = iMax + 1
      If (iMax.gt.30) Call Abend
      v1pos=0
      v2pos=0
      p1pos=0
      p2pos=0
      idum1=index(operator(istart:length),'V')
      idum4=index(operator(istart:length),'X')
      if ((idum4.gt.0.and.idum4.lt.idum1) .or. idum1.eq.0) idum1=idum4
      if (idum1.gt.0) v1pos=idum1+istart-1
      idum2=index(operator(istart:length),'E01')
      idum4=index(operator(istart:length),'CE0')
      if ((idum4.gt.0.and.idum4.lt.idum2) .or. idum2.eq.0) idum2=idum4
      if (idum2.gt.0) then
        if (v1pos.eq.0) then
          v1pos=idum2+istart-1
        else
          if (idum2+istart-1.lt.v1pos) v1pos=idum2+istart-1
        endif
      endif
      idum3=index(operator(istart:length),'S')
      if (idum3.gt.0) then
        if (v1pos.eq.0) then
          v1pos=idum3+istart-1
        else
          if (idum3+istart-1.lt.v1pos) v1pos=idum3+istart-1
        endif
      endif
      if (v1pos.eq.0) then
        goto 200
      endif
c
      idum1=index(operator(istart:length),'P')
      if (idum1.gt.0) p1pos=idum1+istart-1
      if (v1pos.lt.p1pos) then
        istart=v1pos+1
        goto 300
c
      else
c
        idum1=index(operator(v1pos+1:length),'V')
        idum4=index(operator(v1pos+1:length),'X')
        if ((idum4.gt.0.and.idum4.lt.idum1) .or. idum1.eq.0) idum1=idum4
        if (idum1.gt.0) v2pos=v1pos+idum1
        idum2=index(operator(v1pos+1:length),'E01')
        idum4=index(operator(v1pos+1:length),'CE0')
        if ((idum4.gt.0.and.idum4.lt.idum2) .or. idum2.eq.0) idum2=idum4
        if (idum2.gt.0) then
          if (v2pos.eq.0) then
            v2pos=v1pos+idum2
          else
            if (idum2+v1pos.lt.v2pos) v2pos=v1pos+idum2
          endif
        endif
        idum3=index(operator(v1pos+1:length),'S')
        if (idum3.gt.0) then
          if (v2pos.eq.0) then
            v2pos=v1pos+idum3
          else
            if (idum3+v1pos.lt.v2pos) v2pos=v1pos+idum3
          endif
        endif
c
        if (p1pos.eq.0) then
          p2pos=0
        else
          idum1=index(operator(p1pos+1:length),'P')
          if (idum1.gt.0) p2pos=idum1+p1pos
        endif
c
        if (v2pos.lt.p2pos .and. v2pos.ne.0) then
          if (length+3.gt.maxlength) then
            write (stdout,2002) maxlength
2002        format (/2X,'ERROR in insert_ri: maxlength = ',I3,
     *              ' is to small.'//2X,'Increase it in parameters.h.',
     *              //2X,'STOP.',/)
            CALL Abend
          endif
          do 100 k=length,v2pos,-1
             operator(k+3:k+3)=operator(k:k)
 100      continue
          operator(v2pos:v2pos+2)='PQP'
          istart=v2pos+2
          length=length+3
          goto 300
        else
          if (p1pos.eq.0.and.p2pos.eq.0) then
            goto 200
          else
            istart=p2pos+1
            goto 300
          endif
        endif
c
      endif
c
 200  continue
c
      return
      end
