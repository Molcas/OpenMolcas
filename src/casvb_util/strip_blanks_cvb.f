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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine strip_blanks_cvb(line,lenline,blanks,nblank,blankdelim)
c  If BLANKDELIM, a given number of blanks (not leading or trailing)
c  will be subsituted by a single blank. Otherwise all blanks are
c  stripped.
      implicit real*8 (a-h,o-z)
      logical blankdelim
      character*(*) line
      character*1 blanks(nblank)
#include "malloc_cvb.fh"

      do 100 iblank=1,nblank
      if(blanks(iblank).ne.' ')then
        do 200 ich=1,lenline
        if(line(ich:ich).eq.blanks(iblank))line(ich:ich)=' '
200     continue
      endif
100   continue
      ilv=mstacki_cvb(lenline)
      ich2=0
      do 300 ich=1,lenline
      if(line(ich:ich).ne.' ')then
        ich2=ich2+1
        iw(ich2+ilv-1)=ich
c  (Final condition eliminates leading blanks :)
      elseif(blankdelim.and.ich.ge.2.and.line(ich-1:ich-1).ne.' '
     >  .and.ich2.ne.0)then
        ich2=ich2+1
        iw(ich2+ilv-1)=ich
      endif
300   continue
      do 400 ich=1,ich2
      line(ich:ich)=line(iw(ich+ilv-1):iw(ich+ilv-1))
400   continue
      lenline=ich2
      call mfreei_cvb(ilv)
      return
      end
