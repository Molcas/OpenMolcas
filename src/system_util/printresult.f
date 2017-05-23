************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine PrintResult(iUnit,FMT,STR,iCount,STR2,Value,iRank)
c routine to print result in the form:
c      Write(iUnit,FMT) STR,iCount,STR2,(Value(i),i=1,iRank)
c  or, if iCount=0
c      Write(iUnit,FMT) STR, (Value(i),i=1,iRank)
c
      character*(*) STR
      character*(*) STR2
      character*(*) FMT
      character*120 TMP
      character*2 Marker
      real*8 Value(iRank)
      Common /icolorize/icolorize
      if(icolorize.eq.1) then
        Marker='::'
        if(iCount.eq.0) then
          write(TMP,FMT) STR,(Value(i),i=1,iRank)
        else
          Write(TMP,FMT) STR,iCount,STR2,(Value(i),i=1,iRank)
        endif
        init=1
        if(TMP(1:3).eq.'   ') init=3
        Write(iUnit,'(a,a)') Marker,TMP(init:mylen(TMP))
      else
        if(iCount.eq.0) then
          Write(iUnit,FMT) STR,(Value(i),i=1,iRank)
        else
          Write(iUnit,FMT) STR,iCount,STR2,(Value(i),i=1,iRank)
        endif
      endif
      return
      end
