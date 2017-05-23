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
* Copyright (C) 2001, Valera Veryazov                                  *
************************************************************************
      Subroutine SysExpand(strin, strout, iRet)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     general purpose routine for expanding frequently used messages   *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     V.Veryazov University of Lund, 2001                              *
*                                                                      *
************************************************************************
*       character s*256
*       call sysexpand('MSG: open ',s,i)  - correct call
*       call sysexpand('MSG: opens ',s,i) - uncorrect, but could be fixed
*       call sysexpand('MSG: uNit ',s,i)  - mixed case
*       call sysexpand('MSG: blah ',s,i)  - will print BLAH
*             i=real len of s
*       call sysexpand('blah-blah ',s,i)  - will print string as is
*            and return i=0
************************************************************************
      parameter (MAXlabel=8)
      Character*(*) strin, strout
      Character*16 MSGlabel(MAXlabel)
      Character*128 MSGtext(MAXlabel)
      Character*512 sstrin
      Dimension MSGlen(MAXlabel)
      character*27 up,lw
      character*36  printable
      character c
      dimension itab(0:255)
      save up,lw,ifset,itab
      save MSGlabel, MSGtext, MSGlen
      data up /'ABCDEFGHIJKLMNOPQRSTUVWXYZ ' /
      data lw /'abcdefghijklmnopqrstuvwxyz ' /
      data printable /'1234567890-=~!@#$%^&*()_+<>,.?/[]":;'/
      data ifset / 0 /
c FORTRAN hash :-)
      data MSGlabel/'OPEN',
     *  'CLOSE',
     *  'UNIT',
     *  'DELETE',
     *  'SEEK',
     *  'INVALIDOPTION',
     *  'USED',
     *  'NOTOPENED'/
      data MSGtext/'Premature abort while opening file',   !OPEN
     *  'Premature abort while closing the file',          !CLOSE
     *  'Invalid unit number (Lu<=0 or Lu>99)',            !UNIT
     *  'Premature abort while removing the file',         !DELETE
     *  'Premature abort while seeking the file',          !SEEK
     * 'An invalid option or combination of options has been supplied',
     *  'Invalid unit number. The file is already opened', !USED
     *  'File is not Opened'/ !NOT OPENED

c preset of saved data
c this code uses a part of upcase routine
      if (ifset.eq.0) then
        ifset=1
        do i=0,255
          itab(i)=-1
        end do
        do ii=1,26
          i=ichar(up(ii:ii))
          j=ichar(lw(ii:ii))
          itab(j)=i
          itab(i)=i
        end do
        do i=1,MAXlabel
         do j=128,1,-1
         if(MSGtext(i)(j:j).ne.' ') goto 1
         enddo
1        MSGlen(i)=j
        enddo
      end if
c fixation of bug in G77
      sstrin=strin
c no action if it's an ordinary message
      if(sstrin(1:4).ne.'MSG:') then
      do ii=1, len(sstrin)
      c=sstrin(ii:ii)
#ifdef NAGFOR
      if(c.ne.'\') then
#endif
       ip=index(up,c)+index(lw,c)+index(printable,c)
c        print *,'ii',ii
       if(ip.eq.0) sstrin(ii:ii)=' '
#ifdef NAGFOR
       endif
#endif
      enddo
      iRet=0
      return
      endif
c uppercase with removing noncharacters.
      ik=0
      do ii=5,len(sstrin)
        i=ichar(sstrin(ii:ii))
        j=itab(i)
        if(j.ge.0) then
        ik=ik+1
        sstrin(ik:ik)=char(j)
        endif
      end do
      do ii=1,MAXlabel
      if(sstrin(1:ik).eq.MSGlabel(ii)) goto 10
      enddo
c  try again...
      do ii=1,MAXlabel
      if(sstrin(1:4).eq.MSGlabel(ii)(1:4)) goto 10
      enddo
      strout=sstrin(1:ik)
      iRet=ik
      return
10     continue
        strout=MSGtext(ii)(1:MSGlen(ii))
        iRet=MSGlen(ii)
        return
      end
