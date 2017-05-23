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
************************************************************************
*                                                                      *
*     purpose:                                                         *
*       general purpose routine for printing nice messages             *
*       with simple reformatting:                                      *
*           \n and long string converted to new line                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     V.Veryazov University of Lund, 2001                              *
*                                                                      *
************************************************************************
c       call SysPutsStart()
c       call SysPuts('MOL;CAS\n\nmolcas is a quantum '//
c     6  'chemistry software developed by scientists to be '//
c     6   'used by scientists. It is not primarily a commercial '//
c     6   'product and it is not sold in order to produce a '//
c     6   'fortune for its owner (the Lund University).',' ',' ')
c       call SysPutsEnd()
c       end
*
      Subroutine SysPuts(str,str1,str2)
      Character*(*) str,str1,str2
      Character *512 Junk
c because of bug in g77 we can't just concatenate strings and
c  had to have limited length of the string
      iTooLong=60
      iLongEnough=50
      call mycat(Junk,str,str1,str2)
      mlen=mylen(Junk)
      mleni=mlen
      ipos=1
c  check is '\n' .eq. <CR>?
      icr=len('\n')-1
      icr1=0
100   j=100000
      ii=index(Junk(ipos:mleni),'\n')
      if(ii.gt.0) j=min(j,ii)
      iii=index(Junk(ipos:mleni),';')
      if(iii.gt.0) j=min(j,iii)
      if(j.eq.ii.and.ii.gt.0) icr1=icr
      if(j.eq.iii.and.iii.gt.0) icr1=0
      if(j.eq.100000) j=0
      i=j
      if(i.gt.iTooLong.or.(i.eq.0.and.mlen.gt.iTooLong)) then
         ij=index(Junk(ipos+iLongEnough:mleni),' ')
         if(ij.eq.0)then
            call SysDumpStr(Junk(ipos:mleni))
            if(i.eq.0) return
         else
            ij=ij+iLongEnough-1
            call SysDumpStr(Junk(ipos:ipos+ij))
         endif
         ipos=ipos+ij+1
         mlen=mlen-ij-1
         goto 100
      endif
      if(i.eq.0) then
         call SysDumpStr(Junk(ipos:mleni))
         return
      endif
      call SysDumpStr(Junk(ipos:ipos+i-2))
      ipos=ipos+i+icr1
      mlen=mlen-i-icr1
      if(mlen.gt.0) goto 100
      return
      end
************************************************************************
      subroutine  SysDumpStr(str)
      character*(*) str
      character fmt*20
      iTooLong=60
      i=len(str)
      if(i.gt.iTooLong+7) then
c oops! too long
      write (6,'(a,a)')   ' ###    ',str
      return
      endif
      i=iTooLong+8-i
      write(fmt,'(a, i2,a)') '(a,a,',i,'x,a)'
      write (6,fmt) ' ###    ',str,' ###'
      return
      end

************************************************************************
      subroutine mycat(Junk,str0,str1,str2)
c
c  Junk=str0//str1//str2
c
      Character*(*) Junk,str0,str1,str2
      maxlen=len(Junk)

      Junk=' '
      ile=1
      il=mylen(str0)
      if (il.gt.0) then
         ils=1
         ile=il+1
         if(ile.gt.maxlen) goto 100
         Junk(ils:ile)=str0(1:il)
      end if
      il=mylen(str1)
      if (il.gt.0) then
         ils=ile+1
         ile=ile+il
         if(ile.gt.maxlen) goto 100
        Junk(ils:ile)=str1(1:il)
      end if
      il=mylen(str2)
      if (il.gt.0) then
         ils=ile+1
         ile=ile+il
         if(ile.gt.maxlen) goto 100
        Junk(ils:ile)=str2(1:il)
      end if
      return
100     write(6,*) ' too long strings to concatenate: '
        write(6,*) str0,str1,str2
        return
      end
************************************************************************
      function mylen(s)
c
c  return real length of the string without spaces...
c
      Character*(*) s
      il=len(s)
      if(il.eq.0) then
      mylen=0
      return
      endif
      do i=il,1,-1
       if(s(i:i).ne.' ') then
        mylen=i
       return
       endif
      enddo
      mylen=0
      return
      end
************************************************************************
