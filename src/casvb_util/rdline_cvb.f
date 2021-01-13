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
      subroutine rdline_cvb(nfield)
      implicit real*8 (a-h,o-z)
c  BLANKDELIM signifies whether blanks are used to delimit fields :
      logical blankdelim,variat
      logical debug
      logical isitanint_cvb,isitareal_cvb
      external isitanint_cvb,isitareal_cvb
#include "luinp_cvb.fh"
      parameter(nblank=2,ncomeol=3,neol=4,neofield=1,neof=2,nempty=1,
     >  nalias=2)
      character*300 line
      character*(*) string
      character*1 blanks(nblank)
c  COMEOL are comments that comment out the rest of the line
c  (might add 'bracketing' comments (e.g. /* ... */ ) later) :
      character*3 comeol(ncomeol)
      character*1 eol(neol)
      character*1 eofield(neofield)
      character*10 eof(neof)
      character*2 empty(nempty)
      character*5 alias(nalias,2)
      dimension ilv(300)
      save blankdelim
      save blanks,comeol,eof,eol,eofield,empty
      save line,lenline,ilv
      save iline,nline,nlold
      data iline/0/,nline/0/,nlold/0/
      data blankdelim/.true./
      data blanks/' ',','/,comeol/'***','!  ','*  '/
      data eof/'---       ','ENDOFINPUT'/
      data eol/';','=','{','}'/
      data eofield/' '/
      data empty/'--'/
      data alias/'={   ',';    ','}    ',';END;'/
      data debug/.false./

50    continue
      if(nline.eq.-1)then
        nfield=-1
        return
      endif
      if(iline.lt.nline)then
        iline=iline+1
        goto 899
      else
        iline=1
      endif
c  Read full input line from file and make preparations :
100   read(inp,'(a)',end=200)line
      lenline=len_trim_cvb(line)
      call strip_blanks_cvb(line,lenline,blanks,nblank,blankdelim)
      call upper_case_cvb(line,lenline)
c  Check for "end-of-file" character sequence :
      do 300 ieof=1,neof
      ilength=len_trim_cvb(eof(ieof))
      if(line(1:ilength).eq.eof(ieof)(1:ilength))goto 200
300   continue
c  Comment strings (skip rest of line ) :
      indmin=lenline+1
      do 500 icom=1,ncomeol
      ind=index(line(1:lenline),comeol(icom)
     >  (1:len_trim_cvb(comeol(icom))))
      if(ind.ne.0)indmin=min(indmin,ind)
500   continue
      lenline=len_trim_cvb(line(1:indmin-1))
      if(lenline.eq.0)goto 100
c  Aliases :
      do 550 ialias=1,nalias
560   ind=index(line(1:lenline),alias(ialias,1)
     >  (1:len_trim_cvb(alias(ialias,1))))
      if(ind.ne.0)then
        call charinsert_cvb(alias(ialias,2),
     >    len_trim_cvb(alias(ialias,2)),
     >    line,lenline,ind,len_trim_cvb(alias(ialias,1)))
        goto 560
      endif
550   continue
      lenline=len_trim_cvb(line(1:indmin-1))
      if(lenline.eq.0)goto 100
c  Split into lines :
      call izero(ilv,lenline)
      do 600 ieol=1,neol
      if=0
      ilength=len_trim_cvb(eol(ieol))
650   ind=index(line(if+1:lenline),eol(ieol)(1:ilength))
      if(ind.ne.0)then
        ilv(ind+if)=1
        if=ind+if
        goto 650
      endif
600   continue
      nlold=nline
      nline=1
      do 700 ich=1,lenline
      if(ilv(ich).eq.1)nline=nline+1
700   continue
c  Split into fields :
      do 800 ieofield=1,neofield
      if=0
      ilength=max(1,len_trim_cvb(eofield(ieofield)))
850   ind=index(line(if+1:lenline),eofield(ieofield)(1:ilength))
      if(ind.ne.0)then
        ilv(ind+if)=2
        if=ind+if
        goto 850
      endif
800   continue
c  Eliminate field separators at end of lines:
      do 873 ieofield=1,neofield
      ihadchar=0
      do 875 ich=lenline,1,-1
      if(ilv(ich).eq.1)then
        ihadchar=0
      elseif(ilv(ich).eq.2)then
        if(ihadchar.eq.0)ilv(ich)=0
        ilength=len_trim_cvb(eofield(ieofield))
        do 880 i=0,ilength-1
        line(i+ich:i+ich)=' '
880     continue
        ihadchar=1
      else
        ihadchar=1
      endif
875   continue
873   continue
c  Count the number of fields in present card (ILINE).
c  Also make sure line is not all empty :
899   continue
      nfield=1
      jline=1
      ilinebeg=0
      ilineend=-1
      do 900 ich=1,lenline
      if(jline.eq.iline-1)ilinebeg=ich+1
      if(ilv(ich).eq.1)jline=jline+1
      if(jline.eq.iline+1.and.ilineend.eq.-1)ilineend=ich-1
      if(jline.eq.iline.and.ilv(ich).eq.2)nfield=nfield+1
900   continue
      if(iline.eq.1)ilinebeg=1
      if(ilineend.eq.-1)ilineend=lenline
c  Go back and read a new line if this one is empty :
      if(ilinebeg.gt.ilineend)goto 50
      if(len_trim_cvb(line(ilinebeg:ilineend)).eq.0)goto 50
      return
200   nline=-1
      nfield=-1
      return
      entry pushline_cvb()
      if(iline.eq.1.or.nline.eq.-1)then
        backspace(inp)
        iline=nlold
        nline=nlold
      else
        iline=iline-1
      endif
      return
      entry gtany_cvb(string,int,real,ic,ifield,ierr)
c      ic=1
c      goto 1050
c      entry gtint_cvb(int,ifield,ierr)
c      ic=2
c      goto 1050
c      entry gtreal_cvb(real,ifield,ierr)
c      ic=3
c1050  continue
      if(ic.gt.1)ierr=0
      jfield=1
      jline=1
      do 1000 ich=1,lenline
      if(ilv(ich).eq.1)jline=jline+1
      if(jline.eq.iline.and.ilv(ich).eq.2)jfield=jfield+1
      if(jline.eq.iline.and.jfield.eq.ifield)then
        jch=ich
        if(ich.eq.1)jch=0
1100    jch=jch+1
        if(ilv(jch).eq.0.and.jch.ne.lenline+1)goto 1100
        ifirst=ich+1
        if(ich.eq.1)ifirst=1
c  Special character strings to signify empty field ?
        do 1125 iempty=1,nempty
        if(line(ifirst:jch-1).eq.empty(iempty)
     >    (1:len_trim_cvb(empty(iempty))))goto 1130
1125    continue
        if(debug)write(6,*)' Field=',line(ifirst:jch-1)
        if(ic.eq.1)then
          string=line(ifirst:jch-1)
          if(debug)write(6,*)' Field read as string :',string
        elseif(ic.eq.2)then
          if(ifirst.gt.jch-1)then
            ierr=2
            return
          endif
          if(.not.isitanint_cvb(line(ifirst:jch-1)))goto 1150
          read(line(ifirst:jch-1),*,err=1150)int
          if(debug)write(6,*)' Field read as int :',int
        elseif(ic.eq.3)then
          if(ifirst.gt.jch-1)then
            ierr=2
            return
          endif
          if(.not.isitareal_cvb(line(ifirst:jch-1)))goto 1150
          read(line(ifirst:jch-1),*,err=1150)real
          if(debug)write(6,*)' Field read as real :',real
        endif
        return
1130    continue
c  "Empty" field :
        if(ic.eq.1)then
          string=' '
        else
          ierr=2
        endif
        return
      endif
1000  continue
      write(6,*)' Error in input parsing !'
      call abend_cvb()
1150  ierr=1
      if(debug)then
        if(ic.eq.2)then
          write(6,*)' Could not read field as integer.'
        else
          write(6,*)' Could not read field as real.'
        endif
      endif
      return
      entry rdline_init_cvb(variat)
      if(variat)return
      rewind(inp)
2100  read(inp,'(a)',end=2200)line
      lenline=len_trim_cvb(line)
      call strip_blanks_cvb(line,lenline,blanks,nblank,blankdelim)
      call upper_case_cvb(line,lenline)
      if(line(1:6).ne.'&CASVB')goto 2100
C     if(.not.(((.not.blankdelim).and.line(1:10).eq.'&CASVB&END').or.
C    >   (blankdelim.and.line(1:11).eq.'&CASVB &END')))goto 2100
      return
2200  write(6,*)' WARNING: Initiation string not found in input file.'
      return
      end



      logical function isitanint_cvb(a)
      implicit real*8 (a-h,o-z)
      parameter(nallowed=12)
      character*(*) a
      character*(1) allowedchars(nallowed)
      data allowedchars/'+','-','0','1','2','3','4','5','6','7','8','9'/

      do 100 ich=1,len_trim_cvb(a)
      do 200 j=1,nallowed
      if(a(ich:ich).eq.allowedchars(j))goto 100
200   continue
      isitanint_cvb=.false.
      return
100   continue
      isitanint_cvb=.true.
      return
      end
      logical function isitareal_cvb(a)
      implicit real*8 (a-h,o-z)
      parameter(nallowed=17)
      character*(*) a
      character*(1) allowedchars(nallowed)
      data allowedchars/'+','-','0','1','2','3','4','5','6','7','8','9',
     >  '.','E','e','D','d'/

      do 100 ich=1,len_trim_cvb(a)
      do 200 j=1,nallowed
      if(a(ich:ich).eq.allowedchars(j))goto 100
200   continue
      isitareal_cvb=.false.
      return
100   continue
      isitareal_cvb=.true.
      return
      end
