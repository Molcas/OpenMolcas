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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine gsinp_cvb(
     >  orbs,irdorbs,ip_cvb,nvbinp,kbasiscvb_inp,
     >  mxaobf,mxorb,kbasis,strtvb)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      parameter (ngs=7,ncmp=4)
      character*8 guess
      logical firsttime_cvb
      external firsttime_cvb
      dimension guess(ngs)
      dimension orbs(mxaobf,mxorb),irdorbs(mxorb)
      dimension idum(1)
      save guess
      data guess/ 'ORB     ','STRUC   ','READ    ','AOBASIS ',
     >            'MOBASIS ','END     ','ENDGUESS'/

      if(firsttime_cvb())call touch_cvb('INPGS')
      mouse=1
3000  call fstring_cvb(guess,ngs,istr,ncmp,2)
      if(istr.eq.1)then
c 'ORB'
        call int_cvb(idum,1,nread,0)
        iorb=idum(1)
        if(iorb.le.0.or.iorb.gt.mxorb)then
          write(6,*)' Illegal orbital number read :',iorb
          call abend_cvb()
        endif
        if(nread.eq.0)then
          write(6,*)' Orbital label in ORB keyword not found!'
          call abend_cvb()
        endif
        irdorbs(iorb)=mouse
        call fzero(orbs(1,iorb),mxaobf)
        call real_cvb(orbs(1,iorb),mxaobf,nread,0)
      elseif(istr.eq.2)then
c 'STRUC'
c  If previous orb. permutation disable:
        call mhpfreer_cvb(ip_cvb)
        mxread=mavailr_cvb()/2
        ip_cvb=mheapr_cvb(mxread)
        call realz_cvb(w(ip_cvb),mxread,nvbinp,0)
        call mreallocr_cvb(ip_cvb,nvbinp)
        kbasiscvb_inp=kbasis
      elseif(istr.eq.3)then
cc 'READ'
c        call fstring_cvb(readgs,nrdgs,istr2,ncmp,1)
c        if(istr2.eq.1)then
cc 'ORB'
c          iorb1=0
c          iorb2=0
c          call int_cvb(idum,1,nread,0)
c          iorb1=idum(1)
c          if(nread.eq.0)then
c            write(6,*)' No orbital number in READ,ORB keyword!'
c            call abend_cvb()
c          endif
c          call fstring_cvb('TO      ',1,jstr,ncmp,1)
c          if(jstr.ne.0)then
c            call int_cvb(idum,1,nread,0)
c            iorb2=idum(1)
c            if(nread.eq.0)then
c              write(6,*)' No orbital number after READ,...,TO, !'
c              call abend_cvb()
c            endif
c          else
c            iorb2=iorb1
c          endif
c          jorb1=iorb1
c          jorb2=iorb2
c          call setstrtvb_cvb(strtvb)
c          recordnm=strtvb
c3100      call fstring_cvb(readgs2,nrdgs2,istr3,ncmp,1)
c          if(istr3.eq.1)then
c            call real_cvb(recordnm,1,nread,0)
c            if(nread.eq.0)then
c              write(6,*)' No identifier after READ,...,FROM, !'
c              call abend_cvb()
c            endif
c          elseif(istr3.eq.2)then
c            call int_cvb(idum,1,nread,0)
c            jorb1=idum(1)
c            if(nread.eq.0)then
c              write(6,*)' No orbital number after READ,...,AS, !'
c              call abend_cvb()
c            endif
c            call fstring_cvb('TO      ',1,jstr,ncmp,1)
c            if(jstr.ne.0)then
c              call int_cvb(idum,1,nread,0)
c              jorb2=idum(1)
c              if(nread.eq.0)then
c                write(6,*)' No orbital number after READ,...,TO, !'
c                call abend_cvb()
c              endif
c            else
c              jorb2=jorb1
c            endif
c          endif
c          if(istr3.ne.0)goto 3100
c          call othergs_cvb(orbs,w(ip_cvb),recordnm,1,
c     >      iorb1,iorb2,jorb1,jorb2)
c        elseif(istr2.eq.2)then
cc      'STRUC'
c          istruc1=0
c          istruc2=0
c          call int_cvb(idum,1,nread,0)
c          istruc1=idum(1)
c          if(nread.eq.0)then
c            write(6,*)' No structure number in READ,STRUC keyword!'
c            call abend_cvb()
c          endif
c          call fstring_cvb('TO      ',1,jstr,ncmp,1)
c          if(jstr.ne.0)then
c            call int_cvb(idum,1,nread,0)
c            istruc2=idum(1)
c            if(nread.eq.0)then
c              write(6,*)' No structure number after READ,...,TO, !'
c              call abend_cvb()
c            endif
c          else
c            istruc2=istruc1
c          endif
c          jstruc1=istruc1
c          jstruc2=istruc2
c          call setstrtvb_cvb(strtvb)
c          recordnm=strtvb
c3200      call fstring_cvb(readgs2,nrdgs2,istr3,ncmp,1)
c          if(istr3.eq.1)then
c            call real_cvb(recordnm,1,nread,0)
c            if(nread.eq.0)then
c              write(6,*)' No identifier after READ,...,FROM, !'
c              call abend_cvb()
c            endif
c          elseif(istr3.eq.2)then
c            call int_cvb(idum,1,nread,0)
c            jstruc1=idum(1)
c            if(nread.eq.0)then
c              write(6,*)' No structure number after READ,...,AS, !'
c              call abend_cvb()
c            endif
c            call fstring_cvb('TO      ',1,jstr,ncmp,1)
c            if(jstr.ne.0)then
c              call int_cvb(idum,1,nread,0)
c              jstruc2=idum(1)
c              if(nread.eq.0)then
c                write(6,*)' No structure number after READ,...,TO, !'
c                call abend_cvb()
c              endif
c            else
c              jstruc2=jstruc1
c            endif
c          endif
c          if(istr3.ne.0)goto 3200
c          call othergs_cvb(orbs,w(ip_cvb),recordnm,2,
c     >      istruc1,istruc2,jstruc1,jstruc2)
c        elseif(istr2.eq.3)then
cc 'ALL'
c          call setstrtvb_cvb(strtvb)
c          recordnm=strtvb
c          call fstring_cvb('FROM    ',1,jstr,ncmp,1)
c          if(jstr.ne.0)then
c            call real_cvb(recordnm,1,nread,0)
c            if(nread.eq.0)then
c              write(6,*)' No identifier after READ,...,FROM, !'
c              call abend_cvb()
c            endif
c          endif
c          call getguess_cvb(orbs,w(ip_cvb),recordnm,kbasiscvb_inp)
c        endif
      elseif(istr.eq.4)then
c 'AOBASIS'
        mouse=2
      elseif(istr.eq.5)then
c 'MOBASIS'
        mouse=1
      endif
c 'END' , 'ENDGUESS' or unrecognized keyword -- end GUESS input :
      if(.not.(istr.eq.0.or.istr.eq.6.or.istr.eq.7))goto 3000

      return
c Avoid unused argument warnings
      if (.false.) call Unused_real(strtvb)
      end
