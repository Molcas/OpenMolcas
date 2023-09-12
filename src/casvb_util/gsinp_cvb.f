!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************
      subroutine gsinp_cvb(                                             &
     &  orbs,irdorbs,ip_cvb,nvbinp,kbasiscvb_inp,                       &
     &  mxaobf,mxorb,kbasis,strtvb)
      implicit real*8 (a-h,o-z)
#include "WrkSpc.fh"
      parameter (ngs=7,ncmp=4)
      character*8 guess
      logical firsttime_cvb
      external firsttime_cvb
      dimension guess(ngs)
      dimension orbs(mxaobf,mxorb),irdorbs(mxorb)
      dimension idum(1)
      save guess
      data guess/ 'ORB     ','STRUC   ','READ    ','AOBASIS ',          &
     &            'MOBASIS ','END     ','ENDGUESS'/

      if(firsttime_cvb())call touch_cvb('INPGS')
      mouse=1
3000  call fstring_cvb(guess,ngs,istr,ncmp,2)
      if(istr.eq.1)then
! 'ORB'
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
! 'STRUC'
!  If previous orb. permutation disable:
        call mhpfreer_cvb(ip_cvb)
        mxread=mavailr_cvb()/2
        ip_cvb=mheapr_cvb(mxread)
        call realz_cvb(work(ip_cvb),mxread,nvbinp,0)
        call mreallocr_cvb(ip_cvb,nvbinp)
        kbasiscvb_inp=kbasis
      elseif(istr.eq.3)then
!c 'READ'
!        call fstring_cvb(readgs,nrdgs,istr2,ncmp,1)
!        if(istr2.eq.1)then
!c 'ORB'
!          iorb1=0
!          iorb2=0
!          call int_cvb(idum,1,nread,0)
!          iorb1=idum(1)
!          if(nread.eq.0)then
!            write(6,*)' No orbital number in READ,ORB keyword!'
!            call abend_cvb()
!          endif
!          call fstring_cvb('TO      ',1,jstr,ncmp,1)
!          if(jstr.ne.0)then
!            call int_cvb(idum,1,nread,0)
!            iorb2=idum(1)
!            if(nread.eq.0)then
!              write(6,*)' No orbital number after READ,...,TO, !'
!              call abend_cvb()
!            endif
!          else
!            iorb2=iorb1
!          endif
!          jorb1=iorb1
!          jorb2=iorb2
!          call setstrtvb_cvb(strtvb)
!          recordnm=strtvb
!3100      call fstring_cvb(readgs2,nrdgs2,istr3,ncmp,1)
!          if(istr3.eq.1)then
!            call real_cvb(recordnm,1,nread,0)
!            if(nread.eq.0)then
!              write(6,*)' No identifier after READ,...,FROM, !'
!              call abend_cvb()
!            endif
!          elseif(istr3.eq.2)then
!            call int_cvb(idum,1,nread,0)
!            jorb1=idum(1)
!            if(nread.eq.0)then
!              write(6,*)' No orbital number after READ,...,AS, !'
!              call abend_cvb()
!            endif
!            call fstring_cvb('TO      ',1,jstr,ncmp,1)
!            if(jstr.ne.0)then
!              call int_cvb(idum,1,nread,0)
!              jorb2=idum(1)
!              if(nread.eq.0)then
!                write(6,*)' No orbital number after READ,...,TO, !'
!                call abend_cvb()
!              endif
!            else
!              jorb2=jorb1
!            endif
!          endif
!          if(istr3.ne.0)goto 3100
!          call othergs_cvb(orbs,work(ip_cvb),recordnm,1,
!     >      iorb1,iorb2,jorb1,jorb2)
!        elseif(istr2.eq.2)then
!c      'STRUC'
!          istruc1=0
!          istruc2=0
!          call int_cvb(idum,1,nread,0)
!          istruc1=idum(1)
!          if(nread.eq.0)then
!            write(6,*)' No structure number in READ,STRUC keyword!'
!            call abend_cvb()
!          endif
!          call fstring_cvb('TO      ',1,jstr,ncmp,1)
!          if(jstr.ne.0)then
!            call int_cvb(idum,1,nread,0)
!            istruc2=idum(1)
!            if(nread.eq.0)then
!              write(6,*)' No structure number after READ,...,TO, !'
!              call abend_cvb()
!            endif
!          else
!            istruc2=istruc1
!          endif
!          jstruc1=istruc1
!          jstruc2=istruc2
!          call setstrtvb_cvb(strtvb)
!          recordnm=strtvb
!3200      call fstring_cvb(readgs2,nrdgs2,istr3,ncmp,1)
!          if(istr3.eq.1)then
!            call real_cvb(recordnm,1,nread,0)
!            if(nread.eq.0)then
!              write(6,*)' No identifier after READ,...,FROM, !'
!              call abend_cvb()
!            endif
!          elseif(istr3.eq.2)then
!            call int_cvb(idum,1,nread,0)
!            jstruc1=idum(1)
!            if(nread.eq.0)then
!              write(6,*)' No structure number after READ,...,AS, !'
!              call abend_cvb()
!            endif
!            call fstring_cvb('TO      ',1,jstr,ncmp,1)
!            if(jstr.ne.0)then
!              call int_cvb(idum,1,nread,0)
!              jstruc2=idum(1)
!              if(nread.eq.0)then
!                write(6,*)' No structure number after READ,...,TO, !'
!                call abend_cvb()
!              endif
!            else
!              jstruc2=jstruc1
!            endif
!          endif
!          if(istr3.ne.0)goto 3200
!          call othergs_cvb(orbs,work(ip_cvb),recordnm,2,
!     >      istruc1,istruc2,jstruc1,jstruc2)
!        elseif(istr2.eq.3)then
!c 'ALL'
!          call setstrtvb_cvb(strtvb)
!          recordnm=strtvb
!          call fstring_cvb('FROM    ',1,jstr,ncmp,1)
!          if(jstr.ne.0)then
!            call real_cvb(recordnm,1,nread,0)
!            if(nread.eq.0)then
!              write(6,*)' No identifier after READ,...,FROM, !'
!              call abend_cvb()
!            endif
!          endif
!          call getguess_cvb(orbs,work(ip_cvb),recordnm,kbasiscvb_inp)
!        endif
      elseif(istr.eq.4)then
! 'AOBASIS'
        mouse=2
      elseif(istr.eq.5)then
! 'MOBASIS'
        mouse=1
      endif
! 'END' , 'ENDGUESS' or unrecognized keyword -- end GUESS input :
      if(.not.(istr.eq.0.or.istr.eq.6.or.istr.eq.7))goto 3000

      return
! Avoid unused argument warnings
      if (.false.) call Unused_real(strtvb)
      end
