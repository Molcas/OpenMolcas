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
      subroutine loopcntr_init_cvb(inputmode1,initfalse)
      implicit real*8(a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "inpmod_cvb.fh"
#include "loopcntr_cvb.fh"
#include "seth_cvb.fh"
#include "initopt_cvb.fh"
      logical initfalse
      logical guess_available,initial_opts,svbfirst,constrained_opt

      call istkinit_cvb(istackrep,nstackrep)
      inputmode=inputmode1
      ioptim=0
      ioptstep=0
      if(inputmode.eq.2)then
        loopstepmx=loopstep
        noptstep=joptstep

c  Check for "special case" => initial optimizations
        initial_opts=.true.
        guess_available=.false.
        if(nmcscf.ge.2)guess_available=.true.
        if(up2date_cvb('WRITEGS'))guess_available=.true.
        if(.not.up2date_cvb('STRTGS'))guess_available=.true.
        if(.not.up2date_cvb('INPGS'))guess_available=.true.

        if(guess_available)initial_opts=.false.
        if(noptstep.gt.0)initial_opts=.false.
c  Are we doing a constrained optimization ? :
        constrained_opt=.false.
        if(norbrel.gt.0)constrained_opt=.true.
        if(nort.gt.0)constrained_opt=.true.
        if(ndrot.gt.0)constrained_opt=.true.
        if(nfxorb.gt.0)constrained_opt=.true.
        if(ploc)constrained_opt=.true.
        if(nfxvb.gt.0.or.lfxvb.eq.1)constrained_opt=.true.
        if(nzrvb.gt.0.or.lzrvb.eq.1)constrained_opt=.true.
c  If INIT/NOINIT keyword was encountered, override :
        if(initial.eq.0)initial_opts=.false.
        if(initial.eq.1)initial_opts=.true.
c  Finally may be overridden by initfalse :
        if(initfalse)initial_opts=.false.
        if(initial_opts)then
c  IOPTCODE    +1  = REPORT
c              +2  = OPTIM
c              +4  = Svb
c              +8  = freeze structure coefficients
c              +16 = strong-orthogonality constraints

c  Should we do Svb optimization first ? :
          svbfirst=ifcasci_cvb()

          if(.not.constrained_opt)then
c  Two first optimizations are SOPP & PP :
            noptim=0
            if(svbfirst)then
              if(norb.gt.2)then
                noptim=noptim+1
                ioptcode(noptim)=22
              endif
              if(strucopt)then
                noptim=noptim+1
                ioptcode(noptim)=14
                if(noptim.eq.2)ioptcode(1)=ioptcode(1)+8
              endif
            else
              if(norb.gt.2)then
                noptim=noptim+1
                ioptcode(noptim)=18
              endif
              if(strucopt)then
                noptim=noptim+1
                ioptcode(noptim)=10
                if(noptim.eq.2)ioptcode(1)=ioptcode(1)+8
              endif
            endif
          endif
c  Then a third Svb optimization if we are doing Evb :
          if(icrit.ne.1.and.svbfirst)then
            noptim=noptim+1
            ioptcode(noptim)=6
          endif
c  Finally actual optimization :
          noptim=noptim+1
          ioptcode(noptim)=2
c  Add "report" :
          noptim=noptim+1
          ioptcode(noptim)=1

          iopt2step(0)=0
          do 100 i=1,noptim
          iopt2step(i)=1
100       continue
          iopt2step(noptim+1)=noptstep+1
        else
          noptim=noptstep
          call izero(ioptcode,noptim)
          do 200 i=0,noptim
          iopt2step(i)=i
200       continue
c  Append OPTIM keyword if none present
          noptkw=0
          do 300 lll=1,loopstepmx
          if(icode(lll).eq.1)noptkw=noptkw+1
300       continue
          if(noptkw.eq.0)then
            noptim=noptim+1
            ioptcode(noptim)=2
            iopt2step(noptim)=iopt2step(noptim-1)
          endif
c  Append REPORT keyword if none present
          nrepkw=0
          do 400 lll=1,loopstepmx
          if(icode(lll).eq.3)nrepkw=nrepkw+1
400       continue
          if(nrepkw.eq.0)then
            noptim=noptim+1
            ioptcode(noptim)=1
            iopt2step(noptim)=iopt2step(noptim-1)
          endif
          iopt2step(noptim+1)=noptstep+1
        endif
      endif
      return
      end
      logical function loopcntr_iterate_cvb()
      implicit real*8(a-h,o-z)
#include "seth_cvb.fh"
#include "loopcntr_cvb.fh"
#include "initopt_cvb.fh"
      logical unmatched
      external istkprobe_cvb
      logical istkprobe_cvb

      if(iopt2step(ioptim+1).eq.iopt2step(ioptim))then
        ioptim=ioptim+1
        goto 1100
      endif

      joptstep=0
      do 1 ll=1,loopstepmx
      if(joptstep.eq.ioptstep)goto 11
      if(icode(ll).eq.1.or.icode(ll).eq.3)joptstep=joptstep+1
1     continue
11    ll1=ll
10    continue
c  First determine if end of multi-step optimization may have been reached:
      if(istkprobe_cvb(istackrep))then
        call istkpop_cvb(istackrep,nc_zeroed)
        call istkpop_cvb(istackrep,nconvinone)
        call istkpop_cvb(istackrep,italter)
        call istkpop_cvb(istackrep,mxalter)
        call istkpop_cvb(istackrep,kk2)
        call istkpop_cvb(istackrep,ioptstep2)
        call istkpop_cvb(istackrep,ioptstep1)
c  Number of steps is IOPTSTEP2-IOPTSTEP1+1
        if(nconvinone.eq.ioptstep2-ioptstep1+1)then
          ioptstep=ioptstep2
          joptstep=ioptstep2
          ll1=kk2+1
          goto 10
        endif
c  Restore loop information :
        call istkpush_cvb(istackrep,ioptstep1)
        call istkpush_cvb(istackrep,ioptstep2)
        call istkpush_cvb(istackrep,kk2)
        call istkpush_cvb(istackrep,mxalter)
        call istkpush_cvb(istackrep,italter)
        call istkpush_cvb(istackrep,nconvinone)
        call istkpush_cvb(istackrep,nc_zeroed)
      endif
      do 100 ll=ll1,loopstepmx
      if(joptstep.eq.ioptstep)then
c  Looking for next card after previous OPTIM/REPORT:
        if(icode(ll).eq.2.or.icode(ll).eq.4)then
          if(ll.eq.1)then
            unmatched=.true.
          else
            unmatched=(icode(ll)-icode(ll-1).ne.1)
          endif
          if(unmatched)then
            write(6,'(a)')' Unmatched END or closing bracket!'
            call abend_cvb()
          endif
        endif
        if(icode(ll).eq.1.or.icode(ll).eq.3)then
c  OPTIM / REPORT :
          ioptstep=ioptstep+1
          goto 1000
        elseif(icode(ll).eq.5)then
c  ALTERN :
c  Scan rest of file for number of steps in this loop:
          iend=0
          ioptstep2=joptstep
          do 200 kk=ll+1,loopstepmx
          if(icode(kk).eq.1.or.icode(kk).eq.3)then
            ioptstep2=ioptstep2+1
          elseif(icode(kk).eq.5)then
            iend=iend-1
          elseif(icode(kk).eq.6)then
            iend=iend+1
          endif
          if(iend.eq.1)goto 300
200       continue
          write(6,*)' Run-away ENDALTERN or closing bracket!'
          call abend_cvb()
300       continue

          italter=1
          mxalter=ipos(ll)
          nconvinone=-1
          nc_zeroed=0
          call istkpush_cvb(istackrep,ioptstep+1)
          call istkpush_cvb(istackrep,ioptstep2)
          call istkpush_cvb(istackrep,kk)
          call istkpush_cvb(istackrep,mxalter)
          call istkpush_cvb(istackrep,italter)
          call istkpush_cvb(istackrep,nconvinone)
          call istkpush_cvb(istackrep,nc_zeroed)
        elseif(icode(ll).eq.6)then
c  ENDALTERN :
          call istkpop_cvb(istackrep,nc_zeroed)
          call istkpop_cvb(istackrep,nconvinone)
          call istkpop_cvb(istackrep,italter)
          call istkpop_cvb(istackrep,mxalter)
          call istkpop_cvb(istackrep,kk2)
          call istkpop_cvb(istackrep,ioptstep2)
          call istkpop_cvb(istackrep,ioptstep1)
          italter=italter+1
          nstep=ioptstep-ioptstep1+1
          if(nstep.gt.0.and.italter.le.mxalter)then
c  Next loop iteration :
            call istkpush_cvb(istackrep,ioptstep1)
            call istkpush_cvb(istackrep,ioptstep2)
            call istkpush_cvb(istackrep,kk2)
            call istkpush_cvb(istackrep,mxalter)
            call istkpush_cvb(istackrep,italter)
            call istkpush_cvb(istackrep,nconvinone)
            call istkpush_cvb(istackrep,nc_zeroed)
            ioptstep=ioptstep1
            goto 1000
          else
            if(nstep.gt.0)then
              write(6,'(/,a,i4,a)')' Exiting',nstep,
     >          '-step optimization.'
              write(6,'(a,i4)')
     >          ' Maximum number of loop iterations reached :',mxalter
            endif
          endif
        endif
      endif
      if(icode(ll).eq.1.or.icode(ll).eq.3)joptstep=joptstep+1
100   continue
      ioptstep=noptstep+1

1000  continue
      do i=1,noptim
      if(iopt2step(i).eq.ioptstep)goto 1099
      enddo
1099  ioptim=i
1100  loopcntr_iterate_cvb=(ioptim.le.noptim)
      return
      end
