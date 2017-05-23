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
      subroutine setifinish_cvb(icode)
      implicit real*8(a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


      if(icode.eq.1)then
        ifinish=0
      elseif(icode.eq.3)then
        ifinish=1
      endif
      return
      end
      logical function firsttime_cvb()
      implicit real*8(a-h,o-z)
      logical begbracket,second_time_round
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "inpmod_cvb.fh"
#include "seth_cvb.fh"
#include "loopcntr_cvb.fh"
#include "initopt_cvb.fh"
      external istkprobe_cvb
      logical istkprobe_cvb

      if(inputmode.eq.2)then
        if(loopstep.eq.0) then
          begbracket=.false.
        else
          begbracket=
     >      (icode(loopstep).eq.1.and.icode(loopstep+1).eq.2).or.
     >      (icode(loopstep).eq.3.and.icode(loopstep+1).eq.4)
        endif

        firsttime_cvb=(joptstep.eq.ioptstep-1.or.
     >    (ioptstep.eq.0.and.joptstep.eq.0)) .or.
     >    (joptstep.eq.ioptstep.and.begbracket)

        if(ioptim.gt.1.and.iopt2step(ioptim).eq.iopt2step(ioptim-1))
     >    firsttime_cvb=.false.

        if(istkprobe_cvb(istackrep))then
          call istkpop_cvb(istackrep,nc_zeroed_l)
          call istkpop_cvb(istackrep,nconvinone_l)
          call istkpop_cvb(istackrep,italter_l)
          second_time_round=(italter_l.gt.1)
          call istkpush_cvb(istackrep,italter_l)
          call istkpush_cvb(istackrep,nconvinone_l)
          call istkpush_cvb(istackrep,nc_zeroed_l)
        else
          second_time_round=.false.
        endif

        if(second_time_round)firsttime_cvb=.false.

        if(nmcscf.ge.2)firsttime_cvb=.false.
      else
        firsttime_cvb=.false.
      endif
      return
      end
