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
      subroutine loopcntr_cvb(icode1)
      implicit real*8(a-h,o-z)
#include "inpmod_cvb.fh"
#include "seth_cvb.fh"
#include "loopcntr_cvb.fh"
      logical begbracket

      loopstep=loopstep+1
      if(inputmode.eq.2.and.(icode1.eq.5.or.icode1.eq.6))return
      if(loopstep.gt.mxstep)then
        write(6,*)' Loop structure too complicated -- mxstep :',mxstep
        call abend_cvb()
      endif
      if(inputmode.eq.1)then
        icode(loopstep)=icode1
        ipos(loopstep)=ncnt
      endif
      if(icode(loopstep).eq.1.or.icode(loopstep).eq.3)
     >  joptstep=joptstep+1
      if(inputmode.eq.2)then
        if(joptstep.eq.ioptstep)call setifinish_cvb(icode(loopstep))
        begbracket=
     >    (icode(loopstep).eq.1.and.icode(loopstep+1).eq.2).or.
     >    (icode(loopstep).eq.3.and.icode(loopstep+1).eq.4)
        if(joptstep.ge.ioptstep+1.or.(joptstep.eq.ioptstep.
     >    and..not.begbracket))then
c  Goto end of input :
          icnt=ncnt
          loopstep=loopstepmx
          joptstep=noptstep
        elseif(joptstep.lt.ioptstep.and.begbracket)then
c  Goto next closing bracket :
          icnt=ipos(loopstep+1)
          loopstep=loopstep+1
        endif
      endif
      return
      end
