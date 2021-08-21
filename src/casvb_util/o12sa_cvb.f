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
      subroutine o12sa_cvb(nparm1)
      implicit real*8 (a-h,o-z)
c ... Content of CI vectors ...
      logical, external :: tstcnt_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "opt2_cvb.fh"
#include "WrkSpc.fh"

      call ddnewopt_cvb()
      have_solved_it=.false.

c  Find CIVBS :
      ivuse=0
      do 10 iv=1,nv
      if(tstcnt_cvb(work(lc(iv)),4))ivuse=iv
10    continue
      ivuse2=3
      if(ivuse.eq.3)ivuse2=2
      if(ivuse2.gt.nv)ivuse2=1
      if(ivuse.ne.0)then
        i1 = mstackr_cvb(nparm1)
        i2 = mstackr_cvb(nparm1)
        i3 = mstackr_cvb(nvb+nprorb)
        call o12sa2_cvb(work(i1),work(i2),nparm1,
     >    work(lc(ivuse2)),work(lc(ivuse)),
     >    work(lw(9)),work(lv(2)),work(i3))
        call mfreer_cvb(i1)
      else
        if(strucopt)then
          call ddguess_cvb(work(lv(2)),nvb,nprorb)
        else
          call ddguess_cvb([one],1,0)
        endif
      endif

      i1 = mstackr_cvb(nparm1)
      call o12sa3_cvb(work(ix(1)),work(lv(2)),
     >  work(lv(1)),work(lw(4)),work(lw(5)),work(lw(6)),
     >  work(lc(1)),work(lc(2)),work(lc(3)),work(lw(9)),work(i1),
     >  nvb,nprorb,nparm1,strucopt)
      call mfreer_cvb(i1)
      return
      end
