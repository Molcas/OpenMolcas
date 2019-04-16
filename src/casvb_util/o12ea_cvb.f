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
      subroutine o12ea_cvb(nparm1)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "opt2_cvb.fh"
#include "malloc_cvb.fh"

      call ddnewopt_cvb()
      have_solved_it=.false.

c  Find CIVBS & CIVBH :
      ivuse_s=0
      ivuse_h=0
      do 10 iv=1,nv
      if(tstcnt_cvb(w(lc(iv)),4))ivuse_s=iv
10    if(tstcnt_cvb(w(lc(iv)),5))ivuse_h=iv
      ivuse2=3
      if(ivuse_h.eq.3.or.ivuse_s.eq.3)ivuse2=2
      if(ivuse_h.eq.2.or.ivuse_s.eq.2)ivuse2=4
      if(ivuse2.gt.nv)ivuse2=1
      if(ivuse_s.ne.0.and.ivuse_h.ne.0)then
        i1 = mstackr_cvb(nparm1)
        i2 = mstackr_cvb(nparm1)
        i3 = mstackr_cvb(nparm1)
        i4 = mstackr_cvb(nvb+nprorb)
        call o12ea2_cvb(w(i1),w(i2),w(i3),nparm1,
     >    w(lc(ivuse2)),w(lc(ivuse_s)),w(lc(ivuse_h)),
     >    w(lw(9)),w(lv(2)),w(i4))
        call mfreer_cvb(i1)
      else
        if(strucopt)then
          call ddguess_cvb(w(lv(2)),nvb,nprorb)
        else
          call ddguess_cvb([one],1,0)
        endif
      endif
      call str2vbc_cvb(w(lv(2)),w(lw(9)))
      call vb2cic_cvb(w(lw(9)),w(lc(3)))
      return
      end
