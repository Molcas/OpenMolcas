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
      subroutine io_init_cvb()
      implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "io_cvb.fh"

      nrec=0
      thresh_io=1d-5
      call izero(iorder,mxunits)
      call daninit_cvb()

c  Preassign some file names to identifiers :
      call setfn_cvb(recn_jobold,'JOBOLD')
      call setfn_cvb(recn_jobiph,'JOBIPH')
      call setfn_cvb(recn_oneint,'ONEINT')
      call setfn_cvb(recn_vbwfn,'VBWFN')
      call setfn_cvb(recn_tmp01,'TMP01')
      call setfn_cvb(recn_tmp02,'TMP02')
      call setfn_cvb(recn_tmp03,'TMP03')
      call setfn_cvb(recn_tmp04,'TMP04')
      call setfn_cvb(recn_tmp05,'TMP05')
      call setfn_cvb(recn_tmp06,'TMP06')
      call setfn_cvb(recn_tmp07,'TMP07')
      call setfn_cvb(recn_tmp08,'TMP08')
      call setfn_cvb(recn_tmp09,'TMP09')
      return
      end
