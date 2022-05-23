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
      subroutine make_close_cvb(it)
      implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "applyh_cvb.fh"

#include "io_cvb.fh"
#include "idbl_cvb.fh"
#include "fio.fh"
      character*8 vec(11)
      vec(1)='TMP01'
      vec(2)='TMP02'
      vec(3)='TMP03'
      vec(4)='TMP04'
      vec(5)='TMP05'
      vec(6)='TMP06'
      vec(7)='TMP07'
      vec(8)='TMP08'
      vec(9)='TMP09'
      vec(10)='VBWFN'
      vec(11)='JOBIPH'
      il=10
      if(it.eq.0) il=10
      if(it.eq.1) il=11
c  Preassign some file names to identifiers :
      do n=1,MxFile
         do i=1,il
             if(isOpen(n).eq.1) then
              if (LuName(n).eq.vec(i)) then
c         print *,'closing ',LuName(n)
                call daclos(n)
              endif
              endif
      enddo
      enddo
      if(.not.variat) then
         call mkguga_free()
         call getmem('CICTL1','FREE','REAL',lw1_cvb,0)
         call getmem('TUVX','FREE','REAL',ltuvx_cvb,0)
         call getmem('DMAT','FREE','REAL',ldmat_cvb,0)
         call getmem('DSPN','FREE','REAL',ldspn_cvb,0)
         call getmem('PMAT','FREE','REAL',lpmat_cvb,0)
         call getmem('P2AS','FREE','REAL',lpa_cvb,0)
         call getmem('DIAF','FREE','REAL',ldiaf_cvb,0)
         call getmem('FOCC','FREE','REAL',ipfocc_cvb,0)
         call getmem('FI','FREE','REAL',lfi_cvb,0)
         call getmem('FA','FREE','REAL',lfa_cvb,0)
         call getmem('D1I','FREE','REAL',ld1i_cvb,0)
         call getmem('D1A','FREE','REAL',ld1a_cvb,0)
         call getmem('D1tot','FREE','REAL',ld1tot_cvb,0)
         call getmem('OCCN','FREE','REAL',loccn_cvb,0)
         call getmem('LCMO','FREE','REAL',lcmo_cvb,0)
      endif
      return
      end
