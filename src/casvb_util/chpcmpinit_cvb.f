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
      subroutine chpcmpinit_cvb()
      implicit real*8 (a-h,o-z)
#include "lstprm_cvb.fh"

      do 100 i=1,mxprm
100   lstprm(i)=iunset
      iprm=0
      return
      entry chpcmp0_cvb()
      iprm=0
      return
      entry chpcmp2_cvb(itst,iret)
      iprm=iprm+1
      if(iprm.gt.mxprm)then
        write(6,*)' Dimensioning error in CHPCMP2!',iprm,mxprm
        call abend_cvb()
      endif
      iret=lstprm(iprm)
      lstprm(iprm)=itst
      return
      end
      function recinpcmp_cvb(ifield)
      implicit real*8 (a-h,o-z)
      logical recinpcmp_cvb
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"

      if(.not.valid_cvb(recinp_old))then
        recinpcmp_cvb=.true.
      else
        call rdioff_cvb(ifield,recinp,ioff1)
        call rdioff_cvb(ifield+1,recinp,ioff2)
        call rdioff_cvb(ifield,recinp_old,joff1)
        call rdioff_cvb(ifield+1,recinp_old,joff2)
        if(ioff2-ioff1.ne.joff2-joff1)then
          recinpcmp_cvb=.true.
        else
          i1=mstackr_cvb(ioff2-ioff1)
          j1=mstackr_cvb(joff2-joff1)
          call rdr_cvb(w(i1),ioff2-ioff1,recinp,ioff1)
          call rdr_cvb(w(j1),joff2-joff1,recinp_old,joff1)
          do 100 i=0,ioff2-ioff1-1
          if(w(i+i1).ne.w(i+j1))then
            recinpcmp_cvb=.true.
            goto 200
          endif
100       continue
          recinpcmp_cvb=.false.
200       call mfreer_cvb(i1)
        endif
      endif
      return
      end
