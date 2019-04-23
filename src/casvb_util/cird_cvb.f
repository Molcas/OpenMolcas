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
c  ************************************************
c  ** Subroutines to deal with CASSCF CI vectors **
c  ************************************************
c  ********************************
c  ** Routines involving CI only **
c  ********************************
      subroutine cird_cvb(cvec,recn)
      implicit real*8(a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension cvec(*)
      dimension idum(1)
c  *********************************************************************
c  *                                                                   *
c  *  CIRD   := Read CI vector.                                        *
c  *                                                                   *
c  *********************************************************************

      ivec=nint(cvec(1))
      iformat=iform_ci(ivec)
      if(iformat.eq.0)then
        ioffs=0
        call rdis_cvb(idum,1,recn,ioffs)
        iformat=idum(1)
        if(iformat.ne.iform_ci(ivec))then
          write(6,*)' Incompatible vector format on file.'
          write(6,*)' Read :',iformat,' present :',iform_ci(ivec)
          call abend_cvb()
        endif
        call rdis_cvb(icnt_ci(ivec),1,recn,ioffs)
        call rdrs_cvb(w(iaddr_ci(ivec)),ndet,recn,ioffs)
      else
        write(6,*)' Unsupported format in CIRD :',iformat
        call abend_cvb()
      endif
      return
      end
