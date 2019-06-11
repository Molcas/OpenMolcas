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
      subroutine applyhpcx_cvb(civec,c_daxpy)
c  Exact copy if applyh except for c_daxpy in arg list.
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension civec(*)
      save thr2
      data thr2/1.d-20/

      icivec=nint(civec(1))
      n_applyh=n_applyh+1
      call setcnt2_cvb(icivec,0)
      if(iform_ci(icivec).ne.0)then
        write(6,*)' Unsupported format in APPLYH :',iform_ci(icivec)
        call abend_cvb()
      endif

c  (NIRREP may be altered in loop)
      isymmx=nirrep
      do 1000 isyml=1,isymmx
      nci=ncivb(isyml)
      lcim=mstackrz_cvb(nci)
      ibasemx=max(ibasemx,mstackr_cvb(0))
      call vb2mol_cvb(w(iaddr_ci(icivec)),w(lcim),isyml)

c  If only one irrep present keep down memory requirements:
      if(isymmx.gt.1.and.nci.ne.ndet)then
        lcim2=mstackrz_cvb(nci)
        ibasemx=max(ibasemx,mstackr_cvb(0))
        cnrm=ddot_(nci,w(lcim),1,w(lcim),1)
c  If anything there, apply Hamiltonian to vector of this symmetry :
        if(cnrm.gt.thr2)then
          call sigmadet_cvb(w(lcim),w(lcim2),isyml,absym(5),nci)
          if(c_daxpy.ne.zero)
     >      call daxpy_(nci,c_daxpy,w(lcim),1,w(lcim2),1)
        else
          if(c_daxpy.ne.zero)
     >      call daxpy_(nci,c_daxpy,w(lcim),1,w(lcim2),1)
        endif
        call mol2vb_cvb(w(iaddr_ci(icivec)),w(lcim2),isyml)
        call mfreer_cvb(lcim2)
      else
        call fzero(w(iaddr_ci(icivec)),nci)
        cnrm=ddot_(nci,w(lcim),1,w(lcim),1)
c  If anything there, apply Hamiltonian to vector of this symmetry :
        if(cnrm.gt.thr2)then
          call fzero(w(iaddr_ci(icivec)),nci)
          call sigmadet_cvb(w(lcim),w(iaddr_ci(icivec)),isyml,absym(5),
     &      nci)
          if(c_daxpy.ne.zero)
     >      call daxpy_(nci,c_daxpy,w(lcim),1,w(iaddr_ci(icivec)),1)
          call fmove_cvb(w(iaddr_ci(icivec)),w(lcim),nci)
        else
          if(c_daxpy.ne.zero)
     >      call daxpy_(nci,c_daxpy,w(lcim),1,w(iaddr_ci(icivec)),1)
          call fmove_cvb(w(iaddr_ci(icivec)),w(lcim),nci)
        endif
        call mol2vb_cvb(w(iaddr_ci(icivec)),w(lcim),isyml)
      endif
      call mfreer_cvb(lcim)
1000  continue
      return
      end
