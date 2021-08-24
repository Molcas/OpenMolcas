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
      subroutine applyh_cvb(civec)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
#include "casvb.fh"
#include "rasscf_lucia.fh"
      dimension civec(*)
      save thr2
      data thr2/1.d-20/

      kH0_Pointer = lw1_cvb
      icivec=nint(civec(1))
      c_daxpy=zero

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
      call vb2mol_cvb(work(iaddr_ci(icivec)),work(lcim),isyml)

c  If only one irrep present keep down memory requirements:
      if(isymmx.gt.1.and.nci.ne.ndet)then
        lcim2=mstackrz_cvb(nci)
        ibasemx=max(ibasemx,mstackr_cvb(0))
        cnrm=ddot_(nci,work(lcim),1,work(lcim),1)
c  If anything there, apply Hamiltonian to vector of this symmetry :
        if(cnrm.gt.thr2)then
          call sigmadet_cvb(work(lcim),work(lcim2),isyml,nci)
          if(c_daxpy.ne.zero)
     >      call daxpy_(nci,c_daxpy,work(lcim),1,work(lcim2),1)
        else
          if(c_daxpy.ne.zero)
     >      call daxpy_(nci,c_daxpy,work(lcim),1,work(lcim2),1)
        endif
        call mol2vb_cvb(work(iaddr_ci(icivec)),work(lcim2),isyml)
        call mfreer_cvb(lcim2)
      else
        call fzero(work(iaddr_ci(icivec)),nci)
        cnrm=ddot_(nci,work(lcim),1,work(lcim),1)
c  If anything there, apply Hamiltonian to vector of this symmetry :
        if(cnrm.gt.thr2)then
          call fzero(work(iaddr_ci(icivec)),nci)
          call sigmadet_cvb(work(lcim),work(iaddr_ci(icivec)),isyml,nci)
          if(c_daxpy.ne.zero)
     >      call daxpy_(nci,c_daxpy,work(lcim),1,work(iaddr_ci(icivec)),
     >                  1)
          call fmove_cvb(work(iaddr_ci(icivec)),work(lcim),nci)
        else
          if(c_daxpy.ne.zero)
     >      call daxpy_(nci,c_daxpy,work(lcim),1,work(iaddr_ci(icivec)),
     >                  1)
          call fmove_cvb(work(iaddr_ci(icivec)),work(lcim),nci)
        endif
        call mol2vb_cvb(work(iaddr_ci(icivec)),work(lcim),isyml)
      endif
      call mfreer_cvb(lcim)
1000  continue
      return
      end
