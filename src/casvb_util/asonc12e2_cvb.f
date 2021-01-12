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
      subroutine asonc12e2_cvb(c,axc,sxc,nvec,nprm,
     >   civb,civbh,civbs,
     >   orbs,gjorb,gjorb2,gjorb3,cvbdet,
     >   cvb,
     >   vec_all)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension c(nprm,nvec),axc(nprm,nvec),sxc(nprm,nvec)
      dimension civb(ndet),civbh(ndet),civbs(ndet)
      dimension orbs(norb,norb),gjorb(*),gjorb2(*),gjorb3(*)
      dimension cvbdet(ndetvb)
      dimension vec_all(npr)
      dimension cvb(nvb)
      save iter,ipp

      iter=iter+1
      if(ipp.ge.2)then
        write(6,'(/,a,i5,a,f10.3,a)')' Davidson iteration',iter,
     >    ' at',tim_cvb(cpu0),' CPU seconds'
        write(6,'(a)')
     >    ' -----------------------------------------------'
      endif

c  If no optimization of structure coefficients we are doing
c  "Augmented" calc:
      if(strucopt)then
        ic1=1
      else
        ic1=2
      endif

      do 100 ivec=1,nvec
      call free2all_cvb(c(ic1,ivec),vec_all,1)
      if(.not.strucopt)call daxpy_(nvb,c(1,ivec),cvb,1,
     >  vec_all(nprorb+1),1)
c  (CIVB set in O12EA :)
      call cizero_cvb(civbs)
      call oneexc_cvb(civb,civbs,vec_all,.false.,0)
      call str2vbf_cvb(vec_all(nprorb+1),cvbdet)
      call vb2ciaf_cvb(cvbdet,civbs)
      call cicopy_cvb(civbs,civbh)
      call makecivbhs_cvb(civbh,civbs,orbs,gjorb,gjorb2,gjorb3)

      call ci2vbg_cvb(civbh,cvbdet)
      call vb2strg_cvb(cvbdet,vec_all(nprorb+1))
      call fzero(vec_all,nprorb)
      call onedens_cvb(civb,civbh,vec_all,.false.,0)
      call all2free_cvb(vec_all,axc(ic1,ivec),1)
      if(.not.strucopt)axc(1,ivec)=ddot_(nvb,cvb,1,vec_all(nprorb+1),1)

      call ci2vbg_cvb(civbs,cvbdet)
      call vb2strg_cvb(cvbdet,vec_all(nprorb+1))
      call fzero(vec_all,nprorb)
      call onedens_cvb(civb,civbs,vec_all,.false.,0)
      call all2free_cvb(vec_all,sxc(ic1,ivec),1)
      if(.not.strucopt)sxc(1,ivec)=ddot_(nvb,cvb,1,vec_all(nprorb+1),1)
100   continue

      return
      entry asonc12einit_cvb(ippinp)
      iter=0
      ipp=ippinp
      return
      end
