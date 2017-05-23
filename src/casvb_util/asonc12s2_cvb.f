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
      subroutine asonc12s2_cvb(c,sxc,nvec,nprm,
     >   civb,civbs,
     >   orbs,gjorb,gjorb2,gjorb3,cvbdet,
     >   cvb,
     >   vec_all)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension c(nprm,nvec),sxc(nprm,nvec)
      dimension civb(ndet),civbs(ndet)
      dimension orbs(norb,norb),gjorb(*),gjorb2(*),gjorb3(*)
      dimension cvb(nvb),cvbdet(ndetvb)
      dimension vec_all(npr)
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
c  (CIVB set in O12SA :)
      call cizero_cvb(civbs)
      call oneexc_cvb(civb,civbs,vec_all,.false.,0)
      call str2vbf_cvb(vec_all(nprorb+1),cvbdet)
      call vb2ciaf_cvb(cvbdet,civbs)
      call applyts_cvb(civbs,orbs,gjorb,gjorb2,gjorb3)

      call ci2vbg_cvb(civbs,cvbdet)
      call vb2strg_cvb(cvbdet,vec_all(nprorb+1))
      call fzero(vec_all,nprorb)
      call onedens_cvb(civb,civbs,vec_all,.false.,0)
      call all2free_cvb(vec_all,sxc(ic1,ivec),1)
100   if(.not.strucopt)sxc(1,ivec)=ddot_(nvb,cvb,1,vec_all(nprorb+1),1)

      return
      entry asonc12sinit_cvb(ippinp)
      iter=0
      ipp=ippinp
      return
      end
