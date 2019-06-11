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
      subroutine o12ea2_cvb(c,sxc,axc,nprm,
     >   civb,civbs,civbh,
     >   cvbdet,cvb,vec_all)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension c(nprm),sxc(nprm),axc(nprm)
      dimension civb(ndet),civbs(ndet),civbh(ndet)
      dimension cvbdet(ndetvb),cvb(nvb),vec_all(npr)

c  If no optimization of structure coefficients we are doing
c  "Augmented" calc:
      if(strucopt)then
        ic1=1
      else
        ic1=2
      endif

      call str2vbc_cvb(cvb,cvbdet)
      call vb2cic_cvb(cvbdet,civb)

      call ci2vbg_cvb(civbh,cvbdet)
      call vb2strg_cvb(cvbdet,vec_all(nprorb+1))
      call fzero(vec_all,nprorb)
      call onedens_cvb(civb,civbh,vec_all,.false.,0)
      call all2free_cvb(vec_all,axc(ic1),1)
      if(.not.strucopt)axc(1)=ddot_(nvb,cvb,1,vec_all(nprorb+1),1)

      call ci2vbg_cvb(civbs,cvbdet)
      call vb2strg_cvb(cvbdet,vec_all(nprorb+1))
      call fzero(vec_all,nprorb)
      call onedens_cvb(civb,civbs,vec_all,.false.,0)
      call all2free_cvb(vec_all,sxc(ic1),1)
      if(.not.strucopt)sxc(1)=ddot_(nvb,cvb,1,vec_all(nprorb+1),1)

      call fzero(vec_all,nprorb)
      call fmove_cvb(cvb,vec_all(nprorb+1),nvb)
      call all2free_cvb(vec_all,c(ic1),1)
      if(.not.strucopt)c(1)=ddot_(nvb,cvb,1,vec_all(nprorb+1),1)

      cnrm=ddot_(nprm,c,1,sxc,1)
      call dscal_(nprm,one/sqrt(cnrm),c,1)
      call dscal_(nprm,one/sqrt(cnrm),sxc,1)
      call dscal_(nprm,one/sqrt(cnrm),axc,1)
      call ddrestv_cvb(c,axc,sxc,nprm,0,.true.,.true.)

      return
      end
