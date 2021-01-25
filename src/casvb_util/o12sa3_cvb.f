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
      subroutine o12sa3_cvb(vec,cvb,
     >  orbs,gjorb,gjorb2,gjorb3,
     >  civec,civecp,civb,cvbdet,vec_all,
     >  nvb,nprorb,nparm1,strucopt)
      implicit real*8 (a-h,o-z)
      logical strucopt
#include "opt_cvb.fh"
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"

      dimension vec(nparm1)
      dimension cvb(nvb),civec(*),civecp(*),civb(*),cvbdet(*)
      dimension vec_all(nparm1)
      dimension orbs(*),gjorb(*),gjorb2(*),gjorb3(*)

      call makegjorbs_cvb(orbs,gjorb,gjorb2,gjorb3)

      call str2vbc_cvb(cvb,cvbdet)
      call vb2cic_cvb(cvbdet,civb)

      call makecivecp_cvb(civec,civecp,orbs)
      call ci2vbg_cvb(civecp,cvbdet)
      call vb2strg_cvb(cvbdet,vec_all(nprorb+1))
      call fzero(vec_all,nprorb)
      call onedens_cvb(civb,civecp,vec_all,.false.,0)
c  If no optimization of structure coefficients we are doing
c  "Augmented" calc:
      if(strucopt)then
        ic1=1
      else
        ic1=2
      endif
      call all2free_cvb(vec_all,vec(ic1),1)
      if(.not.strucopt)vec(1)=ddot_(nvb,cvb,1,vec_all(nprorb+1),1)
      call ddrhs_cvb(vec,nparm1,0)

      call str2vbc_cvb(cvb,cvbdet)
      call vb2cic_cvb(cvbdet,civb)

      return
      end
