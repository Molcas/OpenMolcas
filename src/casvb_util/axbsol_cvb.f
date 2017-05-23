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
      subroutine axbsol_cvb(dum,rhsp,itdav,maxdav,nfrdim1,
     >  solp,solp_res,eig,eig_res)
      implicit real*8 (a-h,o-z)
#include "direct_cvb.fh"
      dimension rhsp(maxdav),solp(maxdav),solp_res(maxdav)
      save one
      data one/1d0/

      call fmove(rhsp,solp,itdav)
      eig=dnrm2_(itdav,solp,1)
      call dscal_(itdav,one/eig,solp,1)
      eig_res=eig
      call fmove(solp,solp_res,itdav)
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_real(dum)
        call Unused_integer(nfrdim1)
      end if
      end
