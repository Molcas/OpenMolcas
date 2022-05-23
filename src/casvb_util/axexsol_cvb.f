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
      subroutine axexsol_cvb(ap,dum,itdav,maxdav,nfrdim1,
     >  solp,solp_res,eig,eig_res)
c  Diagonalize Hamiltonian in Davidson subspace:
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension ap(maxdav,maxdav),solp(maxdav),solp_res(maxdav)

      i1 = mstackr_cvb(itdav)
      i2 = mstackr_cvb(itdav*itdav)

      call axexsol2_cvb(ap,w(i1),w(i2),dum,itdav,maxdav,
     >  solp,solp_res,eig,eig_res)

      call mfreer_cvb(i1)
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(nfrdim1)
      end
