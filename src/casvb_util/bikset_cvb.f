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
      subroutine bikset_cvb(aikcof,bikcof,
     > nel,nalf,i2s,ndet,ifns,kbasis,share,iprint)
      implicit real*8 (a-h,o-z)
      logical share
      dimension aikcof(ndet,ifns),bikcof(ndet,ifns)
#include "malloc_cvb.fh"

      if(i2s.ne.2*nalf-nel)then
        nalf_use=(i2s+nel)/2
        ndet_use=ndet_cvb(nel,nalf_use)
        i1 = mstackr_cvb(ndet_use*ifns)
        i2 = mstackr_cvb(ndet_use*ifns)
        call biksmain_cvb(w(i1),w(i2),
     >   nel,nalf_use,ndet_use,ifns,kbasis,.false.,iprint)
        call sminus_cvb(w(i1),bikcof,
     >    nel,nalf_use,nalf,ifns)
        if(.not.share)call sminus_cvb(w(i2),aikcof,
     >    nel,nalf_use,nalf,ifns)
        call mfreer_cvb(i1)
      else
        call biksmain_cvb(aikcof,bikcof,
     >   nel,nalf,ndet,ifns,kbasis,share,iprint)
      endif
      return
      end
