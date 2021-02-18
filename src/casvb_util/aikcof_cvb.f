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
      subroutine aikcof_cvb(aikcof,bikcof,
     > ndet,ifns,kbasis,share,
     > sovr)
      implicit real*8 (a-h,o-z)
      logical share
      dimension aikcof(ndet,ifns),bikcof(ndet,ifns)
      dimension sovr(ifns,ifns)

      if(kbasis.eq.6)return
c
c  Generate mapping from determinants to spin functions
c  (If KBASIS<=2 then AIKCOF=BIKCOF and they (probably) share memory)
      if(kbasis.gt.2)then
        call mxattb_cvb(bikcof,bikcof,ifns,ndet,ifns,sovr)
        call mxinv_cvb(sovr,ifns)
        call mxatb_cvb(bikcof,sovr,ndet,ifns,ifns,aikcof)
      elseif(.not.share)then
        call fmove_cvb(bikcof,aikcof,ndet*ifns)
      endif
      return
      end
