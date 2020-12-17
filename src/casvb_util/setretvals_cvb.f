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
      subroutine setretvals_cvb(esym,n_iter)
      implicit real*8 (a-h,o-z)
      dimension esym(*)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"

      if(nac.eq.0) then
        ener(1,iter)=emy
      else
        do 100 jroot=1,lroots
        ener(jroot,iter)=esym(stsym)
100     continue
      endif
      iterci=n_iter
      return
      end
