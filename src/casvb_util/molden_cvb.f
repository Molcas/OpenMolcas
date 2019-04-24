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
      subroutine molden_cvb()
      implicit real*8 (a-h,o-z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "rctfld.fh"

      dimension Dummy(1)

      call daname_cvb(JOBIPH,'JOBIPH')
      idisk=0
      call idafile(JOBIPH,2,iadr15,15,idisk)
      Dummy = 0.0D0
      if(.not.lRf) call interf(0,Dummy,0,1)
      return
      end
