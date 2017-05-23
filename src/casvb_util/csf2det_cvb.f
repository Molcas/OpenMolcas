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
      subroutine csf2det_cvb(vec,detvec,isym_loc,iWay)
      implicit real*8 (a-h,o-z)
#include "csfbas.fh"
#include "ciinfo.fh"
#include "rasdim.fh"
#include "rasscf.fh"
#include "WrkSpc.fh"
      dimension vec(*),detvec(*)

      if(iWay.eq.1)then
        if ( nac.eq.0 ) then
          detvec(1)=vec(1)
          return
        endif

        jCopy = 0
        call csdtvc(vec,detvec,iway,work(kdtoc),
     >              iwork(kicts(1)),isym_loc,jcopy)
      elseif(iWay.eq.2)then
        if ( nac.eq.0 ) then
          vec(1)=detvec(1)
          return
        endif

        jCopy = 0
        call csdtvc(vec,detvec,iway,work(kdtoc),
     >              iwork(kicts(1)),isym_loc,jcopy)
      endif
      return
      end
