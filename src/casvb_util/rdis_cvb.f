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
      subroutine rdis_cvb(ivec,n,file_id,ioffset)
      implicit real*8 (a-h,o-z)
#include "idbl_cvb.fh"
      dimension ivec(n),ibuf(8)

      nreals=n/idbl
      nrem=n-nreals*idbl
      if(nreals.gt.0)call rdlow_cvb(ivec,nreals,file_id,ioffset)
      if(nrem.gt.0)then
        call rdlow_cvb(ibuf,1,file_id,nreals+ioffset)
        call imove_cvb(ibuf,ivec(1+nreals*idbl),nrem)
      endif
      if(nrem.eq.0)then
        ioffset=ioffset+nreals
      else
        ioffset=ioffset+nreals+1
      endif
      return
      end
