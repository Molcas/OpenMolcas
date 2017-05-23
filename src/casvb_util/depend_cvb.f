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
      subroutine depend_cvb(chr1,chr2)
      implicit real*8 (a-h,o-z)
      character*(*) chr1,chr2
#include "make_cvb.fh"

      call mkafter_cvb(chr1,chr2)
      call touchdepend_cvb(chr1,chr2)

      if(iprint.ge.10)then
        write(6,*)' IOFFS :',(ioffs(ii),ii=1,nobj+1)
        write(6,*)' JOFFS :',(joffs(ii),ii=1,nobj+1)
        write(6,*)' I_DEP_ON_J :',(i_dep_on_j(ii),ii=1,ndep_ij)
        write(6,*)' J_DEP_ON_I :',(j_dep_on_i(ii),ii=1,ndep_ji)
      endif
      return
      end
