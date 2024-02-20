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
* Copyright (C) 1992,1994, Per Ake Malmqvist                           *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE PCLOSE()
      use stdalloc, only: mma_deallocate
      use fciqmc_interface, only: DoFCIQMC
      use pt2_guga_data
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  PER-AAKE MALMQUIST 92-12-07
C  Deallocates everything concerned with SGUGA, incl CI array.


#include "rasdim.fh"
#include "caspt2.fh"

      IF(DoCumulant) RETURN
      IF(DoFCIQMC) RETURN
      IF(NACTEL.EQ.0) RETURN
      IF(ISCF.NE.0) RETURN
      Call mma_deallocate(MVL)
      Call mma_deallocate(MVR)
      Call mma_deallocate(NOW1)
      Call mma_deallocate(IOW1)
      Call mma_deallocate(NOCP)
      Call mma_deallocate(IOCP)
      Call mma_deallocate(NOCSF)
      Call mma_deallocate(IOCSF)
      Call mma_deallocate(ICASE)
      Call mma_deallocate(ICOUP)
      Call mma_deallocate(VTAB)

      Call mma_deallocate(SGTMP)

      END SUBROUTINE PCLOSE
