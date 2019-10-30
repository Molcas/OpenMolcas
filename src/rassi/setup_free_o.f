************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Setup_O()
      IMPLICIT REAL*8 (A-H,O-Z)
#include "nq_info.fh"
#include "WrkSpc.fh"
      Call GetMem('ip_O','ALLO','REAL',ip_O,9)
      Call DCOPY_(9,[0.0d0],0,Work(ip_O),1)
      Call DCOPY_(3,[1.0D0],0,Work(ip_O),4)
      Return
      End
*
      Subroutine Free_O()
      IMPLICIT REAL*8 (A-H,O-Z)
#include "nq_info.fh"
#include "WrkSpc.fh"
      Call GetMem('ip_O','FREE','REAL',ip_O,9)
      Return
      End
