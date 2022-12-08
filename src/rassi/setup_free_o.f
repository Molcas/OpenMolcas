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
      use nq_Grid, only: Pax
#include "stdalloc.fh"
      Call mma_allocate(Pax,3,3,Label='Pax')
      Pax(:,:)=0.0D0
      Call DCOPY_(3,[1.0D0],0,Pax,4)
      Return
      End
*
      Subroutine Free_O()
      use nq_Grid, only: Pax
#include "stdalloc.fh"
      Call mma_deallocate(Pax)
      Return
      End
