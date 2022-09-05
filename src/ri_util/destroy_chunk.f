!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Destroy_Chunk()
      use Chunk_Mod
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: Is_Real_Par
#endif
#include "stdalloc.fh"
#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
         Call GA_Destroy(ip_Chunk)
         Call mma_deallocate(iMap)
      Else
         Call mma_deallocate(Chunk)
      End If
#else
      Call mma_deallocate(Chunk)
#endif
      Return
      End
