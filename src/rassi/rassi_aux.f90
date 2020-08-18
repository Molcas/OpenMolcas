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
      Module RASSI_AUX
      Logical :: AO_Mode=.False.
      Integer, Allocatable:: TocM(:), jDisk_TDM(:,:), JOB_INDEX(:)
      Real*8, Allocatable:: CMO1(:), CMO2(:), DMAB(:)
      Integer NASHT_Save, mTRA
      Integer :: JOB1_old=-1, JOB2_old=-1

      Contains

      Integer Function iDisk_TDM(I,J,K)
         Integer I, J, I_Max, J_Min, K
         I_Max=Max(I,J)
         J_Min=Min(I,J)
         iDisk_TDM=jDisk_TDM(K,I_Max*(I_Max-1)/2+J_Min)
      End Function iDisk_TDM

      End Module RASSI_AUX
