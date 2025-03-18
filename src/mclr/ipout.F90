!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************
       Integer Function ipout(ii)
!
!      ipout will page out vector ii to disk and free the memory area
!
       use ipPage
       use stdalloc, only: mma_deallocate
       Implicit Integer (a-h,o-z)
!
       ipout=0
       If (.not.diskbased) Return
!
       If (Status(ii).eq.In_Memory .and. ii.gt.0) Then
          idisk=ida(ii)
          nn=n(ii)
          Call dDafile(Lu_ip,Write,W(ii)%Vec,nn,idisk)
          Status(ii)=On_Disk
          Call mma_deallocate(W(ii)%Vec)
       Else
          ipout=-1
       End if
!
       Return
       End
