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
       Integer Function opout(ii)
!
!      opout will release the memory area of vector ii without updating
!      the disk
!
       use ipPage
       use stdalloc, only: mma_deallocate
       Implicit Integer (a-h,o-z)
!
       If (ii.gt.Max_CI_Vectors) Then
          Write (6,*) 'opout: ii.gt.Max_CI_Vectors'
          Write (6,*) 'ii,Max_CI_Vectors=',ii,Max_CI_Vectors
          Call Abend()
       End If
!
       opout=0
       If (.not.diskbased) Return
!
       If (Status(ii).eq.In_Memory .and. ii.gt.0) Then
          Status(ii)=On_Disk
          Call mma_deallocate(W(ii)%Vec)
       Else
          opout=-1
       End If
!
       Return
       End
