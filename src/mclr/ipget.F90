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
       Integer Function ipget(nn)
!
!      Get the index of a vector with the length nn.
!      Memory or disk space is allocated.
!
       use ipPage
       use stdalloc, only: mma_allocate, mma_deallocate
       use Constants, only: Zero
       Implicit Integer (a-h,o-z)
       Character*4 Label
!
!      Take the next memory slot.
!
       n_CI_Vectors=n_CI_Vectors+1
       ipget=n_CI_Vectors
!
       If (n_CI_Vectors.gt.Max_CI_Vectors) Then
          Write(6,*) 'Number of CI vectors higher than Max_CI_Vectors'
          Write(6,*) 'Max_CI_Vectors=',Max_CI_Vectors
          Call Abend( )
       End If
!
       ida(ipget)=iDisk_Addr_End
       n(ipget)=nn
!
!----- Allocate memory for vector if of non-zero  length
!
       If (nn.gt.0) Then
          Write (Label,'(I3.3)') n_CI_Vectors
          Call mma_allocate(W(ipget)%Vec,nn,Label='ipget'//Label)
          Status(ipget)=In_Memory
          W(ipget)%Vec(:)=Zero
       Else
!         Status(ipget)=Null_Vector
!
!         The calling code doesn't have the logic to handle the
!         case that W(i)%Vec is not allocated. Hence, we have
!         to make a dummy allocation to make sure that the compiler
!         doesn't puke.
          n(ipget)=1
          Write (Label,'(I3.3)') n_CI_Vectors
          Call mma_allocate(W(ipget)%Vec,1,Label='ipget'//Label)
          Status(ipget)=In_Memory
          W(ipget)%Vec(:)=Zero
       End If
!
!      If diskbased mode put vector on disc and release memory
!
       If (DiskBased) then
          If (Status(ipget).ne.Null_Vector) Then
             Call dDafile(Lu_ip,Write,W(ipget)%Vec,nn,iDisk_Addr_End)
             Status(ipget)=On_Disk
             Call mma_deallocate(W(ipget)%Vec)
          End If
       End If
!
       Return
       End
