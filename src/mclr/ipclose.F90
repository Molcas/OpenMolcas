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
       Integer Function ipclose(ia)
       use ipPage
       use stdalloc, only: mma_deallocate
!
!      Object: release all vectors above and including the vector
!              indexed ia.
!
       Real*8 rdum(1)
!
       If (ia.gt.Max_CI_Vectors) Then
          Write (6,*) 'ipclose: ia.gt.Max_CI_Vectors'
          Write (6,*) 'ia,Max_CI_Vectors=',ia,Max_CI_Vectors
          Call Abend()
       End If
!
!
!      Update iDisk_Addr_End
!
       iDisk_Addr_End=0
       If (ia.lt.0) then
!
          n_CI_Vectors=0
!
       Else
!
          n_CI_Vectors=ia-1
          If (DiskBased) Then
             Do ii=1,ia-1
                If (Status(ii).ne.Null_Vector)                          &
     &             Call dDafile(Lu_ip,dWrite,rdum(1),n(ii),             &
     &                          iDisk_Addr_End)
             End do
          End If
!
       End If
!
!      Release memory and flag as a null vector
!
       Do ii=Max(ia,0),Max_CI_Vectors
          If (Status(ii).eq.In_Memory) Then
             Call mma_deallocate(W(ii)%Vec)
             ida(ii)=-1
             n(ii)=0
             Status(ii)=Null_Vector
          End If
       End Do
!
       If (diskbased.and.ia.lt.0)  Then
          Call DACLOS(Lu_ip)
          DiskBased=.False.
       End If
       ipclose=0
!
       Return
       End
