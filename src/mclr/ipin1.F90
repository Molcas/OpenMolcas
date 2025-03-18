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
       Integer Function ipin1(ii,nn)
!
!      Object: return pointer to vector ii with a length of nn and
!              make the vector available in memory as W(ii)%Vec
!
       use ipPage
       use stdalloc, only: mma_allocate, mma_deallocate
       use Constants, only: Zero
       Implicit Integer (a-h,o-z)
       Real*8, Allocatable:: Tmp(:)
!
       If (ii.gt.Max_CI_Vectors) Then
          Write (6,*) 'ipin1: ii.gt.Max_CI_Vectors'
          Write (6,*) 'ii,Max_CI_Vectors=',ii,Max_CI_Vectors
          Call Abend()
       End If
!
!
       If (Status(ii).eq.In_Memory) Then
!
!-------- ii is in memory
!

!         If the size of the vector is larger than was originally set
!         resize the reservation and copy the content

          If (nn.gt.n(ii)) Then
             Call mma_allocate(Tmp,nn,Label='Tmp')
             Tmp(:) = Zero
             Tmp(1:n(ii)) = W(ii)%Vec(:)
             Call mma_deallocate(W(ii)%Vec)
             Call mma_allocate(W(ii)%Vec,nn,Label='ipin1')
             W(ii)%Vec(:) = Tmp(:)
             Call mma_deallocate(Tmp)
             n(ii)=nn
          End If
!
          ip1=ii
!
       Else If (Status(ii).eq.On_Disk) Then
!
!         ii is on disk
!
          Call mma_allocate(W(ii)%Vec,Max(n(ii),nn),Label='ipin1')
          W(ii)%Vec(:)=Zero
!
          nnn=Min(n(ii),nn)
!
!         pick up from disk
!
          idisk=ida(ii)
          Call dDafile(Lu_ip,Read,W(ii)%Vec,nnn,idisk)
          Status(ii)=In_Memory
!
          ip1=ii
!
       Else If (Status(ii).eq.Null_Vector) Then
!
           ip1=-1
!
       Else
!
          ip1=-1
          Write (6,*)
          Write (6,*) 'ipIn1: illegal Status(ii)'
          Write (6,*) 'ii=',ii
          Write (6,*)
          Call Abend()
!
       End If
!
       ipin1 = ip1
!
       Return
       End
