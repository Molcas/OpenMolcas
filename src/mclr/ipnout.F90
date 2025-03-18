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
       Integer Function ipnout(iii)
       use ipPage
       use stdalloc, only: mma_deallocate
!
!      Object: write all vectors in memory on disk but vector iii
!
       Implicit Integer (a-h,o-z)
!
       If (iii.gt.Max_CI_Vectors) Then
          Write (6,*) 'ipout: iii.gt.Max_CI_Vectors'
          Write (6,*) 'iii,Max_CI_Vectors=',iii,Max_CI_Vectors
          Call Abend()
       End If
!
       ipnout=0
       If (.not.DiskBased) Return
!
       Do ii=1,Max_CI_Vectors
!
          If (Status(ii).eq.In_Memory .and. ii.ne.iii) Then
             idisk=ida(ii)
             nn=n(ii)
             Call dDafile(Lu_ip,Write,W(ii)%Vec,nn,idisk)
             Status(ii)=On_Disk
             Call mma_deallocate(W(ii)%Vec)
          End If
!
       End Do
!
       Return
       End
