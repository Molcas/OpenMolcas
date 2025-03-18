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
       Logical Function ipopen(nconf,page)
       use ipPage
       use stdalloc, only: mma_maxDBLE
!
!      Initiate the whole lot.
!
       Implicit Real*8(a-h,o-z)
       Logical page
!
!      Ask how much memory is available
!
       Call mma_maxDBLE(nMax)
       nmax=nmax/2
!
       If (Page) Then
!
!         Initiate for disk based storage.
!
          If (.Not.DiskBased) Then
             Lu_ip=21
             Lu_ip=IsFreeUnit(Lu_ip)
             Call Daname(Lu_ip,'TEMPCIV')
             DiskBased=.True.
          End If
!
!         n  : Length of CI-vector
!         ida: disk address
!
               n(0:Max_CI_Vectors)=0
             ida(0:Max_CI_Vectors)=-1
          Status(0:Max_CI_Vectors)=Null_Vector
!
!         iDisk_Addr_End: next free disk address
!         n_CI_Vectors : number of CI-vectors
!
          iDisk_Addr_End=0
          n_CI_Vectors=0
!
       Else

          If (DiskBased) Then
             Call ipTerm()
             DiskBased=.False.
          End If

       End If
!
       ipopen=DiskBased
!
       Return
! Avoid unused argument warnings
       If (.False.) Call Unused_integer(nconf)
       End
