************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
*                                                                      *
*  A collection of subroutines that makes it possible
*  to page out ci vectors that are not in use at the
*  moment
*
*  These functions and functions are:
*
*      Logical Function ipopen(nconf,page)
*      Integer Function ipclose(ia)
*      Integer Function ipget(nn)
*      Integer Function ipin(ii)
*      Integer Function ipin1(ii,nn)
*      Integer Function ipnout(iii)
*      Integer Function opout(ii)
*      Integer Function ipout(ii)
*      Subroutine ipterm()
*      Subroutine ipinit()
*
*                                                                      *
************************************************************************
*                                                                      *
       Logical Function ipopen(nconf,page)
       use ipPage
*
*      Initiate the whole lot.
*
       Implicit Real*8(a-h,o-z)
       Logical page
#include "stdalloc.fh"
*
*      Ask how much memory is available
*
       Call mma_maxDBLE(nMax)
       nmax=nmax/2
       npp=nconf*10
*
       If (Page) Then
*
*         Initiate for disk based storage.
*
          If (.Not.DiskBased) Then
             Lu_ip=21
             Lu_ip=IsFreeUnit(Lu_ip)
             Call Daname(Lu_ip,'TEMPCIV')
             DiskBased=.True.
          End If
*
*         n  : Length of CI-vector
*         ida: disk address
*
               n(0:Max_CI_Vectors)=0
             ida(0:Max_CI_Vectors)=-1
          Status(0:Max_CI_Vectors)=Null_Vector
*
*         iDisk_Addr_End: next free disk address
*         n_CI_Vectors : number of CI-vectors
*
          iDisk_Addr_End=0
          n_CI_Vectors=0
*
       Else

          If (DiskBased) Then
             Call ipTerm()
             DiskBased=.False.
          End If

       End If
*
       ipopen=DiskBased
*
       Return
       End
*                                                                      *
************************************************************************
*                                                                      *
       Integer Function ipclose(ia)
       use ipPage
*
*      Object: release all vectors above and including the vector
*              indexed ia.
*
#include "stdalloc.fh"
       Real*8 rdum(1)
*
       If (ia.gt.Max_CI_Vectors) Then
          Write (6,*) 'ipclose: ia.gt.Max_CI_Vectors'
          Write (6,*) 'ia,Max_CI_Vectors=',ia,Max_CI_Vectors
          Call Abend()
       End If
*
*
*      Update iDisk_Addr_End
*
       iDisk_Addr_End=0
       If (ia.lt.0) then
*
          n_CI_Vectors=0
*
       Else
*
          n_CI_Vectors=ia-1
          If (DiskBased) Then
             Do ii=1,ia-1
                If (Status(ii).ne.Null_Vector)
     &             Call dDafile(Lu_ip,dWrite,rdum,n(ii),
     &                          iDisk_Addr_End)
             End do
          End If
*
       End If
*
*      Release memory and flag as a null vector
*
       Do ii=Max(ia,0),Max_CI_Vectors
          If (Status(ii).eq.In_Memory) Then
             Call mma_deallocate(W(ii)%Vec)
             ida(ii)=-1
             n(ii)=0
             Status(ii)=Null_Vector
          End If
       End Do
*
       If (diskbased.and.ia.lt.0)  Then
          Call DACLOS(Lu_ip)
          DiskBased=.False.
       End If
       ipclose=0
*
       Return
       End
*                                                                      *
************************************************************************
*                                                                      *
       Integer Function ipget(nn)
*
*      Get the index of a vector with the length nn.
*      Memory or disk space is allocated.
*
       use ipPage
       Implicit Integer (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
       Character*4 Label
*
*      Take the next memory slot.
*
       n_CI_Vectors=n_CI_Vectors+1
       ipget=n_CI_Vectors
*
       If (n_CI_Vectors.gt.Max_CI_Vectors) Then
          Write(6,*) 'Number of CI vectors higher than Max_CI_Vectors'
          Write(6,*) 'Max_CI_Vectors=',Max_CI_Vectors
          Call Abend( )
       End If
*
       ida(ipget)=iDisk_Addr_End
       n(ipget)=nn
*
*----- Allocate memory for vector if of non-zero  length
*
       If (nn.gt.0) Then
          Write (Label,'(I3.3)') n_CI_Vectors
          Call mma_allocate(W(ipget)%Vec,nn,Label='ipget'//Label)
          Status(ipget)=In_Memory
          W(ipget)%Vec(:)=Zero
       Else
          Status(ipget)=Null_Vector
       End If
*
*      If diskbased mode put vector on disc and release memory
*
       If (DiskBased) then
          If (Status(ipget).ne.Null_Vector) Then
             Call dDafile(Lu_ip,Write,W(ipget)%Vec,nn,iDisk_Addr_End)
             Status(ipget)=On_Disk
             Call mma_deallocate(W(ipget)%Vec)
          End If
       End If
*
       Return
       End
*                                                                      *
************************************************************************
*                                                                      *
       Integer Function ipin(ii)
*
*      Object: return pointer to vector ii with a length of n(ii) and
*              make the vector available in memory as W(ii)%Vec
*
*
       use ipPage
       Implicit Integer (a-h,o-z)
*
       nn=n(ii)
       ipin = ipin1(ii,nn)
*
       Return
       End
*                                                                      *
************************************************************************
*                                                                      *
       Integer Function ipin1(ii,nn)
*
*      Object: return pointer to vector ii with a length of nn and
*              make the vector available in memory as W(ii)%Vec
*
       use ipPage
       Implicit Integer (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
       Real*8, Allocatable:: Tmp(:)
*
       If (ii.gt.Max_CI_Vectors) Then
          Write (6,*) 'ipin1: ii.gt.Max_CI_Vectors'
          Write (6,*) 'ii,Max_CI_Vectors=',ii,Max_CI_Vectors
          Call Abend()
       End If
*
*
       If (Status(ii).eq.In_Memory) Then
*
*-------- ii is in memory
*

*         If the size of the vector is larger than was originally set
*         resize the reservation and copy the content

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
*
!         ip1=ii
          ip1=ip_of_Work(W(ii)%Vec(1))
*
       Else If (Status(ii).eq.On_Disk) Then
*
*         ii is on disk
*
          Call mma_allocate(W(ii)%Vec,Max(n(ii),nn),Label='ipin1')
          W(ii)%Vec(:)=Zero
*
          nnn=Min(n(ii),nn)
*
*         pick up from disk
*
          idisk=ida(ii)
          Call dDafile(Lu_ip,Read,W(ii)%Vec,nnn,idisk)
          Status(ii)=In_Memory
*
!         ip1=ii
          ip1=ip_of_Work(W(ii)%Vec(1))
*
       Else If (Status(ii).eq.Null_Vector) Then
*
           ip1=-1
*
       Else
*
          ip1=-1
          Write (6,*)
          Write (6,*) 'ipIn1: illegal Status(ii)'
          Write (6,*) 'ii=',ii
          Write (6,*)
          Call Abend()
*
       End If
*
       ipin1 = ip1
*
       Return
       End
*                                                                      *
************************************************************************
*                                                                      *
       Integer Function ipnout(iii)
       use ipPage
*
*      Object: write all vectors in memory on disk but vector iii
*
       Implicit Integer (a-h,o-z)
#include "stdalloc.fh"
*
       If (iii.gt.Max_CI_Vectors) Then
          Write (6,*) 'ipout: iii.gt.Max_CI_Vectors'
          Write (6,*) 'iii,Max_CI_Vectors=',iii,Max_CI_Vectors
          Call Abend()
       End If
*
       ipnout=0
       If (.not.DiskBased) Return
*
       Do ii=1,Max_CI_Vectors
*
          If (Status(ii).eq.In_Memory .and. ii.ne.iii) Then
             idisk=ida(ii)
             nn=n(ii)
             Call dDafile(Lu_ip,Write,W(ii)%Vec,nn,idisk)
             Status(ii)=On_Disk
             Call mma_deallocate(W(ii)%Vec)
          End If
*
       End Do
*
       Return
       End
*                                                                      *
************************************************************************
*                                                                      *
       Integer Function opout(ii)
*
*      opout will release the memory area of vector ii without updating
*      the disk
*
       use ipPage
       Implicit Integer (a-h,o-z)
#include "stdalloc.fh"
*
       If (ii.gt.Max_CI_Vectors) Then
          Write (6,*) 'opout: ii.gt.Max_CI_Vectors'
          Write (6,*) 'ii,Max_CI_Vectors=',ii,Max_CI_Vectors
          Call Abend()
       End If
*
       opout=0
       If (.not.diskbased) Return
*
       If (Status(ii).eq.In_Memory .and. ii.gt.0) Then
          nn=n(ii)
          Status(ii)=On_Disk
          Call mma_deallocate(W(ii)%Vec)
       Else
          opout=-1
       End If
*
       Return
       End
*                                                                      *
************************************************************************
*                                                                      *
       Integer Function ipout(ii)
*
*      ipout will page out vector ii to disk and free the memory area
*
       use ipPage
       Implicit Integer (a-h,o-z)
#include "stdalloc.fh"
*
       ipout=0
       If (.not.diskbased) Return
*
       If (Status(ii).eq.In_Memory .and. ii.gt.0) Then
          idisk=ida(ii)
          nn=n(ii)
          Call dDafile(Lu_ip,Write,W(ii)%Vec,nn,idisk)
          Status(ii)=On_Disk
          Call mma_deallocate(W(ii)%Vec)
       Else
          ipout=-1
       End if
*
       Return
       End
*                                                                      *
************************************************************************
*                                                                      *
       Subroutine ipterm()
       use ipPage
*
*      Termination
*
       use ipPage
*
       If (DiskBased) Call DaClos(Lu_ip)
*
       Return
       End
*                                                                      *
************************************************************************
*                                                                      *
       Subroutine ipinit()
       use ipPage
*
*      Initialization
*
       use ipPage
*
       DiskBased=.False.
*
       Return
       End
