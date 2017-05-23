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
*                                                                      *
************************************************************************
*                                                                      *
       Logical Function ipopen(nconf,page)
*
*      Initiate the whole lot.
*
       Implicit Real*8(a-h,o-z)
       Logical page
#include "ippage.fh"
#include "WrkSpc.fh"
*
*      Ask how much memory is available
*
       Call GetMem('ipopen','MAX','REAL',ipDum,nmax)
       nmax=nmax/2
       npp=nconf*10
*      If (npp.lt.nmax.or.(.not.page)) Then
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
*         ip_Mem : memory pointer
*         n  : Length of CI-vector
*         ida: disk address
*
          Call ICopy(Max_CI_Vectors+1, 0,0,n,1)
          Call ICopy(Max_CI_Vectors+1,-1,0,ida,1)
          Call ICopy(Max_CI_Vectors+1,ip_Dummy,0,ip_Mem,1)
          Call ICopy(Max_CI_Vectors+1,Null_Vector,0,Status,1)
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
#include "ippage.fh"
#include "WrkSpc.fh"
       Real*8 rdum
*
       If (ia.gt.Max_CI_Vectors) Then
          Write (6,*) 'ipclose: ia.gt.Max_CI_Vectors'
          Write (6,*) 'ia,Max_CI_Vectors=',ia,Max_CI_Vectors
          Call QTrace()
          Call Abend()
       End If
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
       Do ii=Max(ia,0),Max_CI_Vectors
          If (Status(ii).eq.In_Memory) Then
             Call Getmem('ipclose','FREE','REAL',ip_Mem(ii),idum)
             ip_Mem(ii)=ip_Dummy
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
*      Get the index of a vector with the length nn
*      memory or disk space is allocated.
*
       Implicit Integer (a-h,o-z)
#include "ippage.fh"
#include "WrkSpc.fh"
       Character*4 Label
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
*----- Allocate memory for vector
*
       If (nn.gt.0) Then
          Write (Label,'(I3.3)') n_CI_Vectors
          Call GetMem('ipget'//Label,'ALLO','REAL',ipT,nn)
          Status(ipget)=In_Memory
          Call FZero(Work(ipT),nn)
       Else
          ipT=-1
          Status(ipget)=Null_Vector
       End If
*
       If (DiskBased) then
          If (Status(ipget).ne.Null_Vector) Then
             Call dDafile(Lu_ip,Write,Work(ipT),nn,iDisk_Addr_End)
             Status(ipget)=On_Disk
             Call GetMem('ipget','FREE','REAL',ipT,nn)
             ip_Mem(ipget)=ip_Dummy
          End If
       Else
         ip_Mem(ipget)=ipT
       End If
*
       Return
       End
*                                                                      *
************************************************************************
*                                                                      *
       Integer Function ipin(ii)
*
*      Object: return pointer to vector ii
*
       Implicit Integer (a-h,o-z)
#include "ippage.fh"
#include "WrkSpc.fh"
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
       Implicit Integer (a-h,o-z)
#include "ippage.fh"
#include "WrkSpc.fh"
*
       If (ii.gt.Max_CI_Vectors) Then
          Write (6,*) 'ipin1: ii.gt.Max_CI_Vectors'
          Write (6,*) 'ii,Max_CI_Vectors=',ii,Max_CI_Vectors
          Call Abend()
       End If
*
       If (Status(ii).eq.In_Memory) Then
*
*-------- ii is in memory
*
          ip1=ip_Mem(ii)
*
          If (nn.gt.n(ii)) Then
             Call GetMem('ipin1','ALLO','REAL',ip1,nn)
             Call FZero(Work(ip1),nn)
             call dcopy_(n(ii),Work(ip_Mem(ii)),1,Work(ip1),1)
             Call GetMem('ipin1','FREE','REAL',ip_Mem(ii),n(ii))
             n(ii)=nn
             ip_Mem(ii)=ip1
          End If
*
       Else If (Status(ii).eq.On_Disk) Then
*
*         ii is on disk
*
          Call Getmem('ipin1','ALLO','REAL',ip1,Max(n(ii),nn))
          Call FZero(Work(ip1),Max(n(ii),nn))
*
          ip_Mem(ii)=ip1
          nnn=Min(n(ii),nn)

*
*         pick up from disk
*
          idisk=ida(ii)
          Call dDafile(Lu_ip,Read,Work(ip1),nnn,idisk)
          Status(ii)=In_Memory
*
       Else If (Status(ii).eq.Null_Vector) Then
*
           ip1=ip_Dummy ! dummy pointer
*
       Else
*
          ip1=ip_Dummy
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
       Implicit Integer (a-h,o-z)
#include "ippage.fh"
#include "WrkSpc.fh"
*
       If (iii.gt.Max_CI_Vectors) Then
          Write (6,*) 'ipout: iii.gt.Max_CI_Vectors'
          Write (6,*) 'iii,Max_CI_Vectors=',iii,Max_CI_Vectors
          Call QTrace()
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
             ip1=ip_Mem(ii)
             nn=n(ii)
             Call dDafile(Lu_ip,Write,Work(ip1),nn,idisk)
             Status(ii)=On_Disk
             Call Getmem('ipnout','FREE','REAL',ip1,nn)
             ip_Mem(ii)=ip_Dummy
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
*      opout will release the memory area without update the disk
*
       Implicit Integer (a-h,o-z)
#include "ippage.fh"
#include "WrkSpc.fh"
*
       If (ii.gt.Max_CI_Vectors) Then
          Write (6,*) 'opout: ii.gt.Max_CI_Vectors'
          Write (6,*) 'ii,Max_CI_Vectors=',ii,Max_CI_Vectors
          Call QTrace()
          Call Abend()
       End If
*
       opout=0
       If (.not.diskbased) Return
*
       If (Status(ii).eq.In_Memory .and. ii.gt.0) Then
          ip1=ip_Mem(ii)
          nn=n(ii)
          Status(ii)=On_Disk
          Call Getmem('opout','FREE','REAL',ip1,nn)
          ip_Mem(ii)=ip_Dummy
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
*      ipout will page them out to disk and free the memory area
*
       Implicit Integer (a-h,o-z)
#include "ippage.fh"
#include "WrkSpc.fh"
*
       ipout=0
       If (.not.diskbased) Return
*
       If (Status(ii).eq.In_Memory .and. ii.gt.0) Then
          ip1=ip_Mem(ii)
          idisk=ida(ii)
          nn=n(ii)
          Call dDafile(Lu_ip,Write,Work(ip1),nn,idisk)
          Status(ii)=On_Disk
          Call Getmem('ipout','FREE','REAL',ip1,nn)
          ip_Mem(ii)=ip_Dummy
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
*
*      Termination
*
       Implicit Real*8(a-h,o-z)
#include "ippage.fh"
*
       If (DiskBased) Call DaClos(Lu_ip)
*
       Return
       End
*                                                                      *
************************************************************************
*                                                                      *
       Subroutine ipinit()
*
*      Initialization
*
       Implicit Real*8(a-h,o-z)
#include "ippage.fh"
*
       DiskBased=.False.
*
       Return
       End
