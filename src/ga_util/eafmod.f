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
* Copyright (C) Martin Schuetz                                         *
*               Roland Lindh                                           *
*               2015, Steven Vancoillie                                *
*               2017, Ignacio Fdez. Galvan                             *
************************************************************************

************************************************************************
* Wrapper routines agains the GA-ChemIO or MPI-IO package.             *
* (EAF = exclusive access file)                                        *
*                                                                      *
* If not an mpp installation the AIX-IO facility will be used.         *
************************************************************************

#ifdef ADD_
#define eaf_open eaf_open_
#define eaf_close eaf_close_
#define eaf_awrite eaf_awrite_
#define eaf_aread eaf_aread_
#define eaf_write eaf_write_
#define eaf_read eaf_read_
#define eaf_wait eaf_wait_
#endif

************************************************************************
*                                                                      *
      Subroutine EAFOpen(Lu,FName)
      Implicit None
      Character*(*) FName
      Integer Lu
      Integer  isfreeunit
      External isfreeunit
#ifdef _MOLCAS_MPP_
      Character*200 FN
      Integer  iRC
      Logical  Is_Real_Par
      External Is_Real_Par
#ifdef _HAVE_EXTRA_
      Integer  n
#endif
#include "molcas_eaf.fh"
#endif
#ifdef _MOLCAS_MPP_
      If(Is_Real_Par()) Then
      FN=FName
      Call Ext_Rank(FN)
#ifdef _HAVE_EXTRA_
      n=Len_Trim(FName)
      FN(n+1:n+1)=Char(0)
      iRC=molcas_eaf_open(FN,Lu)
#else
      iRC=eaf_open(FN,eaf_rw,Lu)
#endif
      If (iRC.ne.0) Then
         Write (6,*) 'EAFOpen: Abort!'
         Write (6,*) 'iRC=',iRC
         Write (6,*) 'Lu=',Lu
         Write (6,*) 'Fname=',FName
         Write (6,*) 'EAF_Err_Code =',iRC
         Call Abend()
      End If
      Else
#endif
      Lu=isfreeunit(7)
      Call DaName_MF(Lu,FName)
#ifdef _MOLCAS_MPP_
      End If
#endif
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine EAFClose(Lu)
      Implicit None
      Integer Lu
#ifdef _MOLCAS_MPP_
      Integer  iRC
      Logical  Is_Real_Par
      External Is_Real_Par
#include "molcas_eaf.fh"
#endif
#ifdef _MOLCAS_MPP_
      If(Is_Real_Par()) Then
#ifdef _HAVE_EXTRA_
      iRC=molcas_eaf_close(Lu)
#else
      iRC=eaf_close(Lu)
#endif
      If (iRC.ne.0) Then
         Write (6,*) 'EAFClose: Abort!'
         Write (6,*) 'EAF_Err_Code =',iRC
         Call Abend()
      End If
      Else
#endif
      Call DaClos(Lu)
#ifdef _MOLCAS_MPP_
      End If
#endif
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine EAFAWrite(Lu,Buf,nBuf,Disk,id)
      Implicit None
      Integer Lu, nBuf, Buf(nBuf), id
      Real*8 Disk
#ifdef _MOLCAS_MPP_
      Integer  iRC
      Logical  Is_Real_Par
      External Is_Real_Par
#include "molcas_eaf.fh"
#include "SysDef.fh"
#endif
#ifdef _MOLCAS_MPP_
      If(Is_Real_Par()) Then
#ifdef _HAVE_EXTRA_
      iRC=molcas_eaf_awrite(Lu,Disk,Buf,nBuf*ItoB,id)
#else
      iRC=eaf_awrite(Lu,Disk,Buf,nBuf*ItoB,id)
      Disk=Disk+Dble(nBuf*ItoB)
#endif
      If (iRC.ne.0) Then
         Write (6,*) 'EAFAWrite: Abort!'
         Write (6,*) 'EAF_Err_Code =',iRC
         Call Abend()
      End If
      Else
#endif
      id=0
      Call EAFWrite(Lu,Buf,nBuf,Disk)
#ifdef _MOLCAS_MPP_
      End If
#endif
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine EAFARead(Lu,Buf,nBuf,Disk,id)
      Implicit None
      Integer Lu, nBuf, Buf(nBuf), id
      Real*8 Disk
#ifdef _MOLCAS_MPP_
      Integer  iRC
      Logical  Is_Real_Par
      External Is_Real_Par
#include "molcas_eaf.fh"
#include "SysDef.fh"
#endif
#ifdef _MOLCAS_MPP_
      If(Is_Real_Par()) Then
#ifdef _HAVE_EXTRA_
      iRC=molcas_eaf_aread(Lu,Disk,Buf,nBuf*ItoB,id)
#else
      iRC=eaf_aread(Lu,Disk,Buf,nBuf*ItoB,id)
      Disk=Disk+Dble(nBuf*ItoB)
#endif
      If (iRC.ne.0) Then
         Write (6,*) 'EAFARead: Abort!'
         Write (6,*) 'EAF_Err_Code =',iRC
         Call Abend()
      End If
      Else
#endif
      id=0
      Call EAFRead(Lu,Buf,nBuf,Disk)
#ifdef _MOLCAS_MPP_
      End If
#endif
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine EAFWrite(Lu,Buf,nBuf,Disk)
      Implicit None
      Integer Lu, nBuf, Buf(nBuf)
      Real*8 Disk
      Integer iDisk
#ifdef _MOLCAS_MPP_
      Integer  iRC
      Logical  Is_Real_Par
      External Is_Real_Par
#include "molcas_eaf.fh"
#include "SysDef.fh"
#endif
#ifdef _MOLCAS_MPP_
      If(Is_Real_Par()) Then
#ifdef _HAVE_EXTRA_
      iRC=molcas_eaf_write(Lu,Disk,Buf,nBuf*ItoB)
#else
      iRC=eaf_write(Lu,Disk,Buf,nBuf*ItoB)
      Disk=Disk+Dble(nBuf*ItoB)
#endif
      If (iRC.ne.0) Then
         Write (6,*) 'EAFWrite: Abort!'
         Write (6,*) 'EAF_Err_Code =',iRC
         Call Abend()
      End If
      Else
#endif
      iDisk=Int(Disk)
      Call iDaFile(Lu,1,Buf,nBuf,iDisk)
      Disk=Dble(iDisk)
#ifdef _MOLCAS_MPP_
      End If
#endif
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine EAFRead(Lu,Buf,nBuf,Disk)
      Implicit None
      Integer Lu, nBuf, Buf(nBuf)
      Real*8 Disk
      Integer iDisk
#ifdef _MOLCAS_MPP_
      Integer  iRC
      Logical  Is_Real_Par
      External Is_Real_Par
#include "molcas_eaf.fh"
#include "SysDef.fh"
#endif
#ifdef _MOLCAS_MPP_
      If(Is_Real_Par()) Then
#ifdef _HAVE_EXTRA_
      iRC=molcas_eaf_read(Lu,Disk,Buf,nBuf*ItoB)
#else
      iRC=eaf_read(Lu,Disk,Buf,nBuf*ItoB)
      Disk=Disk+Dble(nBuf*ItoB)
#endif
      If (iRC.ne.0) Then
         Write (6,*) 'EAFRead: Abort!'
         Write (6,*) 'EAF_Err_Code =',iRC
         Call Abend()
      End If
      Else
#endif
      iDisk=Int(Disk)
      Call iDaFile(Lu,2,Buf,nBuf,iDisk)
      Disk=Dble(iDisk)
#ifdef _MOLCAS_MPP_
      End If
#endif
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine EAFWait(Lu,id)
      Implicit None
      Integer Lu, id
#ifdef _MOLCAS_MPP_
      Integer iRC
      Logical  Is_Real_Par
      External Is_Real_Par
#include "molcas_eaf.fh"
#endif
#ifdef _MOLCAS_MPP_
      If(Is_Real_Par()) Then
#ifdef _HAVE_EXTRA_
      iRC=molcas_eaf_wait(Lu,id)
#else
      iRC=eaf_wait(Lu,id)
#endif
      If (iRC.ne.0) Then
         Write (6,*) 'EAFWait: Abort!'
         Write (6,*) 'EAF_Err_Code =',iRC
         Call Abend()
      End If
      End If
#endif
      Return
#ifndef _MOLCAS_MPP_
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(Lu)
         Call Unused_integer(id)
      End If
#endif
      End
*                                                                      *
************************************************************************
