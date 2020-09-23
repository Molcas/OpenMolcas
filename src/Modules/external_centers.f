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
* Copyright (C) 2020, Roland Lindh                                     *
************************************************************************
      Module External_Centers
      Implicit None
      Private
#include "stdalloc.fh"
      Public :: nEF, EF_Centers,
     &          OAM_Center, OMQ_Center, nDMS, DMS_Centers, Dxyz,
     &          nWel, Wel_Info, AMP_Center, nRP, RP_Centers,
     &          nData_XF, nXF, nXMolnr, XF, XEle, XMolnr,
     &          nOrdEF, nOrd_XF, iXPolType,
     &          External_Centers_Dmp,
     &          External_Centers_Free,
     &          External_Centers_Get
      Integer :: nEF=0, nOrdEF=-1
      Real*8, Allocatable:: EF_Centers(:,:)
      Real*8, Allocatable:: OAM_Center(:)
      Real*8, Allocatable:: OMQ_Center(:)
      Integer :: nDMS=0
      Real*8, Allocatable:: DMS_Centers(:,:)
      Real*8 :: Dxyz(3)=[0.0D0, 0.0D0, 0.0D0]
      Integer :: nWel=0
      Real*8, Allocatable:: Wel_Info(:,:)
      Real*8, Allocatable:: AMP_Center(:)
      Integer :: nRP=0
      Real*8, Target, Allocatable:: RP_Centers(:,:,:)
      Integer :: nData_XF=0, nXF=0, nXMolnr=0, nOrd_XF=1, iXPolType=0
      Real*8, Allocatable:: XF(:,:)
      Integer, Allocatable:: XEle(:), XMolnr(:,:)
*                                                                      *
************************************************************************
*                                                                      *
      Contains
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine External_Centers_Dmp()
#include "angtp.fh"
#include "info.fh"
      Real*8, Allocatable:: RP_Temp(:,:,:)
      Integer, Allocatable:: iDmp(:)
      Real*8, Allocatable:: DMS_Ext(:,:)
      If (Allocated(EF_Centers)) Then
         Call Put_dArray('EF_Centers',EF_Centers,3*nEF)
      End If
      If (Allocated(OAM_Center)) Then
         Call Put_dArray('OAM_Center',OAM_Center,3)
      End If
      If (Allocated(OMQ_Center)) Then
         Call Put_dArray('OMQ_Center',OMQ_Center,3)
      End If
      If (Allocated(DMS_Centers)) Then
         Call mma_allocate(DMS_Ext,3,nDMS+1,Label='DMS_Ext')
         DMS_ext(1:3,1:nDMS)=DMS_Centers(1:3,1:nDMS)
         DMS_ext(1:3,nDMS+1)=Dxyz(1:3)
         Call Put_dArray('DMS_Centers',DMS_Ext,3*(nDMS+1))
         Call mma_deallocate(DMS_Ext)
      End If
      If (Allocated(Wel_Info)) Then
         Call Put_dArray('Wel_Info',Wel_Info,3*nWel)
      End If
      If (Allocated(AMP_Center)) Then
         Call Put_dArray('AMP_Center',AMP_Center,3)
      End If
      If (Allocated(RP_Centers)) Then
         Call mma_allocate(RP_Temp,3,nRP/3,2)
         RP_Temp(:,:,1)=RP_Centers(:,1:nRP/3,1)
         RP_Temp(:,:,2)=RP_Centers(:,1:nRP/3,2)
         Call Put_dArray('RP_Centers',RP_Temp,nRP*2)
         Call mma_deallocate(RP_Temp)
      End If
      If (Allocated(XF)) Then
         Call Put_dArray('XF',XF,nData_XF*nXF)
      End If
      If (Allocated(XMolnr)) Then
         Call Put_iArray('XMolnr',XMolnr,nXMolnr*nXF)
      End If
      If (Allocated(XEle)) Then
         Call Put_iArray('XEle',XEle,nXF)
      End If
      Call mma_Allocate(iDmp,3,Label='iDmp')
      iDmp(1)=nOrdEF
      iDmp(2)=nOrd_XF
      iDmp(3)=iXPolType
      Call Put_iArray('Misc',iDmp,3)
      Call mma_deallocate(iDmp)
      Return
      End Subroutine External_Centers_Dmp
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine External_Centers_Free()
      If (Allocated(EF_Centers)) Then
         Call mma_deallocate(EF_Centers)
         nEF=0
      End If
      If (Allocated(OAM_Center)) Call mma_deallocate(OAM_Center)
      If (Allocated(OMQ_Center)) Call mma_deallocate(OMQ_Center)
      If (Allocated(DMS_Centers)) Then
         Call mma_deallocate(DMS_Centers)
         nDMS=0
      End If
      If (Allocated(Wel_Info)) Then
         Call mma_deallocate(Wel_Info)
         nWel=0
      End If
      If (Allocated(AMP_Center)) Call mma_deallocate(AMP_Center)
      If (Allocated(RP_Centers)) Then
         Call mma_deallocate(RP_Centers)
         nRP=0
      End If
      If (Allocated(XF)) Then
         Call mma_deallocate(XF)
         Call mma_deallocate(XMolnr)
         Call mma_deallocate(XEle)
         nData_XF=0
         nXF=0
         nXMolnr=0
         nOrdEF=-1
         nOrd_XF=1
         iXPolType=0
      End If
      Return
      End Subroutine External_Centers_Free
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine External_Centers_Get()
      Integer, Allocatable:: iDmp(:)
      Real*8, Allocatable:: DMS_Ext(:,:)
      Logical Found
      Integer Len2
      Call qpg_dArray('EF_Centers',Found,Len2)
      If (Found) Then
         nEF=Len2/3
         If (Allocated(EF_Centers)) Then
            If (SIZE(EF_Centers,2).ne.nEF) Then
               Write (6,*) 'SIZE(EF_Centers,2).ne.nEF'
               Call Abend()
            End If
         Else
            Call mma_allocate(EF_Centers,3,nEF,Label='EF_Centers')
         End If
         Call Get_dArray('EF_Centers',EF_Centers,3*nEF)
      End If
*
      Call qpg_dArray('OAM_Center',Found,Len2)
      If (Found) Then
         If (.Not.Allocated(OAM_Center)) Then
            Call mma_allocate(OAM_Center,3,Label='OAM_Center')
         End If
         Call Get_dArray('OAM_Center',OAM_Center,3)
      End If
*
      Call qpg_dArray('OMQ_Center',Found,Len2)
      If (Found) Then
         If (.Not.Allocated(OMQ_Center)) Then
            Call mma_allocate(OMQ_Center,3,Label='OMQ_Center')
         End If
         Call Get_dArray('OMQ_Center',OMQ_Center,3)
      End If
*
      Call qpg_dArray('DMS_Centers',Found,Len2)
      If (Found) Then
         nDMS=Len2/3-1
         If (Allocated(DMS_Centers)) Then
            If (SIZE(DMS_Centers,2).ne.nDMS) Then
               Write (6,*) 'SIZE(DMS_Centers,2).ne.nDMS'
               Call Abend()
            End If
         Else
            Call mma_allocate(DMS_Centers,3,nDMS,Label='DMS_Centers')
         End If
         call mma_allocate(DMS_Ext,3,nDMS+1,Label='DMS_Ext')
         Call Get_dArray('DMS_Centers',DMS_Ext,3*(nDMS+1))
         DMS_Centers(1:3,1:nDMS)= DMS_Ext(1:3,1:nDMS)
         Dxyz(1:3) =              DMS_Centers(1:3,nDMS+1)
         call mma_deallocate(DMS_Ext)
      End If
*
      Call qpg_dArray('Wel_Info',Found,Len2)
      If (Found) Then
         nWel=Len2/3
         If (Allocated(Wel_Info)) Then
            If (SIZE(Wel_Info,2).ne.nWel) Then
               Write (6,*) 'SIZE(Wel_Info,2).ne.nWel'
               Call Abend()
            End If
         Else
            Call mma_allocate(Wel_Info,3,nWel,Label='Wel_Info')
         End If
         Call Get_dArray('Wel_Info',Wel_Info,3*nWel)
      End If
*
      Call qpg_dArray('AMP_Center',Found,Len2)
      If (Found) Then
         If (.Not.Allocated(AMP_Center)) Then
            Call mma_allocate(AMP_Center,3,Label='AMP_Center')
         End If
         Call Get_dArray('AMP_Center',AMP_Center,3)
      End If
*
      Call qpg_dArray('RP_Centers',Found,Len2)
      If (Found) Then
         nRP=Len2/2
         If (Allocated(RP_Centers)) Then
            If (SIZE(RP_Centers,2).ne.nRP/3) Then
               Write (6,*) 'SIZE(RP_Centers,2).ne.nRP/3'
               Call Abend()
            End If
         Else
            Call mma_allocate(RP_Centers,3,nRP/3,2,Label='RP_Centers')
         End If
         Call Get_dArray('RP_Centers',RP_Centers,2*nRP)
      End If
*
      Call qpg_iArray('XEle',Found,Len2)
      If (Found) Then
         nXF=Len2
         If (.Not.Allocated(XEle)) Then
            Call mma_allocate(XEle,nXF,Label='XEle')
         End If
         Call Get_iArray('XEle',XEle,nXF)
*
         Call qpg_dArray('XMolnr',Found,Len2)
         nXMolnr=Len2/nXF
         If (.Not.Allocated(XMolnr)) Then
            Call mma_allocate(XMolnr,nXMolnr,nXF,Label='XMolnr')
         End If
         Call Get_iArray('XMolnr',XMolnr,nXMolnr*nXF)
*
         Call qpg_dArray('XF',Found,Len2)
         nData_XF=Len2/nXF
         If (.Not.Allocated(XF)) Then
            Call mma_allocate(XF,nData_XF,nXF,Label='XF')
         End If
         Call Get_dArray('XF',XF,nData_XF*nXF)
      End If
      Call mma_Allocate(iDmp,3,Label='iDmp')
      Call Get_iArray('Misc',iDmp,3)
      nOrdEF    = iDmp(1)
      nOrd_XF   = iDmp(2)
      iXPolType = iDmp(3)
      Call mma_deallocate(iDmp)
      End Subroutine External_Centers_Get
*                                                                      *
************************************************************************
*                                                                      *
      End Module External_Centers
