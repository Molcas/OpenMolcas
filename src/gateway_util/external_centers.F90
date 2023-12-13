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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

module External_Centers

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: iXPolType = 0, nData_XF = 0, nDMS = 0, nEF = 0, nOrd_XF = 1, nOrdEF = -1, nRP = 0, nWel = 0, nXF = 0, &
                     nXMolnr = 0
real(kind=wp) :: Dxyz(3) = Zero
integer(kind=iwp), allocatable :: XEle(:), XMolnr(:,:)
real(kind=wp), allocatable :: AMP_Center(:), DMS_Centers(:,:), EF_Centers(:,:), OAM_Center(:), OMQ_Center(:), RP_Centers(:,:,:), &
                              Wel_Info(:,:), XF(:,:)

public :: AMP_Center, DMS_Centers, Dxyz, EF_Centers, External_Centers_Dmp, External_Centers_Free, External_Centers_Get, iXPolType, &
          nData_XF, nDMS, nEF, nOrd_XF, nOrdEF, nRP, nWel, nXF, nXMolnr, OAM_Center, OMQ_Center, RP_Centers, Wel_Info, XEle, XF, &
          XMolnr

!                                                                      *
!***********************************************************************
!                                                                      *
contains
!                                                                      *
!***********************************************************************
!                                                                      *
subroutine External_Centers_Dmp()

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp), allocatable :: iDmp(:)
  real(kind=wp), allocatable :: DMS_Ext(:,:), RP_Temp(:,:,:)

  if (allocated(EF_Centers)) then
    call Put_dArray('EF_Centers',EF_Centers,3*nEF)
  end if
  if (allocated(OAM_Center)) then
    call Put_dArray('OAM_Center',OAM_Center,3)
  end if
  if (allocated(OMQ_Center)) then
    call Put_dArray('OMQ_Center',OMQ_Center,3)
  end if
  if (allocated(DMS_Centers)) then
    call mma_allocate(DMS_Ext,3,nDMS+1,Label='DMS_Ext')
    DMS_ext(1:3,1:nDMS) = DMS_Centers(1:3,1:nDMS)
    DMS_ext(1:3,nDMS+1) = Dxyz(1:3)
    call Put_dArray('DMS_Centers',DMS_Ext,3*(nDMS+1))
    call mma_deallocate(DMS_Ext)
  end if
  if (allocated(Wel_Info)) then
    call Put_dArray('Wel_Info',Wel_Info,3*nWel)
  end if
  if (allocated(AMP_Center)) then
    call Put_dArray('AMP_Center',AMP_Center,3)
  end if
  if (allocated(RP_Centers)) then
    call mma_allocate(RP_Temp,3,nRP/3,2)
    RP_Temp(:,:,1) = RP_Centers(:,1:nRP/3,1)
    RP_Temp(:,:,2) = RP_Centers(:,1:nRP/3,2)
    call Put_dArray('RP_Centers',RP_Temp,nRP*2)
    call mma_deallocate(RP_Temp)
  end if
  if (allocated(XF)) then
    call Put_dArray('XF',XF,nData_XF*nXF)
  end if
  if (allocated(XMolnr)) then
    call Put_iArray('XMolnr',XMolnr,nXMolnr*nXF)
  end if
  if (allocated(XEle)) then
    call Put_iArray('XEle',XEle,nXF)
  end if
  call mma_Allocate(iDmp,3,Label='iDmp')
  iDmp(1) = nOrdEF
  iDmp(2) = nOrd_XF
  iDmp(3) = iXPolType
  call Put_iArray('Misc',iDmp,3)
  call mma_deallocate(iDmp)

  return

end subroutine External_Centers_Dmp
!                                                                      *
!***********************************************************************
!                                                                      *
subroutine External_Centers_Free()

  use stdalloc, only: mma_deallocate

  if (allocated(EF_Centers)) then
    call mma_deallocate(EF_Centers)
    nEF = 0
  end if
  if (allocated(OAM_Center)) call mma_deallocate(OAM_Center)
  if (allocated(OMQ_Center)) call mma_deallocate(OMQ_Center)
  if (allocated(DMS_Centers)) then
    call mma_deallocate(DMS_Centers)
    nDMS = 0
  end if
  if (allocated(Wel_Info)) then
    call mma_deallocate(Wel_Info)
    nWel = 0
  end if
  if (allocated(AMP_Center)) call mma_deallocate(AMP_Center)
  if (allocated(RP_Centers)) then
    call mma_deallocate(RP_Centers)
    nRP = 0
  end if
  if (allocated(XF)) then
    call mma_deallocate(XF)
    call mma_deallocate(XMolnr)
    call mma_deallocate(XEle)
    nData_XF = 0
    nXF = 0
    nXMolnr = 0
    nOrdEF = -1
    nOrd_XF = 1
    iXPolType = 0
  end if

  return

end subroutine External_Centers_Free
!                                                                      *
!***********************************************************************
!                                                                      *
subroutine External_Centers_Get()

  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: u6

  integer(kind=iwp) :: Len2
  logical(kind=iwp) :: Found
  integer(kind=iwp), allocatable :: iDmp(:)
  real(kind=wp), allocatable :: DMS_Ext(:,:)

  call qpg_dArray('EF_Centers',Found,Len2)
  if (Found) then
    nEF = Len2/3
    if (allocated(EF_Centers)) then
      if (size(EF_Centers,2) /= nEF) then
        write(u6,*) 'SIZE(EF_Centers,2) /= nEF'
        call Abend()
      end if
    else
      call mma_allocate(EF_Centers,3,nEF,Label='EF_Centers')
    end if
    call Get_dArray('EF_Centers',EF_Centers,3*nEF)
  end if

  call qpg_dArray('OAM_Center',Found,Len2)
  if (Found) then
    if (.not. allocated(OAM_Center)) then
      call mma_allocate(OAM_Center,3,Label='OAM_Center')
    end if
    call Get_dArray('OAM_Center',OAM_Center,3)
  end if

  call qpg_dArray('OMQ_Center',Found,Len2)
  if (Found) then
    if (.not. allocated(OMQ_Center)) then
      call mma_allocate(OMQ_Center,3,Label='OMQ_Center')
    end if
    call Get_dArray('OMQ_Center',OMQ_Center,3)
  end if

  call qpg_dArray('DMS_Centers',Found,Len2)
  if (Found) then
    nDMS = Len2/3-1
    if (allocated(DMS_Centers)) then
      if (size(DMS_Centers,2) /= nDMS) then
        write(u6,*) 'SIZE(DMS_Centers,2) /= nDMS'
        call Abend()
      end if
    else
      call mma_allocate(DMS_Centers,3,nDMS,Label='DMS_Centers')
    end if
    call mma_allocate(DMS_Ext,3,nDMS+1,Label='DMS_Ext')
    call Get_dArray('DMS_Centers',DMS_Ext,3*(nDMS+1))
    DMS_Centers(1:3,1:nDMS) = DMS_Ext(1:3,1:nDMS)
    Dxyz(1:3) = DMS_Centers(1:3,nDMS+1)
    call mma_deallocate(DMS_Ext)
  end if

  call qpg_dArray('Wel_Info',Found,Len2)
  if (Found) then
    nWel = Len2/3
    if (allocated(Wel_Info)) then
      if (size(Wel_Info,2) /= nWel) then
        write(u6,*) 'SIZE(Wel_Info,2) /= nWel'
        call Abend()
      end if
    else
      call mma_allocate(Wel_Info,3,nWel,Label='Wel_Info')
    end if
    call Get_dArray('Wel_Info',Wel_Info,3*nWel)
  end if

  call qpg_dArray('AMP_Center',Found,Len2)
  if (Found) then
    if (.not. allocated(AMP_Center)) then
      call mma_allocate(AMP_Center,3,Label='AMP_Center')
    end if
    call Get_dArray('AMP_Center',AMP_Center,3)
  end if

  call qpg_dArray('RP_Centers',Found,Len2)
  if (Found) then
    nRP = Len2/2
    if (allocated(RP_Centers)) then
      if (size(RP_Centers,2) /= nRP/3) then
        write(u6,*) 'SIZE(RP_Centers,2) /= nRP/3'
        call Abend()
      end if
    else
      call mma_allocate(RP_Centers,3,nRP/3,2,Label='RP_Centers')
    end if
    call Get_dArray('RP_Centers',RP_Centers,2*nRP)
  end if

  call qpg_iArray('XEle',Found,Len2)
  if (Found) then
    nXF = Len2
    if (.not. allocated(XEle)) then
      call mma_allocate(XEle,nXF,Label='XEle')
    end if
    call Get_iArray('XEle',XEle,nXF)

    call qpg_iArray('XMolnr',Found,Len2)
    nXMolnr = Len2/nXF
    if (.not. allocated(XMolnr)) then
      call mma_allocate(XMolnr,nXMolnr,nXF,Label='XMolnr')
    end if
    call Get_iArray('XMolnr',XMolnr,nXMolnr*nXF)

    call qpg_dArray('XF',Found,Len2)
    nData_XF = Len2/nXF
    if (.not. allocated(XF)) then
      call mma_allocate(XF,nData_XF,nXF,Label='XF')
    end if
    call Get_dArray('XF',XF,nData_XF*nXF)
  end if
  call mma_Allocate(iDmp,3,Label='iDmp')
  call Get_iArray('Misc',iDmp,3)
  nOrdEF = iDmp(1)
  nOrd_XF = iDmp(2)
  iXPolType = iDmp(3)
  call mma_deallocate(iDmp)

end subroutine External_Centers_Get
!                                                                      *
!***********************************************************************
!                                                                      *
end module External_Centers
