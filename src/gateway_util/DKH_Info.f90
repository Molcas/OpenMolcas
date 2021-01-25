!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
Module DKH_Info
Private
Public :: nCtrLD, iCtrLD, radiLD, DKroll, LDKroll, BSS, &
          DKH_Info_Get, DKH_Info_Dmp
#include "stdalloc.fh"
Integer :: i
Integer ::  nCtrLD=0, iCtrLD(10)=[(0,i=1,10)]
Real*8  :: radiLD=0.0D0
Logical :: DKroll=.False.
Logical :: LDKroll=.False.
Logical :: BSS   =.False.
Integer :: iRelae_Info=0

Contains

Subroutine DKH_Info_Dmp()
#include "relae.fh"
  Real*8, Allocatable:: rDmp(:)
  Integer:: Len=1+10+5
  Integer i

  iRelae_Info=iRelae

  Call mma_allocate(rDmp,Len,Label='rDmp:DKH')
  rDmp(1)=DBLE(nCtrLD)
  Do i = 1, 10
    rDMP(i+1)=iCtrLD(i)
  End Do
  rDMP(12)=radiLD
  rDMP(13)=DBLE(0)
  If (DKroll) rDMP(13)=DBLE(1)
  rDMP(14)=DBLE(0)
  If (LDKroll) rDMP(14)=DBLE(1)
  rDMP(15)=DBLE(0)
  If (BSS) rDMP(15)=DBLE(1)
  rDmp(16)=DBLE(iRelae_Info)
  Call Put_dArray('DKH_Info',rDmp,Len)
  Call mma_deallocate(rDmp)
End Subroutine DKH_Info_Dmp


Subroutine DKH_Info_Get()
#include "relae.fh"
  Real*8, Allocatable:: rDmp(:)
  Integer:: Len=1+10+5
  Integer i

  Call mma_allocate(rDmp,Len,Label='rDmp:DKH')
  Call Get_dArray('DKH_Info',rDmp,Len)

  nCtrLD=NINT(rDmp(1))
  Do i = 1, 10
    iCtrLD(i)=NINT(rDMP(i+1))
  End Do
  radiLD=rDMP(12)
  DKroll  = NINT(rDMP(13)).eq.1
  LDKroll = NINT(rDMP(14)).eq.1
  BSS     = NINT(rDMP(15)).eq.1
  iRelae_Info = NINT(rDmp(16))

  iRelae=iRelae_Info

  Call mma_deallocate(rDmp)
End Subroutine DKH_Info_Get

End Module DKH_Info
