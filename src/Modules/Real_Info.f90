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
Module Real_Info
Private
Public :: AccMch, ThrInt, PotNuc, Rtrnc, CutInt, TMass, qNuc, PkAcc, &
          Thrs, RadMax, cdMax, EtMax, E1, E2, RPQMin, SadStep, Shake,&
          ChiI2, &
          Real_Info_Dmp, Real_Info_Get

#include "stdalloc.fh"
Real*8 ::AccMch  =1.d-15
Real*8 ::ThrInt  =1.d-14
Real*8 ::PotNuc  =0.0D0
Real*8 ::Rtrnc   =3.D0
Real*8 ::CutInt  =1.d-16
Real*8 ::TMass   =0.0D0
Real*8 ::qNuc    =0.0D0
Real*8 ::PkAcc   =1.d-14
Real*8 ::Thrs    =1.d-6
Real*8 ::RadMax  = 0.0D0
Real*8 ::cdMax   = 0.0D0
Real*8 ::EtMax   = 0.0D0
Real*8 ::E1      =0.0D0
Real*8 ::E2      =0.0D0
Real*8 ::RPQMin  = 0.4d0
Real*8 ::SadStep =0.1d0
Real*8 ::CLight_Info=0.0D0
Real*8 ::Shake   =-1.0D0
Real*8 ::ChiI2   =0.0D0

Contains

Subroutine Real_Info_Dmp()
#include "RelLight.fh"
  Real*8, Allocatable:: rDmp(:)
  Integer:: Len=19

  CLight_Info=CLightAU
  Call mma_allocate(rDmp,Len,Label='rDmp:Real')

  rDmp(01)=AccMch
  rDmp(02)=ThrInt
  rDmp(03)=PotNuc
  rDmp(04)=Rtrnc
  rDmp(05)=CutInt
  rDmp(06)=TMass
  rDmp(07)=qNuc
  rDmp(08)=PkAcc
  rDmp(09)=Thrs
  rDmp(10)=RadMax
  rDmp(11)=cdMax
  rDmp(12)=EtMax
  rDmp(13)=E1
  rDmp(14)=E2
  rDmp(15)=RPQMin
  rDmp(16)=SadStep
  rDmp(17)=CLight_Info
  rDmp(18)=Shake
  rDmp(19)=ChiI2
  Call Put_dArray('Real_Info',rDmp,Len)
  Call mma_deallocate(rDmp)
End Subroutine Real_Info_Dmp

Subroutine Real_Info_Get()
#include "RelLight.fh"
  Real*8, Allocatable:: rDmp(:)
  Integer:: Len=19

  Call mma_allocate(rDmp,Len,Label='rDmp:Real')
  Call Get_dArray('Real_Info',rDmp,Len)

  AccMch     = rDmp(01)
  ThrInt     = rDmp(02)
  PotNuc     = rDmp(03)
  Rtrnc      = rDmp(04)
  CutInt     = rDmp(05)
  TMass      = rDmp(06)
  qNuc       = rDmp(07)
  PkAcc      = rDmp(08)
  Thrs       = rDmp(09)
  RadMax     = rDmp(10)
  cdMax      = rDmp(11)
  EtMax      = rDmp(12)
  E1         = rDmp(13)
  E2         = rDmp(14)
  RPQMin     = rDmp(15)
  SadStep    = rDmp(16)
  CLight_Info= rDmp(17)
  Shake      = rDmp(18)
  ChiI2      = rDmp(19)

  Call mma_deallocate(rDmp)

  CLightAU = CLight_Info

End Subroutine Real_Info_Get

End Module Real_Info
