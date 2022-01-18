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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************
Module NQ_Structure
Implicit None
Private
Public :: NQ_data, Close_NQ_Data, Info_Ang, Close_Info_Ang, LMax_NQ

#include "stdalloc.fh"

!define declare_ip_dodx     ip_dOdx(iNQ,i)   =ipNQ+(iNQ-1)*l_NQ+15+(iTabMx+1)+(i-1)*9

Type NQ_data_raw
  Sequence
  Real*8, Allocatable:: Coor(:)
  Real*8 :: A_High=-1.0D99
  Real*8 :: A_Low = 1.0D99
  Real*8 :: R_RS  =0.0D0
  Real*8 :: R_max =0.0D0
  Integer :: l_max=-1
  Real*8, Allocatable :: R_Quad(:,:)
  Integer, Allocatable :: Angular(:)
  Integer :: Atom_Nr=-1
  Real*8, Allocatable :: dOdx(:,:,:)
End Type NQ_data_raw

Type (NQ_data_raw), Allocatable:: NQ_data(:)

Type Info_A
  Sequence
  Integer :: L_eff=0
  Integer :: nPoints=0
  Real*8, Allocatable:: R(:,:)
End Type Info_A

Integer, Parameter:: LMax_NQ=62
Type (Info_A) Info_Ang(LMax_NQ)

Contains

Subroutine Close_Info_Ang()

Integer iAngular
Do iAngular = 1, SIZE(Info_Ang)
   Info_Ang(iAngular)%L_eff=0
   Info_Ang(iAngular)%nPoints=0
   If (Allocated(Info_Ang(iAngular)%R)) Call mma_deallocate(Info_Ang(iAngular)%R)
End Do
End Subroutine Close_Info_Ang

Subroutine Close_NQ_Data()
Integer iNQ, nNQ
! Cleanup and close
  nNQ = SIZE(NQ_data)
  Do iNQ = 1, nNQ
     Call mma_deallocate(NQ_data(iNQ)%Coor)
     NQ_data(iNQ)%A_High=-1.0D99
     NQ_data(iNQ)%A_Low = 1.0D99
     NQ_data(iNQ)%R_RS  =0.0D0
     NQ_data(iNQ)%R_max =0.0D0
     NQ_data(iNQ)%l_Max =-1
     If (Allocated(NQ_data(iNQ)%R_Quad))Call mma_deallocate(NQ_data(iNQ)%R_Quad)
     If (Allocated(NQ_data(iNQ)%Angular))Call mma_deallocate(NQ_data(iNQ)%Angular)
     NQ_Data(iNQ)%Atom_Nr=-1
     If (Allocated(NQ_data(iNQ)%dOdx))Call mma_deallocate(NQ_data(iNQ)%dOdx)
   End Do
   Deallocate(NQ_Data)
End Subroutine Close_NQ_Data

End Module NQ_Structure

