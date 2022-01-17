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
Public :: NQ_data, Close_NQ_Data

#include "stdalloc.fh"

!Parameter(l_NQ=15+(iTabMx+1)+27)
!define declare_ip_l_max    ip_lMax(iNQ)     =ipNQ+(iNQ-1)*l_NQ+9
!define declare_ip_r_quad   ip_R_Quad(iNQ)   =ipNQ+(iNQ-1)*l_NQ+11
!define declare_ip_angular  ip_Angular(iNQ)  =ipNQ+(iNQ-1)*l_NQ+12
!define declare_ip_k_max    ip_k_Max(iNQ)    =ipNQ+(iNQ-1)*l_NQ+13
!define declare_ip_atom_nr  ip_Atom_Nr(iNQ)  =ipNQ+(iNQ-1)*l_NQ+14
!define declare_ip_r_low    ip_R_low(iNQ,l)  =ipNQ+(iNQ-1)*l_NQ+15+l
!define declare_ip_dodx     ip_dOdx(iNQ,i)   =ipNQ+(iNQ-1)*l_NQ+15+(iTabMx+1)+(i-1)*9

Type NQ_data_raw
  Sequence
  Real*8, Allocatable:: Coor(:)
  Real*8 :: A_High=0.0D0
  Real*8 :: A_Low =0.0D0
  Real*8 :: R_RS  =0.0D0
  Real*8 :: R_max =0.0D0
End Type NQ_data_raw

Integer :: nNQ=0
Type (NQ_data_raw), Allocatable:: NQ_data(:)

Contains

Subroutine Close_NQ_Data()
Integer iNQ, nNQ
! Cleanup and close
  nNQ = SIZE(NQ_data)
  Do iNQ = 1, nNQ
     Call mma_deallocate(NQ_data(iNQ)%Coor)
     NQ_data(iNQ)%A_High=0.0D0
     NQ_data(iNQ)%A_Low =0.0D0
     NQ_data(iNQ)%R_RS  =0.0D0
     NQ_data(iNQ)%R_max =0.0D0
   End Do
   Deallocate(NQ_Data)
End Subroutine Close_NQ_Data

End Module NQ_Structure

