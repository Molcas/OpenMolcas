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
! Copyright (C) 1992,2020, Roland Lindh                                *
!***********************************************************************
      SubRoutine GetInf(DoRys,nDiff)
!***********************************************************************
!                                                                      *
! Object: to read all input information on the file INFO.              *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January 1992                                             *
!***********************************************************************
      Use Iso_C_Binding
      use Real_Spherical, only: lMax_Internal, Sphere
      use Her_RW, only: nPrp
      use External_Centers, only: nOrdEF
      use Gateway_global, only: Test
      use DKH_Info, only: DKroll
      use Sizes_of_Seward, only: S
      use rctfld_module, only: lMax, cRFStrt,iRFStrt,lRFStrt,rRFStrt,
     &                               cRFEnd ,iRFEnd ,lRFEnd ,rRFEnd
      Implicit None
      Integer nDiff
      Logical DoRys
#include "SysDef.fh"

      Integer Len
      Integer, external:: ip_of_work, ip_of_iWork
!
      Call GetInf_Internal(cRFStrt,iRFStrt,lRFStrt,rRFStrt)
!
!     This is to allow type punning without an explicit interface
      Contains
      SubRoutine GetInf_Internal(cRFStrt,iRFStrt,lRFStrt,rRFStrt)
      Integer, Target :: cRFStrt,iRFStrt,lRFStrt
      Real*8, Target :: rRFStrt
      Integer, Pointer :: p_cRF(:),p_iRF(:),p_lRF(:)
      Real*8, Pointer :: p_rRF(:)
!
!     Load the dynamic input area.
!
      Call Get_Info_Dynamic()
!
!     Load the static input area.
!
      Call Get_Info_Static()
!                                                                      *
!***********************************************************************
!                                                                      *
!     Reaction field parameters
!
      Len = ip_of_iWork(lRFEnd)-ip_of_iWork(lRFStrt)+1
      Call C_F_Pointer(C_Loc(lRFStrt),p_lRF,[Len])
      Call Get_iArray('RFlInfo',p_lRF,Len)
!
      Len = ip_of_Work(rRFEnd)-ip_of_Work(rRFStrt)+1
      Call C_F_Pointer(C_Loc(rRFStrt),p_rRF,[Len])
      Call Get_dArray('RFrInfo',p_rRF,Len)
!
      Len = ip_of_iWork(iRFEnd)-ip_of_iWork(iRFStrt)+1
      Call C_F_Pointer(C_Loc(iRFStrt),p_iRF,[Len])
      Call Get_iArray('RFiInfo',p_iRF,Len)
!
      Len = ip_of_iWork(cRFEnd)-ip_of_iWork(cRFStrt)+1
      Call C_F_Pointer(C_Loc(cRFStrt),p_cRF,[Len])
      Call Get_iArray('RFcInfo',p_cRF,Len)
!
      Nullify(p_lRF,p_rRF,p_iRF,p_cRF)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Generate the transformation matrices
!
      If (S%iAngMx-1.ge.lMax) Then
         Call Sphere(S%iAngMx)
         lmax_internal=S%iAngMx
      Else
         Call Sphere(lMax)
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     Set the highest number of differentiations which will be
!     applied to the basis functions. In this case 2 + 1 ( the
!     kinetic energy operator and a differentiaion with respect
!     to a nuclear coordinate.
!
      nPrp=Max(lMax,3)
!
!     Setup of tables for coefficients of the Rys roots and weights.
!
      If (S%iAngMx.eq.0) nDiff=2
      If (DKroll.and.nOrdEF.gt.0) nDiff=nDiff+nOrdEF
      If (.Not.Test) Call Setup_RW(DoRys,nDiff)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Set up for contracted calculation
!
      Call Flip_Flop(.False.)
!                                                                      *
!***********************************************************************
!                                                                      *
      Call Get_EFP()
!                                                                      *
!***********************************************************************
!
      Return
      End SubRoutine GetInf_Internal
!
      End SubRoutine GetInf
