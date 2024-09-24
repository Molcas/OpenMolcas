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
!#define _DEBUGPRINT_
      Subroutine OptClc_QNR(CInter,nCI,nD,Grd1,Xnp1,mOV,Ind,MxOptm,kOptim,kOV)
      use LnkLst, only: LLGrad,LLx
      Implicit None
      Integer nCI,nD,mOV,MxOptm,kOptim,kOV(2)
      Real*8 CInter(nCI,nD), Grd1(mOV), Xnp1(mOV)
      Integer Ind(MxOptm)
      Interface
         Subroutine OptClc_X(CInter,nCI,nD,Array,mOV,Ind,MxOptm,kOptim,kOV,LL,DD)
         Implicit None
         Integer nCI,nD,mOV,MxOptm,kOptim,kOV(2), LL
         Real*8 CInter(nCI,nD), Array(mOV)
         Integer Ind(MxOptm)
         Real*8, Optional:: DD
         End Subroutine OptClc_X
      End Interface

!
      Call  OptClc_X(CInter,nCI,nD,Grd1,mOV,Ind,MxOptm,kOptim,kOV,LLGrad)
      Call  OptClc_X(CInter,nCI,nD,Xnp1,mOV,Ind,MxOptm,kOptim,kOV,LLx)
      Return
      End Subroutine OptClc_QNR
