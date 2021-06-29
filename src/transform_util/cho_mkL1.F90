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
! Copyright (C) 2005, Giovanni Ghigo                                   *
!***********************************************************************
      Subroutine MkL1(iSymA,iSymI,iI, numV, LyType,iJy, AddLx0, SameLx)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the Cholesky matrix of Inactive(iSymA) for   *
!           occupied iI(iSymI) for numV vectors.                       *
!***********************************************************************
      use Cho_Tra
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Integer iSymA,iSymI,iI, numV, LyType,iJy
      Real*8 AddLx0(*)
      Logical SameLx
#include "rasdim.fh"
#include "SysDef.fh"

!     Build Lx
      If (iI.LE.nIsh(iSymI)) then
        LxType = 1
        iIx = iI
      else
        LxType = 7
        iIx = iI - nIsh(iSymI)
      EndIf

      If (.NOT.SameLx) then
        LyType=LxType
        iJy   =iIx
      else
        If (LyType.EQ.LxType .and. iIx.EQ.iJy) then
          Return
        else
          SameLx=.False.
        EndIf
      EndIf

      iAddTCVX= 1+nIsh(iSymA)*(iIx-1)

      iAddLx  = 1
      Do iV=1,numV
        Call dCopy_(nIsh(iSymA),                                        &
     &              TCVX(LxType,iSymA,iSymI)%A(iAddTCVX,iV),1,          &
     &              AddLx0(iAddLx),1)
        iAddLx  = iAddLx + nIsh(iSymA)
      EndDo

      Return
      End
