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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************
      Subroutine Check_Amp(nSym,nOcc,nVir,iSkip)
      use definitions, only: iwp
      use Symmetry_Info, only: Mul

      Implicit None
      integer(kind=iwp), intent(in)::  nSym, nOcc(nSym), nVir(nSym)
      Integer(kind=iwp), intent(out):: iSkip

      integer(kind=iwp) iSym, iSymi, iSyma
      Integer(kind=iwp) nT1amTot, nT1am(8)

      iSkip=0
      nT1amTot=0
      Do iSym = 1,nSym
         nT1am(iSym) = 0
         Do iSymi = 1,nSym
            iSyma = Mul(iSymi,iSym)
            nT1am(iSym) = nT1am(iSym)                                   &
     &                  + nVir(iSyma)*nOcc(iSymi)
         End Do
         nT1amTot = nT1amTot + nT1am(iSym)
      End Do

      If (nT1amTot .gt. 0) iSkip=1
      End Subroutine Check_Amp
