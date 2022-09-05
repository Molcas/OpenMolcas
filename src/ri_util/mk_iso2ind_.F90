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
      Subroutine Mk_iSO2Ind_(iSO2Sh,iSO2Ind,nSO,nTemp,nShell)
      use Basis_Info, only: nBas_Aux
      use Symmetry_Info, only: nIrrep
      Integer iSO2Sh(nSO), iSO2Ind(nSO), nTemp(nShell)
!
      iSO = 0
      Do iIrrep = 0, nIrrep-1
!
         Call IZero(nTemp,nShell)
         Do iB = 1, nBas_Aux(iIrrep)
            iSO = iSO + 1
            iSh = iSO2Sh(iSO)
            nTemp(iSh) = nTemp(iSh) + 1
            Ind = nTemp(iSh)
!           Write (*,*) 'iSO,iSh,Ind=',iSO,iSh,Ind
            iSO2Ind(iSO)=Ind
         End Do
!
      End Do
!     Call iVcPrt('iSO2Ind','(10I5)',iSO2Ind,nSO)
!
      Return
      End
