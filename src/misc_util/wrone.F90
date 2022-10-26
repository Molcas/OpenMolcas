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
      Subroutine WrOne(rc,Option,InLab,Comp,Data,SymLab)
      Implicit Integer (A-Z)
!
      Character*(*) InLab
      Real*8 Data(*)
!
      Call WrOne_Internal(Data)
!
!     This is to allow type punning without an explicit interface
      Contains
      Subroutine WrOne_Internal(Data)
      Use Iso_C_Binding
      Real*8, Target :: Data(*)
      Integer, Pointer :: iData(:)
      Call C_F_Pointer(C_Loc(Data(1)),iData,[1])
      Call iWrOne(rc,Option,InLab,Comp,iData,SymLab)
      Nullify(iData)
      return
      End Subroutine WrOne_Internal
!
      end
