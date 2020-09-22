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
Module Sizes_of_Seward
Implicit None
Private
Public:: S, Size_Dmp, Size_Get
#include  "itmax.fh"
Integer i
Type Sizes_of_Stuff
     Sequence
     Integer :: Low_Anchor
     Integer :: m2Max=0
     Integer :: mCentr      =0
     Integer :: mCentr_Aux  =0
     Integer :: mCentr_Frag =0
     Integer :: Mx_mdc      =0
     Integer :: Mx_Shll     =0
     Integer :: n2Tot       =0
     Integer :: jMax        =5
     Integer :: MaxPrm(0:iTabMx)=[(0,i=0,iTabMx)]
     Integer :: MaxBas(0:iTabMx)=[(0,i=0,iTabMx)]
     Integer :: High_Anchor
End Type Sizes_of_Stuff
!
Type (Sizes_of_Stuff), Target :: S
Integer, Pointer :: p_ix(:)
Integer iix(2)
Integer nByte_i
Integer Len, Len2
Integer iiLoc
External iiLoc
Logical Found


Interface
   Subroutine Abend()
   End Subroutine Abend
   Subroutine Put_iArray(Label,Data,nData)
   Character*(*) Label
   Integer       nData
   Integer       Data(nData)
   End Subroutine Put_iArray
   Subroutine Get_iArray(Label,Data,nData)
   Character*(*) Label
   Integer       nData
   Integer       Data(nData)
   End Subroutine Get_iArray
   Subroutine Qpg_iArray(Label,Found,nData)
   Character*(*) Label
   Logical       Found
   Integer       nData
   End Subroutine Qpg_iArray
End Interface

Contains

Subroutine Size_Init()
Use Iso_C_Binding
  nbyte_i = iiloc(iix(2)) - iiloc(iix(1))
  Len = iiLoc(S%High_Anchor)-iiLoc(S%Low_Anchor)
  Len = (Len+nByte_i)/nByte_i
  Call C_F_Pointer(C_Loc(S%Low_Anchor),p_ix,[Len])
End Subroutine Size_Init

Subroutine Size_Dmp()
  Call Size_Init()
  Call Put_iArray('Sizes',p_ix,Len)
  Nullify(p_ix)
End Subroutine Size_Dmp

Subroutine Size_Get()
  Call Qpg_iArray('Sizes',Found,Len2)
  If (.NOT.Found) Then
     Write (6,*) 'Size_Get: Sizes not found.'
     Call Abend()
  End If
  Call Size_Init()
  If (Len/=Len2) Then
     Write (6,*) 'Size_Get: Len/=Len2.'
     Call Abend()
  End If
  Call Get_iArray('Sizes',p_ix,Len)
  Nullify(p_ix)
End Subroutine Size_Get

End Module Sizes_of_Seward
