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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUG_
Module Symmetry_Info
Implicit None
Private
Public :: nIrrep, iOper, iChTbl, iChCar, &
          Symmetry_Info_Set, Symmetry_Info_Dmp, Symmetry_Info_Get, Symmetry_Info_Back

#include "stdalloc.fh"
Integer:: nIrrep=1
Integer:: iOper(0:7)=[0,0,0,0,0,0,0,0]
Integer:: iChTbl(0:7,0:7)=Reshape([0,0,0,0,0,0,0,0,      &
                                   0,0,0,0,0,0,0,0,      &
                                   0,0,0,0,0,0,0,0,      &
                                   0,0,0,0,0,0,0,0,      &
                                   0,0,0,0,0,0,0,0,      &
                                   0,0,0,0,0,0,0,0,      &
                                   0,0,0,0,0,0,0,0,      &
                                   0,0,0,0,0,0,0,0],[8,8])
Integer:: iChCar(3)=[0,0,0]

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
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
!
!***********************************************************************
!***********************************************************************
!
Contains
!
!***********************************************************************
!***********************************************************************
!
Subroutine Symmetry_Info_Back(mIrrep,jOper)
Integer:: mIrrep
Integer:: jOper(0:7)
mIrrep=nIrrep
jOper(:)=iOper(0:nIrrep-1)
End Subroutine Symmetry_Info_Back
!
!***********************************************************************
!***********************************************************************
!
Subroutine Symmetry_Info_Set(mIrrep,jOper,jChTab)
Integer:: mIrrep
Integer:: jOper(0:7)
Integer:: jChTab(0:7,0:7)
Integer:: iIrrep, jIrrep
Integer:: iSymX,iSymY,iSymZ, i

nIrrep=mIrrep
iOper(:)=jOper(:)
iChTbl(:,:)=jChTab(:,:)

! Setup characteristics for cartesian basis functions.
! Observe that this is affected by the defined generators.
! In the array we will set the bit corresponding to a symop
! if that symop will alter the sign of the basis function.

iSymX = 0
iSymY = 0
iSymZ = 0
Do i = 0, nIrrep-1
   If (iAnd(iOper(i),1).ne.0) iSymX = 1
   If (iAnd(iOper(i),2).ne.0) iSymY = 2
   If (iAnd(iOper(i),4).ne.0) iSymZ = 4
End Do
iChCar(1) = iSymX
iChCar(2) = iSymY
iChCar(3) = iSymZ

Do iIrrep=0,nIrrep-2
   Do jIrrep=iIrrep+1,nIrrep-1
      If (iOper(iIrrep).eq.iOper(jIrrep)) Then
         Call WarningMessage(2,   &
              ' The generators of the point group are over defined, correct input!;' //' Abend: correct symmetry specifications!')
               Call Quit_OnUserError()
      End If
   End Do
End Do
End Subroutine Symmetry_Info_Set
!
!***********************************************************************
!***********************************************************************
!
Subroutine Symmetry_Info_Dmp()
Integer i, iDmp(1+8+8*8+3)
i=0
iDmp(i+1)=nIrrep
i=i+1
iDmp(i+1:i+8)=iOper(:)
i=i+8
iDmp(i+1:i+8)=iChTbl(:,0)
i=i+8
iDmp(i+1:i+8)=iChTbl(:,1)
i=i+8
iDmp(i+1:i+8)=iChTbl(:,2)
i=i+8
iDmp(i+1:i+8)=iChTbl(:,3)
i=i+8
iDmp(i+1:i+8)=iChTbl(:,4)
i=i+8
iDmp(i+1:i+8)=iChTbl(:,5)
i=i+8
iDmp(i+1:i+8)=iChTbl(:,6)
i=i+8
iDmp(i+1:i+8)=iChTbl(:,7)
i=i+8
iDmp(i+1:i+3)=iChCar(1:3)
i=i+3
Call Put_iArray('Symmetry Info',iDmp,i)
End Subroutine Symmetry_Info_Dmp
!
!***********************************************************************
!***********************************************************************
!
Subroutine Symmetry_Info_Get()
Integer i, iDmp(1+8+8*8+3)
i=1+8+8*8+3
Call Get_iArray('Symmetry Info',iDmp,i)
i=0
nIrrep     =iDmp(i+1)
i=i+1
iOper(:)   =iDmp(i+1:i+8)
i=i+8
iChTbl(:,0)=iDmp(i+1:i+8)
i=i+8
iChTbl(:,1)=iDmp(i+1:i+8)
i=i+8
iChTbl(:,2)=iDmp(i+1:i+8)
i=i+8
iChTbl(:,3)=iDmp(i+1:i+8)
i=i+8
iChTbl(:,4)=iDmp(i+1:i+8)
i=i+8
iChTbl(:,5)=iDmp(i+1:i+8)
i=i+8
iChTbl(:,6)=iDmp(i+1:i+8)
i=i+8
iChTbl(:,7)=iDmp(i+1:i+8)
i=i+8
iChCar(1:3)=iDmp(i+1:i+3)
i=i+3
#ifdef _DEBUG_
Write (6,*)
Write (6,*) 'Symmetry_Info_Get'
Write (6,*)
Write (6,*)
Write (6,*) 'iOper:'
Write (6,'(8I4)') (iOper(i),i=0,nIrrep-1)
#endif
End Subroutine Symmetry_Info_Get
!
!***********************************************************************
!***********************************************************************
!
End Module Symmetry_Info
