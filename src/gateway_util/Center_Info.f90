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
!#define _DEBUGPRINT_

Module Center_Info
Implicit None
Private
Public :: dc, n_dc, Center_Info_Init, Center_Info_Dmp, Center_Info_Get, Center_Info_Free

#include "Molcas.fh"
#include "stdalloc.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    D E C L A R E   D E R I V E D   T Y P E S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
Type Distinct_centers
!   Sequence
    Integer :: iChCnt=0
    Integer :: iStab(0:7)=[0,0,0,0,0,0,0,0]
    Integer :: nStab=0
    Integer :: iCoSet(0:7,0:7)=ReShape([0,0,0,0,0,0,0,0,                    &
                                        0,0,0,0,0,0,0,0,                    &
                                        0,0,0,0,0,0,0,0,                    &
                                        0,0,0,0,0,0,0,0,                    &
                                        0,0,0,0,0,0,0,0,                    &
                                        0,0,0,0,0,0,0,0,                    &
                                        0,0,0,0,0,0,0,0,                    &
                                        0,0,0,0,0,0,0,0],[8,8])
    Character(LEN=LENIN4):: LblCnt=''
End Type Distinct_centers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     E N D   D E C L A R E   D E R I V E D   T Y P E S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Integer, Parameter:: nFields=1+8+1+64
Logical :: Initiated = .FALSE.
Integer :: n_dc=0
Type (Distinct_centers), Allocatable:: dc(:)

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
!     This to make either the initial allocation of dbsc and Shells according to the default sizes
!     as defined by the parameters in Molcas.fh or according to the actual sizes as recorded on the
!     run file.
!
Subroutine Center_Info_Init()
If (Initiated) Then
   Write (6,*) 'Center_Info already initiated!'
   Write (6,*) 'May the is a missing call to Center_Info_Free.'
   Call Abend()
End If
#ifdef _DEBUGPRINT_
Write (6,*)
Write (6,*) 'Enter Center_Info_Init'
#endif
If (n_dc.eq.0) Then
   Allocate(dc(1:MxAtom))
Else
   Allocate(dc(1:n_dc))
End If
Initiated=.True.
#ifdef _DEBUGPRINT_
Write (6,*) 'Exit Center_Info_Init'
#endif
Return
End Subroutine Center_Info_Init
!
!***********************************************************************
!***********************************************************************
!
Subroutine Center_Info_Dmp()
!
Integer i, j, k, licDmp, lcDmp
Integer, Allocatable:: iDmp(:)
Character(LEN=1), Allocatable:: cDmp(:)
!
!     Integer dc stuff
!
#ifdef _DEBUGPRINT_
Write (6,*)
Write (6,*) 'Enter Center_Info_Dmp'
Write (6,*) 'n_dc=',n_dc
#endif
licDmp=n_dc*nFields
Call mma_Allocate(iDmp,licDmp+1,Label='iDmp')
j=0
Do i = 1, n_dc
   iDmp(j+1)=dc(i)%iChCnt
   j=j+1
   iDmp(j+1:j+8)=dc(i)%iStab(0:7)
   j=j+8
   iDmp(j+1)=dc(i)%nStab
   j=j+1
   iDmp(j+1:j+8)=dc(i)%iCoSet(:,0)
   j=j+8
   iDmp(j+1:j+8)=dc(i)%iCoSet(:,1)
   j=j+8
   iDmp(j+1:j+8)=dc(i)%iCoSet(:,2)
   j=j+8
   iDmp(j+1:j+8)=dc(i)%iCoSet(:,3)
   j=j+8
   iDmp(j+1:j+8)=dc(i)%iCoSet(:,4)
   j=j+8
   iDmp(j+1:j+8)=dc(i)%iCoSet(:,5)
   j=j+8
   iDmp(j+1:j+8)=dc(i)%iCoSet(:,6)
   j=j+8
   iDmp(j+1:j+8)=dc(i)%iCoSet(:,7)
   j=j+8
End Do
iDmp(licDmp+1)=n_dc
Call Put_iArray('icDmp',iDmp,licDmp+1)
Call mma_deallocate(iDmp)

lcDmp=n_dc*LENIN4
#ifdef _DEBUGPRINT_
Write (6,*) 'lcDmp=',n_dc
#endif
Call mma_allocate(cDmp,lcDmp,Label='cDmp')
k = 0
Do i = 1, n_dc
   Do j = 1, LENIN4
      cDmp(k+j) = dc(i)%LblCnt(j:j)
   End Do
   k = k + LENIN4
End Do
#ifdef _DEBUGPRINT_
Write (6,*) 'cDmp=',cDmp(1:lcDmp)
#endif
Call Put_cArray('dc: cDmp',cDmp(1),lcDmp)
Call mma_deallocate(cDmp)
#ifdef _DEBUGPRINT_
Write (6,*) 'Exit Center_Info_Dmp'
#endif
Return
End Subroutine Center_Info_Dmp
!
!***********************************************************************
!***********************************************************************
!
Subroutine Center_Info_Get()
!
Integer, Allocatable:: iDmp(:)
Character(LEN=1), Allocatable:: cDmp(:)
Logical Found
Integer i, j, k, Len,lcDmp

#ifdef _DEBUGPRINT_
Write (6,*)
Write (6,*) 'Enter Center_Info_Get'
#endif
Call qpg_iArray('icDmp',Found,Len)
Call mma_Allocate(iDmp,Len,Label='iDmp')
If (Found) Then
   Call Get_iArray('icDmp',iDmp,Len)
Else
   Write (6,*) 'Center_Info_Get: icDmp not found!'
   Call Abend()
End If
lcDmp=Len-1
n_dc=lcDmp/nFields
!
!     Initiate the memory allocation of dc
!
If (.Not.Initiated) Call Center_Info_Init()
#ifdef _DEBUGPRINT_
Write (6,*) 'iDmp(1:Len)=',iDmp(1:Len)
Write (6,*) 'Len=',Len
Write (6,*) 'n_dc=',n_dc
Write (6,*) 'lcDmp=',lcDmp
#endif
j=0
Do i = 1, n_dc
   dc(i)%iChCnt = iDmp(j+1)
   j = j + 1
   dc(i)%iStab(0:7) = iDmp(j+1:j+8)
   j = j + 8
   dc(i)%nStab  = iDmp(j+1)
   j = j + 1
   dc(i)%iCoSet(:,0) = iDmp(j+1:j+8)
   j = j + 8
   dc(i)%iCoSet(:,1) = iDmp(j+1:j+8)
   j = j + 8
   dc(i)%iCoSet(:,2) = iDmp(j+1:j+8)
   j = j + 8
   dc(i)%iCoSet(:,3) = iDmp(j+1:j+8)
   j = j + 8
   dc(i)%iCoSet(:,4) = iDmp(j+1:j+8)
   j = j + 8
   dc(i)%iCoSet(:,5) = iDmp(j+1:j+8)
   j = j + 8
   dc(i)%iCoSet(:,6) = iDmp(j+1:j+8)
   j = j + 8
   dc(i)%iCoSet(:,7) = iDmp(j+1:j+8)
   j = j + 8
End Do
Call mma_deAllocate(iDmp)
!
lcDmp=n_dc*LENIN4
#ifdef _DEBUGPRINT_
Write (6,*) 'lcDmp=',lcDmp
#endif
!
Call qpg_cArray('dc: cDmp',Found,Len)
If (Len/=lcDmp) Then
   Write (6,*) 'Center_Info_Get: Len /= lcDmp'
   Call Abend()
End If
Call mma_Allocate(cDmp,lcDmp,Label='cDmp')
Call Get_cArray('dc: cDmp',cDmp(1),lcDmp)
#ifdef _DEBUGPRINT_
Write (6,*) 'cDmp=',cDmp(1:lcDmp)
#endif
k=0
Do i = 1, n_dc
   Do j = 1, LENIN4
      dc(i)%LblCnt(j:j)  = cDmp(k+j)
   End Do
   k = k + LENIN4
End Do
Call mma_deAllocate(cDmp)
#ifdef _DEBUGPRINT_
Write (6,*) 'Exit Center_Info_Get'
#endif
!
End Subroutine Center_Info_Get
!
!***********************************************************************
!***********************************************************************
!
Subroutine Center_Info_Free()
!
!     Deallocate all allocatable parts of dc.
!
If (.Not.Allocated(dc)) Return
#ifdef _DEBUGPRINT_
Write (6,*)
Write (6,*) 'Enter Center_Info_Free'
#endif
Deallocate(dc)
n_dc=0
Initiated=.False.
!
#ifdef _DEBUGPRINT_
Write (6,*) 'Exit Center_Info_Free'
#endif
Return
End Subroutine Center_Info_Free
!
!***********************************************************************
!***********************************************************************
!
End Module Center_Info
