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
Public :: nIrrep, iOper, iChTbl, iChCar, iChBas, lIrrep, lBsFnc, SymLab, &
          Symmetry_Info_Set, Symmetry_Info_Dmp, Symmetry_Info_Get, Symmetry_Info_Back, Symmetry_Info_Free, &
          Symmetry_Info_Setup

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
Integer:: MxFnc
Integer, Allocatable:: iChBas(:)
Character(LEN=3) :: lIrrep(0:7)=['','','','','','','','']
Character(LEN=80) :: lBsFnc(0:7)=['','','','','','','','']
Character(LEN=3) SymLab

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
! temporary routine!
Subroutine Symmetry_Info_Back(mIrrep)
Integer:: mIrrep
mIrrep=nIrrep
#ifdef _DEBUG_
Write (6,*) 'Call Symmetry_Info_Back'
Write (6,'(A,I4)') 'nIrrep=',nIrrep
#endif
End Subroutine Symmetry_Info_Back
!
!***********************************************************************
!***********************************************************************
!
Subroutine Symmetry_Info_Dmp()
Integer i, j, k, liDmp, lcDmp
Integer, Allocatable:: iDmp(:)
Character(LEN=1), Allocatable:: cDmp(:)

liDmp = 1+8+8*8+3 + MxFnc
Call mma_allocate(iDmp,liDmp,Label='iDmp')

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
iDmp(i+1:i+MxFnc)=iChBas(1:MxFnc)
i=i+MxFnc

#ifdef _DEBUG_
Write (6,*) 'Symmetry_Info_Dmp'
Write (6,*) 'liDmp=',liDmp
Write (6,*) 'MxFnc=',MxFnc
Write (6,*) 'nIrrep=',nIrrep
Write (6,*) 'iOper:'
Write (6,'(8I4)') (iOper(i),i=0,nIrrep-1)
Write (6,*)
Write (6,'(9I4)') (iDmp(i),i=1,liDmp)
Write (6,*)
Write (6,*) 'lIrrep:'
Do i = 0, nIrrep-1
   Write (6,'(A)') lIrrep(i)
End Do
Write (6,*) 'lBsFnc:'
Do i = 0, nIrrep-1
   Write (6,'(A)') lBsFnc(i)
End Do
Write (6,'(2A)') 'SymLab:',SymLab
#endif

Call Put_iArray('Symmetry Info',iDmp,liDmp)
Call mma_deallocate(iDmp)

lcDmp = 3*8 + 80*8 + 3
Call mma_allocate(cDmp,lcDmp,Label='cDmp')
k = 0
Do i = 0, 7
   Do j = 1, 3
      cDmp(j+k)=lIrrep(i)(j:j)
   End Do
   k=k+3
End Do
Do i = 0, 7
   Do j = 1, 80
      cDmp(j+k)=lBsFnc(i)(j:j)
   End Do
   k=k+80
End Do
Do i = 1, 3
   cDmp(i+k)=SymLab(i:i)
End Do
k=k+3
Call put_cArray('SymmetryCInfo',cDmp(1),lcDmp)
Call mma_deallocate(cDmp)

End Subroutine Symmetry_Info_Dmp
!
!***********************************************************************
!***********************************************************************
!
Subroutine Symmetry_Info_Get()
Integer i, j, k, liDmp, lcDmp
Integer, Allocatable:: iDmp(:)
Logical Found
Character(LEN=1), Allocatable:: cDmp(:)

If (Allocated(iChBas)) Return
Call Qpg_iArray('Symmetry Info',Found,liDmp)
Call mma_allocate(iDmp,liDmp,Label='iDmp')
Call Get_iArray('Symmetry Info',iDmp,liDmp)

MxFnc=liDmp - (1+8+8*8+3)
Call mma_allocate(iChBas,MxFnc,Label='iChBas')

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
iChBas(1:MxFnc) = iDmp(i+1:i+MxFnc)
i=i+MxFnc
Call mma_deallocate(iDmp)

lcDmp = 3*8 + 80*8 + 3
Call mma_allocate(cDmp,lcDmp,Label='cDmp')
Call get_carray('SymmetryCInfo',cDmp(1),lcDmp)
k = 0
Do i = 0, 7
   Do j = 1, 3
      lIrrep(i)(j:j)=cDmp(j+k)
   End Do
   k=k+3
End Do
Do i = 0, 7
   Do j = 1, 80
      lBsFnc(i)(j:j)=cDmp(j+k)
   End Do
   k=k+80
End Do
Do i = 1, 3
   SymLab(i:i)=cDmp(i+k)
End Do
Call mma_deallocate(cDmp)
#ifdef _DEBUG_
Write (6,*) 'Symmetry_Info_Get'
Write (6,*) 'liDmp=',liDmp
Write (6,*) 'MxFnc=',MxFnc
Write (6,*) 'nIrrep=',nIrrep
Write (6,*)
Write (6,'(2A)') 'SymLab:',SymLab
Write (6,*) 'iOper:'
Write (6,'(8I4)') (iOper(i),i=0,nIrrep-1)
Write (6,*)
Write (6,*) 'lIrrep:'
Do i = 0, nIrrep-1
   Write (6,'(A)') lIrrep(i)
End Do
Do i = 0, nIrrep-1
   Write (6,'(A)') lBsFnc(i)
End Do
#endif
End Subroutine Symmetry_Info_Get
!
!***********************************************************************
!***********************************************************************
!
Subroutine Symmetry_Info_Free()
If (.not.Allocated(iChBas)) Return
Call mma_deallocate(iChBas)
MxFnc=0
End Subroutine Symmetry_Info_Free
!
!***********************************************************************
!***********************************************************************
!
Subroutine Symmetry_Info_Setup(nOper,Oper,iAng)
Implicit None
Integer :: nOper, iAng, i, j
Character(LEN=3) :: Oper(3)

If (Allocated(iChBas)) Return   ! Return if already initiated.

nIrrep = 2 ** nOper
Call Put_iScalar('NSYM',nIrrep)

iOper(0) = 0
Do i = 1, nOper
   iOper(i) = 0
   Do j = 1, 3
    If(Oper(i)(j:j).eq.'X') iOper(i) = iOper(i) + 1
    If(Oper(i)(j:j).eq.'Y') iOper(i) = iOper(i) + 2
    If(Oper(i)(j:j).eq.'Z') iOper(i) = iOper(i) + 4
   End Do
   If (iOper(i).eq.0) Then
      Call WarningMessage(2,'RdCtl: Illegal symmetry operator!')
      Write (6,*) 'Oper=',Oper(i)
      Write (6,*)
      Call Abend()
   End If
End Do
!                                                                      *
!***********************************************************************
!                                                                      *
!  Generate all operations of the group
!
If (nOper.ge.2) Then
   iOper(4) = iOper(3)
   iOper(3) = iEor(iOper(1),iOper(2))
End If
If (nOper.eq.3) Then
   iOper(5) = iEor(iOper(1),iOper(4))
   iOper(6) = iEor(iOper(2),iOper(4))
   iOper(7) = iEor(iOper(1),iEor(iOper(2),iOper(4)))
End If

Call Put_iArray('Symmetry operations',iOper,nIrrep)

!     Generate the Character table for all Irreps, iChTbl

!     All Irreps are one dimensional, i.e. the Character for the
!     unit operator is 1 in all irreps.
!     The totally symmetric representation will have the character
!     of 1 for any given operation
!     Now, the Irreps are due to classes of operations and will
!     present the character of this class. In case of Abelian groups
!     or other one dimensional groups the classes will have one
!     and only one operation. Hence, the operations themselves can
!     be used to present the character of the Irreps.

! Generate iChTbl, lIrrep, lBsFnc (iSigma)
Call ChTab(iOper,nIrrep,iChTbl)

! Generate iChCar, iChBas, MxFnc
Call Symmetry_Info_Set(iAng)

End Subroutine Symmetry_Info_Setup
!
!***********************************************************************
!***********************************************************************
!
Subroutine Symmetry_Info_Set(iAng)
Integer:: iIrrep, jIrrep
Integer:: iSymX,iSymY,iSymZ, i
Integer:: iAng, lxyz, ixyz, ix, jx, iyMax, iy, jy, iz, jz, jxyz

If (allocated(iChBas)) Return

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

MxFnc=(iAng+1)*(iAng+2)*(iAng+3)/6
Call mma_allocate(iChBas,MxFnc,Label='iChBas')
#ifdef _DEBUG_
Write (6,*) 'Symmetry_Info_Set:'
Write (6,*) 'iAng,MxFnc=',iAng,MxFnc
#endif

lxyz = 0
Do ixyz = 0, iAng
   Do ix = ixyz, 0, -1
      jx = Mod(ix,2)
      iyMax=ixyz-ix
      Do iy = iyMax, 0 , -1
         jy = Mod(iy,2)
         lxyz=lxyz+1
         iz=ixyz-ix-iy
         jz = Mod(iz,2)
         jxyz = jx * iSymX + jy * iSymY + jz * iSymZ
         iChBas(lxyz) = jxyz
      End Do
   End Do
End Do

Do iIrrep=0,nIrrep-2
   Do jIrrep=iIrrep+1,nIrrep-1
      If (iOper(iIrrep).eq.iOper(jIrrep)) Then
         Call WarningMessage(2,   &
              ' The generators of the point group are over defined, correct input!;' //' Abend: correct symmetry specifications!')
               Call Quit_OnUserError()
      End If
   End Do
End Do
#ifdef _DEBUG_
Write (6,*) 'Symmetry_Info_Set:'
Write (6,*) 'MxFnc=',MxFnc
Write (6,*) 'nIrrep=',nIrrep
Write (6,'(A,8I4)') 'iOper:',iOper(0:nIrrep-1)
Write (6,*) 'iChTbl:'
Do i = 0, nIrrep-1
   Write (6,'(8I4)') iChTbl(0:nIrrep-1,i)
End Do
#endif
End Subroutine Symmetry_Info_Set
!
!***********************************************************************
!***********************************************************************
!
SubRoutine ChTab(iOper,nIrrep,iChTbl)
!***********************************************************************
!                                                                      *
! Object: to generate the character table of a point group within      *
!         D2h.                                                         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             September '91                                            *
!***********************************************************************
Implicit None
Integer nIrrep
Integer iOper(nIrrep), iChTbl(1:8,1:8) ! ugly dimensions change to 0:7!
Integer iTest(8)
Integer :: iSigma=1
Character(Len=80) Tmp
Logical Inv, Rot
Character(LEN=6):: xyz(0:7)=['      ','x     ','y     ','xy, Rz', 'z     ','xz, Ry','yz, Rx','I     ']
Integer i, i1, i2, ia, ib, iCh, iFnc, iIrrep, iRot, iSub, iSymX, iSymY, iSymZ, ix, iy, iz, jIrrep
Integer j, jx, jy, jz, Lenlbs, LenlIrr, LenTmp
Integer iclast
External iclast
!                                                                      *
!***********************************************************************
!                                                                      *
If (nIrrep.eq.1) Then
   SymLab='C1 '
   iSigma=1
Else If (nIrrep.eq.2) Then
   If (iOper(2).eq.7) Then
      SymLab='Ci '
      iSigma=1
   Else If (iOper(2).eq.1.or.iOper(2).eq.2.or.iOper(2).eq.4) Then
      SymLab='Cs'
      iSigma=1
   Else
      SymLab='C2'
      iSigma=2
   End If
Else If (nIrrep.eq.4) Then
   If (iOper(2).eq.7.or.iOper(3).eq.7.or.iOper(4).eq.7) Then
      SymLab='C2h'
      iSigma=2
   Else
      Rot = .True.
      Do i = 1, nIrrep
         If (iOper(i).eq.1.or.iOper(i).eq.2.or.iOper(i).eq.4) Rot = .False.
      End Do
      If (Rot) Then
         SymLab='D2 '
         iSigma=2
      Else
         SymLab='C2v'
         iSigma=2
      End If
   End If
Else If (nIrrep.eq.8) Then
   SymLab='D2h'
   iSigma=2
Else
   Call WarningMessage(2,'ChTab: Illegal value of nIrrep')
   Write (6,*) 'nIrrep=',nIrrep
   Call Abend()
End If
ichTbl(:,:)=0
!
!     Go through the functions x, y, and z, and the dyadic functions.
!
iSymX = 0
iSymY = 0
iSymZ = 0
Do i = 1, nIrrep
   If (iAnd(iOper(i),1).ne.0) iSymX = 1
   If (iAnd(iOper(i),2).ne.0) iSymY = 2
   If (iAnd(iOper(i),4).ne.0) iSymZ = 4
End Do
!
!-----Loop over basis functions (a' la Malmqvist)
!
lBsFnc(0:nIrrep-1)='' ! For this to work we need a clean slate.
Do iFnc = 0, 7
   Tmp=xyz(iFnc)

!  Generate a row in the character table of this function

   ix = iAnd(iFnc,iSymX)
   iy = iAnd(iFnc,iSymY)/2
   iz = iAnd(iFnc,iSymZ)/4
!--Loop over all operators
   Do i = 1, nIrrep
      jx = iAnd(iOper(i),iSymX)
      jy = iAnd(iOper(i),iSymY)/2
      jz = iAnd(iOper(i),iSymZ)/4
      iCh = 1
      If (ix.ne.0 .and. jx.ne.0) iCh = -iCh
      If (iy.ne.0 .and. jy.ne.0) iCh = -iCh
      If (iz.ne.0 .and. jz.ne.0) iCh = -iCh
      iTest(i) = iCh
   End Do
!
!--------Compute place of Irrep
!
   If (nIrrep.eq.1) Then
      jIrrep=1
   Else If (nIrrep.eq.2) Then
      jIrrep=1+(1-iTest(2))/2
   Else If (nIrrep.eq.4) Then
      jIrrep=1+( (1-iTest(2))+2*(1-iTest(3)) )/2
   Else If (nIrrep.eq.8) Then
      jIrrep=1+( (1-iTest(2))+2*(1-iTest(3))+4*(1-iTest(5)) )/2
   Else
      jIrrep=-1
      Call WarningMessage(2,'ChTab: Illegal nIrrep value!')
      Write (6,*) 'nIrrep=',nIrrep
      Call Abend()
   End If
   If (lBsFnc(jIrrep-1)(1:1).eq.' ') Then
      lBsFnc(jIrrep-1) = Tmp
      Call ICopy(nIrrep,iTest,1,iChTbl(jIrrep,1),8)
   Else
      LenlBs=Len(lBsFnc(jIrrep-1))
      LenTmp=Len(Tmp)
      i1 = iCLast(lBsFnc(jIrrep-1),LenlBs)
      i2 = iCLast(Tmp,LenTmp)
      lBsFnc(jIrrep-1) = lBsFnc(jIrrep-1)(1:i1)//', '//Tmp(1:i2)
   End If
End Do
!
!     Set up some Mulliken symbols for the irreps
!
Do iIrrep = 1, nIrrep
   lIrrep(iIrrep-1)='a'
   Do i = 1, nIrrep

!     If the character of an rotation in an irreps is -1 then
!     the irreps is assigned the character B, otherwise A.

      If ((iOper(i).eq.3 .or. iOper(i).eq.5 .or. iOper(i).eq.6) .and. iChTbl(iIrrep,i).eq.-1) lIrrep(iIrrep-1)='b'

   End Do
End Do
iSub = 0

!  Subscript according to C2 operations

Rot = .False.
Do i = 1, nIrrep
   If (iOper(i).eq.3 .or. iOper(i).eq.5 .or. iOper(i).eq.6) Rot = .True.
End Do
If (Rot.and.SymLab.ne.'C2v') Then
   iSub = iSub + 1

!  Find the number of A's and B's

   ia = 0
   ib = 0
   Do i = 1, nIrrep
      If (lIrrep(i-1)(1:1).eq.'a') ia = ia + 1
      If (lIrrep(i-1)(1:1).eq.'b') ib = ib + 1
   End Do
   If (nIrrep.eq.8) Then
      ia = ia/2
      ib = ib/2
   End If
   If (SymLab.eq.'C2h') Then
      ia = ia/2
      ib = ib/2
   End If

!  Find the rotations

   iRot = 0
   Do i = 1, nIrrep
      If ( iOper(i).eq.3.or.iOper(i).eq.5.or.iOper(i).eq.6) Then
         iRot = iRot + 1
         Write (Tmp,'(I1)') iRot
         If (ia.gt.1) Then
            Do j = 1, nIrrep
               If (lIrrep(j-1)(1:1).eq.'a'.and.iChTbl(j,i).eq.1) lIrrep(j-1)=lIrrep(j-1)(1:1)//Tmp(1:1)
            End Do
         End If
         If (ib.gt.1) Then
            Do j = 1, nIrrep
               If (lIrrep(j-1)(1:1).eq.'b'.and.iChTbl(j,i).eq.1) lIrrep(j-1)=lIrrep(j-1)(1:1)//Tmp(1:1)
            End Do
         End If
      End If
   End Do
Else If (Rot.and.SymLab.eq.'C2v') Then

!  Find the Rotation

   iRot = -1
   Do i = 1, nIrrep
      If (iOper(i).eq.3.or.iOper(i).eq.5.or.iOper(i).eq.6) iRot = iOper(i)
   End Do

!  Find the first vertical mirror plane to this axis

   Do i = 1, nIrrep
      If (iOper(i).ne.3.and.iOper(i).ne.5.and.iOper(i).ne.6 .and.iOper(i).ne.7.and.iAnd(iOper(i),iRot).ne.1) iRot = i
   End Do
   Do i = 1, nIrrep
      If (iChTbl(i,iRot).eq.1) Then
         j = 1
      Else
         j = 2
      End If
      Write (Tmp,'(I1)') j
      lIrrep(i-1)=lIrrep(i-1)(1:1)//Tmp(1:1)
   End Do
End If

!  Subscript according to inversion if present

Inv=.False.
Do i = 1, nIrrep
   Inv = iOper(i).eq.7 .or. Inv
End Do
If (Inv) Then
   iSub = iSub + 1

!- Loop over each Irrep

   Do iIrrep = 1, nIrrep
      LenlIrr=Len(lIrrep(iIrrep-1))
      i1 = 1 + iCLast(lIrrep(iIrrep-1),LenlIrr)

!---- Loop over operators

      Do i = 1, nIrrep
         If (iOper(i).eq.7) Then
            If (iChTbl(iIrrep,i).eq.1) Then
               lIrrep(iIrrep-1)(i1:i1)='g'
            Else If (iChTbl(iIrrep,i).eq.-1) Then
                lIrrep(iIrrep-1)(i1:i1)='u'
            End If
         End If
      End Do
   End Do
End If

!  Fix labels for Cs

If (SymLab(1:2).eq.'Cs') Then
   lIrrep(0) = 'a'''
   lIrrep(1) = 'a"'
End If
!                                                                      *
!***********************************************************************
!                                                                      *
Call Put_iScalar('Rotational Symmetry Number',iSigma)
Call Put_cArray('Irreps',lIrrep(0),24)
!                                                                      *
!***********************************************************************
!                                                                      *
Return
End Subroutine ChTab
!
!***********************************************************************
!***********************************************************************
!
End Module Symmetry_Info
