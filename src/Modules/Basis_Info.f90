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
      Module Basis_Info
      Implicit None
      Private
      Public :: Basis_Info_Dmp, Basis_Info_Get, Basis_Info_Free,  &
                Distinct_Basis_set_Centers, dbsc, nFrag_LineWords,&
                PAMExp, Shells, Max_Shells
#include "stdalloc.fh"
#include "Molcas.fh"
#include "itmax.fh"
      Integer, Parameter :: Mxdbsc=MxAtom
      Integer, Parameter :: MxShll=iTabMx*MxAtom
!
!     Work in progress
!
!     nCntr : number of centers associated with a dbsc
!     Coor  : the coordinates of a dbsc
!
!     nM1   : number of ECP M1 type terms on the i''th unique center
!     M1xp  : ECP M1-type exponents for i''th unq center
!     M1cf  : ECP M1 type coefficients for i''th unq cntr
!     nM2   : number of ECP M2 type terms on the i''th unique center
!     M2xp  : ECP M2-type exponents for i''th unq center
!     M2cf  : ECP M2 type coefficients for i''th unq cntr
!
!     nFragType:  number of unique centers in a fragment (0=not a frag)
!     nFragCoor:  total number of centers in a fragment
!     nFragEner:  number of orbital energies/occupied orbitals in a given fragment
!     nFragDens:  size of the fragments density matrix
!     FragType : the data of the fragment''s unique centers (associated basis set, size nFragType)
!     FragCoor : the data of all fragment''s centers (atom type / relative coordinates, size nFragCoor)
!     FragEner : the fragment''s orbital''s energies (size nFragEner)
!     FragCoef : the fragment''s MO coefficients (size nFragDens*nFragEner)
!
      Type Distinct_Basis_set_centers
          Sequence
          Real*8, Allocatable:: Coor(:,:)
          Integer:: nCntr=0
          Integer:: nM1=0
          Real*8, Allocatable:: M1xp(:), M1cf(:)
          Integer:: nM2=0
          Real*8, Allocatable:: M2xp(:), M2cf(:)
          Integer:: nFragType=0, nFragCoor=0, nFragEner=0, nFragDens=0
          Real*8, Allocatable:: FragType(:,:), FragCoor(:,:), FragEner(:), FragCoef(:,:)
          Integer:: nPAM2=-1
          Real*8, Allocatable:: PAM2(:)
      End Type Distinct_Basis_set_centers
!
!     Bk  : ECP proj shift parameters for i''th shell.
!           the number of parameters is given by nBasis
!     Occ : Occupation numbers for core ECP orbitals
!
      Type Shell_Info
           Sequence
           Integer :: nBk=0
           Real*8, Allocatable:: Bk(:)
           Real*8, Allocatable:: Occ(:)
      End Type Shell_Info
!
      Real*8, Allocatable:: PAMexp(:,:)
      Integer :: nFrag_LineWords=0, nFields=7, mFields=1
      Type (Distinct_Basis_set_centers) :: dbsc(Mxdbsc)
      Integer :: Max_Shells=0
      Type (Shell_Info) :: Shells(MxShll)
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
         Subroutine Put_dArray(Label,Data,nData)
         Character*(*) Label
         Integer       nData
         Real*8        Data(nData)
         End Subroutine Put_dArray
         Subroutine Get_dArray(Label,Data,nData)
         Character*(*) Label
         Integer       nData
         Real*8        Data(nData)
         End Subroutine Get_dArray
         Subroutine Qpg_dArray(Label,Found,nData)
         Character*(*) Label
         Logical       Found
         Integer       nData
         End Subroutine Qpg_dArray
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
      Subroutine Basis_Info_Dmp()
!
      Integer i, j, nCnttp, nAtoms, nAux, nM1, nM2, nFragCoor, nAux2, nBk
      Integer, Allocatable:: iDmp(:,:)
      Real*8, Allocatable, Target:: rDmp(:,:)
      Real*8, Pointer:: qDmp(:,:)
!     Write (6,*) 'Basis_Info_Dmp()'
!
!     Temporary code until nCnttp has been move over to the Module
!
      i = 0
      Do
          i=i+1
          If (i.gt.Mxdbsc .or. dbsc(i)%nCntr.eq.0) Exit
      End Do
      nCnttp=i-1
#ifdef _DEBUG_
      Write (6,*) 'Basis_Info_Dmp'
      Do i = 1, nCnttp
         Do j = 1, dbsc(i)%nCntr
            Write (6,*) dbsc(i)%Coor(:,j)
         End Do
      End Do
#endif
!
!     Integer dbsc stuff
!
      Call mma_Allocate(iDmp,nFields,nCnttp+1,Label='iDmp')
      nAtoms=0
      nAux   = 0
      Do i = 1, nCnttp
         iDmp(1,i) = dbsc(i)%nCntr
         iDmp(2,i) = dbsc(i)%nM1
         iDmp(3,i) = dbsc(i)%nM2
         iDmp(4,i) = dbsc(i)%nFragType
         iDmp(5,i) = dbsc(i)%nFragCoor
         iDmp(6,i) = dbsc(i)%nFragEner
         iDmp(7,i) = dbsc(i)%nFragDens
         nAtoms=nAtoms+dbsc(i)%nCntr
         nFragCoor=Max(0,dbsc(i)%nFragCoor)  ! Fix the misuse in FragExpand
         nAux = nAux + 2*dbsc(i)%nM1 + 2*dbsc(i)%nM2  &
               +nFrag_LineWords*dbsc(i)%nFragType     &
               +5              *        nFragCoor     &
                               +dbsc(i)%nFragEner     &
             +dbsc(i)%nFragDens*dbsc(i)%nFragEner
#ifdef _DEBUG_
         Write (6,'(A,7I4)') 'iCnttp=',i,                     &
                nFrag_LineWords,dbsc(i)%nFragType,    &
                                dbsc(i)%nFragCoor,    &
                                dbsc(i)%nFragEner,    &
                                dbsc(i)%nFragDens, nAux
#endif
!
         If (dbsc(i)%nPAM2.ne.-1) Then
            Write (6,*) 'Not yet implemented for PAM2 integrals.'
            Call Abend()
         End If
      End Do
      iDmp(1,nCnttp+1)=nFrag_LineWords
      Call Put_iArray('iDmp',iDmp,nFields*(nCnttp+1))
      Call mma_deallocate(iDmp)
!
!     Integer shells stuffs
!
      Call mma_Allocate(iDmp,mFields,Max_Shells-1,Label='iDmp')
      nAux2=0
      Do i = 1, Max_Shells-1
         iDmp(1,i) = Shells(i)%nBK
         nAux2 = nAux2 + 2*Shells(i)%nBK
      End Do
      Call Put_iArray('iDmp:S',iDmp,mFields*(Max_Shells-1))
      Call mma_deallocate(iDmp)
!
!**********************************************************************
!
!**********************************************************************
!
      Call mma_allocate(rDmp,3,nAtoms,Label='rDmp')
      nAtoms = 0
      Do i = 1, nCnttp
!        Call RecPrt('dbsc(i)%Coor',' ',dbsc(i)%Coor(1,1),3,dbsc(i)%nCntr)
         Do j = 1, dbsc(i)%nCntr
            nAtoms=nAtoms+1
            rDmp(1:3,nAtoms)=dbsc(i)%Coor(1:3,j)
         End Do
      End Do
      Call Put_dArray('rDmp',rDmp,3*nAtoms)
      Call mma_deallocate(rDmp)
!
      If (nAux.gt.0) Then
!        Write (*,*) 'nAux=',nAux
         Call mma_allocate(rDmp,nAux,1,Label='rDmp')
         nAux=0
         Do i = 1, nCnttp
            nM1 = dbsc(i)%nM1
            If (nM1.gt.0) Then
!              Call RecPrt('M1xp',' ',dbsc(i)%M1xp,1,nM1)
!              Call RecPrt('M1cf',' ',dbsc(i)%M1cf,1,nM1)
               rDmp(nAux+1:nAux+nM1,1) = dbsc(i)%M1xp(:)
               nAux = nAux + nM1
               rDmp(nAux+1:nAux+nM1,1) = dbsc(i)%M1cf(:)
               nAux = nAux + nM1
            End If
            nM2 = dbsc(i)%nM2
            If (nM2.gt.0) Then
!              Call RecPrt('M2xp',' ',dbsc(i)%M2xp,1,nM2)
!              Call RecPrt('M2cf',' ',dbsc(i)%M2cf,1,nM2)
               rDmp(nAux+1:nAux+nM2,1) = dbsc(i)%M2xp(:)
               nAux = nAux + nM2
               rDmp(nAux+1:nAux+nM2,1) = dbsc(i)%M2cf(:)
               nAux = nAux + nM2
            End If
!
!           Write (*,*) 'iAux=',nAux
!           Write (*,*) nFrag_LineWords, dbsc(i)%nFragType
            If (dbsc(i)%nFragType.gt.0) Then
               qDmp(1:nFrag_LineWords,1:dbsc(i)%nFragType) => rDmp(nAux+1:nAux+nFrag_LineWords*dbsc(i)%nFragType,1)
               qDmp(:,:)=dbsc(i)%FragType(:,:)
               nAux = nAux + nFrag_LineWords*dbsc(i)%nFragType
               Nullify(qDmp)
            End If
!           Write (*,*) dbsc(i)%nFragCoor
            nFragCoor = Max(0,dbsc(i)%nFragCoor)
            If (        nFragCoor.gt.0) Then
               qDmp(1:5,1:        nFragCoor) => rDmp(nAux+1:nAux+5*        nFragCoor,1)
               qDmp(:,:)=dbsc(i)%FragCoor(:,:)
               nAux = nAux + 5*        nFragCoor
               Nullify(qDmp)
            End If
!           Write (*,*) dbsc(i)%nFragEner
            If (dbsc(i)%nFragEner.gt.0) Then
               rDmp(nAux+1:nAux+dbsc(i)%nFragEner,1)=dbsc(i)%FragEner(:)
               nAux = nAux + dbsc(i)%nFragEner
            End If
!           Write (*,*) dbsc(i)%nFragDens
            If (dbsc(i)%nFragDens*dbsc(i)%nFragEner.gt.0) Then
               qDmp(1:dbsc(i)%nFragDens,1:dbsc(i)%nFragEner) => rDmp(nAux+1:nAux+dbsc(i)%nFragDens*dbsc(i)%nFragEner,1)
               qDmp(:,:)=dbsc(i)%FragCoef(:,:)
               nAux = nAux + dbsc(i)%nFragDens*dbsc(i)%nFragEner
               Nullify(qDmp)
            End If
         End Do
!        Call RecPrt('rDmp:A',' ',rDmp,1,nAux)
         Call Put_dArray('rDmp:A',rDmp,nAux)
         Call mma_deallocate(rDmp)
      End If
!
      If (nAux2.gt.0) Then
         Call mma_allocate(rDmp,nAux2,1,Label='rDmp')
         nAux2=0
         Do i = 1, Max_Shells-1
            nBk=Shells(i)%nBk
            If (nBk.gt.0) Then
               rDmp(nAux2+1:nAux2+nBK,1)=Shells(i)%Bk(:)
               nAux2 = nAux2 + nBK
               rDmp(nAux2+1:nAux2+nBK,1)=Shells(i)%Occ(:)
               nAux2 = nAux2 + nBK
            End If
         End Do
         Call Put_dArray('rDmp:S',rDmp,nAux2)
         Call mma_deallocate(rDmp)
      End If
      Return
      End Subroutine Basis_Info_Dmp
!
!***********************************************************************
!***********************************************************************
!
      Subroutine Basis_Info_Get()
!
      Integer, Allocatable:: iDmp(:,:)
      Real*8, Allocatable, Target:: rDmp(:,:)
      Real*8, Pointer:: qDmp(:,:), pDmp(:)
      Logical Found
      Integer Len, i, j, nCnttp, nAtoms, nAux, nM1, nM2, nBK,nAux2
      Integer nFragType, nFragCoor, nFragEner, nFragDens
!     Write (6,*) 'Basis_Info_Get()'
!
      Call qpg_iArray('iDmp',Found,Len)
      nCnttp=Len/nFields-1
      Call mma_Allocate(iDmp,nFields,nCnttp+1,Label='iDmp')
      If (Found) Call Get_iArray('iDmp',iDmp,nFields*(nCnttp+1))
      nAux = 0
      nFrag_LineWords=iDmp(1,nCnttp+1)
      Do i = 1, nCnttp
         dbsc(i)%nCntr     = iDmp(1,i)
         dbsc(i)%nM1       = iDmp(2,i)
         dbsc(i)%nM2       = iDmp(3,i)
         dbsc(i)%nFragType = iDmp(4,i)
         dbsc(i)%nFragCoor = iDmp(5,i)
         dbsc(i)%nFragEner = iDmp(6,i)
         dbsc(i)%nFragDens = iDmp(7,i)
         nFragCoor=Max(0,dbsc(i)%nFragCoor)
         nAux = nAux + 2*dbsc(i)%nM1 + 2*dbsc(i)%nM2  &
               +nFrag_LineWords*dbsc(i)%nFragType     &
               +5              *        nFragCoor     &
                               +dbsc(i)%nFragEner     &
             +dbsc(i)%nFragDens*dbsc(i)%nFragEner
#ifdef _DEBUG_
         Write (6,'(A,7I4)') 'iCnttp=',i,                     &
                nFrag_LineWords,dbsc(i)%nFragType,    &
                                dbsc(i)%nFragCoor,    &
                                dbsc(i)%nFragEner,    &
                                dbsc(i)%nFragDens, nAux
#endif
      End Do
      Call mma_deallocate(iDmp)
!
!
      Call qpg_iArray('iDmp:S',Found,Len)
      Max_Shells = Len/mFields + 1
      Call mma_Allocate(iDmp,mFields,Max_Shells-1,Label='iDmp')
      Call get_iArray('iDmp:S',iDmp,Len)
      nAux2=0
      Do i = 1, Max_Shells-1
         Shells(i)%nBK = iDmp(1,i)
         nAux2 = nAux2 + 2*Shells(i)%nBK
      End Do
      Call mma_deallocate(iDmp)
!
      Call qpg_dArray('rDmp',Found,Len)
      If (.Not.Found) Then
         Write (6,*) 'rDMP not found on the run file.'
         Call Abend()
      End If
      nAtoms=Len/3
      Call mma_allocate(rDmp,3,nAtoms,Label='rDmp')
      Call Get_dArray('rDmp',rDmp,3*nAtoms)
      nAtoms = 0
      Do i = 1, nCnttp
         If (.Not.Allocated(dbsc(i)%Coor)) Then
            Call mma_Allocate(dbsc(i)%Coor,3,dbsc(i)%nCntr,Label='dbsc:C')
         End If
         Do j = 1, dbsc(i)%nCntr
            nAtoms=nAtoms+1
            dbsc(i)%Coor(1:3,j)=rDmp(1:3,nAtoms)
         End Do
      End Do
      Call mma_deallocate(rDmp)
!
      If (nAux.gt.0) Then
         Call qpg_dArray('rDmp:A',Found,Len)
!        Write (*,*) 'nAux=',nAux
         Call mma_allocate(rDmp,nAux,1,Label='rDmp')
         Call Get_dArray('rDmp:A',rDmp,Len)
         nAux=0
         Do i = 1, nCnttp
!
!           ECP stuff
!
            nM1 = dbsc(i)%nM1
            If (nM1.gt.0) Then
               If (.Not.Allocated(dbsc(i)%M1xp)) Call mma_allocate(dbsc(i)%M1xp,nM1,Label='dbsc:M1xp')
               dbsc(i)%M1xp(:)=rDmp(nAux+1:nAux+nM1,1)
               nAux=nAux+nM1
               If (.Not.Allocated(dbsc(i)%M1cf)) Call mma_allocate(dbsc(i)%M1cf,nM1,Label='dbsc:M1cf')
               dbsc(i)%M1cf(:)=rDmp(nAux+1:nAux+nM1,1)
               nAux=nAux+nM1
!              Call RecPrt('M1xp',' ',dbsc(i)%M1xp,1,nM1)
!              Call RecPrt('M1cf',' ',dbsc(i)%M1cf,1,nM1)
            End If
            nM2 = dbsc(i)%nM2
            If (nM2.gt.0) Then
               If (.Not.Allocated(dbsc(i)%M2xp)) Call mma_allocate(dbsc(i)%M2xp,nM2,Label='dbsc:M2xp')
               dbsc(i)%M2xp(:)=rDmp(nAux+1:nAux+nM2,1)
               nAux=nAux+nM2
               If (.Not.Allocated(dbsc(i)%M2cf)) Call mma_allocate(dbsc(i)%M2cf,nM2,Label='dbsc:M2cf')
               dbsc(i)%M2cf(:)=rDmp(nAux+1:nAux+nM2,1)
               nAux=nAux+nM2
!              Call RecPrt('M2xp',' ',dbsc(i)%M2xp,1,nM2)
!              Call RecPrt('M2cf',' ',dbsc(i)%M2cf,1,nM2)
            End If
!
!           Fragment stuff
!
            nFragType  =dbsc(i)%nFragType
            If (nFragType.gt.0) Then
               If (.Not.Allocated(dbsc(i)%FragType)) Call mma_allocate(dbsc(i)%FragType,nFrag_LineWords,nFragType,Label='FragType')
               qDmp(1:nFrag_LineWords,1:nFragType) => rDmp(nAux+1:nAux+nFrag_LineWords*nFragType,1)
               dbsc(i)%FragType(:,:)=qDmp(:,:)
               nAux=nAux+nFrag_LineWords*nFragType
               Nullify(qDmp)
            End If
            nFragCoor  =Max(0,dbsc(i)%nFragCoor)
            If (nFragCoor.gt.0) Then
               If (.Not.Allocated(dbsc(i)%FragCoor)) Call mma_allocate(dbsc(i)%FragCoor,5,nFragCoor,Label='FragCoor')
               qDmp(1:5,1:nFragCoor) => rDmp(nAux+1:nAux+5*nFragCoor,1)
               dbsc(i)%FragCoor(:,:)=qDmp(:,:)
               nAux=nAux+5*nFragCoor
               Nullify(qDmp)
            End If
            nFragEner  =dbsc(i)%nFragEner
            If (nFragEner.gt.0) Then
               If (.Not.Allocated(dbsc(i)%FragEner)) Call mma_allocate(dbsc(i)%FragEner,nFragEner,Label='FragEner')
               pDmp(1:nFragEner) => rDmp(nAux+1:nAux+nFragEner,1)
               dbsc(i)%FragEner(:)=pDmp(:)
               nAux=nAux+nFragEner
               Nullify(pDmp)
            End If
            nFragDens  =dbsc(i)%nFragDens
            If (nFragDens*nFragEner.gt.0) Then
               If (.Not.Allocated(dbsc(i)%FragCoef)) Call mma_allocate(dbsc(i)%FragCoef,nFragDens,nFragEner,Label='FragCoef')
               qDmp(1:nFragDens,1:nFragEner) => rDmp(nAux+1:nAux+nFragDens*nFragEner,1)
               dbsc(i)%FragCoef(:,:)=qDmp(:,:)
               nAux=nAux+nFragDens*nFragEner
               Nullify(qDmp)
            End If
         End Do
         Call mma_deallocate(rDmp)
      End If
!
      If (nAux2.gt.0) Then
         Call qpg_dArray('rDmp:S',Found,Len)
         Call mma_allocate(rDmp,nAux2,1,Label='rDmp')
         Call Get_dArray('rDmp:S',rDmp,Len)
         nAux2=0
         Do i = 1, Max_Shells-1
            nBk=SHells(i)%nBK
            If (nBk.gt.0) Then
               If (.Not.Allocated(Shells(i)%Bk)) Call mma_allocate(Shells(i)%Bk,nBk,Label='Bk')
               Shells(i)%Bk(:)=rDmp(nAux2+1:nAux2+nBk,1)
               nAux2=nAux2+nBk
               If (.Not.Allocated(Shells(i)%Occ)) Call mma_allocate(Shells(i)%Occ,nBk,Label='Occ')
               Shells(i)%Occ(:)=rDmp(nAux2+1:nAux2+nBk,1)
               nAux2=nAux2+nBk
            End If
         End Do
         Call mma_deallocate(rDmp)
      End If
#ifdef _DEBUG_
      Write (6,*) 'Basis_Info_Get'
      Do i = 1, nCnttp
         Do j = 1, dbsc(i)%nCntr
            Write (6,*) dbsc(i)%Coor(:,j)
         End Do
      End Do
#endif
      Return
      End Subroutine Basis_Info_Get
!
!***********************************************************************
!***********************************************************************
!
      Subroutine Basis_Info_Free()
      Integer i
!     Write (6,*) 'Basis_Info_Free()'
!
!     Deallocate all allocatable parts of dbsc.
!
      i = 0
      Do
         i=i+1
         If (i.gt.Mxdbsc .or. dbsc(i)%nCntr.eq.0) Exit
!
!        Molecular Coordinates
!
         If (allocated(dbsc(i)%Coor)) Call mma_deallocate(dbsc(i)%Coor)
         dbsc(i)%nCntr=-1
!
!        ECP stuff
!
         If (allocated(dbsc(i)%M1xp)) Call mma_deallocate(dbsc(i)%M1xp)
         If (allocated(dbsc(i)%M1cf)) Call mma_deallocate(dbsc(i)%M1cf)
         dbsc(i)%nM1=0
         If (allocated(dbsc(i)%M2xp)) Call mma_deallocate(dbsc(i)%M2xp)
         If (allocated(dbsc(i)%M2cf)) Call mma_deallocate(dbsc(i)%M2cf)
         dbsc(i)%nM2=0
!
!        Fragment stuff
!
         If (allocated(dbsc(i)%FragType)) Call mma_deallocate(dbsc(i)%FragType)
         dbsc(i)%nFragType=0
         If (allocated(dbsc(i)%FragCoor)) Call mma_deallocate(dbsc(i)%FragCoor)
         dbsc(i)%nFragCoor=0
         If (allocated(dbsc(i)%FragEner)) Call mma_deallocate(dbsc(i)%FragEner)
         dbsc(i)%nFragEner=0
         If (allocated(dbsc(i)%FragCoef)) Call mma_deallocate(dbsc(i)%FragCoef)
         dbsc(i)%nFragDens=0
!
!        PAM2 stuff
!
         If (allocated(dbsc(i)%PAM2)) Call mma_deallocate(dbsc(i)%PAM2)
         dbsc(i)%nPAM2=-1
      End Do
!
      Do i = 1, Max_Shells-1
         If (Allocated(Shells(i)%Bk)) Call mma_deallocate(Shells(i)%Bk)
         If (Allocated(Shells(i)%Occ)) Call mma_deallocate(Shells(i)%Occ)
      End Do
!
      Return
      End Subroutine Basis_Info_Free
!
!***********************************************************************
!***********************************************************************
!
      End Module Basis_Info
