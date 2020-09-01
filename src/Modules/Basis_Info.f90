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
Public :: Basis_Info_Dmp, Basis_Info_Get, Basis_Info_Free, Distinct_Basis_set_Centers, dbsc, nFrag_LineWords,&
          PAMExp, Shells, Max_Shells, nCnttp, iCnttp_Dummy, Point_Charge, Gaussian_type, mGaussian_Type,     &
          Nuclear_Model, Basis_Info_Init

#include "stdalloc.fh"
#include "Molcas.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    D E C L A R E   D E R I V E D   T Y P E S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!     ECP      : Flag if dbsc is a ECP basis set
!     Frag     : Flag if dbsc is a Fragment basis set
!     Aux      : Flag if dbsc is an auxiliary basis set
!     FOp      : Flag if dbsc has a Fock Operator
!     IsMM     : integer flag to indicate that associated centers are treated as MM centers in QM/MM calculations
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
    Logical:: ECP=.False.
    Logical:: Aux =.False.
    Logical:: Frag=.False.
    Logical:: FOp =.False.
    Integer:: IsMM=0
    Integer:: Parent_iCnttp=0
    Integer:: lOffAO=0
    Integer:: nOpt=0
    Integer:: mdci=0
    Integer:: iVal=0, nVal=0
    Integer:: iPrj=0, nPrj=0
    Integer:: iSRO=0, nSRO=0
    Integer:: iSOC=0, nSOC=0
    Integer:: iPP =0, nPP =0
    Integer:: nShells =0
    Integer:: AtmNr=0
    Real*8:: Charge=0.0D0
    Logical:: NoPair=.False.
End Type Distinct_Basis_set_centers
!
!     nExp  : number of exponents of the i''th shell
!     Exp   : the exponents of the i''th shell
!     nBasis: number of contracted radial functions of the i''th shell
!     Cff_c : Contraction coefficients in processed and raw input form
!     Cff_p : Contraction coefficient in the case of no contraction, processed and raw
!     Cff   : copy of Cff_c or Cff_p
!     Transf: Cartesian transformed to real sphericals.
!     Projct: real sphericals without contaminations (3s, 4d, etc.)
!     Bk    : ECP proj shift parameters for i''th shell.
!             the number of parameters is given by nBasis
!     Occ   : Occupation numbers for core ECP orbitals
!     FockOp: the Fock operator matrix
!     Aux   : Logical flag for auxiliary basis set shells
!     Frag  : Logical flag for fragment shells
!
Type Shell_Info
     Sequence
     Integer :: nExp=0
     Real*8, Allocatable:: Exp(:)
     Integer :: nBasis=0
     Integer :: nBasis_c=0
     Real*8, Allocatable:: pCff(:,:)
     Real*8, Allocatable:: Cff_c(:,:,:), Cff_p(:,:,:)
     Logical :: Transf=.True.
     Logical :: Prjct =.True.
     Integer :: nBk=0
     Real*8, Allocatable:: Bk(:)
     Real*8, Allocatable:: Occ(:)
     Integer :: nAkl=0
     Real*8, Allocatable:: Akl(:,:,:)
     Integer :: nFockOp=0
     Real*8, Allocatable:: FockOp(:,:)
     Logical :: Aux =.False.
     Logical:: Frag=.False.
     Integer:: kOffAO=0
End Type Shell_Info
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     E N D   D E C L A R E   D E R I V E D   T Y P E S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Actual content of Basis_Info
!
Integer, Parameter :: Point_Charge  = 0
Integer, Parameter :: Gaussian_type = 1
Integer, Parameter :: mGaussian_Type= 2

Real*8, Allocatable:: PAMexp(:,:)
Integer :: nFrag_LineWords = 0, nFields =29, mFields = 11
Integer :: nCnttp = 0, iCnttp_Dummy = 0
Integer :: Max_Shells = 0
Logical :: Initiated = .FALSE.
Integer :: Nuclear_Model=Point_Charge

Type (Distinct_Basis_set_centers) , Allocatable:: dbsc(:)
Type (Shell_Info), Allocatable :: Shells(:)

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
!     This to make either the initial allocation of dbsc and Shells according to the default sizes
!     as defined by the parameters in Molcas.fh or according to the actual sizes as recorded on the
!     run file.
!
Subroutine Basis_Info_Init()
If (Initiated) Return
If (nCnttp.eq.0) Then
   Allocate(dbsc(1:Mxdbsc))
Else
   Allocate(dbsc(1:nCnttp))
End If
If (Max_Shells.eq.0) Then
   Allocate(Shells(1:MxShll))
Else
   Allocate(Shells(1:Max_Shells))
End If
Initiated=.True.
Return
End Subroutine Basis_Info_Init
!
!***********************************************************************
!***********************************************************************
!
Subroutine Basis_Info_Dmp()
!
Integer i, j, nAtoms, nAux, nM1, nM2, nFragCoor, nAux2, nBk, nAkl, nFockOp, nExp, nBasis
Integer, Allocatable:: iDmp(:,:)
Real*8, Allocatable, Target:: rDmp(:,:)
Real*8, Pointer:: qDmp(:,:)
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
   iDmp(8,i) = 0
   If (dbsc(i)%ECP) iDmp(8,i)=1
   iDmp(9,i) = 0
   If (dbsc(i)%Frag)iDmp(9,i)=1
   iDmp(10,i) = 0
   If (dbsc(i)%Aux )iDmp(10,i)=1
   iDmp(11,i) = 0
   If (dbsc(i)%FOp )iDmp(11,i)=1
   iDmp(12,i) = dbsc(i)%IsMM
   iDmp(13,i) = dbsc(i)%Parent_iCnttp
   iDmp(14,i) = dbsc(i)%lOffAO
   iDmp(15,i) = dbsc(i)%nOpt
   iDmp(16,i) = dbsc(i)%mdci
   iDmp(17,i) = dbsc(i)%iVal
   iDmp(18,i) = dbsc(i)%nVal
   iDmp(19,i) = dbsc(i)%iPrj
   iDmp(20,i) = dbsc(i)%nPrj
   iDmp(21,i) = dbsc(i)%iSRO
   iDmp(22,i) = dbsc(i)%nSRO
   iDmp(23,i) = dbsc(i)%iSOC
   iDmp(24,i) = dbsc(i)%nSOC
   iDmp(25,i) = dbsc(i)%iPP
   iDmp(26,i) = dbsc(i)%nPP
   iDmp(27,i) = dbsc(i)%nShells
   iDmp(28,i) = dbsc(i)%AtmNr
   iDmp(29,i) = 0
   If (dbsc(i)%NoPair)iDmp(29,i)=1
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
iDmp(2,nCnttp+1)=nCnttp
iDmp(3,nCnttp+1)=iCnttp_Dummy
iDmp(4,nCnttp+1)=Max_Shells
iDmp(5,nCnttp+1)=Nuclear_Model
Call Put_iArray('iDmp',iDmp,nFields*(nCnttp+1))
Call mma_deallocate(iDmp)
!
!     Integer shells stuffs
!
Call mma_Allocate(iDmp,mFields,Max_Shells-1,Label='iDmp')
nAux2=0
Do i = 1, Max_Shells-1
   iDmp(1,i) = Shells(i)%nBK
   iDmp(2,i) = Shells(i)%nAkl
   iDmp(3,i) = Shells(i)%nFockOp
   iDmp(4,i) = Shells(i)%nExp
   iDmp(5,i) = Shells(i)%nBasis
   iDmp(6,i) = Shells(i)%nBasis_c
   iDmp(7,i) = 0
   If (Shells(i)%Transf) iDmp(7,i)=1
   iDmp(8,i) = 0
   If (Shells(i)%Prjct ) iDmp(8,i)=1
   iDmp(9,i) = 0
   If (Shells(i)%Frag  ) iDmp(9,i)=1
   iDmp(10,i) = 0
   If (Shells(i)%Aux   ) iDmp(10,i)=1
   iDmp(11,i) = Shells(i)%kOffAO
   nAux2 = nAux2 + 2*Shells(i)%nBK + 2*Shells(i)%nAkl**2 + Shells(i)%nFockOp**2  &
         + Shells(i)%nExp + 2*Shells(i)%nExp*Shells(i)%nBasis + 2*Shells(i)%nExp**2
#ifdef _DEBUG_
   Write (6,'(A,7I4)') 'iShll=',i,                     &
               Shells(i)%nBK,                          &
               Shells(i)%nAkl,                         &
               Shells(i)%nFockOp,                      &
               Shells(i)%nExp                          &
               Shells(i)%nBasis
#endif

End Do
Call Put_iArray('iDmp:S',iDmp,mFields*(Max_Shells-1))
Call mma_deallocate(iDmp)
!
!**********************************************************************
!
!**********************************************************************
!
Call mma_allocate(rDmp,3,nAtoms+nCnttp,Label='rDmp')
nAtoms = 0
Do i = 1, nCnttp
!  Call RecPrt('dbsc(i)%Coor',' ',dbsc(i)%Coor(1,1),3,dbsc(i)%nCntr)
   Do j = 1, dbsc(i)%nCntr
      nAtoms=nAtoms+1
      rDmp(1:3,nAtoms)=dbsc(i)%Coor(1:3,j)
   End Do
   nAtoms=nAtoms+1
   rDmp(1,nAtoms)=dbsc(i)%Charge
   rDmp(2,nAtoms)=0.0D0
   rDmp(3,nAtoms)=0.0D0
End Do
Call Put_dArray('rDmp',rDmp,3*nAtoms)
Call mma_deallocate(rDmp)
!
If (nAux.gt.0) Then
!  Write (*,*) 'nAux=',nAux
   Call mma_allocate(rDmp,nAux,1,Label='rDmp')
   nAux=0
   Do i = 1, nCnttp
      nM1 = dbsc(i)%nM1
      If (nM1.gt.0) Then
         Call RecPrt('M1xp',' ',dbsc(i)%M1xp,1,nM1)
!        Call RecPrt('M1cf',' ',dbsc(i)%M1cf,1,nM1)
         rDmp(nAux+1:nAux+nM1,1) = dbsc(i)%M1xp(:)
         nAux = nAux + nM1
         rDmp(nAux+1:nAux+nM1,1) = dbsc(i)%M1cf(:)
         nAux = nAux + nM1
      End If
      nM2 = dbsc(i)%nM2
      If (nM2.gt.0) Then
!        Call RecPrt('M2xp',' ',dbsc(i)%M2xp,1,nM2)
!        Call RecPrt('M2cf',' ',dbsc(i)%M2cf,1,nM2)
         rDmp(nAux+1:nAux+nM2,1) = dbsc(i)%M2xp(:)
         nAux = nAux + nM2
         rDmp(nAux+1:nAux+nM2,1) = dbsc(i)%M2cf(:)
         nAux = nAux + nM2
      End If
!
!     Write (*,*) 'iAux=',nAux
!     Write (*,*) nFrag_LineWords, dbsc(i)%nFragType
      If (dbsc(i)%nFragType.gt.0) Then
         qDmp(1:nFrag_LineWords,1:dbsc(i)%nFragType) => rDmp(nAux+1:nAux+nFrag_LineWords*dbsc(i)%nFragType,1)
         qDmp(:,:)=dbsc(i)%FragType(:,:)
         nAux = nAux + nFrag_LineWords*dbsc(i)%nFragType
         Nullify(qDmp)
      End If
!     Write (*,*) dbsc(i)%nFragCoor
      nFragCoor = Max(0,dbsc(i)%nFragCoor)
      If (        nFragCoor.gt.0) Then
         qDmp(1:5,1:        nFragCoor) => rDmp(nAux+1:nAux+5*        nFragCoor,1)
         qDmp(:,:)=dbsc(i)%FragCoor(:,:)
         nAux = nAux + 5*        nFragCoor
         Nullify(qDmp)
      End If
!     Write (*,*) dbsc(i)%nFragEner
      If (dbsc(i)%nFragEner.gt.0) Then
         rDmp(nAux+1:nAux+dbsc(i)%nFragEner,1)=dbsc(i)%FragEner(:)
         nAux = nAux + dbsc(i)%nFragEner
      End If
!     Write (*,*) dbsc(i)%nFragDens
      If (dbsc(i)%nFragDens*dbsc(i)%nFragEner.gt.0) Then
         qDmp(1:dbsc(i)%nFragDens,1:dbsc(i)%nFragEner) => rDmp(nAux+1:nAux+dbsc(i)%nFragDens*dbsc(i)%nFragEner,1)
         qDmp(:,:)=dbsc(i)%FragCoef(:,:)
         nAux = nAux + dbsc(i)%nFragDens*dbsc(i)%nFragEner
         Nullify(qDmp)
      End If
   End Do
!  Call RecPrt('rDmp:A',' ',rDmp,1,nAux)
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
      nAkl=Shells(i)%nAkl
      If (nAkl.gt.0) Then
         Call DCopy_(2*nAkl**2,Shells(i)%Akl,1,rDmp(nAux2+1,1),1)
         nAux2 = nAux2 + 2*nAkl**2
      End If
      nFockOp=Shells(i)%nFockOp
      If (nFockOp.gt.0) Then
         Call DCopy_(nFockOp**2,Shells(i)%FockOp,1,rDmp(nAux2+1,1),1)
         nAux2 = nAux2 + nFockOp**2
      End If
      nExp=Shells(i)%nExp
      If (nExp.gt.0) Then
         Call DCopy_(nExp,Shells(i)%Exp,1,rDmp(nAux2+1,1),1)
         nAux2 = nAux2 + nExp
      End If
      nBasis=Shells(i)%nBasis
!     Note: the contraction coefficients are not always there.
      If (nExp*nBasis.gt.0) Then
         Call DCopy_(2*nExp**2,Shells(i)%Cff_p,1,rDmp(nAux2+1,1),1)
         nAux2 = nAux2 + 2*nExp**2
      End If
      If (nExp*nBasis.gt.0) Then
         Call DCopy_(2*nExp*nBasis,Shells(i)%Cff_c,1,rDmp(nAux2+1,1),1)
         nAux2 = nAux2 + 2*nExp*nBasis
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
Integer Len, i, j, nAtoms, nAux, nM1, nM2, nBK,nAux2, nAkl, nFockOp, nExp, nBasis, Len2
Integer nFragType, nFragCoor, nFragEner, nFragDens
#ifdef _DEBUG_
Write (6,*) 'Basis_Info_Get()'
#endif
Call qpg_iArray('iDmp',Found,Len)
Len2=Len/nFields
Call mma_Allocate(iDmp,nFields,Len2,Label='iDmp')
If (Found) Then
   Call Get_iArray('iDmp',iDmp,Len)
Else
   Write (6,*) 'Basis_Info_Get: iDmp not found!'
   Call Abend()
End If
nFrag_LineWords=iDmp(1,Len2)
nCnttp         =iDmp(2,Len2)
iCnttp_Dummy   =iDmp(3,Len2)
Max_Shells     =iDmp(4,Len2)
Nuclear_Model  =iDmp(5,Len2)
nAux = 0
!
!     Initiate the memory allocation of dsbc and Shells
!
If (.Not.Initiated) Call Basis_Info_Init()
!
Do i = 1, nCnttp
   dbsc(i)%nCntr        = iDmp( 1,i)
   dbsc(i)%nM1          = iDmp( 2,i)
   dbsc(i)%nM2          = iDmp( 3,i)
   dbsc(i)%nFragType    = iDmp( 4,i)
   dbsc(i)%nFragCoor    = iDmp( 5,i)
   dbsc(i)%nFragEner    = iDmp( 6,i)
   dbsc(i)%nFragDens    = iDmp( 7,i)
   dbsc(i)%ECP          = iDmp( 8,i).eq.1
   dbsc(i)%Frag         = iDmp( 9,i).eq.1
   dbsc(i)%Aux          = iDmp(10,i).eq.1
   dbsc(i)%FOp          = iDmp(11,i).eq.1
   dbsc(i)%IsMM         = iDmp(12,i)
   dbsc(i)%Parent_iCnttp= iDmp(13,i)
   dbsc(i)%lOffAO       = iDmp(14,i)
   dbsc(i)%nOpt         = iDmp(15,i)
   dbsc(i)%mdci         = iDmp(16,i)
   dbsc(i)%iVal         = iDmp(17,i)
   dbsc(i)%nVal         = iDmp(18,i)
   dbsc(i)%iPrj         = iDmp(19,i)
   dbsc(i)%nPrj         = iDmp(20,i)
   dbsc(i)%iSRO         = iDmp(21,i)
   dbsc(i)%nSRO         = iDmp(22,i)
   dbsc(i)%iSOC         = iDmp(23,i)
   dbsc(i)%nSOC         = iDmp(24,i)
   dbsc(i)%iPP          = iDmp(25,i)
   dbsc(i)%nPP          = iDmp(26,i)
   dbsc(i)%nShells      = iDmp(27,i)
   dbsc(i)%AtmNr        = iDmp(28,i)
   dbsc(i)%NoPair       = iDmp(29,i).eq.1
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
Call mma_Allocate(iDmp,mFields,Len/mFields,Label='iDmp')
If (Found) Then
   Call get_iArray('iDmp:S',iDmp,Len)
Else
   Write (6,*) 'Basis_Info_Get: iDmp:S not found!'
   Call Abend()
End If
nAux2=0
Do i = 1, Max_Shells-1
   Shells(i)%nBK      = iDmp( 1,i)
   Shells(i)%nAkl     = iDmp( 2,i)
   Shells(i)%nFockOp  = iDmp( 3,i)
   Shells(i)%nExp     = iDmp( 4,i)
   Shells(i)%nBasis   = iDmp( 5,i)
   Shells(i)%nBasis_c = iDmp( 6,i)
   Shells(i)%Transf   = iDmp( 7,i).eq.1
   Shells(i)%Prjct    = iDmp( 8,i).eq.1
   Shells(i)%Frag     = iDmp( 9,i).eq.1
   Shells(i)%Aux      = iDmp(10,i).eq.1
   Shells(i)%kOffAO   = iDmp(11,i)
   nAux2 = nAux2 + 2*Shells(i)%nBK + 2*Shells(i)%nAkl**2 + Shells(i)%nFockOp**2  &
         + Shells(i)%nExp
!  Coefficients only there is nBasis =/=0
   If (Shells(i)%nBasis.gt.0) Then
      nAux2 = nAux2 + 2*Shells(i)%nExp*Shells(i)%nBasis + 2*Shells(i)%nExp**2
   End If
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
   nAtoms=nAtoms+1
   dbsc(i)%Charge    =rDmp(1,nAtoms)
End Do
Call mma_deallocate(rDmp)
!
If (nAux.gt.0) Then
   Call qpg_dArray('rDmp:A',Found,Len)
!  Write (*,*) 'nAux=',nAux
   Call mma_allocate(rDmp,nAux,1,Label='rDmp')
   Call Get_dArray('rDmp:A',rDmp,Len)
   nAux=0
   Do i = 1, nCnttp
!
!     ECP stuff
!
      nM1 = dbsc(i)%nM1
      If (nM1.gt.0) Then
         If (.Not.Allocated(dbsc(i)%M1xp)) Call mma_allocate(dbsc(i)%M1xp,nM1,Label='dbsc:M1xp')
         dbsc(i)%M1xp(:)=rDmp(nAux+1:nAux+nM1,1)
         nAux=nAux+nM1
         If (.Not.Allocated(dbsc(i)%M1cf)) Call mma_allocate(dbsc(i)%M1cf,nM1,Label='dbsc:M1cf')
         dbsc(i)%M1cf(:)=rDmp(nAux+1:nAux+nM1,1)
         nAux=nAux+nM1
!        Call RecPrt('M1xp',' ',dbsc(i)%M1xp,1,nM1)
!        Call RecPrt('M1cf',' ',dbsc(i)%M1cf,1,nM1)
      End If
      nM2 = dbsc(i)%nM2
      If (nM2.gt.0) Then
         If (.Not.Allocated(dbsc(i)%M2xp)) Call mma_allocate(dbsc(i)%M2xp,nM2,Label='dbsc:M2xp')
         dbsc(i)%M2xp(:)=rDmp(nAux+1:nAux+nM2,1)
         nAux=nAux+nM2
         If (.Not.Allocated(dbsc(i)%M2cf)) Call mma_allocate(dbsc(i)%M2cf,nM2,Label='dbsc:M2cf')
         dbsc(i)%M2cf(:)=rDmp(nAux+1:nAux+nM2,1)
         nAux=nAux+nM2
!        Call RecPrt('M2xp',' ',dbsc(i)%M2xp,1,nM2)
!        Call RecPrt('M2cf',' ',dbsc(i)%M2cf,1,nM2)
      End If
!
!     Fragment stuff
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

      nBk=Shells(i)%nBK
      If (nBk.gt.0) Then
         If (.Not.Allocated(Shells(i)%Bk)) Call mma_allocate(Shells(i)%Bk,nBk,Label='Bk')
         Shells(i)%Bk(:)=rDmp(nAux2+1:nAux2+nBk,1)
         nAux2=nAux2+nBk
         If (.Not.Allocated(Shells(i)%Occ)) Call mma_allocate(Shells(i)%Occ,nBk,Label='Occ')
         Shells(i)%Occ(:)=rDmp(nAux2+1:nAux2+nBk,1)
         nAux2=nAux2+nBk
      End If

      nAkl=Shells(i)%nAkl
      If (nAkl.gt.0) Then
         If (.Not.Allocated(Shells(i)%Akl)) Call mma_allocate(Shells(i)%Akl,nAkl,nAkl,2,Label='Akl')
         Call DCopy_(2*nAkl**2,rDmp(nAux2+1,1),1,Shells(i)%Akl,1)
         nAux2=nAux2+2*nAkl**2
      End If

      nFockOp=Shells(i)%nFockOp
      If (nFockOp.gt.0) Then
         If (.Not.Allocated(Shells(i)%FockOp)) Call mma_allocate(Shells(i)%FockOp,nFockOp,nFockOp,Label='FockOp')
         Call DCopy_(nFockOp**2,rDmp(nAux2+1,1),1,Shells(i)%FockOp,1)
         nAux2=nAux2+nFockOp**2
      End If

      nExp=Shells(i)%nExp
      If (nExp.gt.0) Then
         If (.Not.Allocated(Shells(i)%Exp)) Call mma_allocate(Shells(i)%Exp,nExp,Label='Exp')
         Call DCopy_(nExp,rDmp(nAux2+1,1),1,Shells(i)%Exp,1)
         nAux2=nAux2+nExp
      End If

      nBasis=Shells(i)%nBasis
      If (nExp*nBasis.gt.0) Then
         If (.Not.Allocated(Shells(i)%Cff_p)) Call mma_allocate(Shells(i)%Cff_p,nExp,nExp,2,Label='Cff_p')
         Call DCopy_(2*nExp**2,rDmp(nAux2+1,1),1,Shells(i)%Cff_p,1)
         nAux2=nAux2+2*nExp**2
      End If

      If (nExp*nBasis.gt.0) Then
         If (.Not.Allocated(Shells(i)%Cff_c)) Call mma_allocate(Shells(i)%Cff_c,nExp,nBasis,2,Label='Cff_c')
         Call DCopy_(2*nExp*nBasis,rDmp(nAux2+1,1),1,Shells(i)%Cff_c,1)
         nAux2=nAux2+2*nExp*nBasis
!
         If (.Not.Allocated(Shells(i)%pCff)) Call mma_allocate(Shells(i)%pCff,nExp,nBasis,Label='Cff')
         Shells(i)%pCff(:,:)=Shells(i)%Cff_c(:,:,1)
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
#ifdef _DEBUG_
Write (6,*) 'Basis_Info_Free()'
#endif
!
!     Deallocate all allocatable parts of dbsc.
!
Do i = 1, nCnttp
!
!  Molecular Coordinates
!
   If (allocated(dbsc(i)%Coor)) Call mma_deallocate(dbsc(i)%Coor)
   dbsc(i)%nCntr=-1
!
!  ECP stuff
!
   If (allocated(dbsc(i)%M1xp)) Call mma_deallocate(dbsc(i)%M1xp)
   If (allocated(dbsc(i)%M1cf)) Call mma_deallocate(dbsc(i)%M1cf)
   dbsc(i)%nM1=0
   If (allocated(dbsc(i)%M2xp)) Call mma_deallocate(dbsc(i)%M2xp)
   If (allocated(dbsc(i)%M2cf)) Call mma_deallocate(dbsc(i)%M2cf)
   dbsc(i)%nM2=0
!
!  Fragment stuff
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
!  PAM2 stuff
!
   If (allocated(dbsc(i)%PAM2)) Call mma_deallocate(dbsc(i)%PAM2)
   dbsc(i)%nPAM2=-1
End Do
nCnttp=0
iCnttp_Dummy=0
!
!     Stuff on unqiue basis set shells
!
Do i = 1, Max_Shells-1
   If (Allocated(Shells(i)%Bk)) Call mma_deallocate(Shells(i)%Bk)
   If (Allocated(Shells(i)%Occ)) Call mma_deallocate(Shells(i)%Occ)
   Shells(i)%nBk =0
   If (Allocated(Shells(i)%Akl)) Call mma_deallocate(Shells(i)%Akl)
   Shells(i)%nAkl=0
   If (Allocated(Shells(i)%FockOp)) Call mma_deallocate(Shells(i)%FockOp)
   Shells(i)%nFockOp=0
   If (Allocated(Shells(i)%Exp)) Call mma_deallocate(Shells(i)%Exp)
   Shells(i)%nExp=0
   If (Allocated(Shells(i)%pCff)) Call mma_deallocate(Shells(i)%pCff)
   If (Allocated(Shells(i)%Cff_c)) Call mma_deallocate(Shells(i)%Cff_c)
   If (Allocated(Shells(i)%Cff_p)) Call mma_deallocate(Shells(i)%Cff_p)
   Shells(i)%nBasis=0
   Shells(i)%Transf=.True.
End Do
Max_Shells=0
!
If (Allocated(dbsc))   Deallocate(dbsc)
If (Allocated(Shells)) Deallocate(Shells)
Initiated=.False.
!
Return
End Subroutine Basis_Info_Free
!
!***********************************************************************
!***********************************************************************
!
End Module Basis_Info
