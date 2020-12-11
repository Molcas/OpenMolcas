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
Module Slapaf_Info
implicit none
Private
Public:: Cx, Gx, Gx0, NAC, Q_nuclear, dMass, Coor, Grd, ANr, Weights, Shift, GNrm, Lambda, &
         Energy, Energy0, DipM, MF, qInt, dqInt, nSup, Atom, RefGeo, BMx, Degen, jStab, nStab, iCoSet, &
         BM, dBM, iBM, idBM, nqBM, R12, GradRef, KtB, AtomLbl, Smmtrc, &
         Free_Slapaf, Get_Slapaf, Dmp_Slapaf
!
! Arrays always allocated
!
Real*8, Allocatable:: Cx(:,:,:)     ! list of Cartesian coordinates
Real*8, Allocatable:: Gx(:,:,:)     ! list of Cartesian Gradients, State 1
Real*8, Allocatable:: Gx0(:,:,:)    ! list of Cartesian Gradients, State 2 for optimization of conical intersections
Real*8, Allocatable:: Q_nuclear(:)  ! list nuclear charges
Real*8, Allocatable:: dmass(:)      ! list atomic mass in units of (C=12)
Real*8, Allocatable:: Coor(:,:)     ! Cartesian coordinates of the last iteraction
Real*8, Allocatable:: Grd(:,:)      ! gradient of the last iteraction in Cartesian coordinates
Real*8, Allocatable:: Weights(:)    ! list of weights of ALL centers, however, the symmetry unique are first.
Real*8, Allocatable:: Shift(:,:)    ! list of displacements in Cartesian coordinates
Real*8, Allocatable:: GNrm(:)       ! list of the gradient norm for each iteration
Real*8, Allocatable:: Energy(:)     ! list of the energies of each iteration, State 1
Real*8, Allocatable:: Energy0(:)    ! list of the energies of each iteration, State 2 for optimization of conical intersections
Real*8, Allocatable:: MF(:,:)       ! list of Cartesian mode following vectors for each iteration
Real*8, Allocatable:: DipM(:,:)     ! list of dipole moments for each iteration
Real*8, Allocatable:: qInt(:,:)     ! internal coordinates for each iteration
Real*8, Allocatable:: dqInt(:,:)    ! derivatives of internal coordinates for each iteration
Real*8, Allocatable, Target:: RefGeo(:,:)   ! Reference geometry in Cartesian coordinates
Real*8, Allocatable, Target:: R12(:,:)      ! Reference geometry in R-P calculation (not used right now)
Real*8, Allocatable, Target:: GradRef(:,:)  ! Reference gradient
Real*8, Allocatable, Target:: Bmx(:,:)      ! the B matrix
Real*8, Allocatable:: Degen(:,:)    ! list of degeneracy numbers of the unique atoms (three identical entries)
Integer, Allocatable:: jStab(:,:), iCoset(:,:), nStab(:) ! stabilizer and cosets information for the indivudual centers
#include "LenIn.fh"
Character(LEN=LENIN), Allocatable:: AtomLbl(:) ! atomic labels
Logical, Allocatable:: Smmtrc(:,:)    ! Array with logical symmetry information on if a Cartesian is symmetric or not.

! Arrays for automatic internal coordinates
Real*8, Allocatable:: BM(:)         ! ...
Real*8, Allocatable:: dBM(:)        ! ...
Integer, Allocatable:: iBM(:)       ! ...
Integer, Allocatable:: idBM(:)      ! ...
Integer, Allocatable:: nqBM(:)      ! ...

Integer, Allocatable:: ANr(:)       ! list of atomic numbers
!
! Arrays optionally allocated
!
Real*8, Allocatable:: NAC(:,:)      ! list of Cartesian non-adiabatic coupling vector
Real*8, Allocatable:: Lambda(:,:)   ! list of the Lagrange multipiers
!
! Utility arrays with explicit deallocation, i.e. not via Free_Slapaf()
!
Integer, Allocatable:: Atom(:)      ! Temporary arrays for the super symmetry case
Integer, Allocatable:: NSup(:)      ! Temporary arrays for the super symmetry case
Real*8, Allocatable::  KtB(:,:)     ! KtB array for the BMtrx family of subroutines

Logical:: Initiated=.False.
Integer nsAtom
Contains
  Subroutine Free_Slapaf()
#include "stdalloc.fh"
  If (Allocated(Energy)) Call mma_deallocate(Energy)
  If (Allocated(Energy0)) Call mma_deallocate(Energy0)
  If (Allocated(DipM)) Call mma_deallocate(DipM)
  If (Allocated(GNrm)) Call mma_deallocate(GNrm)
  If (Allocated(Cx)) Call mma_deallocate(Cx)
  If (Allocated(Gx)) Call mma_deallocate(Gx)
  If (Allocated(Gx0)) Call mma_deallocate(Gx0)
  If (Allocated(MF)) Call mma_deallocate(MF)
  If (Allocated(Lambda)) Call mma_deallocate(Lambda)
  If (Allocated(Degen)) Call mma_deallocate(Degen)
  If (Allocated(jStab)) Call mma_deallocate(jStab)
  If (Allocated(iCoSet)) Call mma_deallocate(iCoSet)
  If (Allocated(nStab)) Call mma_deallocate(nStab)
  If (Allocated(AtomLbl)) Call mma_deallocate(AtomLbl)
  If (Allocated(Smmtrc)) Call mma_deallocate(Smmtrc)

  If (Allocated(Q_nuclear)) Call mma_deallocate(Q_nuclear)
  If (Allocated(dMass)) Call mma_deallocate(dMass)
  If (Allocated(Coor)) Call mma_deallocate(Coor)
  If (Allocated(Grd)) Call mma_deallocate(Grd)
  If (Allocated(ANr)) Call mma_deallocate(ANr)
  If (Allocated(Weights)) Call mma_deallocate(Weights)
  If (Allocated(Shift)) Call mma_deallocate(Shift)
  If (Allocated(BMx)) Call mma_deallocate(BMx)

  If (Allocated(BM)) Call mma_deallocate(BM)
  If (Allocated(dBM)) Call mma_deallocate(dBM)
  If (Allocated(iBM)) Call mma_deallocate(iBM)
  If (Allocated(idBM)) Call mma_deallocate(idBM)
  If (Allocated(nqBM)) Call mma_deallocate(nqBM)

  If (Allocated(RefGeo)) Call mma_deallocate(RefGeo)
  If (Allocated(R12)) Call mma_deallocate(R12)
  If (Allocated(R12)) Call mma_deallocate(GradRef)

  If (Allocated(NAC)) Call mma_deallocate(NAC)
  If (Allocated(qInt)) Call mma_deallocate(qInt)
  If (Allocated(dqInt)) Call mma_deallocate(dqInt)
  End Subroutine Free_Slapaf



  Subroutine Get_Slapaf(iter,MaxItr,mTROld,lOld_Implicit, nsAtom_In,mLambda)
  Integer iter, MaxItr, mTROld, nsAtom_In, mLambda
  Logical lOld_Implicit
#include "real.fh"
#include "stdalloc.fh"
  Logical Exist
  Integer itmp, iOff, Lngth
  Integer, Allocatable:: Information(:)
  Real*8, Allocatable:: Relax(:)
  Character(LEN=100) SuperName
  Character(LEN=100), External:: Get_SuperName

  Initiated=.True.

  nsAtom=nsAtom_In

  Call mma_allocate(Information,7,Label='Information')

  Call qpg_iArray('Slapaf Info 1',Exist,itmp)
  If (Exist) Call Get_iArray('Slapaf Info 1',Information,7)

  If (.Not.Exist.or.(Exist.and.Information(1).eq.-99)) Then
!    Write (6,*) 'Reinitiate Slapaf fields on runfile'
     Information(:)=0
     Information(3)=-99
     Call Put_iArray('Slapaf Info 1',Information,7)
  End If

  iter  =Information(2)+1
  If (iter.ge.MaxItr+1) Then
     Write (6,*) 'Increase MaxItr in info_slapaf.fh'
     Call WarningMessage(2,'iter.ge.MaxItr+1')
     Call Abend()
  End If
  mTROld=Information(3)
  lOld_Implicit= Information(4).eq.1

  Call mma_deallocate(Information)

  If (.NOT.Allocated(Energy)) Then

  Call mma_allocate(Energy,          MaxItr+1,Label='Energy')
  Energy(:) = Zero
  Call mma_allocate(Energy0,         MaxItr+1,Label='Energy0')
  Energy0(:) = Zero
  Call mma_allocate(DipM,   3,       MaxItr+1,Label='DipM')
  DipM(:,:) = Zero
  Call mma_allocate(GNrm,            MaxItr+1,Label='GNrm')
  GNrm(:) = Zero
  Call mma_allocate(Cx,     3,nsAtom,MaxItr+1,Label='Cx')
  Cx(:,:,:) = Zero
  Call mma_allocate(Gx,     3,nsAtom,MaxItr+1,Label='Gx')
  Gx(:,:,:) = Zero
  Call mma_allocate(Gx0,    3,nsAtom,MaxItr+1,Label='Gx0')
  Gx0(:,:,:) = Zero
  Call mma_allocate(MF,     3,nsAtom,         Label='MF')
  MF(:,:) = Zero
  If (mLambda>0) Then
     Call mma_allocate(Lambda,mLambda,MaxItr+1,Label='Lambda')
     Lambda(:,:)=Zero
  End If

  End If

  If (iter==1) Return

  SuperName=Get_Supername()
  If (SuperName.ne.'numerical_gradient') Then

     Lngth=SIZE(Energy)+SIZE(Energy0)+SIZE(DipM)+SIZE(GNrm)+SIZE(Cx)+SIZE(Gx)  &
          +SIZE(Gx0)+SIZE(MF)+SIZE(Lambda)
     Call mma_allocate(Relax,Lngth,Label='Relax')
     Call Get_dArray('Slapaf Info 2',Relax,Lngth)

     iOff=1
     Call DCopy_(SIZE(Energy ),Relax(iOff),1,Energy ,1)
     iOff=iOff+SIZE(Energy)
     Call DCopy_(SIZE(Energy0),Relax(iOff),1,Energy0,1)
     iOff=iOff+SIZE(Energy0)
     Call DCopy_(SIZE(DipM   ),Relax(iOff),1,DipM   ,1)
     iOff=iOff+SIZE(DipM   )
     Call DCopy_(SIZE(GNrm   ),Relax(iOff),1,GNrm   ,1)
     iOff=iOff+SIZE(GNrm   )
     Call DCopy_(SIZE(Cx     ),Relax(iOff),1,Cx     ,1)
     iOff=iOff+SIZE(Cx     )
     Call DCopy_(SIZE(Gx     ),Relax(iOff),1,Gx     ,1)
     iOff=iOff+SIZE(Gx     )
     Call DCopy_(SIZE(Gx0    ),Relax(iOff),1,Gx0    ,1)
     iOff=iOff+SIZE(Gx0    )
     Call DCopy_(SIZE(MF     ),Relax(iOff),1,MF     ,1)
     iOff=iOff+SIZE(MF     )
     If (Allocated(Lambda)) Then
        Call DCopy_(SIZE(Lambda ),Relax(iOff),1,Lambda ,1)
        iOff=iOff+SIZE(Lambda )
     End If
     Call mma_deallocate(Relax)
   Else
     iter=1
   End If

  End Subroutine Get_Slapaf



  Subroutine Dmp_Slapaf(Stop,Just_Frequencies,Energy_In,Iter,MaxItr,mTROld, &
                        lOld_Implicit,nsAtom)
  Logical Stop, Just_Frequencies, lOld_Implicit
  Real*8 Energy_In
  Integer Iter, MaxItr, mTROld, nsAtom
#include "stdalloc.fh"
  Integer, Allocatable:: Information(:)
  Real*8, Allocatable:: Relax(:)
  Real*8, Allocatable:: GxFix(:,:)
  Integer iOff_Iter, nSlap, iOff, Lngth
  Logical Found
  Character(LEN=100) SuperName
  Character(LEN=100), External:: Get_SuperName

  If (.NOT.Initiated) Then
     Write (6,*) 'Dmp_Slapaf: Slapaf not initiated!'
     Call Abend()
  Else
     Initiated=.False.
  End If

!---  Write information of this iteration to the RLXITR file
  Call mma_allocate(Information,7,Label='Information')
  If (Stop) Then
     Information(1)=-99     ! Deactivate the record
     iOff_Iter=0
     Call Put_iScalar('iOff_Iter',iOff_Iter)
!
!    Restore the runfile data as if the computation was analytic
!    (note the gradient sign must be changed back)
!
     If (Just_Frequencies) Then
        Call Put_dScalar('Last Energy',Energy_In)
        Call mma_allocate(GxFix,3,nsAtom,Label='GxFix')
        call dcopy_(3*nsAtom,Gx,1,GxFix,1)
        GxFix(:,:) = - GxFix(:,:)
        Call Put_Grad(GxFix,3*nsAtom)
        Call mma_deallocate(GxFix)
        Call Put_dArray('Unique Coordinates',Cx,3*nsAtom)
        Call Put_Coord_New(Cx,nsAtom)
     End If
  Else
     Call qpg_iArray('Slapaf Info 1',Found,nSlap)
     If (Found) Then
        Call Get_iArray('Slapaf Info 1',Information,7)
        If (Information(1).ne.-99) Information(1)=MaxItr
     Else
        Information(1)=MaxItr
     End If
  End If

  SuperName=Get_Supername()
  If (SuperName.ne.'numerical_gradient') Then
     Information(2)=Iter
     Information(3)=mTROld ! # symm. transl /rot.
     If (lOld_Implicit) Then
        Information(4)=1
     Else
        Information(4)=0
     End If
     Information(5)=0
     Information(6)=SIZE(Energy)+SIZE(Energy0)+SIZE(DipM)+SIZE(GNrm)
     Information(7)=SIZE(Energy)+SIZE(Energy0)+SIZE(DipM)+SIZE(GNrm)+SIZE(Cx)
     Call Put_iArray('Slapaf Info 1',Information,7)

     Lngth=SIZE(Energy)+SIZE(Energy0)+SIZE(DipM)+SIZE(GNrm)+SIZE(Cx)+SIZE(Gx)  &
          +SIZE(Gx0)+SIZE(MF)+SIZE(Lambda)
     Call mma_allocate(Relax,Lngth,Label='Relax')
     iOff = 1
     Call DCopy_(SIZE(Energy ),Energy ,1,Relax(iOff),1)
     iOff = iOff + SIZE(Energy )
     Call DCopy_(SIZE(Energy0),Energy0,1,Relax(iOff),1)
     iOff = iOff + SIZE(Energy0)
     Call DCopy_(SIZE(DipM   ),DipM   ,1,Relax(iOff),1)
     iOff = iOff + SIZE(DipM   )
     Call DCopy_(SIZE(GNrm   ),GNrm   ,1,Relax(iOff),1)
     iOff = iOff + SIZE(GNrm   )
     Call DCopy_(SIZE(Cx     ),Cx     ,1,Relax(iOff),1)
     iOff = iOff + SIZE(Cx     )
     Call DCopy_(SIZE(Gx     ),Gx     ,1,Relax(iOff),1)
     iOff = iOff + SIZE(Gx     )
     Call DCopy_(SIZE(Gx0    ),Gx0    ,1,Relax(iOff),1)
     iOff = iOff + SIZE(Gx0    )
     Call DCopy_(SIZE(MF     ),MF     ,1,Relax(iOff),1)
     iOff = iOff + SIZE(MF     )
     If (Allocated(Lambda)) Then
        Call DCopy_(SIZE(Lambda ),Lambda ,1,Relax(iOff),1)
        iOff = iOff + SIZE(Lambda )
     End If
     Call Put_dArray('Slapaf Info 2',Relax,Lngth)
     Call mma_deallocate(Relax)
  End If
  Call mma_deallocate(Information)

  End Subroutine Dmp_Slapaf
End Module Slapaf_Info
