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
! Copyright (C) 1999,2023, Roland Lindh                                *
!               2012, Andy May                                         *
!***********************************************************************
!***********************************************************************
!                                                                      *
!                                                                      *
!      These functions should eventually be used globally for any      *
!      seward utility, be it integrals, gradients, direct Fock,        *
!      direct integrals, and second order derivatives.                 *
!                                                                      *
!                          W A R N I N G !                             *
!                                                                      *
!                                                                      *
!      Observe that the index array (IndZ) should always be placed     *
!      after all arrays with real!                                     *
!                                                                      *
!      Real Arrays:                                                    *
!      alpha  : alpha Gaussian exponent                                *
!      beta   : beta  Gaussian exponent                                *
!      Z      : alpha+beta Gaussian exponent                           *
!      Kappa  : Common factor from Gaussian product theorem            *
!      Pcoor  : Coordinates of Gaussian product                        *
!      ab     : Max(|(ab|ab)|^{1/2})    (primitive basis)              *
!               over all angular components                            *
!      abCon  : Max(|(ab|ab)|^{1/2})* C ("contracted basis")           *
!               over all angular components                            *
!                                                                      *
!      Integer Arrays:                                                 *
!      IndZ   : Pair index (rectangular indexation), the last          *
!               entry in the array contains the effective number       *
!               of elements                                            *
!                                                                      *
!      Scalars:                                                        *
!      EstI   : Estimated largest contracted integral |(ab|ab)|^{1/2}  *
!      ZtMax  : Z of the largest abCon                                 *
!      abMax  : largest abCon                                          *
!      ZetaM  : largest Zeta value                                     *
!      ZtMaxD : Z of the largest ab * D                                *
!      abMaxD : largest ab * D                                         *
!                                                                      *
!      Auxiliary arrays:                                               *
!      HrrMtrx: matrices to use for transformation from intermediate   *
!               integrals to real spherical harmonics                  *
!                                                                      *
!      Author: R. Lindh                                                *
!              Dept. of Chemical Physics                               *
!              University of Lund, Sweden                              *
!              April 1999                                              *
!                                                                      *
!      Converted from statement functions by A. May June 2012          *
!      Converted to a user defined type bt R.L. October 2023           *
!***********************************************************************
!                                                                      *
Module k2_structure

use Constants, only: Zero

private

type k2_type
Integer :: nZeta=0
Integer :: ijCmp=0
Integer :: nHm=0
real*8,  Allocatable:: ZZZ(:)
real*8,  Pointer:: Zeta(:)
real*8,  Pointer:: Kappa(:)
real*8,  Pointer:: PCoor(:,:)
real*8,  Pointer:: ZInv(:)
real*8,  Pointer:: ab(:)
real*8,  Pointer:: abG(:,:)
real*8,  Pointer:: abCon(:)
real*8,  Pointer:: Alpha(:)
real*8,  Pointer:: Beta(:)
real*8,  Pointer:: HrrMtrx(:,:)

integer, Allocatable:: IndZ(:)
real*8              :: EstI=Zero
real*8              :: ZtMax=Zero
real*8              :: ZtMaxD=Zero
real*8              :: ZetaM=Zero
real*8              :: abMax=Zero
real*8              :: abMaxD=Zero
End type k2_type

type (k2_type), Allocatable, Target :: k2data(:,:)

Integer, Parameter:: nDArray=11,nDScalar=9

Logical :: k2_processed=.False.
Integer, Allocatable :: IndK2(:,:)
Integer nIndk2

public :: k2_type, k2data, Allocate_k2data, Free_k2data
public :: nDArray, nDScalar, k2_processed, IndK2, nIndk2

contains

Subroutine Allocate_k2data(k2data,nZeta,ijCmp,nHm)
use stdalloc, only: mma_allocate
 use Symmetry_Info, only: nIrrep
Implicit None
Integer nZeta, ijCmp,nHm
type (k2_type), target:: k2data

Integer nk2, iS, iE

k2Data%nZeta=nZeta
k2Data%nHm  =nHm
k2Data%ijCmp=ijCmp

nk2 = nZeta*10 + nHM*nIrrep + nZeta*ijCmp*2
Call mma_allocate(k2Data%ZZZ,nk2,Label='%ZZZ')
iS=1
iE=nZeta
k2Data%Zeta(1:nZeta) => k2Data%ZZZ(iS:iE)
iS=1+iE
iE=iE + nZeta
k2Data%Kappa(1:nZeta) => k2Data%ZZZ(iS:iE)
iS=1+iE
iE=iE + nZeta * 3
k2Data%PCoor(1:nZeta,1:3) => k2Data%ZZZ(iS:iE)
iS=1+iE
iE=iE + nZeta
k2Data%ZInv(1:nZeta) => k2Data%ZZZ(iS:iE)
iS=1+iE
iE=iE + nZeta
k2Data%ab(1:nZeta) => k2Data%ZZZ(iS:iE)
iS=1+iE
iE=iE + nZeta
k2Data%abCon(1:nZeta) => k2Data%ZZZ(iS:iE)
iS=1+iE
iE=iE + nZeta
k2Data%Alpha(1:nZeta) => k2Data%ZZZ(iS:iE)
iS=1+iE
iE=iE + nZeta
k2Data%Beta(1:nZeta) => k2Data%ZZZ(iS:iE)
iS=1+iE
iE=iE + nHm*nIrrep
k2Data%HrrMtrx(1:nHm,1:nIrrep) => k2Data%ZZZ(iS:iE)
If (ijCmp/=0) Then
   iS=1+iE
   iE=iE + nZeta*ijCmp*2
   k2Data%abG(1:nZeta*ijCmp,1:2) => k2Data%ZZZ(iS:iE)
End If

Call mma_allocate(k2Data%IndZ,     nZeta+1,Label='%IndZ')   ! yes +1!
End Subroutine Allocate_k2data

Subroutine Free_k2data()
Implicit None
Integer nIrrep, ik2, iIrrep, i

nIrrep=Size(k2data,1)
ik2   =Size(k2data,2)

Do i = 1, ik2
   Do iIrrep = 1, nIrrep
      Call Free_k2data_Internal(k2data(iIrrep,i))
   End Do
End Do

Deallocate(k2data)

End Subroutine Free_k2data

Subroutine Free_k2data_Internal(k2data_1D)
use stdalloc, only: mma_deallocate
Implicit None
type (k2_type):: k2data_1D

k2Data%nZeta=0
k2Data%nHm  =0
k2Data%ijCmp=0
If (allocated(k2Data_1D%ZZZ))  Call mma_deallocate(k2Data_1D%ZZZ)

Nullify(k2Data_1D%Kappa,k2Data_1D%Pcoor,k2Data_1D%ZInv,k2Data_1D%ab,k2Data_1D%abCon, &
        k2Data_1D%Alpha,k2Data_1D%Beta,k2Data_1D%HRRMtrx,k2Data_1D%abG)

If (allocated(k2Data_1D%IndZ)) Call mma_deallocate(k2Data_1D%IndZ)
End Subroutine Free_k2data_Internal


End Module k2_structure
