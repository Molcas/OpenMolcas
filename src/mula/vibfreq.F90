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
! Copyright (C) 1995,1996, Niclas Forsberg                             *
!***********************************************************************

!module VibMod

!  Contains:
!    VibFreq        (AtCoord,InterVec,Mass,Hess,harmfreq,eigenVec,qMat,
!                    PED,D3,D4,x_anharm,anharmfreq,max_term)
!    CalcS          (AtCoord,InterVec,S)
!    BondStr        (R,i1,i2,j,S)
!    AngBend        (R1,R2,i1,i2,i3,j,S)
!    LinBend        (R1,R2,i1,i2,i3,j,S)
!    Torsion        (R1,R2,R3,i1,i2,i3,i4,j,S)
!    OutOfPl        (R1,R2,R3,i1,i2,i3,i4,j,S)
!    CalcG          (G,Mass,S)
!    Freq           (Hess,G,V,Lambda,B,qMat)
!    CalcGprime     (Gprime,Mass,xvec,InterVec,AtCoord,NumOfAt)
!    CalcGdbleprime (Gdbleprime,Mass,xvec,InterVec,AtCoord,NumOfAt)
!    Anharm         (eigenVec,harmfreq,D3,D4,Gprime,Gdbleprime,x)
!    TransEnergy    (x_anharm,harmfreq,level1,level2)  Result(energy)
!    AnharmonicFreq (x_anharm,harmfreq,anharmfreq)
!    Int_to_Cart1   (InterVec,xvec,AtCoord)
!    Cart_To_Int0   (InterVec,AtCoord,xvec)
!
!  Uses:
!    Constants
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

!use IOTools

!contains

subroutine VibFreq(AtCoord,xvec,InterVec,Mass,Hess,G,Gprime,Gdbleprime,harmfreq,eigenVec,qMat,PED,D3,D4,x_anharm,anharmfreq, &
                   max_term,nOsc,NumOfAt)
!  Purpose:
!    Calculates the vibrational frequencies of a molecule.
!
!  Input:
!    InterVec   : Integer array
!    Mass       : Real array - masses of the atoms.
!    xvec       : Real array - geometry of molecule in internal coordinates.
!    Hess       : Real two dimensional array - force constant matrix.
!    D3         : Real three dimensional array - third derivatives of potential surface.
!    D4         : Real four dimensional array - fourth derivatives of potential surface.
!    max_term   : Integer - highest power of term in polynomial fit.
!
!  Output:
!    AtCoord    : Real two dimensional array - cartesian coordinates of the atoms.
!    harmfreq   : Real array - contains harmonical frequencies.
!    eigenVec   : Real two dimensional array - contains eigenvectors.
!    qMat       : Real two dimensional array - cartesian displacement vectors.
!    PED        : Real three dimensional array - potential energy distribution.
!    x_anharm   : Real two dimensional array - anharmonicity constants.
!    anharmfreq : Real array - contains anharmonical frequencies.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use mula_global, only: ngdim
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), max_term, nOsc, NumOfAt
real(kind=wp), intent(inout) :: AtCoord(3,NumOfAt)
real(kind=wp), intent(out) :: xvec(nosc), G(nosc,nosc), Gprime(ngdim,ngdim,ngdim), Gdbleprime(ngdim,ngdim,ngdim,ngdim), &
                              harmfreq(nosc), eigenVec(nosc,nosc), qMat(3*NumOfAt,nOsc), PED(nosc,nosc,nosc), x_anharm(nosc,nOsc), &
                              anharmfreq(nosc)
real(kind=wp), intent(in) :: Mass(NumOfAt), Hess(nOsc,nOsc), D3(ngdim,ngdim,ngdim), D4(ngdim,ngdim,ngdim,ngdim)
integer(kind=iwp) :: NumInt
real(kind=wp) :: dh
real(kind=wp), allocatable :: B(:,:), Lambda(:), V(:,:)

! Initialize.
!D write(u6,*) ' Entered VIBFREQ.'
NumInt = nOsc
!D write(u6,*) ' NumInt:',NumInt
!D write(u6,*) ' NumOfAt:',NumOfAt
call mma_allocate(V,NumInt,NumInt,label='V')
call mma_allocate(B,3*NumOfAt,NumInt,label='B')

call mma_allocate(Lambda,NumInt,label='Lambda')

! Transform coordinates.
xvec(:) = Zero
!D write(u6,*) ' VIBFREQ, calling Cart_to_Int0.'
call Cart_To_Int0(InterVec,AtCoord,xvec,NumOfAt,NumInt)
!D write(u6,*) ' VIBFREQ, back from Cart_to_Int0.'
!D write(u6,*) ' xvec:'
!D write(u6,'(5f16.8)') xvec

! Calculate the contributions to the B matrix for each internal coordinate.
B(:,:) = Zero
!D write(u6,*) ' VIBFREQ, calling CalcS.'
call CalcS(AtCoord,InterVec,B,NumInt,NumOfAt)
!D write(u6,*) ' VIBFREQ, back from CalcS.'

! Calculate G matrix and first and second derivatives of the G matrix.
!D write(u6,*) ' VIBFREQ, calling CalcG.'
call CalcG(G,Mass,B,NumInt,NumOfAt)
!D Write(u6,*) ' VIBFREQ, back from CalcG.'
Gprime(:,:,:) = Zero
Gdbleprime(:,:,:,:) = Zero
if (max_term > 2) then
  dh = 1.0e-3_wp
  call CalcGprime(Gprime,Mass,xvec,InterVec,AtCoord,NumOfAt,dh,NumInt)
  dh = 1.0e-2_wp
  call CalcGdbleprime(Gdbleprime,Mass,xvec,InterVec,AtCoord,NumOfAt,dh,NumInt)
end if

! Given Hess and G, calculate the eigenvalues and eigenvectors of G*Hess.
!D write(u6,*) ' VIBFREQ, calling Freq.'
call Freq_mula(Hess,G,V,Lambda,B,qMat,nOsc,NumOfAt)
!D write(u6,*) ' VIBFREQ, back from Freq.'
!D write(u6,*) ' Lambda:'
!D write(u6,'(5f16.8)') Lambda

! Calculate harmonic frequencies.

harmfreq(:) = sqrt(abs(Lambda))
!D write(u6,*) ' harmfreq:'
!D write(u6,'(5f16.8)') harmfreq
eigenVec(:,:) = V

! Anharmonicity calculations (if we have third and possibly fourth
! derivatives). First calculation of the anharmonicity constants and
! then calculation of the fundamental frequencies.
x_anharm(:,:) = Zero
if (max_term > 2) then
  call Anharm(eigenVec,harmfreq,D3,D4,Gprime,Gdbleprime,x_anharm,nOsc)
  call AnharmonicFreq(x_anharm,harmfreq,anharmfreq,nOsc)
end if

! Calculate potential energy distribution.
call PotDist(Hess,V,Lambda,PED,NumInt,nOsc)

! Free memory space of B, G and V.
call mma_deallocate(B)
call mma_deallocate(V)
call mma_deallocate(Lambda)

end subroutine VibFreq
