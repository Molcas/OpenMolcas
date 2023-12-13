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
! Copyright (C) 1994-1997, Niclas Forsberg                             *
!               1997,2000, Per Ake Malmqvist                           *
!***********************************************************************

!module VibMod

!  Contains:
!    VibFreq        (AtCoord,InterVec,Mass,Hess,harmfreq,eigenVec,qMat,PED,D3,D4,x_anharm,anharmfreq,max_term)
!    CalcS          (AtCoord,InterVec,S)
!    BondStr        (R,i1,i2,j,S)
!    AngBend        (R1,R2,i1,i2,i3,j,S)
!    LinBend        (R1,R2,i1,i2,i3,j,S)
!    Torsion        (R1,R2,R3,i1,i2,i3,i4,j,S)
!    OutOfPl        (R1,R2,R3,i1,i2,i3,i4,j,S)
!    CalcG          (G,Mass,S)
!    NaNChk         (X,k,i,j)
!    Freq_mula      (Hess,G,V,Lambda,B,qMat)
!    CalcGprime     (Gprime,Mass,xvec,InterVec,AtCoord,NumOfAt)
!    CalcGdbleprime (Gdbleprime,Mass,xvec,InterVec,AtCoord,NumOfAt)
!    Anharm         (eigenVec,harmfreq,D3,D4,Gprime,Gdbleprime,x)
!    TransEnergy    (x_anharm,harmfreq,level1,level2)  Result(energy)
!    AnharmonicFreq (x_anharm,harmfreq,anharmfreq)
!    Int_to_Cart1   (InterVec,xvec,AtCoord)
!    Cart_To_Int0   (InterVec,AtCoord,xvec)
!    Cart_To_Int1   (InterVec,AtCoord,xvec,BMatrix)
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

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

subroutine BondStr(R,i1,i2,j,S,NumOfAt,NumInt)
!  Purpose:
!    Calculate contribution to the S-matrix due to bond stretching.
!
!  Input:
!    R        : Array of Real -  contains the cartesian coordinates of the bond.
!    i1,i2    : Integer - the number of the atom.
!    j        : Integer - the number of the internal coordinate.
!
!  Output:
!    S        : Real three dimensional array - the contributions to S for the parameters specified in the input.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: R(3)
integer(kind=iwp), intent(in) :: i1, i2, j, NumOfAt, NumInt
real(kind=wp), intent(inout) :: S(3,NumOfAt,NumInt)
integer(kind=iwp) :: k
real(kind=wp) :: SR

! Contributions to S.
SR = sqrt(R(1)**2+R(2)**2+R(3)**2)
do k=1,3
  S(k,i1,j) = -R(k)/SR
  S(k,i2,j) = -S(k,i1,j)
  call NaNChk(S(k,i1,j),k,i1,j)
  call NaNChk(S(k,i2,j),k,i1,j)
end do

end subroutine BondStr

subroutine NaNChk(X,k,i,j)

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: X
integer(kind=iwp), intent(in) :: k, i, j
character(len=16) :: str1, str2

write(str2,'(G16.8)') X
call Normalize(str2,str1)
if (str1(1:3) == 'NAN') then
  write(u6,*) ' CalcS subroutine produced Not-a-Number!'
  write(u6,*) ' Internal coordinate nr. j=',j
  write(u6,*) ' Atom nr.                i=',i
  write(u6,*) ' Component               k=',k
  write(u6,*) ' S(k,i,j)=',X
end if

return

end subroutine NaNChk

subroutine AngBend(R1,R2,i1,i2,i3,j,S,NumOfAt,NumInt)
!  Purpose:
!    Calculate contribution to the S-matrix due to valence angle bending.
!
!  Input:
!    R1,R2    : Array of Real - contains the cartesian coordinates of the bond.
!    i1,i2,i3 : Integer - the number of the atom.
!    j        : Integer - the number of the internal coordinate.
!
!  Output:
!    S        : Real three dimensional array - the contributions to S for the parameters specified in the input.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: R1(3), R2(3)
integer(kind=iwp), intent(in) :: i1, i2, i3, j, NumOfAt, NumInt
real(kind=wp), intent(inout) :: S(3,NumOfAt,NumInt)
integer(kind=iwp) :: k
real(kind=wp) :: CosTheta, F1, F2, F3, F4, F5, NR1(3), NR2(3), SinTheta, SR1, SR2, Theta
real(kind=wp), external :: dDot_

! Unit vectors.
SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
do k=1,3
  NR1(k) = R1(k)/SR1
  NR2(k) = R2(k)/SR2
end do

! Angle Theta.
Theta = acos(dDot_(3,NR1,1,NR2,1))
CosTheta = cos(Theta)
SinTheta = sin(Theta)

! Contributions to S.
F1 = One/(SR1*SinTheta)
S(1,i1,j) = (CosTheta*NR1(1)-NR2(1))*F1
S(2,i1,j) = (CosTheta*NR1(2)-NR2(2))*F1
S(3,i1,j) = (CosTheta*NR1(3)-NR2(3))*F1

F2 = SR1-SR2*CosTheta
F3 = SR2-SR1*CosTheta
F4 = F1/SR2
S(1,i2,j) = (F2*NR1(1)+F3*NR2(1))*F4
S(2,i2,j) = (F2*NR1(2)+F3*NR2(2))*F4
S(3,i2,j) = (F2*NR1(3)+F3*NR2(3))*F4

F5 = One/(SR2*SinTheta)
S(1,i3,j) = (CosTheta*NR2(1)-NR1(1))*F5
S(2,i3,j) = (CosTheta*NR2(2)-NR1(2))*F5
S(3,i3,j) = (CosTheta*NR2(3)-NR1(3))*F5

call NaNChk(S(1,i1,j),1,i1,j)
call NaNChk(S(2,i1,j),2,i1,j)
call NaNChk(S(3,i1,j),3,i1,j)
call NaNChk(S(1,i2,j),1,i2,j)
call NaNChk(S(2,i2,j),2,i2,j)
call NaNChk(S(3,i2,j),3,i2,j)
call NaNChk(S(1,i3,j),1,i3,j)
call NaNChk(S(2,i3,j),2,i3,j)
call NaNChk(S(3,i3,j),3,i3,j)

end subroutine AngBend

subroutine LinBend(R1,R2,i1,i2,i3,j,S,NumOfAt,NumInt)
!  Purpose:
!    Calculate contribution to the S-matrix due to valence angle bending for a linear molecule.
!
!  Input:
!    R1,R2    : Array of Real - contains the cartesian coordinates of the bond.
!    i1,i2,i3 : Integer - the number of the atom.
!    j        : Integer - the number of the internal coordinate.
!
!  Output:
!    S        : Real three dimensional array - the contributions to S for the parameters specified in the input.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: R1(3), R2(3)
integer(kind=iwp), intent(in) :: i1, i2, i3, j, NumOfAt, NumInt
real(kind=wp), intent(inout) :: S(3,NumOfAt,NumInt)
real(kind=wp) :: F1, F2, SR1, SR2

! Length of vectors.
SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)

! Contributions to S.
F1 = One/SR1
F2 = One/SR2

S(1,i1,j) = F1
S(2,i1,j) = Zero
S(3,i1,j) = Zero

S(1,i2,j) = -F1-F2
S(2,i2,j) = Zero
S(3,i2,j) = Zero

S(1,i3,j) = F2
S(2,i3,j) = Zero
S(3,i3,j) = Zero

S(1,i1,j+1) = Zero
S(2,i1,j+1) = F1
S(3,i1,j+1) = Zero

S(1,i2,j+1) = Zero
S(2,i2,j+1) = -F1-F2
S(3,i2,j+1) = Zero

S(1,i3,j+1) = Zero
S(2,i3,j+1) = F2
S(3,i3,j+1) = Zero

call NaNChk(S(1,i1,j),1,i1,j)
call NaNChk(S(2,i1,j),2,i1,j)
call NaNChk(S(3,i1,j),3,i1,j)
call NaNChk(S(1,i2,j),1,i2,j)
call NaNChk(S(2,i2,j),2,i2,j)
call NaNChk(S(3,i2,j),3,i2,j)
call NaNChk(S(1,i3,j),1,i3,j)
call NaNChk(S(2,i3,j),2,i3,j)
call NaNChk(S(3,i3,j),3,i3,j)

end subroutine LinBend

subroutine Torsion(R1,R2,R3,i1,i2,i3,i4,j,S,NumOfAt,NumInt)
!  Purpose:
!    Calculate contribution to the S-matrix due to torsion.
!
!  Input:
!    R1,R2,R3 : Array of Real - contains the cartesian coordinates of the bond.
!    i1,i2,
!    i3,i4    : Integer - the number of the atom.
!    j        : Integer - the number of the internal coordinate.
!
!  Output:
!    S        : Real three dimensional array - the contributions to S for the parameters specified in the input.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: R1(3), R2(3), R3(3)
integer(kind=iwp), intent(in) :: i1, i2, i3, i4, j, NumOfAt, NumInt
real(kind=wp), intent(inout) :: S(3,NumOfAt,NumInt)
integer(kind=iwp) :: k
real(kind=wp) :: CosTheta1, CosTheta2, F1, F10, F11, F12, F2, F3, F4, F5, F6, F7, F8, F9, NR1(3), NR2(3), NR3(3), SinTheta1, &
                 SinTheta2, SR1, SR2, SR3, Theta1, Theta2
real(kind=wp), external :: dDot_

! Unit vectors.
SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
do k=1,3
  NR1(k) = R1(k)/SR1
  NR2(k) = R2(k)/SR2
  NR3(k) = R3(k)/SR3
end do

! Angles Theta1 and Theta2.
Theta1 = acos(-dDot_(3,NR1,1,NR2,1))
Theta2 = acos(-dDot_(3,NR2,1,NR3,1))
CosTheta1 = cos(Theta1)
SinTheta1 = sin(Theta1)
CosTheta2 = cos(Theta2)
SinTheta2 = sin(Theta2)

! Contributions to S.
F1 = One/(SR1*(SinTheta1)**2)
F2 = NR1(2)*NR2(3)-NR1(3)*NR2(2)
F3 = NR1(3)*NR2(1)-NR1(1)*NR2(3)
F4 = NR1(1)*NR2(2)-NR1(2)*NR2(1)
S(1,i1,j) = -F2*F1
S(2,i1,j) = -F3*F1
S(3,i1,j) = -F4*F1

F5 = (SR2-SR1*CosTheta1)/(SR1*SR2*(SinTheta1)**2)
F6 = CosTheta2/(SR2*(SinTheta2)**2)
F7 = NR3(2)*NR2(3)-NR3(3)*NR2(2)
F8 = NR3(3)*NR2(1)-NR3(1)*NR2(3)
F9 = NR3(1)*NR2(2)-NR3(2)*NR2(1)
S(1,i2,j) = F5*F2+F6*F7
S(2,i2,j) = F5*F3+F6*F8
S(3,i2,j) = F5*F4+F6*F9

F10 = (SR2-SR3*CosTheta2)/(SR2*SR3*SinTheta2**2)
F11 = CosTheta1/(SR2*(SinTheta1)**2)
S(1,i3,j) = F10*F7+F11*F2
S(2,i3,j) = F10*F8+F11*F3
S(3,i3,j) = F10*F9+F11*F4

F12 = One/(SR3*(SinTheta2)**2)
S(1,i4,j) = -F7*F12
S(2,i4,j) = -F8*F12
S(3,i4,j) = -F9*F12

call NaNChk(S(1,i1,j),1,i1,j)
call NaNChk(S(2,i1,j),2,i1,j)
call NaNChk(S(3,i1,j),3,i1,j)
call NaNChk(S(1,i2,j),1,i2,j)
call NaNChk(S(2,i2,j),2,i2,j)
call NaNChk(S(3,i2,j),3,i2,j)
call NaNChk(S(1,i3,j),1,i3,j)
call NaNChk(S(2,i3,j),2,i3,j)
call NaNChk(S(3,i3,j),3,i3,j)
call NaNChk(S(1,i1,j),1,i4,j)
call NaNChk(S(2,i1,j),2,i4,j)
call NaNChk(S(3,i1,j),3,i4,j)

end subroutine Torsion

subroutine OutOfPl(R1,R2,R3,i1,i2,i3,i4,j,S,NumOfAt,NumInt)
!  Purpose:
!    Calculate contribution to the S-matrix due to the angle between a bond and a plane defined by two bonds.
!
!  Input:
!    R1,R2,R3 : Array of Real - contains the cartesian coordinates of the bond.
!    i1,i2,
!    i3,i4    : Integer - the number of the atom.
!    j        : Integer - the number of the internal coordinate.
!
!  Output:
!    S        : Real three dimensional array - the contributions to S for the parameters specified in the input.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: R1(3), R2(3), R3(3)
integer(kind=iwp), intent(in) :: i1, i2, i3, i4, j, NumOfAt, NumInt
real(kind=wp), intent(inout) :: S(3,NumOfAt,NumInt)
integer(kind=iwp) :: k
real(kind=wp) :: CosPhi, CosTheta, Dot, F0(3), F1, F10, F11, F12, F13, F14, F15, F2, F3, F4, F5, F6, F7, F8, F9, NR1(3), NR2(3), &
                 NR3(3), Phi, SinTheta, SR1, SR2, SR3, TanPhi, Theta
real(kind=wp), external :: Ddot_

! Unit vectors.
SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
do k=1,3
  NR1(k) = R1(k)/SR1
  NR2(k) = R2(k)/SR2
  NR3(k) = R3(k)/SR3
end do

! Calculation of the angle Theta.
Dot = Ddot_(3,NR2,1,NR3,1)
Theta = acos(Dot)
CosTheta = cos(Theta)
SinTheta = sin(Theta)

! Calculation of the out of plane angle Phi.
F1 = NR2(2)*NR3(3)-NR2(3)*NR3(2)
F2 = NR2(3)*NR3(1)-NR2(1)*NR3(3)
F3 = NR2(1)*NR3(2)-NR2(2)*NR3(1)
F0(1) = F1/SinTheta
F0(2) = F2/SinTheta
F0(3) = F3/SinTheta
Phi = asin(Ddot_(3,F0,1,NR1,1))
CosPhi = cos(Phi)
TanPhi = tan(Phi)

! Contributions to S.
F4 = One/(SR1*CosPhi*SinTheta)
F5 = TanPhi/SR1
S(1,i1,j) = F1*F4-NR1(1)*F5
S(2,i1,j) = F2*F4-NR1(2)*F5
S(3,i1,j) = F3*F4-NR1(3)*F5

F6 = NR3(2)*NR1(3)-NR3(3)*NR1(2)
F7 = NR3(3)*NR1(1)-NR3(1)*NR1(3)
F8 = NR3(1)*NR1(2)-NR3(2)*NR1(1)
F9 = One/(SR2*CosPhi*SinTheta)
F10 = TanPhi/(SR2*(SinTheta)**2)
S(1,i2,j) = F6*F9-NR2(1)*F10+CosTheta*F10*NR3(1)
S(2,i2,j) = F7*F9-NR2(2)*F10+CosTheta*F10*NR3(2)
S(3,i2,j) = F8*F9-NR2(3)*F10+CosTheta*F10*NR3(3)

F11 = NR1(2)*NR2(3)-NR1(3)*NR2(2)
F12 = NR1(3)*NR2(1)-NR1(1)*NR2(3)
F13 = NR1(1)*NR2(2)-NR1(2)*NR2(1)
F14 = One/(SR3*CosPhi*SinTheta)
F15 = TanPhi/(SR3*(SinTheta)**2)
S(1,i3,j) = F11*F14-NR3(1)*F15+CosTheta*F15*NR2(1)
S(2,i3,j) = F12*F14-NR3(2)*F15+CosTheta*F15*NR2(2)
S(3,i3,j) = F13*F14-NR3(3)*F15+CosTheta*F15*NR2(3)

S(1,i4,j) = -S(1,i1,j)-S(1,i2,j)-S(1,i3,j)
S(2,i4,j) = -S(2,i1,j)-S(2,i2,j)-S(2,i3,j)
S(3,i4,j) = -S(3,i1,j)-S(3,i2,j)-S(3,i3,j)

call NaNChk(S(1,i1,j),1,i1,j)
call NaNChk(S(2,i1,j),2,i1,j)
call NaNChk(S(3,i1,j),3,i1,j)
call NaNChk(S(1,i2,j),1,i2,j)
call NaNChk(S(2,i2,j),2,i2,j)
call NaNChk(S(3,i2,j),3,i2,j)
call NaNChk(S(1,i3,j),1,i3,j)
call NaNChk(S(2,i3,j),2,i3,j)
call NaNChk(S(3,i3,j),3,i3,j)
call NaNChk(S(1,i1,j),1,i4,j)
call NaNChk(S(2,i1,j),2,i4,j)
call NaNChk(S(3,i1,j),3,i4,j)

end subroutine OutOfPl

subroutine CalcS(AtCoord,InterVec,S,NumInt,NumOfAt)
!  Purpose:
!    Calculate the contribution to the S-matrix for each of the internal coordinates.
!
!  Input:
!    AtCoord  : Two dimensional Real array - contains the cartesian coordinates of the atoms.
!    InterVec : Integer array - contains the atoms that are used in the calculations of each internal coordinate.
!    S        : Real array filled with zeros.
!
!  Output:
!    S        : Real array - contains the contributions of each of the internal coordinates.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), NumInt, NumOfAt
real(kind=wp), intent(in) :: AtCoord(3,NumOfAt)
real(kind=wp), intent(inout) :: S(3,NumOfAt,NumInt)
integer(kind=iwp) :: i1, i2, i3, i4, IntType, k, n_Int
real(kind=wp) :: R(3), R1(3), R2(3), R3(3)

! Calculate the contributions to the S-matrix for each internal
! coordinate specified in the vector InterVec.
k = 1
IntType = InterVec(k)
n_Int = 1
do while (n_Int <= NumInt)
  if (IntType == 1) then
    ! Bond Stretching.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    R(:) = AtCoord(:,i2)-AtCoord(:,i1)
    call BondStr(R,i1,i2,n_Int,S,NumOfAt,NumInt)
    k = k+3
  end if
  if (IntType == 2) then
    ! Valence Angle Bending.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i3)-AtCoord(:,i2)
    call AngBend(R1,R2,i1,i2,i3,n_Int,S,NumOfAt,NumInt)
    k = k+4
  end if
  if (IntType == 3) then
    ! Linear Valence Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i3)-AtCoord(:,i2)
    call LinBend(R1,R2,i1,i2,i3,n_Int-1,S,NumOfAt,NumInt)
    k = k+4
  end if
  if (IntType == 4) then
    ! Torsion.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    R1(:) = AtCoord(:,i2)-AtCoord(:,i1)
    R2(:) = AtCoord(:,i3)-AtCoord(:,i2)
    R3(:) = AtCoord(:,i4)-AtCoord(:,i3)
    call Torsion(R1,R2,R3,i1,i2,i3,i4,n_Int,S,NumOfAt,NumInt)
    k = k+5
  end if
  if (IntType == 5) then
    ! Out of Plane Angle Bending.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i4)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i4)
    R3(:) = AtCoord(:,i3)-AtCoord(:,i4)
    call OutOfPl(R1,R2,R3,i1,i2,i3,i4,n_Int,S,NumOfAt,NumInt)
    k = k+5
  end if
  IntType = InterVec(k)
  n_Int = n_Int+1
end do

end subroutine CalcS

subroutine CalcG(G,Mass,S,NumInt,NumOfAt)
!  Purpose:
!    Calculate G-matrix.
!
!  Input:
!    Mass     : Real array - the mass of the atoms.
!    S        : Real three dimensional array - the contributions from all the internal coordinates.
!
!  Output:
!    G        : Real two dimensional array.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use Constants, only: Zero, One, uToau
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NumInt, NumOfAt
real(kind=wp), intent(out) :: G(NumInt,NumInt)
real(kind=wp), intent(in) :: Mass(NumOfAt), S(3,NumOfAt,NumInt)
integer(kind=iwp) :: i, j, k
real(kind=wp) :: GSum

do i=1,NumInt
  do j=1,NumInt
    GSum = Zero
    do k=1,NumOfAt
      GSum = GSum+(One/(uToau*Mass(k)))*(S(1,k,i)*S(1,k,j)+S(2,k,i)*S(2,k,j)+S(3,k,i)*S(3,k,j))
    end do
    G(i,j) = GSum
  end do
end do

end subroutine CalcG

subroutine Freq_mula(Hess,G,V,Lambda,B,qMat,NumInt,NumOfAt)
!  Purpose:
!    Find eigenvalues and eigenvectors of FG matrix.
!    The eigenvalues are stored in the array Lambda and the eigen-
!    vectors are stored as the columns of the matrix V.
!    The eigenvectors are then used to calculate the cartesian
!    displacements, which are stored in matrix X.
!
!  Input:
!    Hess     : Real two dimensional array -  contains the force constants expressed in internal
!    G        : Real two dimensional array.
!    B        : Real two dimensional array.
!
!  Output:
!    V        : Real two dimensional array - contains the eigenvectors of F*G as columns.
!    Lambda   : Real array - contains the eigenvalues of F*G.
!    qMat     : Real array - contains the cartesian displacements of the atoms.
!
!  Local:
!    U,Tmp    : Real two dimensional arrays.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NumInt, NumOfAt
real(kind=wp), intent(in) :: Hess(NumInt,NumInt), G(NumInt,NumInt), B(3*NumOfAt,NumInt)
real(kind=wp), intent(out) :: V(NumInt,NumInt), Lambda(NumInt), qMat(3*NumOfAt,NumInt)
real(kind=wp) :: Det
real(kind=wp), allocatable :: Temp(:,:), U(:,:)

! Solve secular equation.
call SolveSecEq(Hess,NumInt,V,G,Lambda)

! get memory space for U.
call mma_allocate(U,NumInt,NumInt,label='U')

! Copy matrix containing eigenvectors to U, because this matrix
! will be destroyed when subroutine Dool_MULA is called.
U(:,:) = V

! Calculate the cartesian diplacements, i.e. solve
!     qMat = B * ( B(T) * B )^(-1) * V.

! get memory for matrix Temp.
call mma_allocate(Temp,NumInt,NumInt,label='Temp')

call DGEMM_('T','N',NumInt,NumInt,3*NumOfAt,One,B,3*NumOfAt,B,3*NumOfAt,Zero,Temp,NumInt)
call Dool_MULA(Temp,NumInt,NumInt,U,NumInt,NumInt,Det)
call DGEMM_('N','N',3*NumOfAt,NumInt,NumInt,One,B,3*NumOfAt,U,NumInt,Zero,qMat,3*NumOfAt)

! free memory space of Temp and U.
call mma_deallocate(U)
call mma_deallocate(Temp)

end subroutine Freq_mula

subroutine CalcGprime(Gprime,Mass,xvec,InterVec,AtCoord,NumOfAt,h,NumInt)
!  Purpose:
!    Calculate first derivatives of G.
!
!  Input:
!    Mass     : Real array - the mass of the atoms.
!    xvec     : Real array - the geometry in internal coordinates.
!    InterVec : Integer array.
!    NumOfAt  : Integer - the number of atoms.
!
!  Output:
!    Gprime   : Real two dimensional array - first derivative of the inverse mass tensor G.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use mula_global, only: ngdim
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), NumOfAt, NumInt
real(kind=wp), intent(out) :: Gprime(ngdim,ngdim,ngdim)
real(kind=wp), intent(in) :: Mass(NumOfAt), xvec(NumInt), h
real(kind=wp), intent(inout) :: AtCoord(3,NumOfAt)
integer(kind=iwp) :: icoord, ih, iterm
real(kind=wp) :: xtmp(NumInt)
real(kind=wp), allocatable :: Gtemp(:,:,:), Stemp(:,:,:)

! Initialize.
call mma_allocate(Stemp,3,NumOfAt,NumInt,label='Stemp')
call mma_allocate(Gtemp,NumInt,NumInt,4,label='Gtemp')

do icoord=1,NumInt
  xtmp(:) = xvec
  Gtemp(:,:,:) = Zero
  iterm = 1
  do ih=-3,3,2
    xtmp(icoord) = xvec(icoord)+real(ih,kind=wp)*h
    !call Int_To_Cart(InterVec,xtmp,AtCoord,NumOfAt,NumInt,Mass)
    call Int_to_Cart1(InterVec,xtmp,AtCoord,NumOfAt,NumInt)
    Stemp(:,:,:) = Zero

    call CalcS(AtCoord,InterVec,Stemp,NumInt,NumOfAt)
    call CalcG(Gtemp(:,:,iterm),Mass,Stemp,NumInt,NumOfAt)

    iterm = iterm+1
  end do
  Gprime(1:NumInt,1:NumInt,icoord) = (Gtemp(:,:,1)-27.0_wp*Gtemp(:,:,2)+27.0_wp*Gtemp(:,:,3)-Gtemp(:,:,4))/(48.0_wp*h)
end do
!call Int_To_Cart(InterVec,xvec,AtCoord,NumOfAt,NumInt,Mass)
call Int_to_Cart1(InterVec,xtmp,AtCoord,NumOfAt,NumInt)

call mma_deallocate(Stemp)
call mma_deallocate(Gtemp)

end subroutine CalcGprime

subroutine CalcGdbleprime(Gdbleprime,Mass,xvec,InterVec,AtCoord,NumOfAt,h,NumInt)
!  Purpose:
!    Calculate second derivatives of G.
!
!  Input:
!    Mass       : Real array - the mass of the atoms.
!    xvec       : Real array - the geometry in internal coordinates.
!    InterVec   : Integer array.
!    NumOfAt    : Integer - the number of atoms.
!
!  Output:
!    Gdbleprime : Real two dimensional array - second derivative of the inverse mass tensor G.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use mula_global, only: ngdim
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Three
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), NumOfAt, NumInt
real(kind=wp), intent(out) :: Gdbleprime(ngdim,ngdim,ngdim,ngdim)
real(kind=wp), intent(in) :: Mass(NumOfAt), xvec(NumInt), h
real(kind=wp), intent(inout) :: AtCoord(3,NumOfAt)
integer(kind=iwp) :: icoord, jcoord
real(kind=wp), allocatable :: Gprime1(:,:,:), Gprime2(:,:,:), Gprime3(:,:,:), Gprime4(:,:,:), xtmp(:)

! Initialize.
call mma_allocate(xtmp,NumInt,label='xtmp')
call mma_allocate(Gprime1,NumInt,NumInt,NumInt,label='Gprime1')
call mma_allocate(Gprime2,NumInt,NumInt,NumInt,label='Gprime2')
call mma_allocate(Gprime3,NumInt,NumInt,NumInt,label='Gprime3')
call mma_allocate(Gprime4,NumInt,NumInt,NumInt,label='Gprime4')

do jcoord=1,NumInt
  xtmp(:) = xvec
  xtmp(jcoord) = xvec(jcoord)-Three*h
  call CalcGprime(Gprime1,Mass,xtmp,InterVec,AtCoord,NumOfAt,h,NumInt)
  xtmp(jcoord) = xvec(jcoord)-h
  call CalcGprime(Gprime2,Mass,xtmp,InterVec,AtCoord,NumOfAt,h,NumInt)
  xtmp(jcoord) = xvec(jcoord)+h
  call CalcGprime(Gprime3,Mass,xtmp,InterVec,AtCoord,NumOfAt,h,NumInt)
  xtmp(jcoord) = xvec(jcoord)+Three*h
  call CalcGprime(Gprime4,Mass,xtmp,InterVec,AtCoord,NumOfAt,h,NumInt)
  do icoord=1,NumInt
    Gdbleprime(1:NumInt,1:NumInt,icoord,jcoord) = &
      (Gprime1(:,:,icoord)-27.0_wp*Gprime2(:,:,icoord)+27.0_wp*Gprime3(:,:,icoord)-Gprime4(:,:,icoord))/(48.0_wp*h)
  end do
end do
!call Int_To_Cart(InterVec,xvec,AtCoord,NumOfAt,NumInt,Mass)
call Int_To_Cart1(InterVec,xvec,AtCoord,NumOfAt,NumInt)

call mma_deallocate(xtmp)
call mma_deallocate(Gprime1)
call mma_deallocate(Gprime2)
call mma_deallocate(Gprime3)
call mma_deallocate(Gprime4)

end subroutine CalcGdbleprime

subroutine Anharm(eigenVec,harmfreq,D3,D4,Gprime,Gdbleprime,x,nOsc)
!  Purpose:
!    Calculate the anharmonicity constants.
!    This routine assumes that the curvilinear coordinates are such
!    that they coincide with dimensionless normal coordinates for
!    small displacements.
!
!  Input:
!    eigenVec   : Real two dimensional array - eigenvectors.
!    harmfreq   : Real array - harmonic frequencies.
!    D3         : Real three dimensional array - cubic force constants.
!    D4         : Real four dimensional array - quartic force constants.
!    Gprime     : Real three dimensional array - first derivatives of the inverse mass tensor.
!    Gdbleprime : Real four dimensional array - second derivatives of the inverse mass tensor.
!
!  Output:
!    x          : Real two dimensional array - anharmonicity constants.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use mula_global, only: ngdim
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Four, Five, Six, Eight, Nine, Half, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nOsc
real(kind=wp), intent(in) :: eigenVec(nOsc,nOsc), harmfreq(nOsc), D3(ngdim,ngdim,ngdim), D4(ngdim,ngdim,ngdim,ngdim), &
                             Gprime(ngdim,ngdim,ngdim), Gdbleprime(ngdim,ngdim,ngdim,ngdim)
real(kind=wp), intent(out) :: x(nOsc,nOsc)
integer(kind=iwp) :: i, i1, j, j1, k, k1, l, l1, NumInt
real(kind=wp) :: coef, delta, deltaInv, rsum, sum1, sum2, tmp
real(kind=wp), allocatable :: C(:,:), T3(:,:,:), T4(:,:,:,:), Temp(:,:), V3(:,:,:), V4(:,:,:,:)

NumInt = nOsc

call mma_allocate(C,nOsc,nOsc,label='C')
call mma_allocate(Temp,nOsc,nOsc,label='Temp')
call mma_allocate(V3,nOsc,nOsc,nOsc,label='V3')
call mma_allocate(T3,nOsc,nOsc,nOsc,label='T3')
call mma_allocate(V4,nOsc,nOsc,nOsc,nOsc,label='V4')
call mma_allocate(T4,nOsc,nOsc,nOsc,nOsc,label='T4')

! Calculate the eigenvector matrix, C, in dimensionless normal coordinates.
Temp(:,:) = Zero
do i=1,NumInt
  Temp(i,i) = One/sqrt(harmfreq(i))
end do
call DGEMM_('n','n',NumInt,NumInt,NumInt,One,eigenVec,NumInt,Temp,NumInt,Zero,C,NumInt)

call mma_deallocate(Temp)

! Transform cubic force constants to dimensionless normal coordinates
do i=1,NumInt
  do j=1,NumInt
    do k=1,NumInt
      sum1 = Zero
      sum2 = Zero
      do i1=1,NumInt
        do j1=1,NumInt
          do k1=1,NumInt
            coef = C(i1,i)*C(j1,j)*C(k1,k)
            sum1 = sum1+coef*D3(i1,j1,k1)
            sum2 = sum2+coef*Gprime(i1,j1,k1)
          end do
        end do
      end do
      V3(i,j,k) = sum1
      T3(i,j,k) = sum2
    end do
  end do
end do

! Transform quartic force constants to dimensionless normal coordinates
do i=1,NumInt
  do j=1,NumInt
    do k=1,NumInt
      do l=1,NumInt
        sum1 = Zero
        sum2 = Zero
        do i1=1,NumInt
          do j1=1,NumInt
            do k1=1,NumInt
              do l1=1,NumInt
                coef = C(i1,i)*C(j1,j)*C(k1,k)*C(l1,l)
                sum1 = sum1+coef*D4(i1,j1,k1,l1)
                sum2 = sum2+coef*Gdbleprime(i1,j1,k1,l1)
              end do
            end do
          end do
        end do
        V4(i,j,k,l) = sum1
        T4(i,j,k,l) = sum2
      end do
    end do
  end do
end do

call mma_deallocate(C)

! Calculate diagonal anharmonicity constants.
do i=1,NumInt
  x(i,i) = V4(i,i,i,i)/16.0_wp
  tmp = One/(48.0_wp*harmfreq(i))
  x(i,i) = x(i,i)-tmp*(V3(i,i,i)*(Five*V3(i,i,i)+Six*T3(i,i,i))+Nine*T3(i,i,i)**2)
  do j=1,NumInt
    if (j /= i) then
      tmp = -One/(16.0_wp*harmfreq(j)*(Four*harmfreq(i)**2-harmfreq(j)**2))
      x(i,i) = x(i,i)+tmp*(Eight*harmfreq(i)**2-Three*harmfreq(j)**2)*(V3(i,i,j)**2+T3(i,i,j)**2)
      x(i,i) = x(i,i)+tmp*(Eight*harmfreq(i)**2-harmfreq(j)**2)*(Two*V3(i,i,j)*T3(i,i,j))
      x(i,i) = x(i,i)-tmp*Eight*harmfreq(i)*harmfreq(j)*T3(i,j,i)*(V3(i,i,j)-T3(i,i,j))
      x(i,i) = x(i,i)-tmp*Four*harmfreq(j)**2*T3(i,j,i)**2
    end if
  end do
end do

! Calculate off-diagonal anharmonicity constants.
do i=1,NumInt
  do j=1,NumInt
    if (i /= j) then
      x(i,j) = V4(i,i,j,j)*Quart
      x(i,j) = x(i,j)+T4(i,i,j,j)*Half
      x(i,j) = x(i,j)-(V3(i,i,i)+T3(i,i,i))*(V3(i,j,j)+T3(j,j,i))/(Four*harmfreq(i))
      x(i,j) = x(i,j)-(V3(i,i,j)+T3(i,i,j))*(V3(j,j,j)+T3(j,j,j))/(Four*harmfreq(j))
      do k=1,NumInt
        if ((k /= i) .and. (k /= j)) then
          x(i,j) = x(i,j)-(V3(i,i,k)+T3(i,i,k))*(V3(j,j,k)+T3(j,j,k))/(Two*harmfreq(k))
        end if
      end do
      tmp = Half*harmfreq(i)/(Four*harmfreq(i)**2-harmfreq(j)**2)
      x(i,j) = x(i,j)-tmp*((V3(i,i,j)-T3(i,i,j))**2+Two*(V3(j,j,i)-T3(j,j,i))*T3(i,j,j)+Four*T3(i,j,i)**2)
      tmp = Half*harmfreq(j)/(Four*harmfreq(j)**2-harmfreq(i)**2)
      x(i,j) = x(i,j)-tmp*((V3(j,j,i)-T3(j,j,i))**2+Two*(V3(i,i,j)-T3(i,i,j))*T3(i,j,i)+Four*T3(i,j,j)**2)
      do k=1,NumInt
        if ((k /= i) .and. (k /= j)) then
          delta = One
          rsum = harmfreq(i)+harmfreq(j)+harmfreq(k)
          delta = delta*rsum
          rsum = harmfreq(i)-harmfreq(j)-harmfreq(k)
          delta = delta*rsum
          rsum = -harmfreq(i)+harmfreq(j)-harmfreq(k)
          delta = delta*rsum
          rsum = -harmfreq(i)-harmfreq(j)+harmfreq(k)
          delta = delta*rsum
          deltaInv = One/delta
          x(i,j) = x(i,j)+deltaInv*Half*harmfreq(k)*(harmfreq(i)**2+harmfreq(j)**2-harmfreq(k)**2)* &
                   (V3(i,j,k)**2+T3(i,j,k)**2+T3(j,k,i)**2+T3(k,i,j)**2)
          x(i,j) = x(i,j)+deltaInv*Two*(harmfreq(i)*harmfreq(j)*harmfreq(k))*(V3(i,j,k)*T3(i,j,k)-T3(j,k,i)*T3(k,i,j))
          x(i,j) = x(i,j)-deltaInv*harmfreq(i)*(harmfreq(k)**2+harmfreq(j)**2-harmfreq(i)**2)* &
                   (V3(i,j,k)*T3(k,i,j)-T3(i,j,k)*T3(j,k,i))
          x(i,j) = x(i,j)-deltaInv*harmfreq(j)*(harmfreq(k)**2+harmfreq(i)**2-harmfreq(j)**2)* &
                   (V3(i,j,k)*T3(j,k,i)-T3(i,j,k)*T3(k,i,j))
        end if
      end do
    end if
  end do
end do

call mma_deallocate(V3)
call mma_deallocate(T3)
call mma_deallocate(V4)
call mma_deallocate(T4)

end subroutine Anharm

subroutine TransEnergy(G01,x_anharm1,harmfreq1,level1,G02,x_anharm2,harmfreq2,level2,energy,nDim)
!  Purpose:
!    Calculate transition energy between two states.
!
!  Input:
!    G01,G02    : Real - (G02-G01) = energy difference between the two states.
!    x_anharm1  : Real two dimensional array - anharmonicity constants for ground state.
!    x_anharm2  : Real two dimensional array - anharmonicity constants for excited state.
!    harmfreq1  : Real array - harmonic frequencies for ground state.
!    harmfreq2  : Real array - harmonic frequencies for excited state.
!    level1     : Integer array - quanta for ground state.
!    level2     : Integer array - quanta for excited state.
!
!  Output:
!    energy     : Real variable - energy for transition between level1 and level2.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, level1(nDim), level2(nDim)
real(kind=wp), intent(in) :: G01, x_anharm1(nDim,nDim), harmfreq1(nDim), G02, x_anharm2(nDim,nDim), harmfreq2(nDim)
real(kind=wp), intent(out) :: energy
integer(kind=iwp) :: i, j
real(kind=wp) :: G0, G1

! Calculate energy for level 1.
G0 = G01
do i=1,nDim
  G0 = G0+harmfreq1(i)*(level1(i)+Half)
  do j=i,nDim
    G0 = G0+x_anharm1(i,j)*(level1(i)+Half)*(level1(j)+Half)
  end do
end do

! Calculate energy for level 2.
G1 = G02
do i=1,nDim
  G1 = G1+harmfreq2(i)*(level2(i)+Half)
  do j=i,nDim
    G1 = G1+x_anharm2(i,j)*(level2(i)+Half)*(level2(j)+Half)
  end do
end do

! Calculate energy difference.
energy = G1-G0

end subroutine TransEnergy

subroutine AnharmonicFreq(x_anharm,harmfreq,anharmfreq,nOsc)
!  Purpose:
!    Calculate the anharmonic frequencies.
!
!  Input:
!    x_anharm   : Real two dimensional array - anharmonicity constants.
!    harmfreq   : Real array - harmonic frequencies.
!
!  Output:
!    anharmfreq : Real array - anharmonic frequencies.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nOsc
real(kind=wp), intent(in) :: x_anharm(nOsc,nOsc), harmfreq(nOsc)
real(kind=wp), intent(out) :: anharmfreq(nOsc)
integer(kind=iwp) :: istate
real(kind=wp) :: G1, G2
integer(kind=iwp), allocatable :: level1(:), level2(:)

call mma_allocate(level1,nOsc,label='level1')
call mma_allocate(level2,nOsc,label='level2')

G1 = Zero
G2 = G1
level1(:) = 0
do istate=1,nOsc
  level2(:) = 0
  level2(istate) = 1
  call TransEnergy(G1,x_anharm,harmfreq,level1,G2,x_anharm,harmfreq,level2,anharmfreq(istate),nOsc)
end do

call mma_deallocate(level1)
call mma_deallocate(level2)

end subroutine AnharmonicFreq

subroutine Cart_to_Int1(InterVec,AtCoord,xvec,BMatrix,NumOfAt,NumInt)
!  Purpose:
!    Calculate internal coordinates, and their gradients.
!
!  Input:
!    InterVec : Integer array - contains the atoms that are used in the calculations of each internal coordinate.
!    AtCoord  : Two dimensional Real array - contains the cartesian coordinates of the atoms.
!
!  Output:
!    xvec     : Real array - internal coordinates.
!    BMatrix  : Real array - Their gradients.
!
!  Written by:
!    Niclas Forsberg, Per Ake Malmqvist,
!    Dept. of Theoretical Chemistry, Lund University, 1997.

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), NumOfAt, NumInt
real(kind=wp), intent(in) :: AtCoord(3,NumOfAt)
real(kind=wp), intent(inout) :: xvec(NumInt)
real(kind=wp), intent(inout) :: BMatrix(3,NumOfAt,NumInt)
integer(kind=iwp) :: i1, i2, i3, i4, inttype, k, n_Int
real(kind=wp) :: Area, cos123, cos234, cosdv, cosv, cosv0, cot123, cot234, G1(3), G2(3), G3(3), Normal123(3), Normal234(3), &
                 NR1(3), NR2(3), NR3(3), R(3), R1(3), R2(3), R2R3, R3(3), sin123, sin234, sindv, sinv, sinv0, SR, SR1, SR2, SR3

k = 1
IntType = InterVec(k)
n_Int = 1
do while (n_Int <= NumInt)
  if (IntType == 1) then
    ! Bond.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    R(:) = AtCoord(:,i1)-AtCoord(:,i2)
    SR = sqrt(R(1)**2+R(2)**2+R(3)**2)
    xvec(n_Int) = SR
    BMatrix(:,i1,n_Int) = R/SR
    BMatrix(:,i2,n_Int) = -R/SR
    !write(u6,*) ' Bond. i1,i2=',i1,i2
    !write(u6,'(1x,a,3F16.8)') ' R=',R
    !write(u6,'(1x,a,3F16.8)') 'SR=',SR
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i1,n_Int)=',BMatrix(:,i1,n_Int)
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i2,n_Int)=',BMatrix(:,i2,n_Int)
    k = k+3
    n_Int = n_Int+1
  end if
  if (IntType == 2) then
    ! Valence Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    cosv = -(R1(1)*R2(1)+R1(2)*R2(2)+R1(3)*R2(3))/(SR1*SR2)
    if (cosv > One) cosv = One
    if (cosv < -One) cosv = -One
    sinv = sqrt(One-cosv**2)
    xvec(n_Int) = atan2(sinv,cosv)
    BMatrix(:,i1,n_Int) = (SR1*R2+SR2*cosv*R1)/(SR1**2*SR2*sinv)
    BMatrix(:,i3,n_Int) = -(SR2*R1+SR1*cosv*R2)/(SR2**2*SR1*sinv)
    BMatrix(:,i2,n_Int) = -(BMatrix(:,i1,n_Int)+BMatrix(:,i3,n_Int))
    !write(u6,*) ' Angle. i1,i2,i3=',i1,i2,i3
    !write(u6,'(1x,a,3F16.8)') ' R1=',R1
    !write(u6,'(1x,a,3F16.8)') ' R2=',R2
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i1,n_Int)=',BMatrix(:,i1,n_Int)
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i2,n_Int)=',BMatrix(:,i2,n_Int)
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i3,n_Int)=',BMatrix(:,i3,n_Int)
    k = k+4
    n_Int = n_Int+1
  end if
  if (IntType == 3) then
    ! Linear Valence Angle.
    !vv?? stop
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    cosv = -(R1(1)*R2(1)+R1(2)*R2(2)+R1(3)*R2(3))/(SR1*SR2)
    if (cosv > One) cosv = One
    if (cosv < -One) cosv = -One
    sinv = sqrt(One-cosv**2)
    xvec(n_Int) = atan2(sinv,cosv)
    BMatrix(1,i1,n_Int) = One/SR1
    BMatrix(1,i3,n_Int) = One/SR2
    BMatrix(1,i2,n_Int) = -(BMatrix(1,i1,n_Int)+BMatrix(1,i3,n_Int))
    BMatrix(2,i1,n_Int+1) = One/SR1
    BMatrix(2,i3,n_Int+1) = One/SR2
    BMatrix(2,i2,n_Int+1) = -(BMatrix(1,i1,n_Int)+BMatrix(1,i3,n_Int))
    n_Int = n_Int+2
    k = k+4
  end if
  if (IntType == 4) then
    ! Torsion.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    R3(:) = AtCoord(:,i3)-AtCoord(:,i4)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
    NR1(:) = R1/SR1
    NR2(:) = R2/SR2
    NR3(:) = R3/SR3
    cos123 = -(NR1(1)*NR2(1)+NR1(2)*NR2(2)+NR1(3)*NR2(3))
    if (cos123 > One) cos123 = One
    if (cos123 < -One) cos123 = -One
    sin123 = sqrt(One-cos123**2)
    cot123 = cos123/sin123
    Normal123(1) = (NR1(2)*NR2(3)-NR1(3)*NR2(2))/sin123
    Normal123(2) = (NR1(3)*NR2(1)-NR1(1)*NR2(3))/sin123
    Normal123(3) = (NR1(1)*NR2(2)-NR1(2)*NR2(1))/sin123
    cos234 = -(NR2(1)*NR3(1)+NR2(2)*NR3(2)+NR2(3)*NR3(3))
    if (cos234 > One) cos234 = One
    if (cos234 < -One) cos234 = -One
    sin234 = sqrt(One-cos234**2)
    cot234 = cos234/sin234
    Normal234(1) = (NR2(2)*NR3(3)-NR2(3)*NR3(2))/sin234
    Normal234(2) = (NR2(3)*NR3(1)-NR2(1)*NR3(3))/sin234
    Normal234(3) = (NR2(1)*NR3(2)-NR2(2)*NR3(1))/sin234
    sinv = NR2(1)*(Normal123(2)*Normal234(3)-Normal123(3)*Normal234(2))+ &
           NR2(2)*(Normal123(3)*Normal234(1)-Normal123(1)*Normal234(3))+ &
           NR2(3)*(Normal123(1)*Normal234(2)-Normal123(2)*Normal234(1))
    if (sinv > One) sinv = One
    if (sinv < -One) sinv = -One
    cosv = Normal123(1)*Normal234(1)+Normal123(2)*Normal234(2)+Normal123(3)*Normal234(3)
    if (cosv > One) cosv = One
    if (cosv < -One) cosv = -One
    cosv0 = cos(xvec(n_Int))
    sinv0 = sin(xvec(n_Int))
    sindv = sinv*cosv0-cosv*sinv0
    cosdv = cosv*cosv0+sinv*sinv0
    xvec(n_Int) = xvec(n_Int)+atan2(sindv,cosdv)
    G1(:) = Normal123/(SR1*sin123)
    G3(:) = Normal234/(SR3*sin234)
    G2(:) = (cot123*Normal123+cot234*Normal234)/SR2
    BMatrix(:,i1,n_Int) = G1
    BMatrix(:,i2,n_Int) = -G1+G2
    BMatrix(:,i3,n_Int) = -G2+G3
    BMatrix(:,i4,n_Int) = -G3
    !write(u6,*) ' Torsion. i1,i2,i3,i4=',i1,i2,i3,i4
    !write(u6,'(1x,a,3F16.8)') ' R1=',R1
    !write(u6,'(1x,a,3F16.8)') ' R2=',R2
    !write(u6,'(1x,a,3F16.8)') ' R3=',R3
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i1,n_Int)=',BMatrix(:,i1,n_Int)
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i2,n_Int)=',BMatrix(:,i2,n_Int)
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i3,n_Int)=',BMatrix(:,i3,n_Int)
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i4,n_Int)=',BMatrix(:,i4,n_Int)
    k = k+5
    n_Int = n_Int+1
  end if
  if (IntType == 5) then
    ! Out of Plane Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    R3(:) = AtCoord(:,i3)-AtCoord(:,i4)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
    R2R3 = (R2(1)*R3(1)+R2(2)*R3(2)+R2(3)*R3(3))
    cos234 = -R2R3/(SR2*SR3)
    if (cos234 > One) cos234 = One
    if (cos234 < -One) cos234 = -One
    sin234 = sqrt(One-cos234**2)
    Area = SR2*SR3*sin234
    Normal234(1) = (R2(2)*R3(3)-R2(3)*R3(2))/Area
    Normal234(2) = (R2(3)*R3(1)-R2(1)*R3(3))/Area
    Normal234(3) = (R2(1)*R3(2)-R2(2)*R3(1))/Area
    sinv = (R1(1)*Normal234(1)+R1(2)*Normal234(2)+R1(3)*Normal234(3))/SR1
    if (sinv > One) sinv = One
    if (sinv < -One) sinv = -One
    cosv = sqrt(One-sinv**2)
    xvec(n_Int) = atan2(sinv,cosv)
    G1(1) = (R2(2)*R3(3)-R3(2)*R2(3)-(Area*sinv/SR1)*R1(1))/(SR1*Area*cosv)
    G1(2) = (R2(3)*R3(1)-R3(3)*R2(1)-(Area*sinv/SR1)*R1(2))/(SR1*Area*cosv)
    G1(3) = (R2(1)*R3(2)-R3(1)*R2(2)-(Area*sinv/SR1)*R1(3))/(SR1*Area*cosv)
    G2(1) = (R3(2)*R1(3)-R1(2)*R3(3)-(SR1*sinv/Area)*(SR3*R2(1)-R2R3*R3(1)))/(SR1*Area*cosv)
    G2(2) = (R3(3)*R1(1)-R1(3)*R3(1)-(SR1*sinv/Area)*(SR3*R2(2)-R2R3*R3(2)))/(SR1*Area*cosv)
    G2(3) = (R3(1)*R1(2)-R1(1)*R3(2)-(SR1*sinv/Area)*(SR3*R2(3)-R2R3*R3(3)))/(SR1*Area*cosv)
    G3(1) = (R1(2)*R2(3)-R2(2)*R1(3)-(SR1*sinv/Area)*(SR2*R3(1)-R2R3*R2(1)))/(SR1*Area*cosv)
    G3(2) = (R1(3)*R2(1)-R2(3)*R1(1)-(SR1*sinv/Area)*(SR2*R3(2)-R2R3*R2(2)))/(SR1*Area*cosv)
    G3(3) = (R1(1)*R2(2)-R2(1)*R1(2)-(SR1*sinv/Area)*(SR2*R3(3)-R2R3*R2(3)))/(SR1*Area*cosv)
    BMatrix(:,i1,n_Int) = G1
    BMatrix(:,i2,n_Int) = -G1+G2
    BMatrix(:,i3,n_Int) = -G2+G3
    BMatrix(:,i4,n_Int) = -G3
    k = k+5
    n_Int = n_Int+1
  end if
  IntType = InterVec(k)
end do

end subroutine Cart_to_Int1

subroutine Cart_to_Int0(InterVec,AtCoord,xvec,NumOfAt,NumInt)
!
!  Purpose:
!    Calculate internal coordinates only (not their gradients).
!
!  Input:
!    InterVec : Integer array - contains the atoms that are used in the calculations of each internal coordinate.
!    AtCoord  : Two dimensional Real array - contains the cartesian coordinates of the atoms.
!
!  Output:
!    xvec     : Real array - internal coordinates.
!
!  Written by:
!    Niclas Forsberg, Per Ake Malmqvist,
!    Dept. of Theoretical Chemistry, Lund University, 1997.

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), NumOfAt, NumInt
real(kind=wp), intent(in) :: AtCoord(3,NumOfAt)
real(kind=wp), intent(inout) :: xvec(NumInt)
integer(kind=iwp) :: i1, i2, i3, i4, inttype, k, n_Int
real(kind=wp) :: cos123, cos234, cosdv, cosv, cosv0, Normal123(3), Normal234(3), NR1(3), NR2(3), NR3(3), R(3), R1(3), R2(3), R2R3, &
                 R3(3), sin123, sin234, sindv, sinv, sinv0, SR, SR1, SR2, SR3

k = 1
IntType = InterVec(k)
n_Int = 1
do while (n_Int <= NumInt)
  if (IntType == 1) then
    ! Bond.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    R(:) = AtCoord(:,i1)-AtCoord(:,i2)
    SR = sqrt(R(1)**2+R(2)**2+R(3)**2)
    xvec(n_Int) = SR
    k = k+3
    n_Int = n_Int+1
  end if
  if (IntType == 2) then
    ! Valence Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    cosv = -(R1(1)*R2(1)+R1(2)*R2(2)+R1(3)*R2(3))/(SR1*SR2)
    if (cosv > One) cosv = One
    if (cosv < -One) cosv = -One
    sinv = sqrt(One-cosv**2)
    xvec(n_Int) = atan2(sinv,cosv)
    k = k+4
    n_Int = n_Int+1
  end if
  if (IntType == 3) then
    ! Linear Valence Angle.
    !vv? stop
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i3)-AtCoord(:,i2)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    NR1(:) = R1/SR1
    NR2(:) = R2/SR2
    cosv = NR1(1)*NR2(1)+NR1(2)*NR2(2)+NR1(3)*NR2(3)
    if (cosv > One) cosv = One
    if (cosv < -One) cosv = -One
    xvec(n_Int) = acos(cosv)
    k = k+4
    n_Int = n_Int+2
  end if
  if (IntType == 4) then
    ! Torsion.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    R3(:) = AtCoord(:,i3)-AtCoord(:,i4)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
    NR1(:) = R1/SR1
    NR2(:) = R2/SR2
    NR3(:) = R3/SR3
    cos123 = -(NR1(1)*NR2(1)+NR1(2)*NR2(2)+NR1(3)*NR2(3))
    if (cos123 > One) cos123 = One
    if (cos123 < -One) cos123 = -One
    sin123 = sqrt(One-cos123**2)
    Normal123(1) = (NR1(2)*NR2(3)-NR1(3)*NR2(2))/sin123
    Normal123(2) = (NR1(3)*NR2(1)-NR1(1)*NR2(3))/sin123
    Normal123(3) = (NR1(1)*NR2(2)-NR1(2)*NR2(1))/sin123
    cos234 = -(NR2(1)*NR3(1)+NR2(2)*NR3(2)+NR2(3)*NR3(3))
    if (cos234 > One) cos234 = One
    if (cos234 < -One) cos234 = -One
    sin234 = sqrt(One-cos234**2)
    Normal234(1) = (NR2(2)*NR3(3)-NR2(3)*NR3(2))/sin234
    Normal234(2) = (NR2(3)*NR3(1)-NR2(1)*NR3(3))/sin234
    Normal234(3) = (NR2(1)*NR3(2)-NR2(2)*NR3(1))/sin234
    sinv = NR2(1)*(Normal123(2)*Normal234(3)-Normal123(3)*Normal234(2))+ &
           NR2(2)*(Normal123(3)*Normal234(1)-Normal123(1)*Normal234(3))+ &
           NR2(3)*(Normal123(1)*Normal234(2)-Normal123(2)*Normal234(1))
    if (sinv > One) sinv = One
    if (sinv < -One) sinv = -One
    cosv = Normal123(1)*Normal234(1)+Normal123(2)*Normal234(2)+Normal123(3)*Normal234(3)
    if (cosv > One) cosv = One
    if (cosv < -One) cosv = -One
    cosv0 = cos(xvec(n_Int))
    sinv0 = sin(xvec(n_Int))
    sindv = sinv*cosv0-cosv*sinv0
    cosdv = cosv*cosv0+sinv*sinv0
    xvec(n_Int) = xvec(n_Int)+atan2(sindv,cosdv)
    k = k+5
    n_Int = n_Int+1
  end if
  if (IntType == 5) then
    ! Out of Plane Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    R3(:) = AtCoord(:,i3)-AtCoord(:,i4)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
    R2R3 = (R2(1)*R3(1)+R2(2)*R3(2)+R2(3)*R3(3))
    cos234 = -R2R3/(SR2*SR3)
    if (cos234 > One) cos234 = One
    if (cos234 < -One) cos234 = -One
    sin234 = sqrt(One-cos234**2)
    Normal234(1) = (R2(2)*R3(3)-R2(3)*R3(2))/(SR2*SR3*sin234)
    Normal234(2) = (R2(3)*R3(1)-R2(1)*R3(3))/(SR2*SR3*sin234)
    Normal234(3) = (R2(1)*R3(2)-R2(2)*R3(1))/(SR2*SR3*sin234)
    sinv = (R1(1)*Normal234(1)+R1(2)*Normal234(2)+R1(3)*Normal234(3))/SR1
    if (sinv > One) sinv = One
    if (sinv < -One) sinv = -One
    xvec(n_Int) = asin(sinv)
    k = k+5
    n_Int = n_Int+1
  end if
  IntType = InterVec(k)
end do

end subroutine Cart_to_Int0

subroutine Int_to_Cart1(InterVec,xvec,AtCoord,NumOfAt,NumInt)
!  Purpose:
!    Transform a set of internal coordinates to cartesian coordinates.
!
!  Written by:
!    Per Ake Malmqvist
!    Dept. of Theoretical Chemistry, Lund University, 2000.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), NumOfAt, NumInt
real(kind=wp), intent(in) :: xvec(NumInt)
real(kind=wp), intent(inout) :: AtCoord(3,NumOfAt)
integer(kind=iwp) :: i1, i2, ii, iter, j, j1, j2, jj, k, ncart
real(kind=wp) :: det, RMSErr, rsum
real(kind=wp), allocatable :: BMatrix(:,:,:), EqMat(:,:), EqRHS(:), xvectmp(:)

! Initialize.

call mma_allocate(BMatrix,3,NumOfAt,NumInt,label='BMatrix')

! Purpose: For a given array InterVec, which describes the set of
! internal coordinates, and a given array xvec with values of these
! internal coordinates, find the structure which has the required
! internal coordinates and which is as close as possible to the
! input structure.
! Method: Newton-Raphson iteration, modified by Gauss' approximation
! in each iteration.
ncart = 3*NumOfAt

BMatrix(:,:,:) = Zero
call mma_allocate(EqMat,ncart,ncart,label='EqMat')
call mma_allocate(EqRHS,ncart,label='EqRHS')
call mma_allocate(xvectmp,NumInt,label='xvectmp')
!write(u6,*) ' Int_to_Cart1 has been called with target internal coordinates:'
!write(u6,'(1x,5F16.8)') xvec
iter = 0
RMSErr = huge(RMSErr)
do while (RMSErr > 1.0e-12_wp)
  iter = iter+1
  if (iter > 50) then
    write(u6,*) ' Int_to_Cart1 fails to converge.'
    call abend()
  end if
  ! Find the current internal coordinates and B matrix:
  !xvectmp = xvec
  xvectmp(:) = xvec
  call Cart_to_Int1(InterVec,AtCoord,xvectmp,BMatrix,NumOfAt,NumInt)

  ! Present error:
  rsum = Zero
  do k=1,NumInt
    rsum = rsum+(xvectmp(k)-xvec(k))**2
  end do
  RMSErr = sqrt(rsum)
  ! The B matrix is the matrix of partial derivatives of xvectmp with
  ! respect to AtCoord.
  ! Form the matrix B(transpose)*B -- This will be the matrix of the
  ! linearized equation system.
  jj = 0
  do j1=1,NumOfAt
    do j2=1,3
      jj = jj+1
      ii = 0
      do i1=1,NumOfAt
        do i2=1,3
          ii = ii+1
          rsum = Zero
          do k=1,NumInt
            rsum = rsum+BMatrix(i2,i1,k)*BMatrix(j2,j1,k)
          end do
          EqMat(ii,jj) = rsum
        end do
      end do
    end do
  end do
  ! The right-hand side of the equation system:
  ii = 0
  do i1=1,NumOfAt
    do i2=1,3
      ii = ii+1
      rsum = Zero
      do k=1,NumInt
        rsum = rsum+BMatrix(i2,i1,k)*(xvec(k)-xvectmp(k))
      end do
      EqRHS(ii) = rsum
    end do
  end do
  !D write(u6,*) ' EqRHS:'
  !D write(u6,'(1x,5f16.8)') EqRHS
  ! Automatic step limitation, and in fact ensuring that we go in the
  ! direction of minimization always, is achieved indirectly:
  do j=1,ncart
    EqMat(j,j) = EqMat(j,j)+1.0e-12_wp
  end do
  ! Solve the linear equation system:
  call Dool_MULA(EqMat,ncart,ncart,EqRHS,ncart,1,det)
  !D write(u6,*) ' Solution vector:'
  !D write(u6,'(1x,5f16.8)') EqRHS
  ! After return from dool, EqMat is destroyed and EqRHS is solution.
  ! Update the cartesian coordinates:
  ii = 0
  do i1=1,NumOfAt
    do i2=1,3
      ii = ii+1
      AtCoord(i2,i1) = AtCoord(i2,i1)+EqRHS(ii)
    end do
  end do
  !write(u6,*) ' Iteration iter=',iter
  !write(u6,*) ' AtCoord array:'
  !do i1=1,NumOfAt
  !  write(u6,'(1x,3f16.8)') (AtCoord(i2,i1),i2=1,3)
  !end do
  !write(u6,*) 'RMSErr:',RMSErr
end do
call mma_deallocate(EqMat)
call mma_deallocate(EqRHS)
call mma_deallocate(xvectmp)
call mma_deallocate(BMatrix)

return

end subroutine Int_to_Cart1

!end module VibMod
