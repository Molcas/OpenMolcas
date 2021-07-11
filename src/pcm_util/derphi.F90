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

subroutine DerPhi(IOpt,IC,NESFJ,L1,L2,DP,DA,Vert,Centr,Sphere,IntSph,ISphe)

use PCM_Arrays, only: MxVert
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IOpt, IC, NESFJ, L1, L2, IntSph(MxVert), ISphe
real(kind=wp), intent(in) :: DP(MxVert,3), Vert(3,MxVert), Centr(3,MxVert), Sphere(4,*)
real(kind=wp), intent(out) :: DA
integer(kind=iwp) :: JJ, NS1, NS2
real(kind=wp) :: COSPHI, COSTH, D_COS, Delta, Dist2, DNORM1, DNORM2, DPHI, FACT, PROD, RC2, SENPHI, V1(3), V2(3), T12(3), VEC1(3), &
                 VEC2(3), VEC3(3), VEC4(3)
real(kind=wp), parameter :: Small1 = 1.0e-6_wp, Small = 1.0e-12_wp

! Find the derivative of: Phi(L1) Cos[Theta(L1)]
!
! IOpt = 0 : refers to the L1 side of tessera ITS, knowing the derivatives
!            of the position of the vertices by which is delimited
!
! IOpt = 1 : refers to the L1 side of tessera ITS wrt the radius of the
!            NSJ sphere knowing the derivatives of the positions of the
!            vertices by which is delimited
!
! NS1 is the sphere to which the tessera belongs, NS2 is the sphere that
! creates the side L1 by intersecting NS1

NS1 = ISPHE
NS2 = INTSPH(L1)

! Finds the coordinates of the vertices wrt the center of the circle on
! which is the arc L1 and the radius of the circle

V1(:) = VERT(:,L1)-CENTR(:,L1)
V2(:) = VERT(:,L2)-CENTR(:,L1)
RC2 = V1(1)**2+V1(2)**2+V1(3)**2
PROD = V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
COSPHI = PROD/RC2

! If the tessera is intersected very near its vertex the numerical
! approximation could lead to PROD>RC2 and COSPHI>1. Avoid this by:

Delta = 1.0e-12_wp
if (abs(CosPhi) > One) CosPhi = sign(One-Delta,CosPhi)
SENPHI = sqrt(One-COSPHI*COSPHI)

! Compute the derivative of Phi(L1)

VEC1(:) = V1(:)-COSPHI*V2(:)
VEC2(:) = DP(L2,:)
VEC3(:) = V2(:)-COSPHI*V1(:)
VEC4(:) = DP(L1,:)

! If the side is just that created by the moving sphere some components
! are corrected
! IOpt = 0: due to the derivative of the distance between sphere centers
! IOpt = 1: due to the derivative of the radius of sphere NESFJ

if (NS2 == NESFJ) then
  T12(:) = Sphere(1:3,NESFJ)-Sphere(1:3,NS1)
  DIST2 = T12(1)*T12(1)+T12(2)*T12(2)+T12(3)*T12(3)
  if (IOpt == 0) then
    FACT = (Sphere(4,NS1)**2-Sphere(4,NESFJ)**2+DIST2)/(Two*DIST2)
    VEC2(IC) = VEC2(IC)-FACT
    VEC4(IC) = VEC4(IC)-FACT
  else if (IOpt == 1) then
    VEC2(:) = VEC2(:)+Sphere(4,NESFJ)*T12(:)/DIST2
    VEC4(:) = VEC4(:)+Sphere(4,NESFJ)*T12(:)/DIST2
  else
    write(u6,'(a)') 'Illegal IOpt in DerPhi.'
    call Abend()
  end if
end if

DPHI = Zero
do JJ=1,3
  DPHI = DPHI-(VEC1(JJ)*VEC2(JJ)+VEC3(JJ)*VEC4(JJ))
end do
if (abs(SenPhi) < Small) then
  if (abs(DPhi) > Small1) then
    write(u6,'(a)') 'SenPhi small but not DPhi in DerPhi.'
    call Abend()
  end if
  DPhi = Zero
else
  DPHI = DPHI/(RC2*SENPHI)
end if

! Compute the cosine of the polar angle

DNORM1 = Zero
DNORM2 = Zero
V1(:) = VERT(:,L1)-Sphere(1:3,NS1)
V2(:) = Sphere(1:3,NS2)-Sphere(1:3,NS1)
do JJ=1,3
  DNORM1 = DNORM1+V1(JJ)*V1(JJ)
  DNORM2 = DNORM2+V2(JJ)*V2(JJ)
end do
DNORM1 = sqrt(DNORM1)
DNORM2 = sqrt(DNORM2)
COSTH = (V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3))/(DNORM1*DNORM2)

! If the side is not that formed by the moving sphere the derivative of
! the polar angle vanishes

D_COS = Zero

! Otherwise:

if (NS2 == NESFJ) then
  do JJ=1,3
    D_COS = D_COS+V2(JJ)*DP(L1,JJ)
  end do
  if (IOpt == 0) D_COS = D_COS+V1(IC)-Sphere(4,NS1)*COSTH*V2(IC)/DNORM2
  D_COS = D_COS/(Sphere(4,NS1)*DNORM2)
end if
DA = COSTH*DPHI+acos(COSPHI)*D_COS
DA = Sphere(4,NS1)*Sphere(4,NS1)*DA

return

end subroutine DerPhi
