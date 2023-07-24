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

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine AF3X3LEV(RDIST,DELTAE,C3val,C6val,C8val,De,ULR)
!=======================================================================
!*** Simplified version of AF3x3potRet which does not return derivatives

use Constants, only: Zero, One, Two, Three, Six, Eight, Twelve
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: RDIST, DELTAE, C3val, C6val, C8val, De
real(kind=wp), intent(out) :: ULR
integer(kind=iwp) :: I, J, L
real(kind=wp) :: DDe(3,3), DM1(3,3), DM3(3,3), DM5(3,3), DR(3,3), EIGVEC(3,1), H(3,3), M1, M3, M5, Q(3,3), RDIST2, RDIST3, RET, &
                 RETPi, RETSig, W(3)

M1 = C3val
M3 = C6val
M5 = C8val
! what is this number?
RET = 9.36423830e-4_wp*RDIST
RETSig = cos(RET)+(RET)*sin(RET)
RETPi = RETSig-RET**2*cos(RET)
RDIST2 = RDIST**2
RDIST3 = RDIST*RDIST2
!write(25,*) 'Variables = "r", "U(r)","U(r)-U(r)^2/(4De)"'
!write(25,*) 'zone T = "U(r)"'
! Initialize interaction matrix to 0.0
do I=1,3
  H(I,I) = Zero
end do
! Prepare interation matrix  H
H(1,1) = -(M1*RETSig+M3/(RDIST3)+M5/(RDIST3*RDIST2))/(Three*RDIST3)
H(1,2) = -sqrt(Two)*H(1,1)
H(2,1) = H(1,2)
H(1,3) = M1*RETPi/(sqrt(Six)*RDIST3)
H(3,1) = H(1,3)
H(2,2) = 2*H(1,1)+DELTAE
H(2,3) = H(1,3)/sqrt(Two)
H(3,2) = H(2,3)
H(3,3) = DELTAE
! Prepare radial derivative of interaction matrix (? is it needed ?)
DR(1,1) = (THree*M1*RETSig+Six*M3/RDIST3+Eight*M5/(RDIST3*RDIST2))/(Three*RDIST3*RDIST)
DR(1,2) = -sqrt(Two)*DR(1,1)
DR(2,1) = DR(1,2)
DR(2,2) = Two*DR(1,1)
DR(1,3) = -Three*H(1,3)/RDIST
DR(3,1) = DR(1,3)
DR(2,3) = -Three*H(2,3)/RDIST
DR(3,2) = DR(2,3)
DR(3,3) = Zero
! Partial derivative of interaction matric  H  w.r.t.  C3
DM1(1,1) = -(RETSig+M1/(Two*De*RDIST3))/(Three*RDIST3)
DM1(1,2) = -sqrt(Two)*DM1(1,1)
DM1(2,1) = DM1(1,2)
DM1(2,2) = Two*DM1(1,1)
DM1(1,3) = RETPi/(sqrt(Six)*RDIST3)
DM1(3,1) = DM1(1,3)
DM1(2,3) = DM1(1,3)/sqrt(Two)
DM1(3,2) = DM1(2,3)
DM1(3,3) = Zero
! Partial derivative of interaction matric  H  w.r.t.  C6
DM3(1,1) = -One/(Three*RDIST3**2)
DM3(1,2) = -sqrt(Two)*DM3(1,1)
DM3(1,3) = Zero
DM3(2,1) = DM3(1,2)
DM3(2,2) = Two*DM3(1,1)
DM3(2,3) = Zero
DM3(3,1) = DM3(1,3)
DM3(3,2) = DM3(2,3)
DM3(3,3) = Zero
! Partial derivative of interaction matric  H  w.r.t.  C8
DM5(1,1) = DM3(1,1)/(RDIST2)
DM5(1,2) = DM3(1,2)/(RDIST2)
DM5(1,3) = Zero
DM5(2,1) = DM3(1,2)
DM5(2,2) = DM3(2,2)/(RDIST2)
DM5(2,3) = Zero
DM5(3,1) = DM5(1,3)
DM5(3,2) = DM5(2,3)
DM5(3,3) = Zero
! Partial derivative of interaction matric  H  w.r.t.  De
DDe(1,1) = M1**2/(Twelve*(RDIST3*De)**2)
DDe(1,2) = -sqrt(Two)*DDe(1,1)
DDe(1,3) = Zero
DDe(2,1) = DDe(1,2)
DDe(2,2) = Two*DDe(1,1)
DDe(2,3) = Zero
DDe(3,1) = DDe(1,3)
DDe(3,2) = DDe(2,3)
DDe(3,3) = Zero
! Call subroutine to prepare and invert interaction matrix  H
call ZHEEVJ3(H,Q,W)
L = 1
! Nor - identify the lowest eigenvalue of  H  and label it  L
do J=2,3
  if (W(J) < W(L)) L = J
end do
ULR = -W(L)
EIGVEC(:,1) = Q(:,L)
write(u6,*) EIGVEC
!write(25,600) RDIST,ULR
!600 format(2D16.7)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Modulus = SQRABS(Z)
!Modulus = real(Z,kind=wp)**2

return

end subroutine AF3X3LEV
