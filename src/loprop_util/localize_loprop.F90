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

subroutine Localize_LoProp(Ttot,Ttot_Inv,nBas,SMatrix,iCenter,iType)
!                                                                      *
!***********************************************************************
!                                                                      *
! H2 molecule with 2s basis functions on each H: A and B are the two
! centers
! Original S matrix  (nz=non zero; 10-8 should be zero but it is not for
! numerical imprecision
!       1     10-8  nz   nz
! S0=   10-8  1     nz   nz
!       nz    nz    1    10-8
!       nz    nz    10-8 1
!
! Step 1. LO S0 in order to generate S1 (orthogonalization of all the
! blocks for the same center
!       1     0     nz   nz
! S1=   0     1     nz   nz
!       nz    nz    1    0
!       nz    nz    0    1
!
! Step 2. LO S1 -> S2 (orthogonalization of the OAOB and VAVB blocks
!       1     0     0    nz
! S2=   0     1     nz   0
!       0     nz    1    0
!       nz    0     0    1
!
! Step 3. GS S2 -> S3 (orthogonalization of the VBOA and VAOB blocks
!       1     0     0    0
! S3=   0     1     0    0
!       0     0     1    0
!       0     0     0    1
!                                                                      *
!***********************************************************************
!                                                                      *

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBas, iCenter(nBas), iType(nBas)
real(kind=wp), intent(out) :: Ttot(nBas,nBas), Ttot_Inv(nBas,nBas)
real(kind=wp), intent(inout) :: SMatrix(nBas,nBas)
integer(kind=iwp) :: iBas, IndType(7), iUHF, LuOut, nOcc, nSym
character(len=128) :: OrbName
character(len=80) :: Note
character(len=6) :: Filename
real(kind=wp), allocatable :: E(:), S(:,:), SSave(:,:), T1(:,:), T2(:,:), T3(:,:), T4(:,:), tmp(:,:)
integer(kind=iwp), parameter :: Occ = 1, Vir = 0

!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate temporary memory

call mma_allocate(T1,nBas,nBas,label='T1')
call mma_allocate(T2,nBas,nBas,label='T2')
call mma_allocate(T3,nBas,nBas,label='T3')
call mma_allocate(T4,nBas,nBas,label='T4')
call mma_allocate(tmp,nBas,nBas,label='tmp')
call mma_allocate(S,nBas,nBas,label='S')
call mma_allocate(SSave,nBas,nBas,label='SSave')

! Save S because GS will destroy it!

S(:,:) = SMatrix(:,:)
!                                                                      *
!***********************************************************************
!                                                                      *
! Step 1. GS S0 ->S1

call unitmat(T1,nBas)
T2(:,:) = T1

call Step1(iCenter,S,nBas,T1,iType,SMatrix,tmp)
!                                                                      *
!***********************************************************************
!                                                                      *
! Step 2. LO S1 ->S2

SSave(:,:) = S(:,:)

call Step2(iCenter,S,nBas,T2,iType,SSave,tmp)
!                                                                      *
!***********************************************************************
!                                                                      *
! Now set up things for the GS orthogonalization
!
! Step 3. GS S2 ->S3

call Step3(iCenter,S,nBas,T3,iType)

!                                                                      *
!***********************************************************************
!                                                                      *
! Now do a final Lowdin to remove virtual-virtual
!
! Step 4. LO S3 ->S4

call Step4(S,nBas,T4,iType)

!                                                                      *
!***********************************************************************
!                                                                      *
! Now I have LO orthog twice and GS orthog
! The corresponding transformation matrices are T1 T2 and TGS
! Now I generate the total transformation matrix TTOT=T1*T2*TGS
!
! ...now T is T1*T2*T3*T4

call Ttotal(T1,T2,T3,T4,Ttot,Ttot_Inv,nBas)
!                                                                      *
!***********************************************************************
!                                                                      *
! Check that all this works by applaying TTot to original overlap
! matrix and you should get unit matrix
!
!call RecPrt('Old S',' ',SMatrix,nBas,nBas)
!call xxDGeMul(SMatrix,nBas,'N',Ttot,nBas,'N',tmp,nBas,nBas,nBas,nBas)
!call xxDGeMul(Ttot,nBas,'T',tmp,nBas,'N',SMatrix,nBas,nBas,nBas,nBas)
!call RecPrt('New S',' ',SMatrix,nBas,nBas)
!                                                                      *
!***********************************************************************
!                                                                      *
! Dealloctate memory

call mma_deallocate(T1)
call mma_deallocate(T2)
call mma_deallocate(T3)
call mma_deallocate(T4)
call mma_deallocate(tmp)
call mma_deallocate(S)
call mma_deallocate(SSave)
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out the transformation matrix as molecular orbitals.

nOcc = 0
do iBas=1,nBas
  if (iType(iBas) == Occ) nOcc = nOcc+1
end do

IndType(1) = 0
IndType(2) = nOcc
IndType(3) = 0
IndType(4) = 0
IndType(5) = 0
IndType(6) = nBas-nOcc
IndType(7) = 0

call mma_allocate(E,nBas,label='E')
E(:) = Zero
OrbName = 'LPRORB'
LuOut = 20
iUHF = 0
nSym = 1
Note = 'LoProp localized orbitals'
call WrVec_(OrbName,LuOut,'COEI',iUHF,nSym,[nBas],[nBas],TTot,[Zero],E,[Zero],E,[Zero],IndType,Note,0)
call mma_deallocate(E)
!                                                                      *
!***********************************************************************
!                                                                      *
! Not implemented for symmetry!

Filename = 'MD_LPR'
call Get_iScalar('nSym',nSym)
if (nSym == 1) call Molden_Interface(iUHF,OrbName,Filename)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Localize_LoProp
