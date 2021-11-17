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

subroutine Move_EC(rMP,EC,lMax,nij,nElem,Coor,nAtoms,Q_Nuc,C_o_C,nPert,Bond_Threshold,iANr,T_Values,iT_Sets,iWarnings, &
                   Num_Warnings,Opt_Method,iPlot)

use stdalloc, only: mma_allocate, mma_deallocate
use Real_Spherical, only: Sphere, Sphere_Free
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lMax, nij, nElem, nAtoms, nPert, iANr(nAtoms), iPlot
real(kind=wp), intent(inout) :: rMP(nij,0:nElem-1,0:nPert-1), EC(3,nij)
real(kind=wp), intent(in) :: Coor(3,nAtoms), Q_Nuc(nAtoms), C_o_C(3), Bond_Threshold
real(kind=wp), intent(out) :: T_Values(nij)
integer(kind=iwp), intent(inout) :: iT_Sets(nij), Num_Warnings
integer(kind=iwp), intent(out) :: iWarnings(nij)
character(len=12) :: Opt_Method
integer(kind=iwp) :: iAtom, ii, ij, iPert, jAtom, jj
real(kind=wp) :: CX(3), Dipole_Rot_A, Dipole_Rot_B, Dipole_Rot_AB, R_A, R_B
logical(kind=iwp) :: Bond_OK, Check_Bond
real(kind=wp), allocatable :: A(:,:), B(:,:), Scratch_New(:), Scratch_Org(:), xnrMP(:,:), xrMP(:,:), xxrMP(:,:)

call Sphere(lMax)

call mma_allocate(B,3,nij,label='B')
call mma_allocate(xnrMP,nij,nElem,label='xnrMP')

B(:,:) = Zero

xnrMP(:,:) = Zero
ij = 0
do iAtom=1,nAtoms
  do jAtom=1,iAtom
    ij = ij+1
    if (iAtom == jAtom) then
      xnrMP(ij,1) = Q_nuc(iAtom)
      CX(1) = Coor(1,iAtom)
      CX(2) = Coor(2,iAtom)
      CX(3) = Coor(3,iAtom)
      call ReExpand(xnrMP,nij,nElem,CX,C_o_C,ij,lMax)
    end if
  end do
end do

call mma_allocate(A,3,nij,label='A')
call mma_allocate(xrMP,nij,nElem,label='xrMP')
call mma_allocate(xxrMP,nij,nElem,label='xxrMP')
call mma_allocate(Scratch_Org,nij*(2+lMax+1),label='Scratch_Org')
call mma_allocate(Scratch_New,nij*(2+lMax+1),label='Scratch_New')

do iAtom=1,nAtoms
  ii = iAtom*(iAtom+1)/2
  do jAtom=iAtom,1,-1
    jj = jAtom*(jAtom+1)/2
    ij = iAtom*(iAtom-1)/2+jAtom
    A(:,:) = EC(:,:)

    if (iAtom == jAtom) then
      ! Atomic domain
      ! Do nothing for now - put appropriate expresion in later.
      iT_sets(ij) = 1
      T_values(ij) = Zero
      B(:,ij) = EC(:,ij)
      !?? end if
    else
      ! Bond domain
      !
      ! Check whether bond is too long - if so skip it.

      Bond_OK = Check_Bond(EC(:,ii),EC(:,jj),iANr(iAtom),iANr(jAtom),Bond_Threshold)
      if (.not. Bond_OK) then
        T_values(ij) = Zero
        B(:,ij) = EC(:,ij)
      else

        ! Bond is real bond - proceed

        iT_sets(ij) = 1
        if (Opt_Method(1:9) == 'Multipole') then
          call Rotate_Dipole(rMP,EC,nij,nElem,nPert,ij,ii,jj,Dipole_Rot_A,Dipole_Rot_B,Dipole_Rot_AB,R_A,R_B)

          call Find_Dipole_Center(rMP(ii,0,0),rMP(jj,0,0),Dipole_Rot_A,Dipole_Rot_B,xnrMP(ii,1),xnrMP(jj,1),R_A,R_B,T_Values(ij), &
                                  iPlot)

          B(:,ij) = EC(:,ij)+T_Values(ij)*(EC(:,jj)-EC(:,ii))
        else
          call Min_Mult_Error(EC,A,B,EC(:,ii),EC(:,jj),rMP,xrMP,xxrMP,xnrMP,lMax,nij,nElem,iAtom,jAtom,nAtoms,nPert,C_o_C, &
                              Scratch_New,Scratch_Org,iPlot,T_Values(ij),iWarnings(ij),Num_Warnings)
        end if
      end if

    end if
  end do
end do

call mma_deallocate(A)
call mma_deallocate(xrMP)
call mma_deallocate(xxrMP)
call mma_deallocate(xnrMP)
call mma_deallocate(Scratch_Org)
call mma_deallocate(Scratch_New)

do ij=1,nij
  do iPert=0,nPert-1
    call ReExpand(rMP(:,:,iPert),nij,nElem,EC(:,ij),B(:,ij),ij,lMax)
  end do
end do
EC(:,:) = B(:,:)

call mma_deallocate(B)

call Sphere_Free()

return

end subroutine Move_EC
