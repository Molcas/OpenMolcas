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
! Copyright (C) 2026, Lila Zapp                                        *
!***********************************************************************

subroutine RotateNxN(CMO,Ovlp,nOrb2Loc,nBasis,Ovlp_sqrt, Gradient, Hdiag, BName,nAtoms,nBas_per_Atom,nBas_Start,PA)

use definitions, only: wp,iwp,u6
use stdalloc, only: mma_allocate, mma_deallocate
use constants, only: Zero,One, Pi
use Molcas, only: LenIn8
use Localisation_globals, only: Debug, OptMeth

implicit none

integer(kind=iwp), intent(in) :: nAtoms, nBas_per_Atom(nAtoms), nBas_Start(nAtoms), nBasis, nOrb2Loc
real(kind=wp), intent(inout) :: CMO(nBasis,nOrb2Loc)
real(kind=wp), intent(in) :: Ovlp_sqrt(nBasis,nBasis),Gradient(nOrb2Loc,nOrb2Loc), Hdiag(nOrb2Loc,nOrb2Loc),Ovlp(nBasis,nBasis)
real(kind=wp), intent(out) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
character(len=LenIn8), intent(in) :: BName(nBasis)

real(kind=wp), allocatable :: kappa(:,:),kappa_cnt(:,:),xkappa_cnt(:,:), unitary_mat(:,:), rotated_CMO(:,:)
real(kind=wp), parameter :: thrsh_taylor = 1.0e-16_wp, alpha = 0.3
real(kind=wp) :: factor, ithrsh, DD, Thr
integer(kind=iwp) :: cnt,i,k, iBas
logical(kind=iwp), parameter :: debug_exp = .false.
real(kind=wp), External :: DDot_

call mma_Allocate(kappa,nOrb2Loc,nOrb2Loc,Label='kappa')
call mma_Allocate(kappa_cnt,nOrb2Loc,nOrb2Loc,Label='kappa_cnt') != kappa^cnt
call mma_Allocate(xkappa_cnt,nOrb2Loc,nOrb2Loc,Label='xkappa_cnt') !saves the previous kappa_cnt
call mma_Allocate(unitary_mat,nOrb2Loc,nOrb2Loc,Label='unitary_mat')
call mma_Allocate(rotated_cmo,nBasis,nOrb2Loc,Label='rotated_cmo')

! define the transformation matrix
kappa(:,:) = Zero
kappa_cnt(:,:) = Zero
xkappa_cnt(:,:) = Zero

if (OptMeth == 2) then
    kappa(:,:) = -Gradient(:,:)/Hdiag(:,:)
else if (OptMeth == 3) then
    kappa(:,:) = alpha*Gradient(:,:)
end if
DD=Sqrt(DDot_(nOrb2Loc**2,Kappa,1,Kappa,1))
Thr= 0.5E0_wp * Pi
If (DD>=Thr)Then
!           Write(6,*) 'Rescale Kappa(:,:)'
    Kappa(:,:) = (Thr/DD)*Kappa(:,:)
End If

kappa_cnt(:,:) = kappa !kappa^cnt = kappa since cnt=1
xkappa_cnt(:,:) = kappa_cnt

unitary_mat(:,:) = Zero
call unitmat(unitary_mat,nOrb2Loc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! analogous to exp_series in scf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cnt = 1
factor = One
ithrsh = 2.0e-16_wp

unitary_mat(:,:) =  unitary_mat(:,:) - kappa(:,:)

if (debug_exp) then
    write(u6,*) 'Taylor expansion: n=1'
    call RecPrt('unitary_mat = I - kappa^1',' ',unitary_mat(:,:), nOrb2Loc, nOrb2Loc)
    call RecPrt('kappa',' ',kappa(:,:), nOrb2Loc, nOrb2Loc)
    write(u6,*) 'Taylor expansion: more terms'
end if


do while (ithrsh > thrsh_taylor)

    !the number of the term = the exponent for kappa in that term
    cnt = cnt+1

    !the faculty value that the matrix will be divided by
    factor = factor*DBLE(cnt)

    !calculate the cnt'th exponent of the kappa matrix
    ! this works by multiplying the matrix from the previous term (kappa^cnt) by the
    ! initial kappa matrix (just kappa^1)
    ! C <= alpha*A*B + beta*C
    ! kappa_cnt <= 1*kappa_cnt*kappa + 0*kappa_cnt
    call dgemm_('N','N',nOrb2Loc,nOrb2Loc,nOrb2Loc,(One/DBLE(cnt)),xkappa_cnt,nOrb2Loc,kappa,nOrb2Loc,Zero,&
                kappa_cnt,norb2Loc)
    xkappa_cnt(:,:) = kappa_cnt

    ! differentiation of odd and even cases, because this expands exp(-kappa)
    ! all terms starting at n=2
    if (mod(cnt,2) == 0) then
        unitary_mat(:,:) =  unitary_mat + kappa_cnt(:,:)
        if (debug_exp) then
            write(u6,'(A,F10.1,A,I2,A,ES12.4)') 'term: + 1/',factor,' * kappa^',cnt, &
            ', current ithrsh = ', ithrsh
        end if
    else
        unitary_mat(:,:) =  unitary_mat - kappa_cnt(:,:)
        if (debug_exp) then
            write(u6,'(A,F10.1,A,I2,A,ES12.4)') 'term: - 1/',factor,' * kappa^',cnt, &
            ', current ithrsh = ', ithrsh
        end if
    end if

    ithrsh = maxval(abs(Kappa_Cnt(:,:))/(abs(unitary_mat)+thrsh_taylor))

    if (debug_exp) then
        write(u6,'(A,F10.1,A,I2,A,ES12.4)') 'term: + 1/',factor,' * kappa^',cnt, &
            ', current ithrsh = ', ithrsh
        call RecPrt('kappa^cnt',' ',kappa_cnt(:,:), nOrb2Loc, nOrb2Loc)
        call RecPrt('unitary_mat',' ',unitary_mat(:,:), nOrb2Loc, nOrb2Loc)
    end if
end do

if (debug) then
    write(u6,"(//A)") "rotating the orbitals with:"
    call RecPrt('kappa',' ',kappa(:,:), nOrb2Loc, nOrb2Loc)
    call RecPrt('unitary transformation matrix (exp(-kappa))',' ',unitary_mat(:,:), nOrb2Loc, nOrb2Loc)
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! transform the orbitals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rotated_CMO(:,:) = Zero

do iBas = 1, nBasis
    do k = 1,nOrb2Loc
        do i = 1,nOrb2Loc
            rotated_CMO(iBas,k) = rotated_CMO(iBas,k) + CMO(iBas,i) * unitary_mat(i,k)
        end do
    end do
end do

!reset CMO to be updated
CMO(:,:) = rotated_CMO(:,:)
call GenerateP(Ovlp,CMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Ovlp_sqrt)

call mma_Deallocate(kappa)
call mma_Deallocate(kappa_cnt)
call mma_Deallocate(xkappa_cnt)
call mma_Deallocate(unitary_mat)
call mma_Deallocate(rotated_CMO)

end subroutine RotateNxN
