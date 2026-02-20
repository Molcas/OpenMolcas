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

subroutine RotateNxN(CMO,kappa,nOrb2Loc,nBasis,nAtoms,kappa_cnt,xkappa_cnt,unitary_mat,rotated_CMO)
! this subroutine rotates the orbitals as CMO = exp(-kappa) * CMO
!
! for technical reasons, auxiliary matrices (kappa_cnt, xkappa_cnt, unitary_mat, rotated_CMO) are allocated outside of the loop
! that calls this routine. However it would work perfectly fine if these are allocated and deallocated within this routine to reduce
! the number of arguments (not recommended for speed)

use definitions, only: wp,iwp,u6
use stdalloc, only: mma_allocate, mma_deallocate
use constants, only: Zero,One
use Localisation_globals, only: Debug, OptMeth

implicit none

integer(kind=iwp), intent(in) :: nAtoms, nBasis, nOrb2Loc
real(kind=wp), intent(inout) :: CMO(nBasis,nOrb2Loc)
real(kind=wp), intent(inout) :: kappa(nOrb2Loc,nOrb2Loc),kappa_cnt(nOrb2Loc,nOrb2Loc),xkappa_cnt(nOrb2Loc,nOrb2Loc),&
                             unitary_mat(nOrb2Loc,nOrb2Loc), rotated_CMO(nBasis,nOrb2Loc)
real(kind=wp), parameter :: thrsh_taylor = 1.0e-16_wp
real(kind=wp) :: factor, ithrsh
integer(kind=iwp) :: cnt,i,k, iBas
logical(kind=iwp), parameter :: debug_exp = .false.

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

end subroutine RotateNxN
