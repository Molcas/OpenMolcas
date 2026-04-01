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

subroutine expkap_localisation(kappa,nOrb2Loc,kappa_cnt,xkappa_cnt,unitary_mat)
! analogous to exp_series in scf, only for the specified orbital subspace
!
! purpose: returns exp(-kappa) as unitary_mat
!
! kappa_cnt and xkappa_cnt are auxiliary matrices of dim (norb2loc,norb2loc) required for the Taylor expansion.
! they need to be allocated outside of this subroutine

use definitions, only: wp,iwp,u6
use constants, only: Zero,One
use Localisation_globals, only: Debug

implicit none

integer(kind=iwp), intent(in) :: nOrb2Loc
real(kind=wp), intent(inout) :: kappa(nOrb2Loc,nOrb2Loc),kappa_cnt(nOrb2Loc,nOrb2Loc),xkappa_cnt(nOrb2Loc,nOrb2Loc),&
                             unitary_mat(nOrb2Loc,nOrb2Loc)
real(kind=wp), parameter :: thrsh_taylor = 1.0e-18_wp
real(kind=wp) :: factor, ithrsh
integer(kind=iwp) :: cnt,maxel(2)
logical(kind=iwp), parameter :: debug_exp = .false.
real(kind=wp),External :: DDot_

kappa_cnt(:,:) = kappa !kappa^cnt = kappa since cnt=1
xkappa_cnt(:,:) = kappa_cnt

unitary_mat(:,:) = Zero
call unitmat(unitary_mat,nOrb2Loc)

cnt = 1
factor = One
ithrsh = One

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
    xkappa_cnt(:,:) = kappa_cnt(:,:)

    ! differentiation of odd and even cases, because this expands exp(-kappa)
    ! all terms starting at n=2
    if (mod(cnt,2) == 0) then
        unitary_mat(:,:) =  unitary_mat + kappa_cnt(:,:)
        if (debug_exp) then
            write(u6,'(A,ES10.1,A,I2,A,ES12.4)') 'term: +',1/factor,' * kappa^',cnt, &
            ', current ithrsh = ', ithrsh
        end if
    else
        unitary_mat(:,:) =  unitary_mat - kappa_cnt(:,:)
        if (debug_exp) then
            write(u6,'(A,ES10.1,A,I2,A,ES12.4)') 'term: -',1/factor,' * kappa^',cnt, &
            ', current ithrsh = ', ithrsh
        end if
    end if

    ithrsh = sqrt(DDot_(nOrb2Loc**2,Kappa_Cnt(:,:),1,Kappa_Cnt(:,:),1))

    ! sanity check for divergence
    if (ithrsh/720 > One) then
        write(u6,*) "Bug: elements of the kappa matrix fed to expkap_localisation() are too large",&
                    "- the Taylor expansion diverges"
        write(u6,*) "Stopping Taylor expansion at ",cnt,"-th term. Rescale kappa before feeding it this subroutine."
        call Abend()
    end if

    if (debug_exp) then
        write(u6,'(A,ES10.1,A,I2,A,ES12.4)') 'term:  ',1/factor,' * kappa^',cnt, &
            ', new ithrsh     = ', ithrsh
        call RecPrt('kappa^cnt',' ',kappa_cnt(:,:), nOrb2Loc, nOrb2Loc)
        call RecPrt('unitary_mat',' ',unitary_mat(:,:), nOrb2Loc, nOrb2Loc)
    end if

    maxel(:) = maxloc(abs(unitary_mat(:,:)),2)
    if (unitary_mat(maxel(1),maxel(2)) > One) then
        write(u6,*) "element of U bigger cannot be larger than 1, as U is supposed to be orthogonal/unitary:", &
                     unitary_mat(maxel(1),maxel(2))
        call Abend()
    end if
end do
if (debug) then
    write(u6,"(//A)") "rotating the orbitals with:"
    call RecPrt('kappa',' ',kappa(:,:), nOrb2Loc, nOrb2Loc)
    call RecPrt('unitary transformation matrix (exp(-kappa))',' ',unitary_mat(:,:), nOrb2Loc, nOrb2Loc)
end if
end subroutine expkap_localisation
