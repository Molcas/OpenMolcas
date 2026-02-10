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
! Copyright (C) Yannick Carissan                                       *
!***********************************************************************

subroutine ComputeFunc(nAtoms,nOrb2Loc,PA,Functional, eval_func)
! Author: Y. Carissan

use Constants, only: Zero,One
use Definitions, only: wp, iwp, u6
use Localisation_globals, only: Debug

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nOrb2Loc
real(kind=wp), intent(in) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(out) :: Functional
logical(kind=iwp), intent(in) :: eval_func
integer(kind=iwp) :: iAt, iMO_s
real(kind=wp) :: thr, d_s

thr = 0.01

if (Debug) then
    write(u6,*) "In ComputeFunc: "
    write(u6,*) "-----------------"
end if

Functional = Zero
do iMO_s=1,nOrb2Loc
    ! calculating the Pipek measure of localization, that tells over how many atoms, MO s extends,
    d_s = Zero
    do iAt=1,nAtoms
        d_s = d_s + PA(iMO_s,iMO_s,iAt)**2
        !write(u6,"(A,I5,A,F10.4)") "Atom ", iAt," PA_ss = ", PA(iMO_s,iMO_s,iAt)
        !if (eval_func) then
        !    write(u6,"(A,I5,A,F10.4)") "Contribution from Atom ", iAt," Q_A = ", PA(iMO_s,iMO_s,iAt)
        !end if
    end do
    Functional = Functional + d_s

    if (d_s>1e-12) then
        d_s = 1/d_s
    else
        if (eval_func) then
            write(u6,*) "WARNING: d_s^{-1} is zero for MO ", iMO_s
        end if
        d_s = huge(d_s)
    end if
    if (eval_func) then
        write(u6,"(A,I4,A,F8.3,1X,A)") "MO ",iMO_s," extends over ",d_s, " atoms"
    end if
end do

if (eval_func) then
    !write(u6,"(A,F18.10)")  'Sum (d_s^-1)         = ',Functional
    !write(u6,"(//A,F8.4,A,F6.2,A)") "Functional / nOrb2Loc = ", Functional/nOrb2Loc, ", so a mean localisation of ", &
    !                Functional/nOrb2Loc*100, "% was reached"
    write(u6,"(/A, F6.2,A,/)") "Mean localisation (Functional/nOrb2Loc) = ", Functional/nOrb2Loc*100, " %"
end if

return

end subroutine ComputeFunc
