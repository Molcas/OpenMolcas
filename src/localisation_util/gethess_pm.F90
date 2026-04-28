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
subroutine GetHess_PM(nAtoms,nOrb2Loc,PA,Hessian)
!
! Purpose: compute the full Hessian of the Pipek-Mezey functional w.r.t. elements of the kappa matrix

use Constants, only: Zero,One,Two, Eight
use Definitions, only: wp, iwp, u6
use Localisation_globals, only: Debug

implicit none

integer(kind=iwp), intent(in) :: nAtoms, nOrb2Loc
real(kind=wp), intent(in) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(out) :: Hessian(nOrb2Loc*(nOrb2Loc-1)/2,nOrb2Loc*(nOrb2Loc-1)/2)
integer(kind=iwp) :: iAtom, a,b,c,d,ab,cd
real(kind=wp) :: Q_aa, Q_bb, Q_cc, Q_dd, Q_ab, Q_ac, Q_bc, Q_ad, Q_bd, d_ac, d_ad, d_bc, d_bd

Q_aa = Zero
Q_bb = Zero
Q_cc = Zero
Q_dd = Zero

Q_ab = Zero
Q_ac = Zero
Q_ad = Zero
Q_bd = Zero
Q_bc = Zero

d_ac = One
d_ad = One
d_bc = One
d_bd = One


Hessian(:,:) = Zero
ab = 0
cd = 0

do a=1,nOrb2Loc-1
    do b=a+1,nOrb2Loc

        cd = 0
        ab = ab + 1 ! compound index rows

        ! get off diagonal elements
        do c = 1, nOrb2Loc-1
            do d = c+1,nOrb2Loc

                !reset Kronecker deltas
                d_ac = One
                d_ad = One
                d_bc = One
                d_bd = One

                ! evaluate Kronecker deltas
                if (a /=c ) d_ac = Zero
                if (a /=d ) d_ad = Zero
                if (b /=c ) d_bc = Zero
                if (b /=d ) d_bd = Zero

                cd = cd + 1 ! compound index columns

                do iAtom =1,nAtoms
                    Q_aa=PA(a,a,iAtom)
                    Q_bb=PA(b,b,iAtom)
                    Q_cc=PA(c,c,iAtom)
                    Q_dd=PA(d,d,iAtom)

                    Q_ab=PA(a,b,iAtom)
                    Q_ac=PA(a,c,iAtom)
                    Q_ad=PA(a,d,iAtom)
                    Q_bc=PA(b,c,iAtom)
                    Q_bd=PA(b,d,iAtom)

                    Hessian(ab,cd) = Hessian(ab,cd) + Eight * Q_ab * (d_ac*Q_ad - d_ad*Q_ac -d_bc*Q_bd + d_bd*Q_bc)
                    Hessian(ab,cd) = Hessian(ab,cd) + Two * d_ac * Q_bd * (Two*Q_aa - Q_bb - Q_dd)
                    Hessian(ab,cd) = Hessian(ab,cd) + Two * d_ad * Q_bc * (-Two*Q_aa + Q_bb + Q_cc)
                    Hessian(ab,cd) = Hessian(ab,cd) + Two * d_bc * Q_ad * (-Two*Q_bb + Q_aa + Q_dd)
                    Hessian(ab,cd) = Hessian(ab,cd) + Two * d_bd * Q_ac * (Two*Q_bb - Q_aa - Q_cc)
                end do

                !write(u6,"(A,4(I2),4X,A,2(I4))") "a,b,c,d=",a,b,c,d, "ab,cd=",ab,cd
            end do
        end do

    end do
end do

if (Debug) then
    write(u6,*) ' '
    write(u6,*) 'In GetHess_PM'
    write(u6,*) '-------------'
    call RecPrt("Hessian","",Hessian,nOrb2Loc*(nOrb2Loc-1)/2,nOrb2Loc*(nOrb2Loc-1)/2)
    write(u6,*) ' '
end if

end subroutine GetHess_PM
