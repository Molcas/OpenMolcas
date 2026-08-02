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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine RotateOrbB(CMO,Col,Lbl,nComp,nBas,nOrb2Loc,Maximisation,ThrRot,PctSkp,Debug)
! Author: T.B. Pedersen
!
! Purpose: rotate orbitals (Jacobi Sweeps) for Boys localisation.

use Constants, only: Zero, One, Half, Quart, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nComp, nBas, nOrb2Loc
real(kind=wp), intent(inout) :: CMO(nBas,*), Lbl(nOrb2Loc,nOrb2Loc,nComp)
real(kind=wp), intent(out) :: Col(nOrb2Loc,2), PctSkp
logical(kind=iwp), intent(in) :: Maximisation, Debug
real(kind=wp), intent(in) :: ThrRot
integer(kind=iwp) :: iComp, iCouple, iMO1, iMO2, iMO_s, iMO_t
real(kind=wp) :: Alpha, Alpha1, Alpha2, Ast, Bst, cos4alpha, Gamma_rot, sin4alpha, Tst, Tstc, Tsts, xDone, xOrb2Loc, xTotal
character(len=80) :: Txt
character(len=*), parameter :: SecNam = 'RotateOrbB'

xDone = Zero
if (Debug) iCouple = 0
do iMO1=1,nOrb2Loc-1
  do iMO2=iMO1+1,nOrb2Loc

    if (Debug) then
      iCouple = iCouple+1
      write(u6,'(A9,I5)') 'Couple n:',iCouple
      write(u6,'(A9,I5)') '    MO1 :',iMO1
      write(u6,'(A9,I5)') '    MO2 :',iMO2
    end if

    iMO_s = iMO1
    iMO_t = iMO2

    Ast = Zero
    Bst = Zero
    do iComp=1,nComp
      Ast = Ast+Lbl(iMO_s,iMO_t,iComp)**2-Quart*(Lbl(iMO_s,iMO_s,iComp)-Lbl(iMO_t,iMO_t,iComp))**2
      Bst = Bst+Lbl(iMO_s,iMO_t,iComp)*(Lbl(iMO_s,iMO_s,iComp)-Lbl(iMO_t,iMO_t,iComp))
    end do

    if ((Ast == Zero) .and. (Bst == Zero)) then
      cos4alpha = -One
      sin4alpha = Zero
    else
      cos4alpha = -Ast/sqrt(Ast**2+Bst**2)
      sin4alpha = Bst/sqrt(Ast**2+Bst**2)
    end if
    Tst = abs(cos4alpha)-One
    if (Tst > Zero) then
      if (Tst > 1.0e-10_wp) then
        write(Txt,'(A,ES18.10)') 'Actual: cos4alpha = ',cos4alpha
        call SysAbendMsg(SecNam,'-1.0 < cos4alpha < 1.0',Txt)
      else
        if (cos4alpha < Zero) then
          cos4alpha = -One
        else
          cos4alpha = One
        end if
      end if
    end if

    Alpha1 = acos(cos4alpha)*Quart
    Alpha2 = asin(sin4alpha)*Quart
    if (Alpha2 < Zero) Alpha1 = Alpha2+Pi
    Alpha = Alpha1
    if (.not. Maximisation) then
      Gamma_rot = Alpha-Pi*Quart
    else
      Gamma_rot = Alpha
    end if
    if (Debug) then
      write(u6,'(A9,F10.5)') '   Ast :',Ast
      write(u6,'(A9,F10.5)') '   Bst :',Bst
      write(u6,'(A9,F10.5)') 'Alpha1 :',Alpha1
      write(u6,'(A9,F10.5)') 'Alpha2 :',Alpha2
      write(u6,'(A9,F10.5)') ' Gamma :',Gamma_rot
    end if

    Tsts = abs(sin(Gamma_rot))
    Tstc = One-abs(cos(Gamma_rot))
    if ((Tsts > ThrRot) .or. (Tstc > ThrRot)) then
      call Rot_st(CMO(1,iMO_s),CMO(1,iMO_t),nBas,Gamma_rot,Debug)
      call UpdateB(Col,nOrb2Loc,Lbl,nComp,Gamma_rot,iMO_s,iMO_t,Debug)
      xDone = xDone+One
    end if

  end do
end do

if (nOrb2Loc > 1) then
  xOrb2Loc = real(nOrb2Loc,kind=wp)
  xTotal = xOrb2Loc*(xOrb2Loc-One)*Half
  PctSkp = 1.0e2_wp*(xTotal-xDone)/xTotal
else
  PctSkp = Zero
end if

end subroutine RotateOrbB
