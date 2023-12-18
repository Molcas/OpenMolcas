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
!               2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RotateOrb(cMO,PACol,nBasis,nAtoms,PA,Maximisation,nOrb2loc,BName,nBas_per_Atom,nBas_Start,ThrRot,PctSkp,Debug)
! Author: Yannick Carissan.
!
! Modifications:
!    - October 6, 2005 (Thomas Bondo Pedersen):
!      Array PACol introduced in argument list.

use Constants, only: Zero, One, Half, Quart, Pi
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: nBasis, nAtoms, nOrb2Loc, nBas_per_Atom(nAtoms), nBas_Start(nAtoms)
real(kind=wp), intent(inout) :: cMO(nBasis,*), PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(out) :: PACol(nOrb2Loc,2), PctSkp
logical(kind=iwp), intent(in) :: Maximisation, Debug
character(len=LenIn8), intent(in) :: BName(*)
real(kind=wp), intent(in) :: ThrRot
integer(kind=iwp) :: iAt, iCouple, iMO1, iMO2, iMO_s, iMO_t
real(kind=wp) :: Alpha, Alpha1, Alpha2, Ast, Bst, cos4alpha, Gamma_rot, PA_ss, PA_st, PA_tt, sin4alpha, SumA, SumB, Tst, Tstc, &
                 Tsts, xDone, xOrb2Loc, xTotal
character(len=LenIn8) :: PALbl
character(len=80) :: Txt

xDone = Zero
if (Debug) then
  write(u6,*) 'RotateOrb[Debug]: nBas_per_Atom: ',nBas_per_Atom(:)
  iCouple = 0
end if
do iMO1=1,nOrb2Loc-1
  do iMO2=iMO1+1,nOrb2Loc

    if (Debug) then
      iCouple = iCouple+1
      write(u6,'(a9,i5)') 'Couple n:',iCouple
      write(u6,'(a9,i5)') '    MO1 :',iMO1
      write(u6,'(a9,i5)') '    MO2 :',iMO2
    end if

    iMO_s = iMO1
    iMO_t = iMO2
    SumA = Zero
    SumB = Zero
    do iAt=1,nAtoms
      PA_st = PA(iMO_t,iMO_s,iAt)
      PA_ss = PA(iMO_s,iMO_s,iAt)
      PA_tt = PA(iMO_t,iMO_t,iAt)
      if (Debug) then
        write(u6,*) 'In RotateOrb'
        write(u6,*) '------------'
        PALbl = 'PA__'//BName(nBas_Start(iAt))(1:LenIn)
        call RecPrt(PALbl,' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
        write(u6,*) '**************************'
        write(u6,*) 'A :',iAt
        write(u6,*) '<',iMO_s,'|PA|',iMO_t,'> = ',PA_st
        write(u6,*) '<',iMO_s,'|PA|',iMO_s,'> = ',PA_ss
        write(u6,*) '<',iMO_t,'|PA|',iMO_t,'> = ',PA_tt
        write(u6,*) '**************************'
      end if
      SumA = SumA+PA_st**2-Quart*(PA_ss-PA_tt)**2
      SumB = SumB+PA_st*(PA_ss-PA_tt)
    end do
    Ast = SumA
    Bst = SumB

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
        call SysAbendMsg('RotateOrb','-1.0 < cos4alpha < 1.0',Txt)
      else
        if (cos4alpha < Zero) then
          cos4alpha = -One
        else
          cos4alpha = One
        end if
      end if
    end if

    ! On choisit le cos car Alpha IN [0,Pi/2]

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
      write(u6,'(a9,f10.5)') '   Ast :',Ast
      write(u6,'(a9,f10.5)') '   Bst :',Bst
      write(u6,'(a9,f10.5)') 'Alpha1 :',Alpha1
      write(u6,'(a9,f10.5)') 'Alpha2 :',Alpha2
      write(u6,'(a9,f10.5)') ' Gamma :',Gamma_rot
    end if

    Tsts = abs(sin(Gamma_rot))
    Tstc = One-abs(cos(Gamma_rot))
    if ((Tsts > ThrRot) .or. (Tstc > ThrRot)) then
      call Rot_st(cMO(1,iMO_s),cMO(1,iMO_t),nBasis,Gamma_rot,Debug)
      call UpdateP(PACol,BName,nBas_Start,nOrb2Loc,nAtoms,PA,Gamma_rot,iMO_s,iMO_t,Debug)
      xDone = xDone+One
    end if

    if (Debug) then
      call RecPrt('MO after rotation',' ',cMO,nBasis,nBasis)
    end if

  end do !iMO2
end do !iMO1

if (nOrb2Loc > 1) then
  xOrb2Loc = real(nOrb2Loc,kind=wp)
  xTotal = xOrb2Loc*(xOrb2Loc-One)*Half
  PctSkp = 1.0e2_wp*(xTotal-xDone)/xTotal
else
  PctSkp = Zero
end if

return

end subroutine RotateOrb
