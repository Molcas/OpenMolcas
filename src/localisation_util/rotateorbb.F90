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

subroutine RotateOrbB(CMO,Col,ipLbl,nComp,nBas,nOrb2Loc,Maximisation,ThrRot,PctSkp,Debug)
! Author: T.B. Pedersen
!
! Purpose: rotate orbitals (Jacobi Sweeps) for Boys localisation.

implicit real*8(a-h,o-z)
real*8 CMO(nBas,*), Col(nOrb2Loc,2)
integer ipLbl(nComp)
logical Maximisation, Debug
#include "WrkSpc.fh"
#include "real.fh"
character*10 SecNam
parameter(SecNam='RotateOrbB')
character*80 Txt

xDone = 0.0d0
if (Debug) iCouple = 0
do iMO1=1,nOrb2Loc-1
  do iMO2=iMO1+1,nOrb2Loc

    if (Debug) then
      iCouple = iCouple+1
      write(6,'(A9,I5)') 'Couple n:',iCouple
      write(6,'(A9,I5)') '    MO1 :',iMO1
      write(6,'(A9,I5)') '    MO2 :',iMO2
    end if

    iMO_s = iMO1
    iMO_t = iMO2

    Ast = 0.0d0
    Bst = 0.0d0
    do iComp=1,nComp
      ip0 = ipLbl(iComp)-1
      iss = ip0+nOrb2Loc*(iMO_s-1)+iMO_s
      itt = ip0+nOrb2Loc*(iMO_t-1)+iMO_t
      ist = ip0+nOrb2Loc*(iMO_t-1)+iMO_s
      Ast = Ast+Work(ist)**2-2.5d-1*(Work(iss)-Work(itt))**2
      Bst = Bst+Work(ist)*(Work(iss)-Work(itt))
    end do

    if ((Ast == 0.0d0) .and. (Bst == 0.0d0)) then
      cos4alpha = -1.0d0
      sin4alpha = 0.0d0
    else
      cos4alpha = -Ast/sqrt(Ast**2+Bst**2)
      sin4alpha = Bst/sqrt(Ast**2+Bst**2)
    end if
    Tst = abs(cos4alpha)-1.0d0
    if (Tst > 0.0d0) then
      if (Tst > 1.0d-10) then
        write(Txt,'(A,D18.10)') 'Actual: cos4alpha = ',cos4alpha
        call SysAbendMsg(SecNam,'-1.0d0 < cos4alpha < 1.0d0',Txt)
      else
        if (cos4alpha < 0.0d0) then
          cos4alpha = -1.0d0
        else
          cos4alpha = 1.0d0
        end if
      end if
    end if

    Alpha1 = acos(cos4alpha)/4.0d0
    Alpha2 = asin(sin4alpha)/4.0d0
    if (Alpha2 < 0.0d0) Alpha1 = Alpha2+PI
    Alpha = Alpha1
    if (.not. Maximisation) then
      Gamma_rot = Alpha-PI/4.0d0
    else
      Gamma_rot = Alpha
    end if
    if (Debug) then
      write(6,'(A9,F10.5)') '   Ast :',Ast
      write(6,'(A9,F10.5)') '   Bst :',Bst
      write(6,'(A9,F10.5)') 'Alpha1 :',Alpha1
      write(6,'(A9,F10.5)') 'Alpha2 :',Alpha2
      write(6,'(A9,F10.5)') ' Gamma :',Gamma_rot
    end if

    Tsts = sin(Gamma_rot)
    Tstc = 1.0d0-cos(Gamma_rot)
    if ((abs(Tsts) > ThrRot) .or. (abs(Tstc) > ThrRot)) then
      call Rot_st(CMO(1,iMO_s),CMO(1,iMO_t),nBas,Gamma_rot,Debug)
      call UpdateB(Col,nOrb2Loc,ipLbl,nComp,Gamma_rot,iMO_s,iMO_t,Debug)
      xDone = xDone+1.0d0
    end if

  end do
end do

if (nOrb2Loc > 1) then
  xOrb2Loc = dble(nOrb2Loc)
  xTotal = xOrb2Loc*(xOrb2Loc-1.0d0)/2.0d0
  PctSkp = 1.0d2*(xTotal-xDone)/xTotal
else
  PctSkp = 0.0d0
end if

end subroutine RotateOrbB
