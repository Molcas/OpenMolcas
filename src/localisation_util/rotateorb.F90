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

subroutine RotateOrb(cMO,PACol,nBasis,nAtoms,PA,Maximisation,nOrb2loc,Name,nBas_per_Atom,nBas_Start,ThrRot,PctSkp,Debug)
! Author: Yannick Carissan.
!
! Modifications:
!    - October 6, 2005 (Thomas Bondo Pedersen):
!      Array PACol introduced in argument list.

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
#include "real.fh"
#include "Molcas.fh"
real*8 cMO(nBasis,*), PACol(nOrb2Loc,2)
integer nBas_per_Atom(*), nBas_Start(*)
real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
logical Maximisation, Debug
character*80 Txt
character*(LENIN8) Name(*), PALbl

xDone = 0.0d0
if (Debug) then
  write(6,*) 'RotateOrb[Debug]: nBas_per_Atom: ',(nBas_per_Atom(i),i=1,nAtoms)
  iCouple = 0
end if
do iMO1=1,nOrb2Loc-1
  do iMO2=iMO1+1,nOrb2Loc

    if (Debug) then
      iCouple = iCouple+1
      write(6,'(a9,i5)') 'Couple n:',iCouple
      write(6,'(a9,i5)') '    MO1 :',iMO1
      write(6,'(a9,i5)') '    MO2 :',iMO2
    end if

    iMO_s = iMO1
    iMO_t = iMO2
    SumA = Zero
    SumB = Zero
    do iAt=1,nAtoms
      PAst = PA(iMO_t,iMO_s,iAt)
      pass = PA(iMO_s,iMO_s,iAt)
      PAtt = PA(iMO_t,iMO_t,iAt)
      if (Debug) then
        write(6,*) 'In RotateOrb'
        write(6,*) '------------'
        PALbl = 'PA__'//Name(nBas_Start(iAt))(1:LENIN)
        call RecPrt(PALbl,' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
        write(6,*) '**************************'
        write(6,*) 'A :',iAt
        write(6,*) '<',iMO_s,'|PA|',iMO_t,'> = ',PAst
        write(6,*) '<',iMO_s,'|PA|',iMO_s,'> = ',pass
        write(6,*) '<',iMO_t,'|PA|',iMO_t,'> = ',PAtt
        write(6,*) '**************************'
      end if
      SumA = SumA+PAst**2-0.25d0*(pass-PAtt)**2
      SumB = SumB+PAst*(pass-PAtt)
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
    Tst = abs(cos4alpha)-1.0d0
    if (Tst > 0.0d0) then
      if (Tst > 1.0d-10) then
        write(Txt,'(A,D18.10)') 'Actual: cos4alpha = ',cos4alpha
        call SysAbendMsg('RotateOrb','-1.0d0 < cos4alpha < 1.0d0',Txt)
      else
        if (cos4alpha < 0.0d0) then
          cos4alpha = -1.0d0
        else
          cos4alpha = 1.0d0
        end if
      end if
    end if

    ! On choisit le cos car Alpha IN [0;PI/2]

    Alpha1 = acos(cos4alpha)/4.0d0
    Alpha2 = asin(sin4alpha)/4.0d0
    if (Alpha2 < Zero) Alpha1 = Alpha2+PI
    Alpha = Alpha1
    if (.not. Maximisation) then
      Gamma_rot = Alpha-PI/4.0d0
    else
      Gamma_rot = Alpha
    end if
    if (Debug) then
      write(6,'(a9,f10.5)') '   Ast :',Ast
      write(6,'(a9,f10.5)') '   Bst :',Bst
      write(6,'(a9,f10.5)') 'Alpha1 :',Alpha1
      write(6,'(a9,f10.5)') 'Alpha2 :',Alpha2
      write(6,'(a9,f10.5)') ' Gamma :',Gamma_rot
    end if

    Tsts = sin(Gamma_rot)
    Tstc = 1.0d0-cos(Gamma_rot)
    if ((abs(Tsts) > ThrRot) .or. (abs(Tstc) > ThrRot)) then
      call Rot_st(cMO(1,iMO_s),cMO(1,iMO_t),nBasis,Gamma_rot,Debug)
      call UpdateP(PACol,Name,nBas_Start,nOrb2Loc,nAtoms,PA,Gamma_rot,iMO_s,iMO_t,Debug)
      xDone = xDone+1.0d0
    end if

    if (Debug) then
      call RecPrt('MO after rotation',' ',cMO,nBasis,nBasis)
    end if

  end do !iMO2
end do !iMO1

if (nOrb2Loc > 1) then
  xOrb2Loc = dble(nOrb2Loc)
  xTotal = xOrb2Loc*(xOrb2Loc-1.0d0)/2.0d0
  PctSkp = 1.0d2*(xTotal-xDone)/xTotal
else
  PctSkp = 0.0d0
end if

return

end subroutine RotateOrb
