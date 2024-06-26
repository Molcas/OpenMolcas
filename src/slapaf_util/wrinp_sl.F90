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

subroutine WrInp_sl()

use Slapaf_Info, only: Analytic_Hessian, AtomLbl, Baker, Beta, Beta_Disp, Coor, Cubic, Curvilinear, ddV_Schlegel, Delta, eMEPTest, &
                       FindTS, GNrm_Threshold, Header, HWRS, iOptC, iOptH, IRC, iRow, Line_Search, lNmHss, lOld, MEP, MEP_Algo, &
                       MEP_Type, Mode, MxItr, nMEP, nWndw, Redundant, rHidden, rMEP, ThrEne, ThrGrd
use kriging_mod, only: blaAI, blavAI, blAI, blvAI, Kriging, Max_Microiterations, mblAI, nD_In, set_l
use Constants, only: Two, auTokJmol
use Definitions, only: wp, iwp, u6

implicit none
#include "print.fh"
integer(kind=iwp) :: iPrint, iRout, nsAtom
real(kind=wp) :: Value_l

iRout = 3
iPrint = nPrint(iRout)

if (lNmHss) lOld = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 5) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  write(u6,*)
  write(u6,*)
  call CollapseOutput(1,'      Slapaf input parameters:')
  write(u6,'(3X,A)') '      ------------------------'
  write(u6,*)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  write(u6,'(A,I5)') ' Maximum number of iterations:             ',MxItr
  if (Baker) then
    write(u6,'(A)') ' Convergence test a la Baker.'
  else
    write(u6,'(A)') ' Convergence test a la Schlegel.'
  end if
  write(u6,'(A,ES8.1)') ' Convergence criterion on gradient/para.<=:',ThrGrd
  write(u6,'(A,ES8.1)') ' Convergence criterion on step/parameter<=:',ThrGrd
  write(u6,'(A,ES8.1)') ' Convergence criterion on energy change <=:',ThrEne
  write(u6,'(A)') ' Parameters for step-restricted optimization'
  if (.not. Kriging) then
    write(u6,'(A,ES9.2)') ' Maximum step length (initial seed):      ',Beta
  else
    write(u6,'(A,ES9.2)') ' Maximum step length (micro iterations):  ',Beta
  end if
  write(u6,*)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (Kriging) then
    write(u6,*) '-RVO activated with parameters:'
    !write(u6,'(A,I6)') '   GEK starts at iteration:                   ',nspAI
    write(u6,'(A,I6)') '   Maximum number of sample points (energies) used in GEK: ',nWndw/2
    write(u6,'(A,I6)') '   Maximum number of sample points (gradients) used in GEK: ',nWndw/2-nD_In
    !write(u6,'(A,I6)') '   Parameter of diff. for Matern (p):         ',pAI
    write(u6,'(A,I6)') '   Maximum number of micro iterations:        ',Max_Microiterations
    if (set_l) then
      call Get_dScalar('Value_l',Value_l)
      write(u6,*) '  Global characteristic length scale, l:     ',Value_l
    else
      write(u6,*) '  Individual characteristic length scales set to reproduce HMF Hessian.'
    end if

    if (blaAI) then
      write(u6,'(A,F10.5,A)') '   Baseline is highest energy plus: ',blavAI,' a.u'
    else
      if (mblAI) then
        write(u6,*) '  Baseline set to maximum value of the energy'
      else if (blAI) then
        write(u6,'(A,F9.5,A,/,A,F9.5,A)') '  Baseline (trend function) changed to value:',blvAI,'a.u.', &
                                          '                                             ',blvAI*auTokJmol,' kJ/mol'
      end if
    end if
    write(u6,'(A,F10.5,A)') '   Maximum dispersion accepted:     ',Beta_disp,' * abs(g.max.comp)'
  else
    write(u6,*) '-RFO activated with parameters:'
    write(u6,'(A,I6)') '   Maximum number of data points used in RFO: ',nWndw
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (Line_Search) then
    write(u6,'(A)') ' Line search is performed'
    write(u6,*)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (btest(iOptC,8)) then

    write(u6,'(1X,A)') '-Constrained optimization.'

    if (MEP) then
      if (IRC == 0) then
        write(u6,'(1X,A)') ' Minimum Energy Path (MEP) search'
      else if (IRC == 1) then
        write(u6,'(1X,A)') ' IRC forward search'
      else
        write(u6,'(1X,A)') ' IRC backward search'
      end if
      write(u6,'(1X,A,I5)') ' Maximum number of points:',nMEP
      if (eMEPtest) write(u6,'(1X,A)') ' Stop when energy increases'
      if (MEP_Algo == 'GS') then
        write(u6,'(1X,A)') ' MEP optimization algorithm: Gonzalez-Schlegel'
      else if (MEP_Algo == 'MB') then
        write(u6,'(1X,A)') ' MEP optimization algorithm: Mueller-Brown'
      end if
      if (MEP_Type == 'SPHERE    ') then
        write(u6,'(1X,A)') ' Type of constraint: Hypersphere'
      else if (MEP_Type == 'TRANSVERSE') then
        write(u6,'(1X,A)') ' Type of constraint: Hyperplane'
      end if
    end if

    if (rMEP) then
      write(u6,'(1X,A)') ' Reverse Minimum Energy Path (rMEP) search'
      write(u6,'(1X,A,I3)') ' Maximum number of points:',nMEP
      if (eMEPtest) write(u6,'(1X,A)') ' Stop when energy decreases'
      if (MEP_Type == 'SPHERE    ') then
        write(u6,'(1X,A)') ' Type of constraint: Hypersphere'
      else if (MEP_Type == 'TRANSVERSE') then
        write(u6,'(1X,A)') ' Type of constraint: Hyperplane'
      end if
    end if

    if (FindTS) then
      write(u6,'(1X,A)') '-The optimization will home in on a transition state if:'
      write(u6,'(A)') '  a) Negative curvature is encountered, and'
      write(u6,'(A,F10.4)') '  b) the norm of the gradient is below:',GNrm_Threshold
      if (btest(iOptC,9)) then
        write(u6,'(A)') '  TS-search by RS-I-RFO.'
        !if (Kriging) then
        !  write(u6,'(A)') '  TS-search by RV-I-RFO.'
        !else
        !  write(u6,'(A)') '  TS-search by RS-I-RFO.'
        !end if
      else
        write(u6,'(A)') '  TS-search by RS-P-RFO.'
        !if (Kriging) then
        !  write(u6,'(A)') '  TS-search by RV-P-RFO.'
        !else
        !  write(u6,'(A)') '  TS-search by RS-P-RFO.'
        !end if
      end if
    end if

  end if
  write(u6,*)

  if (btest(iOptC,7)) then
    write(u6,'(1X,A)') '-Optimization for minimum.'
    if (btest(iOptC,0)) then
      write(u6,'(A)') '  Optimization method: quasi-NR.'
    else if (btest(iOptC,1)) then
      write(u6,'(A)') '  Optimization method: C1-DIIS.'
    else if (btest(iOptC,2)) then
      write(u6,'(A)') '  Optimization method: C2-DIIS.'

    else if (btest(iOptC,3)) then
      if (Kriging) then
        write(u6,'(A)') '  Optimization method: RVO.'
      else
        write(u6,'(A)') '  Optimization method: RS-RFO.'
      end if
    else
      call WarningMessage(2,' WrInp: Wrong iOptC setting!')
      write(u6,*) ' iOptC=',iOptC
      call Abend()
    end if
  else
    write(u6,'(1X,A)') '-Optimization for transition state.'
    if (btest(iOptC,9)) then
      write(u6,'(A)') '  Optimization method: RS-I-RFO'
      !if (Kriging) then
      !  write(u6,'(A)') '  Optimization method: RV-I-RFO'
      !else
      !  write(u6,'(A)') '  Optimization method: RS-I-RFO'
      !end if
    else
      write(u6,'(A)') '  Optimization method: RS-P-RFO'
      !if (Kriging) then
      !  write(u6,'(A)') '  Optimization method: RV-P-RFO'
      !else
      !  write(u6,'(A)') '  Optimization method: RS-P-RFO'
      !end if
    end if
    if (Mode > 0) then
      write(u6,'(A,I2)') '  Original mode to follow:',Mode
    else
      write(u6,'(A)') '  No mode to follow is specified!'
      write(u6,'(A)') '  Optimization will follow mode with the lowest eigenvalue.'
    end if
  end if
  write(u6,*)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Special options on DIIS

  if (btest(iOptC,1) .or. btest(iOptC,2)) then
    if (btest(iOptC,4)) then
      write(u6,'(1X,A)') '-DIIS based on <dx|dx>.'
    else if (btest(iOptC,5)) then
      write(u6,'(1X,A)') '-DIIS based on <g|dx>.'
    else if (btest(iOptC,6)) then
      write(u6,'(1X,A)') '-DIIS based on <g|g>.'
    else
      call WarningMessage(2,' WrInp: Wrong iOptC setting!')
      write(u6,*) ' iOptC=',iOptC
      call Abend()
    end if
    write(u6,*)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Hessian information

  if (Analytic_Hessian) then
    write(u6,'(1X,A)') '-The Hessian is analytic.'
    write(u6,'(1X,A)') ' Hessian from either input or runfile.'
  else
    if ((.not. lOld) .and. (.not. lNmHss)) then
      if (DDV_Schlegel) then
        write(u6,'(1X,A)') '-Initial Hessian guessed a la Schlegel.'
      else
        if (Kriging) then
          write(u6,'(1X,A)') '-Hessian guessed by Kriging surrogate surface.'
        else
          write(u6,'(1X,A)') '-Initial Hessian guessed by Hessian Model Function (HMF).'
          if (btest(iOptC,10)) write(u6,'(A)') '  HMF augmented with weak interactions.'
        end if
      end if
    else if (lOld .and. (.not. lNmHss)) then
      write(u6,'(1X,A)') '-Initial Hessian guess was read from a RUNFILE file.'
    else
      write(u6,'(1X,A,/,A,ES9.2)') '-Initial Hessian guess is estimated with finite differences.', &
                                   '    Two point symmetric formula, Delta=',Delta
      if (Cubic) write(u6,'(1X,A)') '-Cubic force constants evaluated numerically.'
    end if
  end if
  write(u6,*)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Hessian update method

  if (.not. Kriging) then
    if (btest(iOptH,0)) then
      write(u6,'(1X,A)') '-Hessian update method: Fletcher-Meyer'
    else if (btest(iOptH,1)) then
      write(u6,'(1X,A)') '-Hessian update method: Broyden-Powell'
    else if (btest(iOptH,2)) then
      write(u6,'(1X,A)') '-Hessian update method: Broyden-Fletcher-Goldfarb-Shanno'
    else if (btest(iOptH,3)) then
      write(u6,'(1X,A)') '-Hessian update method: none'
    else if (btest(iOptH,4)) then
      write(u6,'(1X,A)') '-Hessian update method: Murtagh-Sargent-Powell'
    else if (btest(iOptH,5)) then
      write(u6,'(1X,A)') '-Hessian update method: EU update by Bofill'
    else if (btest(iOptH,6)) then
      write(u6,'(1X,A)') '-Hessian update method: TS-BFGS update by Bofill'
    else
      call WarningMessage(2,' WrInp: Wrong iOptH setting!')
      write(u6,*) ' Nonrecognizable iOptH setting:',iOptH
      call Abend()
    end if
    if (.not. btest(iOptH,3)) write(u6,'(A,I3)') '  Maximum number of points in Hessian update:',nWndw
    write(u6,*)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (rHidden >= Two) then
    write(u6,'(1X,A,/,1X,A,F6.2,A)') '-Improved QM/MM Hessian.',' Hidden atoms until ',rHidden,' bohrs are included.'
    write(u6,*)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (iRow > 0) then
    if (Redundant) then
      write(u6,'(1X,A)') '-Relaxation will be done in user supplied redundant internal coordinates.'
    else
      write(u6,'(1X,A)') '-Relaxation will be done in user supplied non-redundant internal coordinates.'
    end if
  else
    if (CurviLinear) then
      if (HWRS) then
        if (Redundant) then
          write(u6,'(1X,A)') '-Relaxation will be done on redundant internal coordinates, based on'
          write(u6,*) ' force constant weighted redundant internal coordinates.'
        else
          write(u6,'(1X,A)') '-Relaxation will be done on non-redundant internal coordinates, based on'
          write(u6,*) ' force constant weighted redundant internal coordinates.'
        end if
      else
        if (Redundant) then
          write(u6,'(1X,A)') '-Relaxation will be done in redundant delocalized internal coordinates.'
        else
          write(u6,'(1X,A)') '-Relaxation will be done in non-redundant delocalized internal coordinates.'
        end if
      end if
    else
      if (Redundant) then
        write(u6,'(1X,A)') '-Relaxation will be done in redundant Cartesian coordinates.'
      else
        write(u6,'(1X,A)') '-Relaxation will be done in approximate non-redundant Cartesian normal mode coordinates.'
      end if
    end if
  end if
  write(u6,*)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 6) then
  nsAtom = size(Coor,2)
  write(u6,*)
  write(u6,'(A)') ' Header from ONEINT:'
  call Banner(Header,2,len(Header(1))+12)
  write(u6,*)
  call PrList('Symmetry Distinct Nuclear Coordinates / bohr',AtomLbl,nsAtom,Coor,3,nsAtom)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 5) call CollapseOutput(0,'      Slapaf input parameters:')

return

end subroutine WrInp_sl
