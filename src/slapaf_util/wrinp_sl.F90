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

use kriging_mod
use Slapaf_Info, only: Coor, AtomLbl
use Slapaf_Parameters, only: iRow, ddV_Schlegel, HWRS, iOptH, IRC, Curvilinear, Redundant, FindTS, Analytic_Hessian, iOptC, &
                             rHidden, lOld, Beta, Beta_Disp, Line_Search, GNrm_Threshold, Mode, ThrEne, ThrGrd, Baker, eMEPTest, &
                             rMEP, MEP, nMEP, MEP_Type, MEP_Algo, Header, Delta, lNmHss, Cubic, MxItr, nWndw

implicit real*8(a-h,o-z)
#include "real.fh"
#include "print.fh"
#include "constants.fh"

iRout = 3
iPrint = nPrint(iRout)

Lu = 6

if (lNmHss) then
  lOld = .false.
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 5) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  write(Lu,*)
  write(Lu,*)
  call CollapseOutput(1,'      Slapaf input parameters:')
  write(Lu,'(3X,A)') '      ------------------------'
  write(Lu,*)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  write(Lu,'(A,I5)') ' Maximum number of iterations:             ',MxItr
  if (Baker) then
    write(Lu,'(A)') ' Convergence test a la Baker.'
  else
    write(Lu,'(A)') ' Convergence test a la Schlegel.'
  end if
  write(Lu,'(A,ES8.1)') ' Convergence criterion on gradient/para.<=:',ThrGrd
  write(Lu,'(A,ES8.1)') ' Convergence criterion on step/parameter<=:',ThrGrd
  write(Lu,'(A,ES8.1)') ' Convergence criterion on energy change <=:',ThrEne
  write(Lu,'(A)') ' Parameters for step-restricted optimization'
  if (.not. Kriging) then
    write(Lu,'(A,ES9.2)') ' Maximum step length (initial seed):      ',Beta
  else
    write(Lu,'(A,ES9.2)') ' Maximum step length (micro iterations):  ',Beta
  end if
  write(Lu,*)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (Kriging) then
    write(Lu,*) '-RVO activated with parameters:'
    !write(Lu,'(A,I6)') '   GEK starts at iteration:                   ',nspAI
    write(Lu,'(A,I6)') '   Maximum number of sample points (energies) used in GEK: ',nWndw/2
    write(Lu,'(A,I6)') '   Maximum number of sample points (gradients) used in GEK: ',nWndw/2-nD_In
    !write(Lu,'(A,I6)') '   Parameter of diff. for Matern (p):         ',pAI
    write(Lu,'(A,I6)') '   Maximum number of micro iterations:        ',Max_Microiterations
    if (set_l) then
      call Get_dScalar('Value_l',Value_l)
      write(Lu,*) '  Global characteristic length scale, l:     ',Value_l
    else
      write(Lu,*) '  Individual characteristic length scales set to reproduce HMF Hessian.'
    end if

    if (blaAI) then
      write(6,'(A,F10.5,A)') '   Baseline is highest energy plus: ',blavAI,' a.u'
    else
      if (mblAI) then
        write(6,*) '  Baseline set to maximum value of the energy'
      else if (blAI) then
        write(6,'(A,F9.5,A,/,A,F9.5,A)') '  Baseline (trend function) changed to value:',blvAI,'a.u.', &
                                         '                                             ', &
                                         blvAI*CONV_AU_TO_KJ_PER_MOLE_,' kJ/mol'
      end if
    end if
    write(6,'(A,F10.5,A)') '   Maximum dispersion accepted:     ',Beta_disp,' * abs(g.max.comp)'
  else
    write(Lu,*) '-RFO activated with parameters:'
    write(Lu,'(A,I6)') '   Maximum number of data points used in RFO: ',nWndw
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (Line_Search) then
    write(Lu,'(A)') ' Line search is performed'
    write(Lu,*)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (iand(iOptC,256) == 256) then

    write(Lu,'(1X,A)') '-Constrained optimization.'

    if (MEP) then
      if (IRC == 0) then
        write(Lu,'(1X,A)') ' Minimum Energy Path (MEP) search'
      else if (IRC == 1) then
        write(Lu,'(1X,A)') ' IRC forward search'
      else
        write(Lu,'(1X,A)') ' IRC backward search'
      end if
      write(Lu,'(1X,A,I5)') ' Maximum number of points:',nMEP
      if (eMEPtest) write(Lu,'(1X,A)') ' Stop when energy increases'
      if (MEP_Algo == 'GS') then
        write(Lu,'(1X,A)') ' MEP optimization algorithm: Gonzalez-Schlegel'
      else if (MEP_Algo == 'MB') then
        write(Lu,'(1X,A)') ' MEP optimization algorithm: Mueller-Brown'
      end if
      if (MEP_Type == 'SPHERE    ') then
        write(Lu,'(1X,A)') ' Type of constraint: Hypersphere'
      else if (MEP_Type == 'TRANSVERSE') then
        write(Lu,'(1X,A)') ' Type of constraint: Hyperplane'
      end if
    end if

    if (rMEP) then
      write(Lu,'(1X,A)') ' Reverse Minimum Energy Path (rMEP) search'
      write(Lu,'(1X,A,I3)') ' Maximum number of points:',nMEP
      if (eMEPtest) write(Lu,'(1X,A)') ' Stop when energy decreases'
      if (MEP_Type == 'SPHERE    ') then
        write(Lu,'(1X,A)') ' Type of constraint: Hypersphere'
      else if (MEP_Type == 'TRANSVERSE') then
        write(Lu,'(1X,A)') ' Type of constraint: Hyperplane'
      end if
    end if

    if (FindTS) then
      write(Lu,'(1X,A)') '-The optimization will home in on a transition state if:'
      write(Lu,'(A)') '  a) Negative curvature is encountered, and'
      write(Lu,'(A,F10.4)') '  b) the norm of the gradient is below:',GNrm_Threshold
      if (iand(iOptC,512) == 512) then
        write(Lu,'(A)') '  TS-search by RS-I-RFO.'
        !if (Kriging) then
        !  write(Lu,'(A)') '  TS-search by RV-I-RFO.'
        !else
        !  write(Lu,'(A)') '  TS-search by RS-I-RFO.'
        !end if
      else
        write(Lu,'(A)') '  TS-search by RS-P-RFO.'
        !if (Kriging) then
        !  write(Lu,'(A)') '  TS-search by RV-P-RFO.'
        !else
        !  write(Lu,'(A)') '  TS-search by RS-P-RFO.'
        !end if
      end if
    end if

  end if
  write(Lu,*)

  if (iand(iOptC,128) == 128) then
    write(Lu,'(1X,A)') '-Optimization for minimum.'
    if (iand(iOptC,1) == 1) then
      write(Lu,'(A)') '  Optimization method: quasi-NR.'
    else if (iand(iOptC,2) == 2) then
      write(Lu,'(A)') '  Optimization method: C1-DIIS.'
    else if (iand(iOptC,4) == 4) then
      write(Lu,'(A)') '  Optimization method: C2-DIIS.'

    else if (iand(iOptC,8) == 8) then
      if (Kriging) then
        write(Lu,'(A)') '  Optimization method: RVO.'
      else
        write(Lu,'(A)') '  Optimization method: RS-RFO.'
      end if
    else
      call WarningMessage(2,' WrInp: Wrong iOptC setting!')
      write(Lu,*) ' iOptC=',iOptC
      call Abend()
    end if
  else
    write(Lu,'(1X,A)') '-Optimization for transition state.'
    if (iand(iOptC,512) == 512) then
      write(Lu,'(A)') '  Optimization method: RS-I-RFO'
      !if (Kriging) then
      !  write(Lu,'(A)') '  Optimization method: RV-I-RFO'
      !else
      !  write(Lu,'(A)') '  Optimization method: RS-I-RFO'
      !end if
    else
      write(Lu,'(A)') '  Optimization method: RS-P-RFO'
      !if (Kriging) then
      !  write(Lu,'(A)') '  Optimization method: RV-P-RFO'
      !else
      !  write(Lu,'(A)') '  Optimization method: RS-P-RFO'
      !end if
    end if
    if (Mode > 0) then
      write(Lu,'(A,I2)') '  Original mode to follow:',Mode
    else
      write(Lu,'(A)') '  No mode to follow is specified!'
      write(Lu,'(A)') '  Optimization will follow mode with the lowest eigenvalue.'
    end if
  end if
  write(Lu,*)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Special options on DIIS

  if ((iand(iOptC,2) == 2) .or. (iand(iOptC,4) == 4)) then
    if (iand(iOptC,16) == 16) then
      write(Lu,'(1X,A)') '-DIIS based on <dx|dx>.'
    else if (iand(iOptC,32) == 32) then
      write(Lu,'(1X,A)') '-DIIS based on <g|dx>.'
    else if (iand(iOptC,64) == 64) then
      write(Lu,'(1X,A)') '-DIIS based on <g|g>.'
    else
      call WarningMessage(2,' WrInp: Wrong iOptC setting!')
      write(Lu,*) ' iOptC=',iOptC
      call Abend()
    end if
    write(Lu,*)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Hessian information

  if (Analytic_Hessian) then
    write(Lu,'(1X,A)') '-The Hessian is analytic.'
    write(Lu,'(1X,A)') ' Hessian from either input or runfile.'
  else
    if ((.not. lOld) .and. (.not. lNmHss)) then
      if (DDV_Schlegel) then
        write(Lu,'(1X,A)') '-Initial Hessian guessed a la Schlegel.'
      else
        if (Kriging) then
          write(Lu,'(1X,A)') '-Hessian guessed by Kriging surrogate surface.'
        else
          write(Lu,'(1X,A)') '-Initial Hessian guessed by Hessian Model Function (HMF).'
          if (iand(iOptC,1024) == 1024) then
            write(Lu,'(A)') '  HMF augmented with weak interactions.'
          end if
        end if
      end if
    else if (lOld .and. (.not. lNmHss)) then
      write(Lu,'(1X,A)') '-Initial Hessian guess was read from a RUNFILE file.'
    else
      write(Lu,'(1X,A,/,A,E9.2)') '-Initial Hessian guess is estimated with finite differences.', &
                                  '    Two point symmetric formula, Delta=',Delta
      if (Cubic) then
        write(Lu,'(1X,A)') '-Cubic force constants evaluated numerically.'
      end if
    end if
  end if
  write(Lu,*)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Hessian update method

  if (.not. Kriging) then
    if (iand(iOptH,1) == 1) then
      write(Lu,'(1X,A)') '-Hessian update method: Fletcher-Meyer'
    else if (iand(iOptH,2) == 2) then
      write(Lu,'(1X,A)') '-Hessian update method: Broyden-Powell'
    else if (iand(iOptH,4) == 4) then
      write(Lu,'(1X,A)') '-Hessian update method: Broyden-Fletcher-Goldfarb-Shanno'
    else if (iand(iOptH,8) == 8) then
      write(Lu,'(1X,A)') '-Hessian update method: none'
    else if (iand(iOptH,16) == 16) then
      write(Lu,'(1X,A)') '-Hessian update method: Murtagh-Sargent-Powell'
    else if (iand(iOptH,64) == 64) then
      write(Lu,'(1X,A)') '-Hessian update method: EU update by Bofill'
    else if (iand(iOptH,128) == 128) then
      write(Lu,'(1X,A)') '-Hessian update method: TS-BFGS update by Bofill'
    else
      call WarningMessage(2,' WrInp: Wrong iOptH setting!')
      write(Lu,*) ' Nonrecognizable iOptH setting:',iOptH
      call Abend()
    end if
    if (.not. (iand(iOptH,8) == 8)) then
      write(Lu,'(A,I3)') '  Maximum number of points in Hessian update:',nWndw
    end if
    if (iand(iOptH,32) == 32) then
      write(Lu,'(A)') '  Hessian update order according to Schlegel'
    end if
    write(Lu,*)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (rHidden >= Two) then
    write(Lu,'(1X,A,/,1X,A,F6.2,A)') '-Improved QM/MM Hessian.',' Hidden atoms until ',rHidden,' bohrs are included.'
    write(Lu,*)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (iRow > 0) then
    if (Redundant) then
      write(Lu,'(1X,A)') '-Relaxation will be done in user supplied redundant internal coordinates.'
    else
      write(Lu,'(1X,A)') '-Relaxation will be done in user supplied non-redundant internal coordinates.'
    end if
  else
    if (CurviLinear) then
      if (HWRS) then
        if (Redundant) then
          write(Lu,'(1X,A)') '-Relaxation will be done on redundant internal coordinates, based on'
          write(Lu,*) ' force constant weighted redundant internal coordinates.'
        else
          write(Lu,'(1X,A)') '-Relaxation will be done on non-redundant internal coordinates, based on'
          write(Lu,*) ' force constant weighted redundant internal coordinates.'
        end if
      else
        if (Redundant) then
          write(Lu,'(1X,A)') '-Relaxation will be done in redundant delocalized internal coordinates.'
        else
          write(Lu,'(1X,A)') '-Relaxation will be done in non-redundant delocalized internal coordinates.'
        end if
      end if
    else
      if (Redundant) then
        write(Lu,'(1X,A)') '-Relaxation will be done in redundant Cartesian coordinates.'
      else
        write(Lu,'(1X,A)') '-Relaxation will be done in approximate non-redundant Cartesian normal mode coordinates.'
      end if
    end if
  end if
  write(Lu,*)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 6) then
  nsAtom = size(Coor,2)
  write(Lu,*)
  write(Lu,'(A)') ' Header from ONEINT:'
  call Banner(Header,2,len(Header(1))+12)
  write(Lu,*)
  call PrList('Symmetry Distinct Nuclear Coordinates / bohr',AtomLbl,nsAtom,Coor,3,nsAtom)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 5) call CollapseOutput(0,'      Slapaf input parameters:')

return

end subroutine WrInp_sl
