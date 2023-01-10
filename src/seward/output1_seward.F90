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
! Copyright (C) 2006, Roland Lindh                                     *
!***********************************************************************

subroutine Output1_Seward(lOPTO)
!***********************************************************************
!                                                                      *
!     Object: to write the output of seward                            *
!                                                                      *
!     Author: Roland Lindh, Dept Chem. Phys., Lund University, Sweden  *
!             September '06                                            *
!***********************************************************************

use Basis_Info, only: dbsc, Gaussian_Type, MolWgh, nCnttp, Nuclear_Model, Point_Charge
use Center_Info, only: dc
use Period, only: Cell_l, ispread
use GeoList, only: Centr, Mass
use MpmC, only: Coor_MPM
use EFP_Module, only: lEFP
#ifdef _EFP_
use EFP_Module, only: nEFP_Fragments
#endif
use External_centers, only: AMP_Center, DMS_Centers, nDMS, nEF, nOrdEF, nWel, nXF, OAM_Center, OMQ_Center, XF
use DKH_Info, only: BSS, DKroll, iCtrLD, iRELAE, iRELMP, LDKroll, nCtrLD, radiLD
use Sizes_of_Seward, only: S
use Gateway_Info, only: CutInt, DoFMM, EMFR, FNMC, GIAO, kVector, lAMFI, lMXTC, lRel, RPQMin, ThrInt, Vlct
use RICD_Info, only: iRI_Type, LDF, Do_RI, Cholesky, Do_acCD_Basis, Skip_High_AC, Cho_OneCenter, LocalDF, Thrshld_CD
use Symmetry_Info, only: nIrrep
use Gateway_global, only: GS_Mode, Onenly, Run_Mode, Prprt, Test
use Constants, only: Zero, One, Two, Ten, Pi, Angstrom
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: lOPTO
#include "Molcas.fh"
#include "rmat.fh"
#include "rctfld.fh"
#include "print.fh"
#include "localdf.fh"
integer(kind=iwp) :: i, iCnttp, iDKH_H_Order, iDKH_X_Order, iParam, iPrint, iRout, iTtl, LuWr, nTtl
real(kind=wp) :: temp
logical(kind=iwp) :: l_aCD_Thr, Found, lNoPair, lPam2, lECP, lPP
character(len=80) :: Title(10)

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)
LuWr = u6
!                                                                      *
!***********************************************************************
!                                                                      *
lNoPair = .false.
lPam2 = .false.
lECP = .false.
lPP = .false.
do i=1,nCnttp
  lNoPair = lNoPair .or. dbsc(i)%NoPair
  lPam2 = lPam2 .or. dbsc(i)%lPam2
  lECP = lECP .or. dbsc(i)%ECP
  lPP = lPP .or. (dbsc(i)%nPP /= 0)
  lECP = lECP .or. dbsc(i)%ECP
  lPP = lECP .or. dbsc(i)%ECP
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Start of output

if (Test) then
  write(LuWr,*)
  write(LuWr,'(15X,88A)') ('*',i=1,45)
  write(LuWr,'(15X,A)') '* TEST: SEWARD will only process the input! *'
  write(LuWr,'(15X,88A)') ('*',i=1,45)
else

  iDKH_X_Order = iRELAE/10000
  iParam = iRELAE-1000-iDKH_X_Order*10000
  iDKH_H_Order = iParam/10
  iParam = iParam-iDKH_H_Order*10

  write(LuWr,'(15X,A)') 'SEWARD will generate:'
  write(LuWr,'(15X,A,I2)') '   Multipole Moment integrals up to order ',S%nMltpl
  if (.not. Prprt) then
    write(LuWr,'(15X,A)') '   Kinetic Energy integrals'
    if (Nuclear_Model == Gaussian_Type) then
      write(LuWr,'(15X,A)') '   Nuclear Attraction integrals (finite nuclei - Gaussian type)'
    else if (Nuclear_Model == Point_Charge) then
      write(LuWr,'(15X,A)') '   Nuclear Attraction integrals (point charge)'
    else
      write(LuWr,'(15X,A)') '   Nuclear Attraction integrals (finite nuclei -  Modified Gaussian type)'
    end if
    if (lECP) then
      if (lNoPair) then
        select case (IRELMP)
          case (0)
            write(LuWr,'(15X,A)') '   One-Electron Hamiltonian integrals modified with ECP and No-Pair contributions'
          case (1)
            write(LuWr,'(15X,A)') '   One-Electron Hamiltonian integrals modified with ECP and No-Pair (DK1) contributions'
          case (2)
            write(LuWr,'(15X,A)') '   One-Electron Hamiltonian integrals modified with ECP and No-Pair (DK2) contributions'
          case (3)
            write(LuWr,'(15X,A)') '   One-Electron Hamiltonian integrals modified with ECP and No-Pair (DK3) contributions'
          case (11)
            write(LuWr,'(15X,A)') '   One-Electron Hamiltonian integrals modified with ECP and RESC contributions'
          case (21)
            write(LuWr,'(15X,A)') '   One-Electron Hamiltonian integrals modified with ECP and ZORA contributions'
          case (22)
            write(LuWr,'(15X,A)') '   One-Electron Hamiltonian integrals modified with ECP and ZORA-FP contributions'
          case (23)
            write(LuWr,'(15X,A)') '   One-Electron Hamiltonian integrals modified with ECP and IORA contributions'
        end select
      else
        write(LuWr,'(15X,A)') '   One-Electron Hamiltonian integrals modified with ECP contributions'
      end if
    else
      write(LuWr,'(15X,A)') '   One-Electron Hamiltonian integrals'
    end if
    if (FNMC) then
      write(LuWr,'(15X,A)') '   Finite nuclear mass correction added'
    end if
    if (lPAM2) then
      write(LuWr,'(15X,A)') '   Include potentials for DMFT calculation'
    end if
    if (lRel) then
      write(LuWr,'(15X,A)') '   Mass-Velocity integrals'
      write(LuWr,'(15X,A)') '   Darwin One-Electron Contact Term integrals'
    end if
    if (Vlct) write(LuWr,'(15X,A)') '   Velocity integrals'
    if (DKroll) then
      if (iRELAE < 1000) then
        select case (IRELAE)
          case (0)
            write(LuWr,'(15X,A)') '   Relativistic Douglas-Kroll integrals'
          case (1)
            write(LuWr,'(15X,A)') '   Relativistic Douglas-Kroll (DK1) integrals'
          case (2)
            write(LuWr,'(15X,A)') '   Relativistic Douglas-Kroll (DK2) integrals'
          case (3)
            write(LuWr,'(15X,A)') '   Relativistic Douglas-Kroll (DK3) integrals'
          case (4)
            write(LuWr,'(15X,A)') '   full Relativistic Douglas-Kroll (DK3) integrals'
          case (11)
            write(LuWr,'(15X,A)') '   Relativistic RESC integrals'
          case (21)
            write(LuWr,'(15X,A)') '   Relativistic ZORA integrals'
          case (22)
            write(LuWr,'(15X,A)') '   Relativistic ZORA-FP integrals'
          case (23)
            write(LuWr,'(15X,A)') '   Relativistic IORA integrals'
          case (101)
            write(LuWr,'(15X,A)') '   Relativistic X2C integrals'
          case (102)
            write(LuWr,'(15X,A)') '   Relativistic BSS integrals'
          case default
            if (BSS) write(LuWr,'(15X,A)') '   Relativistic Barysz-Sadlej-Snijders integrals'
        end select
      else
        if (LDKroll) then
          write(LuWr,'(17X,A)') ' Relativistic Local-Douglas-Kroll-Hess integrals:'
          if (nCtrLD == 0) then
            if (radiLD == Zero) then
              write(LuWr,'(17X,A)') '   - Atomic approximation'
            else
              write(LuWr,'(17X,A)') '   - Full local approximation'
            end if
          else
            write(LuWr,'(17X,A)') '   - Partial local approximation:'
            write(LuWr,'(17X,A,10(A4))') '     - Centers: ',(dc(iCtrLD(i))%LblCnt,i=1,nCtrLD)
            write(LuWr,'(17X,A,F6.2,A)') '     - Cutoff radius: ',radiLD,' Bohr'
          end if
        else
          write(LuWr,'(17X,A)') ' Relativistic Douglas-Kroll-Hess integrals:'
        end if
        if (iParam == 1) then
          write(LuWr,'(17X,A)') '   - Parametrization         : OPT'
        else if (iParam == 2) then
          write(LuWr,'(17X,A)') '   - Parametrization         : EXP'
        else if (iParam == 3) then
          write(LuWr,'(17X,A)') '   - Parametrization         : SQR'
        else if (iParam == 4) then
          write(LuWr,'(17X,A)') '   - Parametrization         : MCW'
        else if (iParam == 5) then
          write(LuWr,'(17X,A)') '   - Parametrization         : CAY'
        end if
        write(LuWr,'(17X,A,I2)') '   - DKH order of Hamiltonian:',iDKH_H_order
        write(LuWr,'(17X,A,I2)') '   - DKH order of Properties :',iDKH_X_order
        write(LuWr,'(17X,A)') '        - multipole moment operators'
        write(LuWr,'(17X,A)') '        - electric potential operators'
        write(LuWr,'(17X,A)') '        - contact operators'
      end if
    end if
    if (lRF) then
      if (PCM) then
        write(LuWr,'(15X,A)') '   Reaction Field integrals (PCM)'
      else if (lLangevin) then
        write(LuWr,'(15X,A)') '   Reaction Field integrals (Langevin)'
      else
        write(LuWr,'(15X,A)') '   Reaction Field integrals (KirkWood-Onsager)'
      end if
    end if
    if (lPP) write(LuWr,'(15X,A)') '   Pseudo Potential integrals'
  end if
  if (allocated(XF)) write(LuWr,'(15X,A,I6,A)') '   External field from',nXF,' point(s) added to the one-electron Hamiltonian'
  if ((nEF > 0) .and. (nOrdEF >= 0)) write(LuWr,'(15X,A,I6,A)') '   Electric potential for',nEF,' points'
  if ((nEF > 0) .and. (nOrdEF >= 1)) write(LuWr,'(15X,A,I6,A)') '   Electric field integrals for',nEF,' points'
  if ((nEF > 0) .and. (nOrdEF >= 2)) write(LuWr,'(15X,A,I6,A)') '   Electric field gradient integrals for',nEF,' points'
  if ((nEF > 0) .and. (nOrdEF >= 2)) write(LuWr,'(15X,A,I6,A)') '   Contact term integrals for',nEF,' points'
  if (allocated(DMS_Centers)) write(LuWr,'(15X,A,I6,A)') '   Diamagnetic shielding integrals for',nDMS,' points'
  if (allocated(OAM_Center)) write(LuWr,'(15X,A,3(F7.4,1X),A)') '   Orbital angular momentum around (',(OAM_Center(i),i=1,3),')'
  if (allocated(OMQ_Center)) write(LuWr,'(15X,A,3(F7.4,1X),A)') '   Orbital magnetic quadrupole around (',(OMQ_Center(i),i=1,3),')'
  if (Vlct .and. (S%nMltpl >= 2)) write(LuWr,'(15X,A,3(F7.4,1X),A)') '   Velocity quadrupole around (',(Coor_MPM(i,3),i=1,3),')'
  if (allocated(AMP_Center)) write(LuWr,'(15X,A,3(F7.4,1X),A)') &
                             '   Products of Orbital angular momentum operators around (',(AMP_Center(i),i=1,3),')'
  if (nWel /= 0) write(LuWr,'(15X,A,I4,A)') '   Spherical well for',nWel,' exponent(s) added to the one-electron Hamiltonian'
  if (lAMFI) write(LuWr,'(15X,A)') '   Atomic mean-field integrals'
  if (lMXTC) write(LuWr,'(15X,A)') '   Hyperfine Magnetic integrals(MAG) calculated from Gen1Int F90 library'
  if (DoFMM) then
    write(LuWr,'(15X,A)') '   Integral environment set up for FMM option'
    write(LuWr,'(15X,A,F10.5)') '    - RPQMin: ',RPQMin
  end if
  if (lEFP) then
#   ifdef _EFP_
    write(LuWr,'(15X,A)') '   Effective Fragment potentials added       '
    write(LuWr,'(15X,A,I4)') '    - # of fragments: ',nEFP_fragments
#   else
    write(LuWr,'(15X,A)') '   EFP input specified but code not enabled for the option'
    call Abend()
#   endif
  end if
  if (.not. Onenly) then
    if (Cholesky) then
      write(LuWr,'(15X,A)') '   Cholesky decomposed two-electron repulsion integrals'
      if (Cho_OneCenter) then
        write(LuWr,'(17X,A,G10.2)') '  - 1C-CD Threshold: ',Thrshld_CD
      else
        write(LuWr,'(17X,A,G10.2)') '  - CD Threshold: ',Thrshld_CD
      end if
    else if (Do_RI) then
      if (LocalDF) then
        if (LDF_Constraint == -1) then
          write(LuWr,'(15X,A)') '   Local Density Fitting coefficients'
        else
          write(LuWr,'(15X,A)') '   Constrained Local Density Fitting coefficients'
          if (LDF_Constraint == 0) then
            write(LuWr,'(17X,A)') '  - constraint type: charge'
          else
            call WarningMessage(2,'Unknown constraint!')
            write(LuWr,'(A,I10)') 'LDF_Constraint=',LDF_Constraint
            call LDF_Quit(-1)
          end if
        end if
        if (LDF2) then
          write(LuWr,'(17X,A,G10.2)') '  - two-center auxiliary functions included (when needed); target accuracy: ',Thr_Accuracy
        else
          write(LuWr,'(17X,A)') '  - two-center auxiliary functions not included'
        end if
      else if (LDF) then
        write(LuWr,'(15X,A)') '   LDF decomposed two-electron repulsion integrals stored Cholesky style'
        write(LuWr,'(15X,A)') '    Concept demonstration only!'
      else
        write(LuWr,'(15X,A)') '   RI decomposed two-electron repulsion integrals stored Cholesky style'
      end if
      if (iRI_Type == 1) then
        write(LuWr,'(17X,A)') '  - RIJ auxiliary basis'
      else if (iRI_Type == 2) then
        write(LuWr,'(17X,A)') '  - RIJK auxiliary basis'
      else if (iRI_Type == 3) then
        write(LuWr,'(17X,A)') '  - RIC auxiliary basis'
      else if (iRI_Type == 5) then
        write(LuWr,'(17X,A)') '  - External RICD auxiliary basis'
      else
        if (Do_acCD_Basis) then
          write(LuWr,'(17X,A)') '  - acCD auxiliary basis'
        else
          write(LuWr,'(17X,A)') '  - aCD auxiliary basis'
        end if
        write(LuWr,'(17X,A,G10.2)') '  - CD Threshold: ',Thrshld_CD
        l_aCD_Thr = .false.
        do iCnttp=1,nCnttp
          l_aCD_Thr = l_aCD_Thr .or. (dbsc(iCnttp)%aCD_Thr /= One)
        end do
        if (l_aCD_Thr) then
          write(LuWr,'(17X,A)') '     Note that the threshold for individual basis sets might be modified!'
        end if
        if (Skip_High_AC) then
          write(LuWr,'(17X,A)') '  - Skip high angular momentum combinations'
        end if
      end if
    else
      write(LuWr,'(15X,A)') '   Two-Electron Repulsion integrals'
    end if
  end if
  if (RMat_On) then
    write(LuWr,*)
    write(LuWr,'(15X,A)') '   OBSERVE that some integrals are modified to enable variational R-matrix calculations!'
  end if
  if (GIAO) then
    write(LuWr,'(15X,A)') '   GIAO integrals differentiated with respect to B'
    write(LuWr,'(15X,A)') '     dS/dB                                        '
    write(LuWr,'(15X,A)') '     dT/dB                                        '
    write(LuWr,'(15X,A)') '     dV/dB                                        '
  end if

  ! Transition moment integrals for oscillator strengths of
  ! electronic transitions.

  if (EMFR) then
    write(LuWr,'(15X,A)') '   Transition moment integrals'
    write(LuWr,'(15X,A,3(F7.4,1X),A)') '   The wavevector k: (',(kVector(i),i=1,3),')'
    temp = sqrt(KVector(1)**2+KVector(2)**2+kVector(3)**2)
    temp = (Two*Pi)/temp
    write(LuWr,'(15X,A,(F10.4,1X),A)') '   Wavelength:        ',Temp,'a.u.'
    write(LuWr,'(15X,A,(F10.4,1X),A)') '                      ',Temp*Angstrom,'Angstrom'
    write(LuWr,'(15X,A,(F10.4,1X),A)') '                      ',Temp*Angstrom/Ten,'nm'
  end if

end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Qpg_cArray('SewardXTitle',Found,nTtl)
if (Found) then
  nTtl = nTtl/80
  call Get_cArray('SewardXTitle',Title(1),nTtl*80)
  if (iPrint >= 6) then
    write(LuWr,*)
    write(LuWr,'(15X,88A)') ('*',i=1,88)
    write(LuWr,'(15X,88A)') '*',(' ',i=1,86),'*'
    do iTtl=1,nTtl
      write(LuWr,'(15X,A,A,A)') '*   ',Title(iTtl),'   *'
    end do
    write(LuWr,'(15X,88A)') '*',(' ',i=1,86),'*'
    write(LuWr,'(15X,88A)') ('*',i=1,88)
  else
    write(LuWr,*)
    write(LuWr,'(A)') ' Title:'
    do iTtl=1,nTtl
      write(LuWr,'(8X,A)') Title(iTtl)
    end do
    write(LuWr,*)
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
write(LuWr,*)
write(LuWr,'(19X,A,E9.2)') 'Integrals are discarded if absolute value <:',ThrInt
write(LuWr,'(19X,A,E9.2)') 'Integral cutoff threshold is set to       <:',CutInt
!                                                                      *
!***********************************************************************
!                                                                      *
if (Run_Mode == GS_Mode) call Print_Symmetry()
!                                                                      *
!***********************************************************************
!                                                                      *
if (nIrrep > 1) then
  if (MolWgh == 0) then
    write(LuWr,*)
    write(LuWr,'(19X,A)') ' Symmetry adaptation a la DCR.'
    write(LuWr,*)
  end if
  if (MolWgh == 1) then
    write(LuWr,*)
    write(LuWr,'(19X,A)') ' Symmetry adaptation a la MOLECULE.'
    write(LuWr,*)
  end if
  if (MolWgh == 2) then
    write(LuWr,*)
    write(LuWr,'(19X,A)') ' Unitary symmetry adaptation'
    write(LuWr,*)
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Print Cell unit information

if (Cell_l) then
  write(LuWr,'(6X,30(''-''))')
  write(LuWr,'(6X,A)') '* - the centers of the Unit Cell'
  write(LuWr,'(A)') ' '
  write(LuWr,'(A,3I3)') 'Spread of the unit cell:',(ispread(i),i=1,3)
  write(LuWr,'(A)') ' '
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out basis set information

if (Run_Mode == GS_Mode) then
  call Print_Basis(lOPTO)
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out coordinates, bond, angles and torsional angles

  if (lOPTO) then
    call Print_Geometry(1)
  else
    call Print_Geometry(0)
  end if
  call Print_Isotopes()
!                                                                      *
!***********************************************************************
!                                                                      *
! Rigid Rotor analysis etc.

  call RigRot(Centr,Mass,S%kCentr)
!                                                                      *
!***********************************************************************
!                                                                      *
  call Print_Basis2()
!                                                                      *
!***********************************************************************
!                                                                      *
  call Print_OpInfo()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Output1_Seward
