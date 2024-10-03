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
! Copyright (C) 2016,2017, Roland Lindh                                *
!***********************************************************************

subroutine PrFin0(Dens,Dens_ab,nDT,EOrb,nEO,CMO,nCMO,KntE)
!***********************************************************************
!                                                                      *
!     purpose: Final printout                                          *
!                                                                      *
!     input:                                                           *
!       Dens    : the last total density                               *
!       EOrb    : orbital energies of length nEO                       *
!       CMO     : molecular orbital coefficients of length nCMO        *
!                                                                      *
!***********************************************************************

#ifdef _HDF5_
use mh5, only: mh5_put_dset
use SCFWfn, only: wfn_energy
#endif
#ifdef _FDE_
use Embedding_Global, only: Eemb, embPot
#endif
use KSDFT_Info, only: CoefR, CoefX
use OFembed, only: Do_OFemb
use SpinAV, only: Do_SpinAV
use InfSCF, only: Addc_KSDFT, DE_KSDFT_c, DMOMax, Do_Addc, Do_Tw, E1V, E2V, E_nondyn, EKin, EneV, Erest_xc, FMOMax, iPrint, &
                  jPrint, KSDFT, lPaper, MxConstr, nBas, nBT, nD, nIter, nIterP, nOcc, nOrb, nSym, PotNuc, s2CNO, s2UHF, WarnCfg, &
                  WarnPocc, WarnSlow
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDT, nEO, nCMO
real(kind=wp), intent(inout) :: Dens(nDT)
real(kind=wp), intent(in) :: Dens_ab(nDT), EOrb(nEO), CMO(nCMO), KntE(nDT)
integer(kind=iwp) :: iDumm, iMult, iPL, iSpn, iTol
real(kind=wp) :: Dumm1(1), E_Tw, ECNO, Ecorr, sUHF, Virial
character(len=80) :: Lines(6)
character(len=60) :: Frmt
integer(kind=iwp), external :: Cho_X_GetTol, iPrintLevel
real(kind=wp), external :: DDot_
logical(kind=iwp), external :: Reduce_Prt

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

jPrint = iPrint
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0
if (iPL <= 1) jPrint = 1

!----------------------------------------------------------------------*

! Calculate kinetic energy
if (nD == 2) Dens(:) = Dens(:)+Dens_ab(:)
EKin = DDot_(nBT,KntE,1,Dens,1)

! Print out header to final results
Lines(1) = 'SCF/KS-DFT Program, Final results'
Lines(2) = ' '
Lines(3) = ' '
!if (.false.) Lines(3) = Molcas_revision
Lines(4) = ' '
Lines(5) = 'Final Results'
if (jPrint >= 2) then
  call Banner(Lines,5,lPaper-7)
  write(6,*)
end if

! Compute energies from aufbau density matrix
if (abs(EKin) > 1.0D-6) then
  Virial = -EneV/EKin
else
  Virial = Zero
end if

! Compute estimate of FNO correlation energy from Delta_Tw
if (Do_Tw) call Tw_corr_drv(EOrb,nEO,CMO,nCMO,Ecorr)

! Print out final results
if (WarnCfg) &
  call WarningMessage(1,'Warning:; The program may have converged to a solution;that does not correspond to the lowest energy!')
if (WarnPocc) call WarningMessage(1,'Warning:; The program may have converged to a solution;with partial occupation numbers!')
if (WarnSlow) call WarningMessage(1,'Warning:; The program had convergence problems;and terminated with looser convergence')
Frmt = '(6X,A,T50,F19.10)'
suhf = -Half+sqrt(0.25d0+s2uhf)
call put_dscalar('UHFSPIN',SUHF)
iTol = min(Cho_X_GetTol(8),8)
if (jPrint >= 2) then
  if (MxConstr > 0) then
    DE_KSDFT_c = Zero
    if (Do_Addc) then
      call SetUp_iSD()
      call Get_DEcorr(nBT,Dumm1,iDumm,'SCF ')
      call Free_iSD()
    end if
    ECNO = EneV+E_nondyn+DE_KSDFT_c
    if (KSDFT /= 'SCF') ECNO = ECNO+Erest_xc
    call PrintResult(6,Frmt,'Total energy',0,' ',[ECNO],1)
    call PrintResult(6,Frmt,'Nondynamical correlation energy',0,' ',[E_nondyn],1)
    if (KSDFT /= 'SCF') call PrintResult(6,Frmt,'Energy-restoring term',0,' ',[Erest_xc],1)
    if (Do_Addc) call PrintResult(6,Frmt,'Added correlation energy ('//ADDC_KSDFT(1:4)//') ',0,' ',[DE_KSDFT_c],1)
    call Add_Info('E_CNO',[ECNO],1,iTol)
  end if
  if (Do_Tw) then
    E_Tw = EneV+Ecorr
    call PrintResult(6,Frmt,'Total energy',0,' ',[E_Tw],1)
    call PrintResult(6,Frmt,'Delta_Tw correlation energy',0,' ',[Ecorr],1)
    call Add_Info('E_Tw',[E_Tw],1,iTol)
  end if
  if (KSDFT == 'SCF') then
    call PrintResult(6,Frmt,'Total SCF energy',0,' ',[EneV],1)
    !write(6,Frmt) 'Total SCF energy',EneV
  else
    call PrintResult(6,Frmt,'Total KS-DFT energy',0,' ',[EneV],1)
    !write(6,Frmt) 'Total KS-DFT energy',EneV
  end if
  write(6,Frmt) 'One-electron energy',E1V
# ifdef _FDE_
  ! Embedding
  if (embPot) write(6,Frmt) 'E from embedding potential(<Psi|v_emb|Psi>)',Eemb
# endif
  write(6,Frmt) 'Two-electron energy',E2V
  write(6,Frmt) 'Nuclear repulsion energy',PotNuc
  write(6,Frmt) 'Kinetic energy (interpolated)',EKin
  write(6,Frmt) 'Virial theorem',Virial
  if (.not. Do_SpinAV) then
    write(6,Frmt) 'Total spin, S(S+1)',s2uhf
    write(6,Frmt) 'Total spin, S',suhf
  end if
  if (MxConstr > 0) write(6,Frmt) 'Spin deviation',s2uhf-s2CNO
end if
iSpn = int(suhf+Half)
iMult = 2*iSpn+1
call Put_iScalar('Multiplicity',iMult)
call Add_Info('E_SCF',[EneV],1,iTol)
#ifdef _HDF5_
call mh5_put_dset(wfn_energy,EneV)
#endif

if ((nIter(nIterP) > 0) .and. (jPrint >= 2)) then
  write(6,Frmt) 'Max non-diagonal density matrix element',DMOMax
  write(6,Frmt) 'Max non-diagonal Fock matrix element',FMOMax
end if
if ((CoefX /= 1.0) .or. (CoefR /= 1.0)) then
  write(6,Frmt) 'Exchange scaling factor',CoefX
  write(6,Frmt) 'Correlation scaling factor',CoefR
end if
if (jPrint >= 2) write(6,*)
!if ((jPrint >= 2) .and. Do_OFemb) Call OFE_print(EneV)
if (Do_OFemb) call OFE_print(EneV)

! xml tagging

if (KSDFT == 'SCF') then
  call xml_dDump('energy','Total SCF energy','a.u.',1,[EneV],1,1)
else
  call xml_dDump('energy','Total KS-DFT energy','a.u.',1,[EneV],1,1)
end if
call xml_dDump('kinetic','Kinetic energy','a.u.',2,[Ekin],1,1)
call xml_dDump('virial','Virial coefficient','a.u.',2,[Virial],1,1)
call xml_dDump('spin','UHF spin','',1,[suhf],1,1)
call xml_dDump('potnuc','Nuclear repulsion energy','a.u.',1,[potnuc],1,1)
call xml_dDump('energy1el','One electron energy','a.u.',1,[E1V],1,1)
call xml_dDump('energy2el','Two electron energy','a.u.',1,[E2V],1,1)
call xml_iDump('nsym','Number of irreps','',1,[nSym],1,1)
call xml_iDump('nbas','Number of basis functions','',1,nBas,nSym,1)
call xml_iDump('norb','Number of orbitals','',1,nOrb,nSym,1)
if (nD == 1) then
  call xml_iDump('nocc','Number of occupied orbitals','',1,nOcc(:,1),nSym,1)
else
  call xml_iDump('nocc_a','Number of occupied alpha orbitals','',1,nOcc(:,1),nSym,1)
  call xml_iDump('nocc_b','Number of occupied beta orbitals','',1,nOcc(:,2),nSym,1)
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine PrFin0
