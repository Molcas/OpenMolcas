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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               1995, Martin Schuetz                                   *
!***********************************************************************

subroutine ReadIn_SCF(SIntTh)
!***********************************************************************
!                                                                      *
!     purpose: Read input to SCF: one-electron integrals, informations *
!              about basis set and molecule and options to SCF.        *
!                                                                      *
!***********************************************************************

use Gateway_Info, only: PkAcc
use InfSCF, only: DSCF, EThr, KSDFT, nCore, nDisc, TimFld
use Definitions, only: wp

implicit none
real(kind=wp), intent(inout) :: SIntTh
real(kind=wp) :: CPU1, CPU2, Tim1, Tim2, Tim3

call Timing(Cpu1,Tim1,Tim2,Tim3)
!                                                                      *
!***********************************************************************
!                                                                      *
! read one electron integral file header

call R1IBas()
!                                                                      *
!***********************************************************************
!                                                                      *
! read input

call RdInp_SCF()
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate memory for SCF procedure

call MemAlo()
!                                                                      *
!***********************************************************************
!                                                                      *
! read one-electron integrals

call R1IntA()
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize seward

call IniSew_scf(DSCF,EThr,SIntTh,KSDFT)
!                                                                      *
!***********************************************************************
!                                                                      *
! setup for direct or conventional integral calculations

if (DSCF) then
  call Set_Basis_Mode('Valence')
  call Setup_iSD()
  call AlloK2()
  call Free_iSD()

  ! Initiate integral packing for semi-direct implementation

  if (nDisc /= 0) call inipkr8(PkAcc,.true.)

  ! Allocate buffers for semi-direct SCF

  call IniBuf(nDisc,nCore)

else

  ! Read the header of the two-electron integral file

  call Rd2Int_SCF()

end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(1) = TimFld(1)+(Cpu2-Cpu1)

end subroutine ReadIn_SCF
