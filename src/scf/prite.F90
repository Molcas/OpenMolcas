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
!***********************************************************************

subroutine PrIte(QNR,CMO,mBB,nD,Ovrlp,mBT,OccNo,mmB)
!***********************************************************************
!                                                                      *
!     purpose: Print out informations in every iteration               *
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!***********************************************************************

use InfSCF, only: AccCon, CPUItr, DltNrm, DltNTh, DMOMax, DNorm, DThr, E1V, E2V, EDiff, EneV, EThr, FMOMax, FThr, Iter, IterPrLv, &
                  jPrint, TNorm
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: QNR
integer(kind=iwp), intent(in) :: mBB, nD, mBT, mmB
real(kind=wp), intent(in) :: CMO(mBB,nD), Ovrlp(mBT), OccNo(mmB,nD)
real(kind=wp) :: Shift = Zero
logical(kind=iwp) :: Set_Shift = .false.
character :: cDltNrm, cDMOMax, cEDiff, cFMOMax

if (iterprlv > 0) then
  write(u6,*)
  write(u6,'(a)') '*******************'
  write(u6,'(a,i3,a)') '** Iteration ',iter,' **'
  write(u6,'(a)') '*******************'
  write(u6,*)
  write(u6,'(a,f10.2)') 'Cpu time [sec]        ',CpuItr

  select case (AccCon)
    case ('None','NoneDa')
      write(u6,'(a)') 'No convergence acceleration'
    case ('EDIIS','ADIIS')
      write(u6,'(a)') 'Convergence is accelerated by damping'
    case ('QNRc1D')
      write(u6,'(2a)') 'Convergence is accelerated by QNR with ','c1-DIIS'
    case ('QNRc2D')
      write(u6,'(2a)') 'Convergence is accelerated by QNR with ','c2-DIIS'
    case default
      write(u6,'(2a)') 'Convergence accelerations is ',AccCon
  end select

  write(u6,*)
  write(u6,'(a,f16.8)') 'Total energy          ',EneV
  write(u6,'(a,f16.8)') 'One electron energy   ',E1V
  write(u6,'(a,f16.8)') 'Two electron energy   ',E2V

  if ((abs(Ediff) > Ethr) .or. (iter <= 1)) then
    write(u6,'(a,f16.8)') 'Energy difference     ',Ediff
  else
    write(u6,'(a,f16.8,a)') 'Energy difference     ',Ediff,' is converged'
  end if

  if (QNR) then
    if (DltNrm > DltNth) then
      write(u6,'(a,f16.8)') 'Delta norm            ',DltNrm
    else
      write(u6,'(a,f16.8,a)') 'Delta norm            ',DltNrm,' is converged'
    end if
  else
    if (abs(DMOMax) > Dthr) then
      write(u6,'(a,f16.8)') 'Max offdiagonal Dij   ',DMOmax
    else
      write(u6,'(a,f16.8,a)') 'Max offdiagonal Dij   ',DMOmax,' is converged'
    end if
  end if

  if (abs(FMOMax) > Fthr) then
    write(u6,'(a,f16.8)') 'Max offdiagonal Fij   ',FMOmax
  else
    write(u6,'(a,f16.8,a)') 'Max offdiagonal Fij   ',FMOmax,' is converged'
  end if

  write(u6,'(a,f16.8)') 'D-norm                ',sqrt(Dnorm)
  write(u6,'(a,f16.8)') 'T-norm                ',sqrt(Tnorm)
  if (iterprlv >= 2) call MulPop(CMO,mBB,nD,Ovrlp,mBT,OccNo,mmB)

else if (jPrint >= 2) then

  if (.not. Set_Shift) then
    if (abs(EneV) > 1.0e3_wp) then
      Shift = abs(EneV)*1.0e-3_wp
      Shift = real(int(Shift),kind=wp)*1.0e3_wp
      write(u6,*)
      write(u6,'(1X,A,f10.0,A)') 'The total and one-electron energies are shifted by a value of ',Shift,' a.u.'
      write(u6,*)
    end if
    Set_Shift = .true.
  end if
  cEDiff = ' '
  if ((abs(Ediff) > Ethr) .or. (EDiff > Zero)) cEDiff = '*'
  cFMOMax = ' '
  if (abs(FMOMax) > Fthr) cFMOMax = '*'

  if (QNR) then
    cDltNrm = ' '
    if (DltNrm > DltNth) cDltNrm = '*'
    write(u6,'(1X,i3,3f16.9,1x,3(es10.2,a1,1x),2es11.2,3x,A,f6.0)') &
      Iter,EneV+Shift,E1V+Shift,E2V,EDiff,cEDiff,DltNrm,cDltNrm,FMOMax,cFMOMax,sqrt(DNorm),sqrt(TNorm),AccCon,CpuItr
  else
    cDMOMax = ' '
    if (abs(DMOMax) > Dthr) cDMOMax = '*'
    write(u6,'(1X,i3,3f16.9,1x,3(es10.2,a1,1x),2es11.2,3x,A,f6.0)') &
      Iter,EneV+Shift,E1V+Shift,E2V,EDiff,cEDiff,DMOMax,cDMOMax,FMOMax,cFMOMax,sqrt(DNorm),sqrt(TNorm),AccCon,CpuItr

  end if
end if
call XFlush(u6)

end subroutine prite
