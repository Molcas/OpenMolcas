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
! Copyright (C) 1989-1992, Roland Lindh                                *
!               1990, IBM                                              *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine McKinley(ireturn)
!***********************************************************************
!                                                                      *
!  Object: Driver for the one and two electron integral second order   *
!          derivative program McKinley.                                *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA     *
!          July '89 - May '90                                          *
!                                                                      *
!          Roland Lindh, Dept. of Theoretical Chemistry, University of *
!          Lund, SWEDEN. Modified to gradient calculations September   *
!          1991 - February 1992.                                       *
!                                                                      *
!          Anders Bernhardsson, Dept. of Theoretical Chemistry,        *
!          University of  Lund, SWEDEN.                                *
!          Modified to  second order derivatives October '94 -         *
!          '95                                                         *
!***********************************************************************

use McKinley_global, only: CPUStat, lGrd, lHss, Nona, nOneel, nTotal
use Index_Functions, only: nTri_Elem
use Basis_Info, only: dbsc, nCnttp
use Gateway_global, only: Onenly, Test
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
#include "Molcas.fh"
#include "disp.fh"
#include "print.fh"
#include "etwas.fh"
#include "rctfld.fh"
#include "warnings.h"
integer(kind=iwp) :: i, iCnttp, iDummer, iopt, iPrint, irc, iRout, lLine, nDiff, nGrad, nHess, nsAtom
real(kind=wp) :: dum1, dum2, dum3, TCpu1, TCpu2, Time, TWall1, TWall2
character(len=120) :: Lines
logical(kind=iwp) :: DoRys, Run_MCLR
real(kind=wp), allocatable :: GradN(:), Hess(:), Temp(:)
!integer(kind=iwp), parameter :: nLines = 12

!                                                                      *
!***********************************************************************
!                                                                      *
!call McKinley_banner()
call CWTime(TCpu1,TWall1)
iRout = 1
CpuStat(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Print program header
!                                                                      *
!***********************************************************************
!                                                                      *
!Lines(1) = _MOLCAS_VERSION_
!#ifdef _DEMO_
!Lines(2) = 'DEMO VERSION'
!#else
!Lines(2) = ' '
!#endif
!Lines(3) = ' '
!Lines(4) = Vrsn
!Lines(5) = 'A Vectorized Direct Integral Program for derivatives'
!Lines(6) = 'of Cartesian and Spherical Harmonic Gaussians'
!Lines(7) = 'Written by Anders Bernhardsson and Roland Lindh '
!Lines(8) = 'Backtransformation of the 2nd order density matrix from MO to SO by Per-AAke Malmqvist'
!Lines(9) = 'Dept. of Theoretical Chemistry, Chemical Centre, Lund (Sweden)'
!Lines(10) = ' '
!Lines(11) = ' '
!Lines(12) = 'Compiled at '//_BUILD_DATE_
!lLine = Len(Lines(1))
!call Banner(Lines,nLines,lLine)
!                                                                      *
!***********************************************************************
!                                                                      *
! Set error conditions
!
!call XuFlow()
!Call ErrSet(209,1,1,2,1,209)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Check if a numerical procedure will be used
!
call SuperMac()
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the input information as Seward dumped on INFO.
! Set up some info
! Read input

nDiff = 2
DoRys = .true.
call IniSew(DoRys,nDiff)
! pcm_solvent
! check if there is a reaction field
!write(u6,*) 'In mckinley PCM',pcm
call Init_RctFld(.false.,iCharge_ref)
! pcm_solvent end
nsAtom = 0
do iCnttp=1,nCnttp
  nsAtom = nsAtom+dbsc(iCnttp)%nCntr
end do
call Inputh(Run_MCLR)
iPrint = nPrint(iRout)
nGrad = 0
do i=0,nIrrep-1
  nGrad = nGrad+lDisp(i)
end do
call OpnFls_Mckinley()
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate area for hessian etc

nHess = nTri_Elem(nGrad)

call mma_allocate(Hess,nHess,Label='Hess')
Hess(:) = Zero
call mma_allocate(Temp,nHess,Label='Temp')
Temp(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
!    Calculate the second order derivatives of the one electron        *
!    integrals and contract with D.                                    *
!                                                                      *
!***********************************************************************
!                                                                      *
if (lHss) then
  if (iPrint >= 6) then
    write(u6,*)
    write(u6,'(A)') 'The 2nd order derivatives of the one-electron integrals are calculated and contracted with the '// &
                    'one-electron density matrix. '
    write(u6,*)
  end if
  call Timing(dum1,Time,dum2,dum3)
  call Drvh2(Hess,Temp,nHess,show)
  call DrvEtc(nGrad)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!        Compute contribution from the nuclear repulsion.              *
!                                                                      *
!***********************************************************************
!                                                                      *
if (lHss) then
  call DrvN2(Temp,nGrad)
  if (SHOW) call HssPrt(Temp,nHess)
  Hess(:) = Hess+Temp
  if (Show) call HssPrt(Hess,nHess)
end if
if (lGrd) then
  call mma_allocate(GradN,nGrad,Label='GradN')
  call DrvN1_mck(GradN,nGrad)
  iopt = 0
  irc = -1
  call dWrMCK(iRC,iOpt,'NUCGRAD',1,GradN,1)
  if (irc /= 0) call SysAbendMsg('mckinley','Error in writing','Option=NUCGRAD')
  call mma_deallocate(GradN)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!       Calculate the first order derivatives of the one electron      *
!       integrals and store on disc in file MCKINT                     *
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 6) then
  write(u6,*)
  write(u6,'(A)') 'The 1st order derivatives of the one-electron integrals are calculated and stored on disk'
  write(u6,*)
end if
call Drvh1_mck(Nona)
!                                                                      *
!***********************************************************************
!                                                                      *
!      Calculate two electron integrals. First order is contracted     *
!      to Fock matrixes and MO (IJKl) on the fly. Second order         *
!      derivatives are contracted with P.                              *
!      Derivatives are stored in MCKINT.                               *
!                                                                      *
!***********************************************************************
!                                                                      *
nhess = nTri_Elem(ngrad)
call Timing(dum1,Time,dum2,dum3)
CPUStat(nOneel) = CPUStat(nOneel)+Time
if (.not. Onenly) then

  nIsh(:) = 0
  nAsh(:) = 0

  call PrepP()

  iOpt = 0
  iRC = -1
  call WrMck(iRC,iOpt,'NISH',1,nIsh,iDummer)
  if (iRC /= 0) then
    write(u6,*) 'Mckinley: Error writing to MckInt!'
    call Abend()
  end if
  iOpt = 0
  iRC = -1
  call WrMck(iRC,iOpt,'NASH',1,nAsh,iDummer)
  if (iRC /= 0) then
    write(u6,*) 'Mckinley: Error writing to MckInt!'
    call Abend()
  end if

  call Drvg2(Temp,nhess,lGrd,lHss)

  call CloseP()

  if (lHss) then
    call GADSum(Temp,nHess)
    Temp(:) = Half*Temp
    if (Show) call HssPrt(Temp,nHess)

    ! Accumulate contribution to the hessian!

    Hess(:) = Hess+Temp

    if (Show) then
      call Banner('Complete static Hessian',1,23)
      call HssPrt(Hess,nHess)
    end if
    call WrHDsk(Hess,ngrad)
  end if

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Close 'MCKINT' file
iRc = -1
iOpt = 0
call ClsMck(iRC,iOpt)
if (iRc /= 0) then
  write(u6,*) 'McKinley: Error closing MCKINT!'
  call Abend()
end if
call mma_deallocate(Temp)
call mma_deallocate(Hess)

call ClsSew()
!                                                                      *
!***********************************************************************
!                                                                      *
! Epilogue

Lines = 'All data is written to disk, and could be accessed through the MCLR or RASSI program.'
lLine = len(Lines)
call Banner(Lines,1,lLine)

call CWTime(TCpu2,TWall2)

call Timing(Time,dum1,dum2,dum3)
CPUStat(nTotal) = Time
if (iPrint >= 6) call Sttstc()
if (Test) then
  ireturn = _RC_INPUT_ERROR_
else
  call Request_MCLR_Run(Run_MCLR,ireturn,iPrint)
end if

return

end subroutine McKinley
