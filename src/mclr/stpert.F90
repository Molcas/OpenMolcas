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

subroutine StPert()

use MckDat, only: sNew
use ipPage, only: ipin, W
use MCLR_Data, only: Hss, FAMO_SpinP, FAMO_SpinM, G2mm, G2mp, G2pp, Fp, Fm, G1p, G1m
use MCLR_Data, only: ipCI, nDens2
use MCLR_Data, only: RMS, rBetaA, rBetaS
use MCLR_Data, only: lDisp, SwLbl
use MCLR_Data, only: MS2
use MCLR_Data, only: FnMck, LuMck
use input_mclr, only: nSym, McKinley, PT2, nDisp, SpinPol, nAsh, nBas, nIsh, nTPert, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero, Half
use Definitions, only: wp, u6

implicit none
character(len=16) Label
character(len=8) MckLbl
character(len=288) Header
integer idum(1)
real*8, allocatable :: Tmp1(:), Tmp2(:)
integer nHss, iS, iRC, iOpt, nAct, iSym, nG, nG2, iType, nMax, iDummer
real*8 rAlphas

nHss = 0
do iS=1,nSym
  nHss = nHss+lDisp(is)*(lDisp(is)+1)/2
end do
call mma_allocate(Hss,nHss,Label='Hss')
Hss(:) = Zero

if (.not. Mckinley) then
  irc = -1
  iopt = ibset(0,sNew)
  call OPNMCK(irc,iopt,FNMCK,LUMCK)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error opening MCKINT'
    call Abend()
  end if
  irc = -1
  iopt = 0
  LABEL = 'SEWARD'
  if (PT2) LABEL = 'PT2LAG'
  MckLbl = 'PERT'
  call cWrMck(iRC,iOpt,MckLbl,1,LABEL,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'NDISP'
  idum(1) = ndisp
  call WrMck(iRC,iOpt,MckLbl,1,idum,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'TDISP'
  call WrMck(iRC,iOpt,MckLbl,1,ntpert,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'Title'
  call cWrMck(iRC,iOpt,MckLbl,1,Header,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'nSym'
  idum(1) = nSym
  call WrMck(iRC,iOpt,MckLbl,1,idum,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'nBas'
  call WrMck(iRC,iOpt,MckLbl,1,nBas,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'ldisp'
  call WrMck(iRC,iOpt,MckLbl,1,ldisp,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'chdisp'
  call cWrMck(iRC,iOpt,MckLbl,1,swlbl(1),iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'NISH'
  call WrMck(iRC,iOpt,MckLbl,1,nish,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'NASH'
  call WrMck(iRC,iOpt,MckLbl,1,nash,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
end if

if (SPINPOL) then
  call coeff(ralphas,rbetaa,rbetas)
  rms = real(ms2,kind=wp)*Half
  nAct = 0
  do iSym=1,nSym
    nAct = nAct+nAsh(iSym)
  end do
  nG = nAct**2
  nG2 = nAct**4
  call mma_allocate(famo_spinp,ndens2,Label='famo_spinp')
  call mma_allocate(famo_spinm,ndens2,Label='famo_spinm')
  call mma_allocate(G2mp,nG2,Label='G2mp')
  call mma_allocate(G2pp,nG2,Label='G2pp')
  call mma_allocate(G2mm,nG2,Label='G2mm')
  call mma_allocate(Fm,nG2,Label='Fm')
  call mma_allocate(Fp,nG2,Label='Fp')
  call mma_allocate(G1p,nG,Label='G1p')
  call mma_allocate(G1m,nG,Label='G1m')
  itype = 2
  call ipin(ipCI)
  call SpinDens(W(ipCI)%A,W(ipCI)%A,STATE_SYM,STATE_SYM,G2mm,G2mp,G2pp,Fm,Fp,G1m,G1p,iType)

  call mma_allocate(Tmp2,ndens2,Label='Tmp2')
  call mma_MaxDBLE(nMax)
  call mma_allocate(Tmp1,nMax/2,Label='Tmp1')

  call Ex_spin(G1p,FAMO_Spinp,Tmp1,nMax/2,Tmp2)
  call Ex_spin(G1m,FAMO_Spinm,Tmp1,nMax/2,Tmp2)

  call mma_deallocate(Tmp1)
  call mma_deallocate(Tmp2)
end if

end subroutine StPert
