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
! Copyright (C) 1989, Bjorn O. Roos                                    *
!               1989, Per Ake Malmqvist                                *
!               1991,1993,1996, Markus P. Fuelscher                    *
!***********************************************************************

subroutine Fmat(CMO,PUVX,D,D1A,FI,FA)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Update the Fock matrix for the active orbitals and transform     *
!     it to MO basis as well as the matrix FI (Fock matrix) for        *
!     frozen and inactive orbitals).                                   *
!                                                                      *
!     calling arguments:                                               *
!     CMO     : array of real*8                                        *
!               MO-coefficients                                        *
!     PUVX    : array of real*8                                        *
!               two-electron integrals (pu!vx)                         *
!     D       : array of real*8                                        *
!               averaged one-body density matrix                       *
!     D1A     : array of real*8                                        *
!               active one body density matrix in AO-basis             *
!     FI      : array of real*8                                        *
!               inactive Fock matrix. In input is in AO basis.         *
!               In output in MO basis.                                 *
!     FA      : array of real*8                                        *
!               active Fock matrix. In input in AO Basis.              *
!               In output in MO basis. It is also modified by scaling  *
!               exchange part in ExFac /= 1.0d0                        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     B.O. Roos and P.Aa. Malmqvist                                    *
!     University of Lund, Sweden, 1989                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history:                                                         *
!     - updated for MOLCAS version 2                                   *
!       M.P. Fuelscher, University of Lund, Sweden, 1991               *
!     - updated for MOLCAS version 3                                   *
!       M.P. Fuelscher, University of Lund, Sweden, 1993               *
!     - updated for integral direct and reaction field calculations    *
!       M.P. Fuelscher, University of Lund, Sweden, 1996               *
!                                                                      *
!***********************************************************************

use RunFile_procedures, only: Get_dExcdRa
use rasscf_global, only: KSDFT, DFTFOCK, ECAS, EMY, ExFac, NAC, NewFock, nFint, VIA, VIA_DFT, l_casdft
use PrintLevel, only: DEBUG
use output_ras, only: LF, IPRLOC
use general_data, only: NSYM, NTOT1, NASH, NBAS, NFRO, NISH, NORB
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
real*8 CMO(*), PUVX(*), D(*), D1A(*), FI(*), FA(*)
character(len=16), parameter :: ROUTINE = 'FMAT    '
real*8, allocatable :: TmpFck(:), Tmp1(:), Tmp2(:), TmpD1A(:)
integer iPrLev, iOff, iSym, iBas, i, iFro, ij, iOff1, iOff2, iOff3, iOrb, ipTmpFck, ipTmpFckA, ipTmpFckI, j, jOrb, nTmpFck
real*8, external :: DDot_

! Local print level (if any)
IPRLEV = IPRLOC(4)
if (IPRLEV >= DEBUG) then
  write(LF,*) ' Entering ',ROUTINE

  write(6,*) repeat('*',65)
  write(6,*) 'Entering FMAT routine called by SXCTL!'
  write(6,*) repeat('*',65)
  write(6,*) 'printing input matrices :'
  write(6,*) repeat('*',65)
  write(LF,*)
  write(LF,*) ' CMOs in FMAT'
  write(LF,*) ' ---------------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    if (iBas /= 0) then
      write(6,*) 'Sym =',iSym
      do i=1,iBas
        write(6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
      end do
      iOff = iOff+(iBas*iBas)
    end if
  end do

  write(LF,*)
  write(LF,*) ' PUVX in FMAT'
  write(LF,*) ' ---------------------'
  write(LF,*)
  call wrtmat(PUVX,1,nFint,1,nFint)

  write(LF,*)
  write(LF,*) ' ---------------------'
  call TRIPRT('Averaged one-body density matrix D, in MO in FMAT',' ',D,NAC)

  write(LF,*)
  write(LF,*) ' D1A in AO basis in FMAT'
  write(LF,*) ' ---------------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call wrtmat(D1A(iOff),iBas,iBas,iBas,iBas)
    iOff = iOff+iBas*iBas
  end do

  write(LF,*)
  write(LF,*) ' FI in AO-basis in FMAT'
  write(LF,*) ' --------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iOrb = nOrb(iSym)
    call TriPrt(' ',' ',FI(iOff),iOrb)
    iOff = iOff+(iOrb*iOrb+iOrb)/2
  end do

  write(LF,*)
  write(LF,*) ' FA in AO-basis in FMAT'
  write(LF,*) ' --------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iOrb = nOrb(iSym)
    call TriPrt(' ',' ',FA(iOff),iOrb)
    iOff = iOff+(iOrb*iOrb+iOrb)/2
  end do
end if

! create FA in AO basis
call mma_allocate(Tmp1,nTot1,Label='Tmp1')
call Fold(nSym,nBas,D1A,Tmp1)
if (KSDFT /= 'SCF') NewFock = 0
!if (NewFock == 0) then
!  nBMX = 0
!  do iSym=1,nSym
!    nBMX = max(nBMX,nBas(iSym))
!  end do
!  FA(1:nTot1) = Zero
!  call FTwo_Drv(nSym,nBas,nAsh,nSkipX,Tmp1,D1A,FA,nTot1,ExFac,nBMX,CMO)
!end if

! Inactive-active contribution to ECAS
VIA = dDot_(nTot1,FI,1,Tmp1,1)
ECAS = EMY+VIA
if (iPrLev >= DEBUG) then
  write(LF,*) ' Total core energy fmat:       ',EMY
  write(LF,*) ' inactive-active interaction:  ',VIA
  write(LF,*) ' CAS energy (core+interaction):',ECAS
end if
call mma_deallocate(Tmp1)

! print FI and FA
if (iPrLev >= DEBUG) then
  write(LF,*)
  write(LF,*) ' FI in AO-basis in fmat'
  write(LF,*) ' --------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iOrb = nOrb(iSym)
    call TriPrt(' ',' ',FI(iOff),iOrb)
    iOff = iOff+(iOrb*iOrb+iOrb)/2
  end do
  write(LF,*)
  write(LF,*) ' FA in AO-basis in fmat'
  write(LF,*) ' --------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iOrb = nOrb(iSym)
    call TriPrt(' ',' ',FA(iOff),iOrb)
    iOff = iOff+(iOrb*iOrb+iOrb)/2
  end do
end if

! transform FI from AO to MO basis
iOff1 = 1
iOff2 = 1
iOff3 = 1
do iSym=1,nSym
  iBas = nBas(iSym)
  if (iBas == 0) cycle
  iOrb = nOrb(iSym)
  if (iOrb == 0) cycle
  iFro = nFro(iSym)
  call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
  call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
  call Square(FI(iOff1),Tmp1,1,iBas,iBas)
  call DGEMM_('N','N',iBas,iOrb,iBas,1.0d0,Tmp1,iBas,CMO(iOff2+(iFro*iBas)),iBas,0.0d0,Tmp2,iBas)
  call DGEMM_Tri('T','N',iOrb,iOrb,iBas,1.0d0,Tmp2,iBas,CMO(iOff2+(iFro*iBas)),iBas,0.0d0,FI(iOff3),iOrb)
  call mma_deallocate(Tmp2)
  call mma_deallocate(Tmp1)
  iOff1 = iOff1+(iBas*iBas+iBas)/2
  iOff2 = iOff2+iBas*iBas
  iOff3 = iOff3+(iOrb*iOrb+iOrb)/2
end do

! transform FA from AO to MO basis
iOff1 = 1
iOff2 = 1
iOff3 = 1
do iSym=1,nSym
  iBas = nBas(iSym)
  if (iBas == 0) cycle
  iOrb = nOrb(iSym)
  if (iOrb == 0) cycle
  iFro = nFro(iSym)
  call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
  call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
  call Square(FA(iOff1),Tmp1,1,iBas,iBas)
  call DGEMM_('N','N',iBas,iOrb,iBas,1.0d0,Tmp1,iBas,CMO(iOff2+(iFro*iBas)),iBas,0.0d0,Tmp2,iBas)
  call DGEMM_Tri('T','N',iOrb,iOrb,iBas,1.0d0,Tmp2,iBas,CMO(iOff2+(iFro*iBas)),iBas,0.0d0,FA(iOff3),iOrb)
  call mma_deallocate(Tmp2)
  call mma_deallocate(Tmp1)
  iOff1 = iOff1+(iBas*iBas+iBas)/2
  iOff2 = iOff2+iBas*iBas
  iOff3 = iOff3+(iOrb*iOrb+iOrb)/2
end do

!***********************************************************************
!              Add DFT part to Fock matrix:                            *
!***********************************************************************
if ((KSDFT(1:3) /= 'SCF') .and. (KSDFT(1:3) /= 'PAM') .and. (.not. l_casdft)) then
  ipTmpFckI = -99999
  ipTmpFckA = -99999
  call Get_dExcdRa(TmpFck,nTmpFck)
  ipTmpFck = 1
  if (nTmpFck == NTOT1) then
    ipTmpFckI = ipTmpFck
  else if (nTmpFck == 2*NTOT1) then
    ipTmpFckI = ipTmpFck
    ipTmpFckA = ipTmpFck+nTot1
  else
    write(LF,*) ' Somethings wrong in dim. DFT',nTmpFck
    call Abend()
  end if
  call mma_allocate(TmpD1A,nTot1,Label='TmpD1A')
  call Fold(nSym,nBas,D1A,TmpD1A)
  VIA_DFT = dDot_(nTot1,TmpFck(ipTmpFckI),1,TmpD1A,1)
  call mma_deallocate(TmpD1A)

  ! Transform alpha density from AO to MO

  iOff1 = 1
  iOff2 = 1
  iOff3 = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    if (iBas == 0) cycle
    iOrb = nOrb(iSym)
    if (iOrb == 0) cycle
    iFro = nFro(iSym)
    call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
    call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
    call Square(TmpFck(ipTmpFckI+iOff1-1),Tmp1,1,iBas,iBas)
    call DGEMM_('N','N',iBas,iOrb,iBas,1.0d0,Tmp1,iBas,CMO(iOff2+(iFro*iBas)),iBas,0.0d0,Tmp2,iBas)
    call DGEMM_Tri('T','N',iOrb,iOrb,iBas,1.0d0,Tmp2,iBas,CMO(iOff2+(iFro*iBas)),iBas,0.0d0,TmpFck(ipTmpFckI+iOff3-1),iOrb)
    call mma_deallocate(Tmp2)
    call mma_deallocate(Tmp1)
    iOff1 = iOff1+(iBas*iBas+iBas)/2
    iOff2 = iOff2+iBas*iBas
    iOff3 = iOff3+(iOrb*iOrb+iOrb)/2
  end do

  ! Transform Active DFT Fock from AO to MO

  if (ipTmpFckA /= -99999) then
    iOff1 = 1
    iOff2 = 1
    iOff3 = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      if (iBas == 0) cycle
      iOrb = nOrb(iSym)
      if (iOrb == 0) cycle
      iFro = nFro(iSym)
      call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
      call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
      call Square(TmpFck(ipTmpFckA+iOff1-1),Tmp1,1,iBas,iBas)
      call DGEMM_('N','N',iBas,iOrb,iBas,1.0d0,Tmp1,iBas,CMO(iOff2+(iFro*iBas)),iBas,0.0d0,Tmp2,iBas)
      call DGEMM_Tri('T','N',iOrb,iOrb,iBas,1.0d0,Tmp2,iBas,CMO(iOff2+(iFro*iBas)),iBas,0.0d0,TmpFck(ipTmpFckA+iOff3-1),iOrb)
      call mma_deallocate(Tmp2)
      call mma_deallocate(Tmp1)
      iOff1 = iOff1+(iBas*iBas+iBas)/2
      iOff2 = iOff2+iBas*iBas
      iOff3 = iOff3+(iOrb*iOrb+iOrb)/2
    end do
  end if

  !if (DFTFOCK(1:4) /= 'ROKS') then
  !  write(LF,*) ' Just add a,b to FA,FI',DFTFOCK(1:4)
  !else
  !  write(LF,*) ' ROKS formula',DFTFOCK(1:4)
  !end if

  if (DFTFOCK(1:4) /= 'ROKS') then
    call daxpy_(NTOT1,1.0d0,TmpFck(ipTmpFckI),1,FI,1)
    if (ipTmpFckA /= -99999) call daxpy_(NTOT1,1.0d0,TmpFck(ipTmpFckA),1,FA,1)
  else if (DFTFOCK(1:4) == 'ROKS') then
    iOff1 = 0
    do iSym=1,nSym
      do iOrb=1,nOrb(iSym)
        do jOrb=1,iOrb
          ij = iOff1+iOrb*(iOrb-1)/2+jOrb
          if (iOrb <= nIsh(iSym)) FI(ij) = FI(ij)+0.5d0*(TmpFck(ipTmpFckI+ij-1)+TmpFck(ipTmpFckA+ij-1))
          if ((iOrb > nIsh(iSym)) .and. (iOrb <= nIsh(iSym)+nAsh(iSym))) then
            if (jOrb <= nIsh(iSym)) then
              FI(ij) = FI(ij)+TmpFck(ipTmpFckA+ij-1)
            else
              FI(ij) = FI(ij)+0.5d0*(TmpFck(ipTmpFckI+ij-1)+TmpFck(ipTmpFckA+ij-1))
            end if
          end if
          if (iOrb > nIsh(iSym)+nAsh(iSym)) then
            if ((jOrb > nIsh(iSym)) .and. (jOrb <= nIsh(iSym)+nAsh(iSym))) then
              FI(ij) = FI(ij)+TmpFck(ipTmpFckI+ij-1)
            else
              FI(ij) = FI(ij)+0.5d0*(TmpFck(ipTmpFckI+ij-1)+TmpFck(ipTmpFckA+ij-1))
            end if

          end if
        end do
      end do
      iOff1 = iOff1+(nOrb(iSym)*nOrb(iSym)+nOrb(iSym))/2
    end do
  else
    write(LF,*) ' Not implemented yet'
  end if
  call mma_deallocate(TmpFck)
end if
!***********************************************************************
if (iPrLev >= DEBUG) then
  write(LF,*)
  write(LF,*) ' FA in MO-basis in fmat'
  write(LF,*) ' --------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iOrb = nOrb(iSym)
    call TriPrt(' ',' ',FA(iOff),iOrb)
    iOff = iOff+(iOrb*iOrb+iOrb)/2
  end do
end if
! update Fock matrix by rescaling exchange term...
if (NewFock == 1) call Upd_FA(PUVX,FA,D,ExFac)

! print FI and FA
if (iPrLev >= DEBUG) then
  write(LF,*)
  write(LF,*) ' FI in MO-basis in fmat'
  write(LF,*) ' --------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iOrb = nOrb(iSym)
    call TriPrt(' ',' ',FI(iOff),iOrb)
    iOff = iOff+(iOrb*iOrb+iOrb)/2
  end do
end if
if (iPrLev >= DEBUG) then
  write(LF,*)
  write(LF,*) ' FA in MO-basis in fmat after upd_FA'
  write(LF,*) ' --------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iOrb = nOrb(iSym)
    call TriPrt(' ',' ',FA(iOff),iOrb)
    iOff = iOff+(iOrb*iOrb+iOrb)/2
  end do
end if

end subroutine Fmat
