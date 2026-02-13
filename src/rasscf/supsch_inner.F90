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
! Copyright (C) 1997, Luis Serrano-Andres                              *
!***********************************************************************

subroutine SUPSCH_INNER(SMAT,CMOO,CMON,Temp1,Temp2,nOrbMX,IxSym2,nOrb_tot)
! Program RASSCF
!
! Objective: To check the order of the input of natural orbitals
!            to obtain the right labels for the Supersymmetry matrix.
!
! Called from ortho, neworb, fckpt2, and natorb.
!
! Luis Serrano-Andres
! University of Lund, Sweden, 1997
! **** Molcas-4 *** Release 97 04 01 **********

use OneDat, only: sNoNuc, sNoOri
use rasscf_global, only: FDIAG, iSupSM, Iter, ixsym
use PrintLevel, only: DEBUG
use output_ras, only: IPRLOC
use general_data, only: NSYM, NBAS
use Constants, only: Zero, One
use Definitions, only: u6

implicit none
integer nOrbMX, nOrb_Tot
real*8 CMOO(*), CMON(*), SMAT(*)
real*8 Temp1(nOrbMX*nOrbMX), Temp2(nOrbMX*nOrbMX)
integer IxSym2(nOrb_tot)
character(len=16), parameter :: ROUTINE = 'SUPSCH_INNER'
integer pSij
real*8 DUM(1)
character(len=8) :: Label
integer :: i_Component, i_Opt, i_RC, i_SymLbl, iGroup, iLabel, iOrb, iOrder, iPrLev, iSafe, jOrb, kCof, kGroup, kOrb, nBs, nnOrb, &
           nOGr1, nOGr2, iSym
real*8 :: OldOvlp, Ovlp1, Ovlp2, xOvlp
#include "warnings.h"

! Local print level (if any)
IPRLEV = IPRLOC(4)
if (IPRLEV >= DEBUG) write(u6,*) ' Entering ',ROUTINE

! Read overlap matrix SMAT:

i_Rc = 0
i_Opt = ibset(ibset(0,sNoOri),sNoNuc)
i_Component = 1
i_SymLbl = 1
Label = 'Mltpl  0'
call RdOne(i_Rc,i_Opt,Label,i_Component,Smat,i_SymLbl)
if (i_Rc /= 0) then
  write(u6,*)
  write(u6,*) ' ********************* ERROR **********************'
  write(u6,*) ' SUPSCH: Failed to read overlap from ONEINT.       '
  write(u6,*) ' RASSCF is using overlaps to compare old and new   '
  write(u6,*) ' orbitals, but could not read overlaps from ONEINT.'
  write(u6,*) ' Something is wrong with the file, or possibly with'
  write(u6,*) ' the program. Please check.                        '
  write(u6,*) ' **************************************************'
  call Quit(_RC_IO_ERROR_READ_)
end if

! Check order of the orbitals for supersymmetry option

if ((ISupsm == 1) .and. (Iter >= 1)) then

  if (IPRLEV >= DEBUG) then
    call PRIMO_RASSCF('Testing old orb for supersymmetry',FDIAG,DUM,CMOO)
    call PRIMO_RASSCF('Testing new orb for supersymmetry',FDIAG,DUM,CMON)
  end if

  kOrb = 0
  kCof = 0
  pSij = 0
  do iOrb=1,nOrb_tot
    IxSym2(iOrb) = 0
  end do
  do iSym=1,nSym
    iSafe = 0
    nBs = nBas(iSym)
    if (nBs <= 0) goto 1966

    ! Computing orbital overlaping Sum(p,q) C1kp* C2lq Spq

    call Square(SMat(pSij+1),Temp1,1,nBs,nBs)
    call DGEMM_('N','N',nBs,nBs,nBs,One,Temp1,nBs,CMON(kCof+1),nBs,Zero,Temp2,nBs)
    call DGEMM_('T','N',nBs,nBs,nBs,One,CMOO(kCof+1),nBs,Temp2,nBs,Zero,Temp1,nBs)

    ! Checking the maximum overlap between the orbitals
    ! and building the new supersymmetry matrix

    nnOrb = 1
    do iOrb=1,nBs
      iLabel = IxSym(kOrb+iOrb)
      iOrder = 1
      xOvlp = Zero
      do jOrb=0,nBs-2
        Ovlp1 = abs(Temp1(nnOrb+(jOrb*nBs)))
        Ovlp2 = abs(Temp1(nnOrb+((jOrb+1)*nBs)))
        OldOvlp = xOvlp
        xOvlp = max(Ovlp1,Ovlp2,xOvlp)
        if (jOrb == 0) then
          if (Ovlp1 < Ovlp2) iOrder = 2
        else
          if (xOvlp /= OldOvlp) iOrder = jOrb+2
        end if
      end do
      IxSym2(kOrb+iOrder) = iLabel
      nnOrb = nnOrb+1
    end do

    ! Number of orbital groups on the symmetry

    kGroup = 0
    do iOrb=1,nBs-1
      kGroup = max(IxSym(kOrb+iOrb),IxSym(kOrb+iOrb+1),kGroup)
    end do
    do iGroup=0,kGroup
      nOGr1 = 0
      nOGr2 = 0
      do iOrb=1,nBs
        if (IxSym(kOrb+iOrb) == iGroup) nOGr1 = nOGr1+1
        if (IxSym2(kOrb+iOrb) == iGroup) nOGr2 = nOGr2+1
      end do

      ! Checking if we have the same number of orbitals per group

      if ((nOGr2 /= nOGr1) .and. (iGroup /= 0)) then
        iSafe = 1
        call WarningMessage(1,'Supersymmetry may have failed.')
        write(u6,*) ' Check orbital order or try cleaning orbitals.'
      end if
    end do

    ! New matrix replaces the old one

    if (iSafe == 0) then
      do iOrb=1,nBs
        IxSym(kOrb+iOrb) = IxSym2(kOrb+iOrb)
      end do
    end if
    kCof = kCof+(nBs*nBs)
    pSij = pSij+(nBs*nBs+nBs)/2
    kOrb = kOrb+nBs
1966 continue
  end do
end if

end subroutine SUPSCH_INNER
