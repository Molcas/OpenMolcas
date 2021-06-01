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
! Copyright (C) 2007, Francesco Aquilante                              *
!***********************************************************************

subroutine Cho_SOSmp2_Setup(irc)
! Francesco Aquilante   May 2007.
!
! Purpose: setup of SOS-MP2 program.

#include "implicit.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
!***********************************************************************
MulD2h(i,j) = ieor(i-1,j-1)+1
!***********************************************************************

irc = 0

! Setup index arrays and counters.
! --------------------------------

if (DecoMP2 .and. (ThrMP2 <= 0.0d0)) then
  call Get_dScalar('Cholesky Threshold',ThrMP2)
end if

call ChoMP2_GetInf(nOrb,nOcc,nFro,nDel,nVir)
iOcc(1) = 0
iVir(1) = 0
nOccT = nOcc(1)
nVirT = nVir(1)
do iSym=2,nSym
  iOcc(iSym) = nOccT
  iVir(iSym) = nVirT
  nOccT = nOccT+nOcc(iSym)
  nVirT = nVirT+nVir(iSym)
end do

do iSym=1,nSym
  nT1am(iSym) = 0
  do iSymi=1,nSym
    iSyma = MulD2h(iSymi,iSym)
    iT1am(iSyma,iSymi) = nT1am(iSym)
    nT1am(iSym) = nT1am(iSym)+nVir(iSyma)*nOcc(iSymi)
  end do
end do

do iSym=1,nSym
  nT1AOT(iSym) = 0
  do iSymAl=1,nSym
    iSymi = MulD2h(iSymAl,iSym)
    iT1AOT(iSymi,iSymAl) = nT1AOT(iSym)
    nT1AOT(iSym) = nT1AOT(iSym)+nOcc(iSymi)*nBas(iSymAl)
  end do
end do

do iSym=1,nSym
  nAOVir(iSym) = 0
  do iSyma=1,nSym
    iSymAl = MulD2h(iSyma,iSym)
    iAOVir(iSymAl,iSyma) = nAOVir(iSym)
    nAOVir(iSym) = nAOVir(iSym)+nBas(iSymAl)*nVir(iSyma)
  end do
end do

if (ChoAlg == 2) then
  do iSym=1,nSym
    nMatab(iSym) = 0
    do iSymb=1,nSym
      iSyma = MulD2h(iSymb,iSym)
      iMatab(iSyma,iSymb) = nMatab(iSym)
      nMatab(iSym) = nMatab(iSym)+nVir(iSyma)*nVir(iSymb)
    end do
  end do
else
  call Cho_iZero(nMatab,8)
  call Cho_iZero(iMatab,64)
end if

! If batching over occuped orbitals is forced by user, then turn it Off!
! ----------------------------------------------------------------------

ForceBatch = .false.

nBatch = 1

! Initialize file units.
! ----------------------

do iSym=1,nSym
  do iTyp=1,nTypF
    call ChoMP2_OpenF(0,iTyp,iSym)
  end do
end do

end subroutine Cho_SOSmp2_Setup
