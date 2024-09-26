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

subroutine SetUp_SCF()
!***********************************************************************
!                                                                      *
!     purpose: Set up needed parameters                                *
!                                                                      *
!***********************************************************************

use InfSCF, only: DSCF, kOV, MaxBas, MaxBOF, MaxBOO, MaxBXO, MaxOrb, MaxORF, MaxORO, mOV, nBas, nBB, nBO, nBT, nD, nFro, nnB, &
                  nnFr, nnO, nnOc, nOCC, nOFS, nOO, nOrb, nOV, nSym
use Definitions, only: iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: iSym, maxnOcc(MxSym), minnOcc(MxSym)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Set up global parameters used later
nnOc = 0
nnFr = 0
nnB = 0
nnO = 0
nBT = 0
nBO = 0
nBB = 0
nOO = 0
nOV = 0
mOV = 0
kOV(:) = 0
nOFS = 0
MaxBas = 0
MaxOrb = 0
MaxOrF = 0
MaxOrO = 0
MaxBxO = 0
MaxBOF = 0
MaxBOO = 0
do iSym=1,nSym
  if (nD == 1) then
    maxnOcc(iSym) = nOcc(iSym,1)
    minnOcc(iSym) = nOcc(iSym,1)
  else
    maxnOcc(iSym) = max(nOcc(iSym,1),nOcc(iSym,2))
    minnOcc(iSym) = min(nOcc(iSym,1),nOcc(iSym,2))
  end if
  if (nBas(iSym) > MxBas) then
    write(u6,*) 'SetUp: nBas(iSym) > MxBas'
    write(u6,*) 'nBas(iSym),MxBas=',nBas(iSym),MxBas
    call Abend()
  end if
  if (nOrb(iSym) > nBas(iSym)) then
    write(u6,*) 'SetUp: nOrb(iSym) > nBas(iSym)'
    write(u6,*) 'nOrb(iSym),nBas(iSym)=',nOrb(iSym),nBas(iSym)
    call Abend()
  end if
  if (maxnOcc(iSym) > nOrb(iSym)) then
    write(u6,*) 'iSym=',iSym
    write(u6,*) 'SetUp: nOcc(iSym) > nOrb(iSym)'
    write(u6,*) 'nOcc(iSym),nOrb(iSym)=',maxnOcc(iSym),nOrb(iSym)
    call Abend()
  end if
  if (nFro(iSym) > minnOcc(iSym)) then
    write(u6,*) 'SetUp: nFro(iSym) > nOcc(iSym)'
    write(u6,*) 'nFro(iSym),nOcc(iSym)=',nFro(iSym),minnOcc(iSym)
    call Abend()
  end if
  nnOc = nnOc+nOcc(iSym,1)
  if (nD == 2) nnOc = nnOc+nOcc(iSym,2)
  nnFr = nnFr+nFro(iSym)
  nnB = nnB+nBas(iSym)
  nnO = nnO+nOrb(iSym)
  nBT = nBT+nBas(iSym)*(nBas(iSym)+1)/2
  nBO = nBO+nBas(iSym)*nOrb(iSym)
  nBB = nBB+nBas(iSym)*nBas(iSym)
  nOO = nOO+nOrb(iSym)*nOrb(iSym)
  kOV(:) = kOV(:)+(nOcc(iSym,:)-nFro(iSym))*(nOrb(iSym)-nOcc(iSym,:))
  nOV = nOV+(maxnOcc(iSym)-nFro(iSym))*(nOrb(iSym)-minnOcc(iSym))
  nOFS = nOFS+(nOrb(iSym)-nFro(iSym))**2
  MaxBas = max(MaxBas,nBas(iSym))
  MaxOrb = max(MaxOrb,nOrb(iSym))
  MaxOrF = max(MaxOrF,nOrb(iSym)-nFro(iSym))
  MaxOrO = max(MaxOrO,nOrb(iSym)-minnOcc(iSym))
  MaxBxO = max(MaxBxO,nBas(iSym)*nOrb(iSym))
  MaxBOF = max(MaxBOF,nBas(iSym)*(nOrb(iSym)-nFro(iSym)))
  MaxBOO = max(MaxBOO,nBas(iSym)*(nOrb(iSym)-minnOcc(iSym)))
end do
mOV = kOV(1)+kOV(2)

if ((nnB > 2*MxBas) .and. (.not. DSCF)) then
  write(u6,*) 'SetUp: nnB > 2*MxBas .and. .not.DSCF'
  write(u6,*) 'nnB,MxBas=',nnB,MxBas
  call Abend()
else if ((nnB > 4*MxBas) .and. DSCF) then
  write(u6,*) 'SetUp: nnB > 4*MxBas .and. DSCF'
  write(u6,*) 'nnB,MxBas=',nnB,MxBas
  call Abend()
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine SetUp_SCF
