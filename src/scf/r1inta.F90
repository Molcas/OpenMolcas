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

!#define _DEBUGPRINT_
subroutine R1IntA()
!***********************************************************************
!     purpose: Read one-electron hamiltonian and overlap matrix        *
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!***********************************************************************

use OneDat, only: sNoNuc, sNoOri
use InfSCF, only: nBT, OneHam, Ovrlp, Tot_Nuc_Charge
#ifdef _DEBUGPRINT_
use Index_Functions, only: nTri_Elem
use InfSCF, only: nBas, nSym
#endif
#ifdef _FDE_
use Embedding_Global, only: embInt, embPot, embPotInBasis, embPotPath
#endif
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iComp, iOpt, iRC, iSyLbl
character(len=8) :: Label
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ist1Hm, istOvl, iSym
#endif
#ifdef _FDE_
integer(kind=iwp) :: iDummyEmb, iEmb, iUnit
integer(kind=iwp), external :: isFreeUnit
#endif

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

#ifdef _FDE_
! Embedding
iDummyEmb = 0
call Get_iScalar('embpot',iDummyEmb)
if (iDummyEmb == 1) embPot = .true.
if (embPot) then
  call mma_allocate(embInt,nBT,Label='Emb')
  call EmbPotRdRun()
end if
#endif
! Allocate memory for one-electron integrals
call mma_allocate(OneHam,nBT,Label='OneHam')
call mma_allocate(Ovrlp,nBT+4,Label='Ovrlp')
OneHam(:) = Zero
Ovrlp(:) = Zero

! Read core Hamiltonian
iRc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iComp = 1
iSyLbl = 1
Label = 'OneHam  '
call RdOne(iRc,iOpt,Label,iComp,OneHam,iSyLbl)
if (iRc /= 0) then
  write(u6,*) 'R1Inta: Error readin ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
#ifdef _FDE_
! Embedding
if (embPot) then
  if (embPotInBasis) then
    ! If the potential is given in basis set representation it
    ! has not been calculated with a OneEl call and is just read
    ! from file here.
    iunit = isFreeUnit(1)
    call molcas_open(iunit,embPotPath)
    do iEmb=1,nBT
      read(iunit,*) embInt(iEmb)
    end do
    close(iunit)
  else
    ! Read one-electron integrals due to embedding potential
    iRc = -1
    Label = 'embpot  '
    call RdOne(iRc,iOpt,Label,iComp,embInt,iSyLbl)
    if (iRc /= 0) then
      write(u6,*) 'R1Inta: Error readin ONEINT'
      write(u6,'(A,A)') 'Label=',Label
      call Abend()
    end if
  end if
end if
#endif
#ifdef _DEBUGPRINT_
ist1Hm = 1
write(u6,*)
write(u6,*) ' One electron Hamiltonian at start'
write(u6,*) ' ---------------------------------'
do iSym=1,nSym
  write(u6,*) ' symmetry block',iSym
  call TriPrt(' ',' ',OneHam(ist1Hm),nBas(iSym))
  ist1Hm = ist1Hm+nTri_Elem(nBas(iSym))
end do
#endif

! Read overlap integrals and total effective nuclear charge
iRc = -1
iOpt = ibset(0,sNoOri)
iComp = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(iRc,iOpt,Label,iComp,Ovrlp,iSyLbl)
if (iRc /= 0) then
  write(u6,*) 'R1Inta: Error readin ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
Tot_Nuc_Charge = Ovrlp(nBT+4)
#ifdef _DEBUGPRINT_
istOvl = 1
write(u6,*)
write(u6,*) ' Overlap matrix at start'
write(u6,*) ' -----------------------'
do iSym=1,nSym
  write(u6,*) ' symmetry block',iSym
  call TriPrt(' ',' ',Ovrlp(istOvl),nBas(iSym))
  istOvl = istOvl+nTri_Elem(nBas(iSym))
end do
#endif

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine R1IntA
