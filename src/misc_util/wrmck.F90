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
! Copyright (C) 1993, Per-Olof Widmark                                 *
!               1993, Markus P. Fuelscher                              *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine WrMCK(rc,Option,InLab,iComp,iData,iSymLab)
!***********************************************************************
!                                                                      *
!     Purpose: write data to one-electron integral file                *
!                                                                      *
!     Calling parameters:                                              *
!     rc      : Return code                                            *
!     Option  : Switch to set options                                  *
!     InLab   : Identifier for the data to write                       *
!     Comp    : Composite identifier to select components              *
!     iData   : contains on input the data to store on disk            *
!     SymLab  : symmetry label of the provided data                    *
!                                                                      *
!     Local data declarations:                                         *
!     Label   : character*8, used to covert incoming names             *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark and M.P. Fuelscher                                  *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!     Modified by:                                                     *
!     Anders Bernhardsson                                              *
!     University of Lund, Sweden, 1995                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use MckDat, only: AuxMck, LenOp, lTocMck, MxOp, NaN, NotNaN, nTitle, oAddr, oComp, oLabel, oSymLb, pASh, pBas, pchdisp, pdegdisp, &
                  pEnd, pFID, pish, pldisp, pndisp, pNext, pnrdisp, pOp, pPert, pSym, pSymOp, ptdisp, pTitle, pVersN, rcMck, sDbg, &
                  sLength, TocMck
use Definitions, only: iwp, u6, RtoI, ItoB

#include "intent.fh"

implicit none
integer(kind=iwp), intent(inout) :: rc
integer(kind=iwp), intent(in) :: Option, iComp, iSymLab
character(len=*), intent(in) :: InLab
integer(kind=iwp), intent(_IN_) :: iData(*)
#include "Molcas.fh"
integer(kind=iwp) :: Comp, i, iBas, icpi, iDisk, ij, ijS, iS, iSym, j, jBas, jS, k, Len_, Length, LuMCK, nA, SymLab, &
                     LabTmp(2)
logical(kind=iwp) :: Debug
!character(len=11) :: Label_Add
character(len=8) :: Label
character(len=*), parameter :: TheName = 'WrMck'
!real(kind=wp), external :: Check, DDot_

!----------------------------------------------------------------------*
! Start procedure:                                                     *
! Pick up the file definitions                                         *
!----------------------------------------------------------------------*
SymLab = iSymLab
icpi = ItoB
Comp = iComp
Len_ = rc
rc = rcMck%good
LuMCK = AuxMck%Lu
!----------------------------------------------------------------------*
! Check the file status                                                *
!----------------------------------------------------------------------*
if (.not. AuxMck%Opn) call SysFileMsg(TheName,'MSG: open',LuMck,' ')
!----------------------------------------------------------------------*
! Truncate the label to 8 characters and convert it to upper case      *
!----------------------------------------------------------------------*
!call StdFmt(InLab,Label)
Label = InLab
call UpCase(Label)
Length = len(Label)/ItoB
LabTmp(:Length) = transfer(Label,LabTmp,Length)
!----------------------------------------------------------------------*
! Print debugging information                                          *
!----------------------------------------------------------------------*
Debug = btest(option,sDbg)
if (Debug) then
  write(u6,*) '<<< Entering WrMck >>>'
  write(u6,'(a,z8)') ' rc on entry:     ',rc
  write(u6,'(a,a)') ' Label on entry:  ',Label
  write(u6,'(a,z8)') ' Comp on entry:   ',Comp
  write(u6,'(a,z8)') ' SymLab on entry: ',SymLab
  write(u6,'(a,z8)') ' Option on entry: ',Option
  write(u6,*) ' Contents of the Toc'
  write(u6,*) ' ==================='
  write(u6,'(i6,z8)') 'pFID,TocMck(pFID)=',pFID,TocMck(pFID)
  write(u6,'(i6,z8)') pVersN,TocMck(pVersN)
  write(u6,'(i6,z8)') pTitle,TocMck(pTitle)
  write(u6,'(i6,z8)') pOp,TocMck(pOp)
  write(u6,'(i6,z8)') pSym,TocMck(pSym)
  write(u6,'(i6,z8)') pSymOp,TocMck(pSymOp)
  write(u6,'(i6,z8)') pBas,TocMck(pBas)
  write(u6,'(i6,z8)') pNext,TocMck(pNext)
  write(u6,'(i6,z8)') pEnd,TocMck(pEnd)
end if
!----------------------------------------------------------------------*
! Store data in TocMck                                                 *
!----------------------------------------------------------------------*

if (Label == 'TITLE') then

  TocMck(pTitle) = NotNaN
  TocMck(pTitle+1:pTitle+nTitle) = iData(1:nTitle)
  !--------------------------------------------------------------------*

else if (Label == 'NSYM') then

  if ((iData(1) > MxSym) .or. (iData(1) < 1)) then
    call SysWarnMsg(TheName,'Label=',Label)
    call SysValueMsg('iData(1)=',iData(1))
  end if
  TocMck(pSym) = iData(1)
  !--------------------------------------------------------------------*

else if (Label == 'NBAS') then

  if (TocMck(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocMck(pSym)
  TocMck(pbas:pbas+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'NISH') then

  if (TocMck(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocMck(pSym)
  TocMck(pISH:pISH+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'NASH') then

  if (TocMck(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocMck(pSym)
  TocMck(pASH:pASH+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'LDISP') then

  if (TocMck(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocMck(pSym)
  TocMck(pldisp:pldisp+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'TDISP') then

  if (TocMck(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocMck(pndisp)
  TocMck(ptdisp:ptdisp+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'NDISP') then

  TocMck(pnDisp) = iData(1)
  !--------------------------------------------------------------------*

else if (label == 'CHDISP') then

  if (TocMck(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocMck(pndisp)*30/icpi+1
  TocMck(pchdisp:pchdisp+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'NRCTDISP') then

  if (TocMck(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocMck(pndisp)
  TocMck(pnrdisp:pnrdisp+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'DEGDISP ') then

  if (TocMck(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocMck(pndisp)
  TocMck(pdegdisp:pdegdisp+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (Label == 'SYMOP') then

  if (TocMck(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = (3*TocMck(pSym)+ItoB-1)/ItoB
  TocMck(pSymOp:pSymOp+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'PERT') then

  Length = 16/icpi
  TocMck(pPert:pPert+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else
  !--------------------------------------------------------------------*
  ! Store operators as individual records on disk                      *
  ! Note: If the incoming operator has already been stored             *
  ! previously (label, component and symmetry labels are identical)    *
  ! it will replace the existing one.                                  *
  !--------------------------------------------------------------------*
  if ((label == 'STATHESS') .or. (label == 'RESPHESS') .or. (label == 'CONNHESS') .or. (label == 'HESS') .or. (label == 'NUCGRAD') &
      .or. (label == 'TWOGRAD')) then
    Comp = 1
    SymLab = 1
  end if

  if (TocMck(pBas) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)

  k = 0
  do i=MxOp,1,-1
#   ifdef _I8_
    if ((TocMck(pOp+LenOp*(i-1)+oLabel) == LabTmp(1)) .and. (TocMck(pOp+LenOp*(i-1)+oComp) == Comp)) k = i
    !   .and. (TocMck(pOp+LenOp*(i-1)+oSymLb) == SymLab)) k = i
#   else
    if ((TocMck(pOp+LenOp*(i-1)+oLabel) == LabTmp(1)) .and. (TocMck(pOp+LenOp*(i-1)+oLabel+1) == LabTmp(2)) .and. &
        (TocMck(pOp+LenOp*(i-1)+oComp) == Comp)) k = i
    !   .and. (TocMck(pOp+LenOp*(i-1)+oSymLb) == SymLab)) k = i
#   endif
  end do

  iDisk = TocMck(pOp+LenOp*(k-1)+oAddr)
  if (k == 0) then
    do i=MxOp,1,-1
      if (TocMck(pOp+LenOp*(i-1)+oLabel) == NaN) k = i
    end do
    iDisk = TocMck(pNext)
    if (Debug) then
      write(u6,*) ' This is a new field!'
      write(u6,*) ' iDisk=',iDisk
      write(u6,*) ' FldNo=',k
      write(u6,*) ' pNext=',pNext
    end if
  else
    if (Debug) then
      write(u6,*) ' This is an old field!'
      write(u6,*) ' iDisk=',iDisk
      write(u6,*) ' FldNo=',k
      write(u6,*) ' pNext=',pNext
    end if
  end if
  if (k == 0) call SysAbendMsg(TheName,'Undefined Label:',Label)
  !write(u6,*) isymlab,label
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (Label == 'MOPERT') then

    nA = 0
    do i=0,TocMck(psym)-1
      nA = TocMck(pAsh+i)+nA
    end do
    Length = nTri_Elem(nA)
    Length = nTri_Elem(Length)
    !write(u6,*) Length
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (label == 'NUCGRAD') then
    Comp = 1
    SymLab = 1
    if (TocMck(pldisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
    Length = TocMck(pldisp)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (label == 'TWOGRAD') then
    Comp = 1
    SymLab = 1
    if (TocMck(pldisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
    Length = TocMck(pldisp)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if ((label == 'STATHESS') .or. (label == 'RESPHESS') .or. (label == 'CONNHESS') .or. (label == 'HESS')) then
    Comp = 1
    SymLab = 1
    if (TocMck(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
    Length = 0
    do iSym=0,TocMck(pSym)-1
      Length = Length+nTri_Elem(TocMck(pldisp+isym))
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (label == 'INACTIVE') then
    Length = 0
    do iS=1,TocMck(pSym)
      do jS=1,TocMck(pSym)
        ijS = Mul(iS,jS)-1
        if (btest(iSymLab,ijS)) then
          jBas = TocMck(pbas+jS-1)
          iBas = TocMck(pbas+iS-1)
          if (jBas == NaN) call SysAbendMsg(TheName,'jBas == NaN at label',Label)
          if (iBas == NaN) call SysAbendMsg(TheName,'iBas == NaN at label',Label)
          Length = Length+iBas*jBas
        end if
      end do
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (label == 'TOTAL') then
    Length = 0
    do iS=1,TocMck(pSym)
      do jS=1,TocMck(pSym)
        ijS = Mul(iS,jS)-1
        if (btest(iSymLab,ijS)) then
          jBas = TocMck(pbas+jS-1)
          iBas = TocMck(pbas+iS-1)
          if (jBas == NaN) call SysAbendMsg(TheName,'jBas == NaN at label',Label)
          if (iBas == NaN) call SysAbendMsg(TheName,'iBas == NaN at label',Label)
          Length = Length+iBas*jBas
        end if
      end do
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else
    Length = 0
    do i=1,TocMck(pSym)
      do j=1,i
        ij = Mul(i,j)-1
        if (btest(iSymLab,ij)) then
          if (i == j) then
            Length = Length+nTri_Elem(TocMck(pBas-1+i))
          else
            Length = Length+TocMck(pBas-1+i)*TocMck(pBas-1+j)
          end if
        end if
      end do
    end do
    if (btest(Option,sLength)) Length = Len_
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  Length = RtoI*Length
  TocMck(pOp+LenOp*(k-1)+oLabel) = LabTmp(1)
# ifndef _I8_
  TocMck(pOp+LenOp*(k-1)+oLabel+1) = LabTmp(2)
# endif
  TocMck(pOp+LenOp*(k-1)+oComp) = Comp
  TocMck(pOp+LenOp*(k-1)+oSymLb) = iSymLab
  TocMck(pOp+LenOp*(k-1)+oAddr) = iDisk
  !write(u6,*) Length,idisk
  call iDAFile(LuMCK,1,iData,Length,iDisk)
  !if ((Label == 'TOTAL') .or. (Label == 'INACTIVE') .or. (Label == 'MOPERT')) then
  !  !write(u6,*) 'iComp=',iComp
  !  !call RecPrt(Label,' ',iData,1,Length/RtoI)
  !  Check = DDot_(Length/RtoI,iData,1,iData,1)
  !  Label_Add = ' '
  !  Label_Add = Label
  !  write(Label_Add(9:11),'(I3.3)') iComp
  !  !write(u6,*) Label_Add,Label,iComp,Check(1)
  !  if (Label == 'TOTAL') then
  !    call Add_Info(Label_Add,Check,1,6)
  !  else if (Label == 'INACTIVE') then
  !    call Add_Info(Label_Add,Check,1,7)
  !  else
  !    call Add_Info(Label_Add,Check,1,1)
  !  end if
  !end if
  !TocMck(pNext) = iDisk
  TocMck(pNext) = max(TocMck(pNext),iDisk)
end if
!----------------------------------------------------------------------*
! Finally copy the TocMck back to disk                                 *
!----------------------------------------------------------------------*
iDisk = 0
call iDaFile(LuMCK,1,TocMck,lTocMck,iDisk)
!----------------------------------------------------------------------*
! Print debugging information                                          *
!----------------------------------------------------------------------*
if (Debug) then
  write(u6,*) '<<< Exiting WrMck >>>'
  write(u6,'(a,z8)') ' rc on exit:     ',rc
  write(u6,'(a,a)') ' Label on exit:  ',Label
  write(u6,'(a,z8)') ' Comp on exit:   ',Comp
  write(u6,'(a,z8)') ' SymLab on exit: ',SymLab
  write(u6,'(a,z8)') ' Option on exit: ',Option
end if

!----------------------------------------------------------------------*
! Terminate procedure                                                  *
!----------------------------------------------------------------------*
return

end subroutine WrMCK
