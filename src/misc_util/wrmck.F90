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
!     Global data declarations (Include file):                         *
!     Parm    : Table of paramaters                                    *
!     rcParm  : Table of return codes                                  *
!     Switch  : Table of options                                       *
!     Common  : Common block containing TocOne                         *
!     Data    : Data definitions                                       *
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
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: rc, Option, iComp, iData(*), iSymLab
character(len=*) :: InLab
#include "MckDat.fh"
integer(kind=iwp) :: Comp, i, iBas, icpi, iDisk, ij, ijS, iS, isopen, iSym, j, jBas, jS, k, Len_, Length, LuMCK, nA, SymLab, &
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
icpi = itob
Comp = iComp
Len_ = rc
rc = rc0000
LuMCK = AuxMCK(pLu)
isopen = AuxMCK(pOpen)
!----------------------------------------------------------------------*
! Check the file status                                                *
!----------------------------------------------------------------------*
if (isopen /= 1) call SysFileMsg(TheName,'MSG: open',LuMck,' ')
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
Debug = btest(option,10)
if (Debug) then
  write(u6,*) '<<< Entering WrMck >>>'
  write(u6,'(a,z8)') ' rc on entry:     ',rc
  write(u6,'(a,a)') ' Label on entry:  ',Label
  write(u6,'(a,z8)') ' Comp on entry:   ',Comp
  write(u6,'(a,z8)') ' SymLab on entry: ',SymLab
  write(u6,'(a,z8)') ' Option on entry: ',Option
  write(u6,*) ' Contents of the Toc'
  write(u6,*) ' ==================='
  write(u6,'(i6,z8)') 'pFID,TocOne(pFID)=',pFID,TocOne(pFID)
  write(u6,'(i6,z8)') pVersN,TocOne(pVersN)
  write(u6,'(i6,z8)') pTitle,TocOne(pTitle)
  write(u6,'(i6,z8)') pOp,TocOne(pOp)
  write(u6,'(i6,z8)') pSym,TocOne(pSym)
  write(u6,'(i6,z8)') pSymOp,TocOne(pSymOp)
  write(u6,'(i6,z8)') pBas,TocOne(pBas)
  write(u6,'(i6,z8)') pNext,TocOne(pNext)
  write(u6,'(i6,z8)') pEnd,TocOne(pEnd)
end if
!----------------------------------------------------------------------*
! Store data in TocOne                                                 *
!----------------------------------------------------------------------*

if (Label == 'TITLE') then

  TocOne(pTitle) = NotNaN
  TocOne(pTitle+1:pTitle+nTitle) = iData(1:nTitle)
  !--------------------------------------------------------------------*

else if (Label == 'NSYM') then

  if ((iData(1) > MxSym) .or. (iData(1) < 1)) then
    call SysWarnMsg(TheName,'Label=',Label)
    call SysValueMsg('iData(1)=',iData(1))
  end if
  TocOne(pSym) = iData(1)
  !--------------------------------------------------------------------*

else if (Label == 'NBAS') then

  if (TocOne(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocOne(pSym)
  TocOne(pbas:pbas+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'NISH') then

  if (TocOne(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocOne(pSym)
  TocOne(pISH:pISH+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'NASH') then

  if (TocOne(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocOne(pSym)
  TocOne(pASH:pASH+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'LDISP') then

  if (TocOne(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocOne(pSym)
  TocOne(pldisp:pldisp+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'TDISP') then

  if (TocOne(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocOne(pndisp)
  TocOne(ptdisp:ptdisp+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'NDISP') then

  ToCOne(pnDisp) = iData(1)
  !--------------------------------------------------------------------*

else if (label == 'CHDISP') then

  if (TocOne(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocOne(pndisp)*30/icpi+1
  TocOne(pchdisp:pchdisp+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'NRCTDISP') then

  if (TocOne(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocOne(pndisp)
  TocOne(pnrdisp:pnrdisp+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'DEGDISP ') then

  if (TocOne(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocOne(pndisp)
  TocOne(pdegdisp:pdegdisp+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (Label == 'SYMOP') then

  if (TocOne(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = (3*TocOne(pSym)+ItoB-1)/ItoB
  TocOne(pSymOp:pSymOp+Length-1) = iData(1:Length)
  !--------------------------------------------------------------------*

else if (label == 'PERT') then

  Length = 16/icpi
  TocOne(pPert:pPert+Length-1) = iData(1:Length)
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

  if (TocOne(pBas) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)

  k = 0
  do i=MxOp,1,-1
#   ifdef _I8_
    if ((TocOne(pOp+LenOp*(i-1)+oLabel) == LabTmp(1)) .and. (TocOne(pOp+LenOp*(i-1)+oComp) == Comp)) k = i
    !   .and. (TocOne(pOp+LenOp*(i-1)+oSymLb) == SymLab)) k = i
#   else
    if ((TocOne(pOp+LenOp*(i-1)+oLabel) == LabTmp(1)) .and. (TocOne(pOp+LenOp*(i-1)+oLabel+1) == LabTmp(2)) .and. &
        (TocOne(pOp+LenOp*(i-1)+oComp) == Comp)) k = i
    !   .and. (TocOne(pOp+LenOp*(i-1)+oSymLb) == SymLab)) k = i
#   endif
  end do

  iDisk = TocOne(pOp+LenOp*(k-1)+oAddr)
  if (k == 0) then
    do i=MxOp,1,-1
      if (TocOne(pOp+LenOp*(i-1)+oLabel) == NaN) k = i
    end do
    iDisk = TocOne(pNext)
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
    do i=0,TocOne(psym)-1
      nA = TocOne(pAsh+i)+nA
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
    if (TocOne(pldisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
    Length = TocOne(pldisp)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (label == 'TWOGRAD') then
    Comp = 1
    SymLab = 1
    if (TocOne(pldisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
    Length = TocOne(pldisp)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if ((label == 'STATHESS') .or. (label == 'RESPHESS') .or. (label == 'CONNHESS') .or. (label == 'HESS')) then
    Comp = 1
    SymLab = 1
    if (TocOne(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
    Length = 0
    do iSym=0,TocOne(pSym)-1
      Length = Length+nTri_Elem(TocOne(pldisp+isym))
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (label == 'INACTIVE') then
    Length = 0
    do iS=1,TocOne(pSym)
      do jS=1,TocOne(pSym)
        ijS = Mul(iS,jS)-1
        if (btest(iSymLab,ijS)) then
          jBas = TocOne(pbas+jS-1)
          iBas = TocOne(pbas+iS-1)
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
    do iS=1,TocOne(pSym)
      do jS=1,TocOne(pSym)
        ijS = Mul(iS,jS)-1
        if (btest(iSymLab,ijS)) then
          jBas = TocOne(pbas+jS-1)
          iBas = TocOne(pbas+iS-1)
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
    do i=1,TocOne(pSym)
      do j=1,i
        ij = Mul(i,j)-1
        if (btest(iSymLab,ij)) then
          if (i == j) then
            Length = Length+nTri_Elem(TocOne(pBas-1+i))
          else
            Length = Length+TocOne(pBas-1+i)*TocOne(pBas-1+j)
          end if
        end if
      end do
    end do
    if (iand(Option,slength) /= 0) Length = Len_
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  Length = rtoi*Length
  TocOne(pOp+LenOp*(k-1)+oLabel) = LabTmp(1)
# ifndef _I8_
  TocOne(pOp+LenOp*(k-1)+oLabel+1) = LabTmp(2)
# endif
  TocOne(pOp+LenOp*(k-1)+oComp) = Comp
  TocOne(pOp+LenOp*(k-1)+oSymLb) = iSymLab
  TocOne(pOp+LenOp*(k-1)+oAddr) = iDisk
  !write(u6,*) Length,idisk,nauxdt
  call iDAFile(LuMCK,1,iData,Length+nAuxDt,iDisk)
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
  !TocOne(pNext) = iDisk
  TocOne(pNext) = max(TocOne(pNext),iDisk)
end if
!----------------------------------------------------------------------*
! Finally copy the TocOne back to disk                                 *
!----------------------------------------------------------------------*
iDisk = 0
call iDaFile(LuMCK,1,TocOne,lToc,iDisk)
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
