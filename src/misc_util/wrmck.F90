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

subroutine WrMCK(rc,Option,InLab,iComp,data,iSymLab)
!***********************************************************************
!                                                                      *
!     Purpose: write data to one-electron integral file                *
!                                                                      *
!     Calling parameters:                                              *
!     rc      : Return code                                            *
!     Option  : Switch to set options                                  *
!     InLab   : Identifier for the data to write                       *
!     Comp    : Composite identifier to select components              *
!     Data    : contains on input the data to store on disk            *
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

implicit integer(A-Z)
#include "MckDat.fh"
character*(*) InLab
dimension data(*)
!real*8 DDot_, Check
!external ddot_
logical Debug
!character*8 Label, Label_Add*11
character*8 Label
dimension LabTmp(2)
character*16 TheName
data TheName/'WrMck'/
data Debug/.false./
! Statement function
MulTab(i,j) = ieor(i-1,j-1)+1

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
open = AuxMCK(pOpen)
!----------------------------------------------------------------------*
! Check the file status                                                *
!----------------------------------------------------------------------*
if (open /= 1) call SysFileMsg(TheName,'MSG: open',LuMck,' ')
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
if (iand(option,1024) /= 0) debug = .true.
if (Debug) then
  write(6,*) '<<< Entering WrMck >>>'
  write(6,'(a,z8)') ' rc on entry:     ',rc
  write(6,'(a,a)') ' Label on entry:  ',Label
  write(6,'(a,z8)') ' Comp on entry:   ',Comp
  write(6,'(a,z8)') ' SymLab on entry: ',SymLab
  write(6,'(a,z8)') ' Option on entry: ',Option
  write(6,*) ' Contents of the Toc'
  write(6,*) ' ==================='
  write(6,'(i6,z8)') 'pFID,TocOne(pFID)=',pFID,TocOne(pFID)
  write(6,'(i6,z8)') pVersN,TocOne(pVersN)
  write(6,'(i6,z8)') pTitle,TocOne(pTitle)
  write(6,'(i6,z8)') pOp,TocOne(pOp)
  write(6,'(i6,z8)') pSym,TocOne(pSym)
  write(6,'(i6,z8)') pSymOp,TocOne(pSymOp)
  write(6,'(i6,z8)') pBas,TocOne(pBas)
  write(6,'(i6,z8)') pNext,TocOne(pNext)
  write(6,'(i6,z8)') pEnd,TocOne(pEnd)
end if
!----------------------------------------------------------------------*
! Store data in TocOne                                                 *
!----------------------------------------------------------------------*

if (Label == 'TITLE') then

  TocOne(pTitle) = NotNaN
  call iCopy(nTitle,data(1),1,TocOne(pTitle+1),1)
  !--------------------------------------------------------------------*

else if (Label == 'NSYM') then

  if ((data(1) > MxSym) .or. (data(1) < 1)) then
    call SysWarnMsg(TheName,'Label=',Label)
    call SysValueMsg('Data(1)=',data(1))
  end if
  TocOne(pSym) = data(1)
  !--------------------------------------------------------------------*

else if (Label == 'NBAS') then

  if (TocOne(pSym) == NaN) then
    call SysAbendMsg(TheName,'Undefined Label:',Label)
  end if
  Length = TocOne(pSym)
  call iCOPY(Length,data,1,TocOne(pbas),1)
  !--------------------------------------------------------------------*

else if (label == 'NISH') then

  if (TocOne(pSym) == NaN) then
    call SysAbendMsg(TheName,'Undefined Label:',Label)
  end if
  Length = TocOne(pSym)
  call iCOPY(Length,data,1,TocOne(pISH),1)
  !--------------------------------------------------------------------*

else if (label == 'NASH') then

  if (TocOne(pSym) == NaN) then
    call SysAbendMsg(TheName,'Undefined Label:',Label)
  end if
  Length = TocOne(pSym)
  call iCOPY(Length,data,1,TocOne(pASH),1)
  !--------------------------------------------------------------------*

else if (label == 'LDISP') then

  if (TocOne(pSym) == NaN) then
    call SysAbendMsg(TheName,'Undefined Label:',Label)
  end if
  Length = TocOne(pSym)
  call iCOPY(Length,data,1,TocOne(pldisp),1)
  !--------------------------------------------------------------------*

else if (label == 'TDISP') then

  if (TocOne(pndisp) == NaN) then
    call SysAbendMsg(TheName,'Undefined Label:',Label)
  end if
  Length = TocOne(pndisp)
  call iCOPY(Length,data,1,TocOne(ptdisp),1)
  !--------------------------------------------------------------------*

else if (label == 'NDISP') then

  ToCOne(pnDisp) = data(1)
  !--------------------------------------------------------------------*

else if (label == 'CHDISP') then

  if (TocOne(pndisp) == NaN) then
    call SysAbendMsg(TheName,'Undefined Label:',Label)
  end if
  Length = TocOne(pndisp)*30/icpi+1
  call iCOPY(Length,data,1,TocOne(pchdisp),1)
  !--------------------------------------------------------------------*

else if (label == 'NRCTDISP') then

  if (TocOne(pndisp) == NaN) then
    call SysAbendMsg(TheName,'Undefined Label:',Label)
  end if
  Length = TocOne(pndisp)
  call iCOPY(Length,data,1,TocOne(pnrdisp),1)
  !--------------------------------------------------------------------*

else if (label == 'DEGDISP ') then

  if (TocOne(pndisp) == NaN) then
    call SysAbendMsg(TheName,'Undefined Label:',Label)
  end if
  Length = TocOne(pndisp)
  call iCOPY(Length,data,1,TocOne(pdegdisp),1)
  !--------------------------------------------------------------------*

else if (Label == 'SYMOP') then

  if (TocOne(pSym) == NaN) then
    call SysAbendMsg(TheName,'Undefined Label:',Label)
  end if
  Length = (3*TocOne(pSym)+ItoB-1)/ItoB
  call iCopy(Length,data,1,TocOne(pSymOp),1)
  !--------------------------------------------------------------------*

else if (label == 'PERT') then

  call icopy(16/icpi,data,1,TocOne(pPert),1)
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

  if (TocOne(pBas) == NaN) then
    call SysAbendMsg(TheName,'Undefined Label:',Label)
  end if

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
      write(6,*) ' This is a new field!'
      write(6,*) ' iDisk=',iDisk
      write(6,*) ' FldNo=',k
      write(6,*) ' pNext=',pNext
    end if
  else
    if (Debug) then
      write(6,*) ' This is an old field!'
      write(6,*) ' iDisk=',iDisk
      write(6,*) ' FldNo=',k
      write(6,*) ' pNext=',pNext
    end if
  end if
  if (k == 0) call SysAbendMsg(TheName,'Undefined Label:',Label)
  !write(6,*) isymlab,label
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (Label == 'MOPERT') then

    nA = 0
    do i=0,TocOne(psym)-1
      nA = TocOne(pAsh+i)+nA
    end do
    Length = nA*(na+1)/2
    Length = Length*(Length+1)/2
    !write(6,*) Length
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (label == 'NUCGRAD') then
    Comp = 1
    SymLab = 1
    if (TocOne(pldisp) == NaN) then
      call SysAbendMsg(TheName,'Undefined Label:',Label)
    end if
    Length = TocOne(pldisp)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (label == 'TWOGRAD') then
    Comp = 1
    SymLab = 1
    if (TocOne(pldisp) == NaN) then
      call SysAbendMsg(TheName,'Undefined Label:',Label)
    end if
    Length = TocOne(pldisp)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if ((label == 'STATHESS') .or. (label == 'RESPHESS') .or. (label == 'CONNHESS') .or. (label == 'HESS')) then
    Comp = 1
    SymLab = 1
    if (TocOne(pndisp) == NaN) then
      call SysAbendMsg(TheName,'Undefined Label:',Label)
    end if
    Length = 0
    do iSym=0,TocOne(pSym)-1
      Length = Length+TocOne(pldisp+isym)*(TocOne(pldisp+isym)+1)/2
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (label == 'INACTIVE') then
    Length = 0
    do iS=1,TocOne(pSym)
      do jS=1,TocOne(pSym)
        ijS = MulTab(iS,jS)
        if (iand(2**(ijS-1),iSymLab) /= 0) then
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
        ijS = MulTab(iS,jS)
        if (iand(2**(ijS-1),iSymLab) /= 0) then
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
        ij = MulTab(i,j)-1
        if (iand(2**ij,iSymLab) /= 0) then
          if (i == j) then
            Length = Length+TocOne(pBas-1+i)*(TocOne(pBas-1+i)+1)/2
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
  !write(6,*) Length,idisk,nauxdt
  call iDAFile(LuMCK,1,data,Length+nAuxDt,iDisk)
  !if ((Label == 'TOTAL') .or. (Label == 'INACTIVE') .or. (Label == 'MOPERT')) then
  !  !write(6,*) 'iComp=',iComp
  !  !call RecPrt(Label,' ',Data,1,Length/RtoI)
  !  Check = DDot_(Length/RtoI,Data,1,Data,1)
  !  Label_Add = ' '
  !  Label_Add = Label
  !  write(Label_Add(9:11),'(I3.3)') iComp
  !  !write(6,*) Label_Add,Label,iComp,Check(1)
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
  write(6,*) '<<< Exiting WrMck >>>'
  write(6,'(a,z8)') ' rc on exit:     ',rc
  write(6,'(a,a)') ' Label on exit:  ',Label
  write(6,'(a,z8)') ' Comp on exit:   ',Comp
  write(6,'(a,z8)') ' SymLab on exit: ',SymLab
  write(6,'(a,z8)') ' Option on exit: ',Option
end if

!----------------------------------------------------------------------*
! Terminate procedure                                                  *
!----------------------------------------------------------------------*
return

end subroutine WrMCK
