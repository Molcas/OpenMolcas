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

subroutine RdMCK(rc,Option,InLab,iComp,iData,iSymLab)
!***********************************************************************
!                                                                      *
!     Purpose: Read data from one-electron integral file               *
!                                                                      *
!     Calling parameters:                                              *
!     rc      : Return code                                            *
!     Option  : Switch to set options                                  *
!     InLab   : Identifier for the data to read                        *
!     Comp    : Composite identifier to select components              *
!     iData   : contains on output the requested data                  *
!     SymLab  : symmetry label of the requested data                   *
!                                                                      *
!     Local data declarations:                                         *
!     Label   : character*8, used to covert incoming names             *
!     HldBuf  : I/O buffer                                             *
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
use MckDat, only: AuxMck, LenOp, MxOp, NaN, nBuf, nTitle, oAddr, oComp, oLabel, oSymLb, pASh, pBas, pChdisp, pdegdisp, pldisp, &
                  pndisp, pnrdisp, pOp, pPert, pSym, pSymOp, ptdisp, pTitle, rcMck, sLength, sDbg, sOpSiz, sRdCur, sRdFst, sRdNxt, &
                  TmpBuf, TocMck
use Definitions, only: iwp, u6, RtoI, ItoB

#include "intent.fh"

implicit none
integer(kind=iwp), intent(inout) :: rc
integer(kind=iwp), intent(in) :: Option, iComp, iSymLab
character(len=*), intent(inout) :: InLab
integer(kind=iwp), intent(_OUT_) :: iData(*)
integer(kind=iwp) :: CmpTmp, Comp, CurrOp = 1, i, iBas, icpi, iDisk, idx, iIrr, ij, ijS, iLen, IndDta, iS, iSym, iTmp, j, jBas, &
                     jS, k, Len_, Length, LuMck, na, SymLab, tBuf, TmpCmp
logical(kind=iwp) :: Debug, NoGo, NoOpSiz
character(len=16), parameter :: TheName = 'RdMck'
character(len=8) :: TmpLab, Label

!----------------------------------------------------------------------*
! Start procedure:                                                     *
! Pick up the file definitions                                         *
!----------------------------------------------------------------------*
icpi = ItoB
Len_ = rc
rc = rcMck%good
LuMck = AuxMck%Lu
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
TmpLab = Label
iLen = len(TmpLab)/ItoB
Comp = icomp
Symlab = isymlab
if ((label == 'STATHESS') .or. (label == 'RESPHESS') .or. (label == 'CONNHESS') .or. (label == 'HESS') .or. &
    (label == 'NUCGRAD') .or. (label == 'TWOGRAD')) then
  Comp = 1
  Symlab = 1
end if
!----------------------------------------------------------------------*
! Print debugging information                                          *
!----------------------------------------------------------------------*
Debug = btest(option,sDbg)
if (Debug) then
  write(u6,*) '<<< Entering RdMck  >>>'
  write(u6,'(a,z8)') ' rc on entry:     ',rc
  write(u6,'(a,a)') ' Label on entry:  ',Label
  write(u6,'(a,z8)') ' Comp on entry:   ',Comp
  write(u6,'(a,z8)') ' SymLab on entry: ',SymLab
  write(u6,'(a,z8)') ' Option on entry: ',Option
end if
!----------------------------------------------------------------------*
! Check reading mode                                                   *
!----------------------------------------------------------------------*
if (btest(option,sRdFst) .and. btest(option,sRdNxt)) then
  call SysWarnMsg(TheName,'Invalid value','sRdFst and sRdNxt')
else if (btest(option,sRdFst) .and. btest(option,sRdCur)) then
  call SysWarnMsg(TheName,'Invalid value','sRdFst and sRdNxt')
else if (btest(option,sRdNxt) .and. btest(option,sRdCur)) then
  call SysWarnMsg(TheName,'Invalid value','sRdNxt and sRdCur')
end if
!----------------------------------------------------------------------*
! Read data from ToC                                                   *
!----------------------------------------------------------------------*
NoOpSiz = .not. btest(option,sOpSiz)
NoGo = .not. (btest(option,sRdFst) .or. btest(option,sRdNxt) .or. btest(option,sRdCur))
if ((Label == 'TITLE') .and. NoGo) then
  if (TocMck(pTitle) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (NoOpSiz) then
    iData(1:nTitle) = TocMck(pTitle+1:pTitle+nTitle)
    if (debug) then
      write(u6,'(a,z8)') ' Reading Title:'
      write(u6,'(8(1x,z8))') iData(1:nTitle)
    end if
  else
    iData(1) = nTitle
    if (debug) write(u6,'(a,z8)') ' Reading Title:',iData(1)
  end if
else if ((Label == 'CHDISP') .and. NoGo) then
  if (TocMck(pCHDISP) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (NoOpSiz) then
    Length = TocMck(pnDisp)*30/icpi+1
    iData(1:Length) = TocMck(pchdisp:pchdisp+Length-1)
    if (debug) then
      write(u6,'(a,z8)') ' Reading perturbations:'
      write(u6,'(8(1x,z8))') iData(1:Length)
    end if
  else
    iData(1) = TocMck(pnDisp)*30/icpi+1
    if (debug) write(u6,'(a,z8)') ' Reading perturbations:',iData(1)
  end if
else if ((Label == 'NDISP') .and. NoGo) then
  if (TocMck(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (NoOpSiz) then
    iData(1) = TocMck(pndisp)
  else
    iData(1) = 1
  end if
  if (debug) write(u6,'(a,z8)') ' Reading nSym: ',iData(1)
else if ((Label == 'NSYM') .and. NoGo) then
  if (TocMck(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (NoOpSiz) then
    iData(1) = TocMck(pSym)
  else
    iData(1) = 1
  end if
  if (debug) write(u6,'(a,z8)') ' Reading nSym: ',iData(1)
else if ((Label == 'NBAS') .and. NoGo) then
  if (TocMck(pBas) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (NoOpSiz) then
    Length = TocMck(pSym)
    iData(1:Length) = TocMck(pBas:pBas+Length-1)
    if (debug) write(u6,'(a,8z8)') ' Reading nBas: ',iData(1:Length)
  else
    iData(1) = TocMck(pSym)
    if (debug) write(u6,'(a,z8)') ' Reading nBas: ',iData(1)
  end if
else if ((Label == 'LDISP') .and. NoGo) then
  if (TocMck(pldisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (NoOpSiz) then
    Length = TocMck(psym)
    iData(1:Length) = TocMck(pldisp:pldisp+Length-1)
    if (debug) write(u6,'(a,8z8)') ' Reading ldisp: ',iData(1:Length)
  else
    iData(1) = TocMck(pSym)
    if (debug) write(u6,'(a,z8)') ' Reading ldisp: ',iData(1)
  end if
else if ((Label == 'TDISP') .and. NoGo) then
  if (TocMck(ptdisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (NoOpSiz) then
    Length = TocMck(pndisp)
    iData(1:Length) = TocMck(ptdisp:ptdisp+Length-1)
    if (debug) write(u6,'(a,8z8)') ' Reading nBas: ',iData(1:Length)
  else
    iData(1) = TocMck(pSym)
    if (debug) write(u6,'(a,z8)') ' Reading nBas: ',iData(1)
  end if
else if ((Label == 'NASH') .and. NoGo) then
  if (TocMck(pASH) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (NoOpSiz) then
    Length = TocMck(pSYM)
    iData(1:Length) = TocMck(pASH:pASH+length-1)
    if (debug) write(u6,'(a,8z8)') ' Reading nASH: ',iData(1:Length)
  else
    iData(1) = TocMck(pSym)
    if (debug) write(u6,'(a,z8)') ' Reading nASH: ',iData(1)
  end if
else if (label == 'PERT') then
  Length = 16/icpi
  iData(1:Length) = TocMck(pPert:pPert+Length-1)
else if (label == 'NRCTDISP') then
  if (TocMck(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocMck(pndisp)
  iData(1:Length) = TocMck(pnrdisp:pnrdisp+Length-1)
else if ((Label == 'SYMOP') .and. NoGo) then
  if (TocMck(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (NoOpSiz) then
    !Length = (3*TocMck(pSym)-1)/icpi+1
    Length = (3*TocMck(pSym)+ItoB-1)/ItoB
    iData(1:Length) = TocMck(pSymOp:pSymOp+Length-1)
    if (debug) then
      write(u6,'(a)') ' Reading symmetry operators:'
      write(u6,'(8(1x,z8))') iData(1:Length)
    end if
  else
    iData(1) = TocMck(pSym)
    if (debug) write(u6,'(a,z8)') ' Reading symmetry operators:',iData(1)
  end if
else if (label == 'DEGDISP ') then
  if (TocMck(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocMck(pndisp)
  iData(1:Length) = TocMck(pdegdisp:pdegdisp+Length-1)
else
!----------------------------------------------------------------------*
! Read operators from integral records                                 *
!----------------------------------------------------------------------*
  if (btest(option,sRdNxt)) then
    if (debug) write(u6,'(a)') ' Reading next item'
    CurrOp = CurrOp+1
    if (CurrOp > MxOp) then
      CurrOp = 0
    else if (TocMck(pOp+LenOp*(CurrOp-1)+oLabel) == Nan) then
      CurrOp = 0
    else
      i = CurrOp
      !LabTmp(1) = TocMck(pOp+LenOp*(i-1)+oLabel)
      !#ifndef _I8_
      !LabTmp(2) = TocMck(pOp+LenOp*(i-1)+oLabel+1)
      !#endif
      idx = pOp+LenOp*(i-1)+oLabel
      TmpLab = transfer(TocMck(idx:idx+iLen-1),TmpLab)
      Label = TmpLab
      InLab = Label
      SymLab = TocMck(pOp+LenOp*(i-1)+oSymLb)
      Comp = TocMck(pOp+LenOp*(i-1)+oComp)
    end if
  else if (btest(option,sRdFst)) then
    if (debug) write(u6,'(a)') ' Reading first item'
    CurrOp = 1
    if (TocMck(pOp+LenOp*(CurrOp-1)+oLabel) == NaN) then
      CurrOp = 0
    else
      i = CurrOp
      !LabTmp(1) = TocMck(pOp+LenOp*(i-1)+oLabel)
      !#ifndef _I8_
      !LabTmp(2) = TocMck(pOp+LenOp*(i-1)+oLabel+1)
      !#endif
      idx = pOp+LenOp*(i-1)+oLabel
      TmpLab = transfer(TocMck(idx:idx+iLen-1),TmpLab)
      Label = TmpLab
      InLab = Label
      SymLab = TocMck(pOp+LenOp*(i-1)+oSymLb)
      Comp = TocMck(pOp+LenOp*(i-1)+oComp)
    end if
  else if (btest(option,sRdCur)) then
    if (debug) write(u6,'(a)') ' Reading current item'
    if ((CurrOp < 1) .or. (CurrOp > MxOp)) then
      CurrOp = 0
    else if (TocMck(pOp+LenOp*(CurrOp-1)+oLabel) == NaN) then
      CurrOp = 0
    else
      i = CurrOp
      !LabTmp(1) = TocMck(pOp+LenOp*(i-1)+oLabel)
      !#ifndef _I8_
      !LabTmp(2) = TocMck(pOp+LenOp*(i-1)+oLabel+1)
      !#endif
      idx = pOp+LenOp*(i-1)+oLabel
      TmpLab = transfer(TocMck(idx:idx+iLen-1),TmpLab)
      Label = TmpLab
      InLab = Label
      SymLab = TocMck(pOp+LenOp*(i-1)+oSymLb)
      Comp = TocMck(pOp+LenOp*(i-1)+oComp)
    end if
  else
    CurrOp = 0
    do i=MxOp,1,-1
      !LabTmp(1) = TocMck(pOp+LenOp*(i-1)+oLabel)
      !#ifndef _I8_
      !LabTmp(2) = TocMck(pOp+LenOp*(i-1)+oLabel+1)
      !#endif
      idx = pOp+LenOp*(i-1)+oLabel
      TmpLab = transfer(TocMck(idx:idx+iLen-1),TmpLab)
      CmpTmp = TocMck(pOp+LenOp*(i-1)+oComp)
      TmpCmp = Comp
      if ((TmpLab == Label) .and. (CmpTmp == TmpCmp)) CurrOp = i
    end do
  end if
  if (CurrOp == 0) call SysAbendMsg(TheName,'Current Operation == 0',' ')
  SymLab = TocMck(pOp+LenOp*(CurrOp-1)+oSymLb)
  if (Label == 'MOPERT') then
    iTmp = SymLab
    do iIrr=1,TocMck(pSym)
      iTmp = iTmp/2
    end do
    na = 0
    do i=0,TocMck(psym)-1
      na = TocMck(pAsh+i)+na
    end do
    Length = nTri_Elem(na)
    Length = nTri_Elem(Length)
  else if (label == 'NUCGRAD') then
    Comp = 1
    SymLab = 1
    if (TocMck(pldisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
    Length = TocMck(pldisp)
  else if (label == 'TWOGRAD') then
    Comp = 1
    SymLab = 1
    if (TocMck(pldisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
    Length = TocMck(pldisp)
  else if ((label == 'STATHESS') .or. (label == 'RESPHESS') .or. (label == 'CONNHESS') .or. (label == 'HESS')) then
    Comp = 1
    SymLab = 1
    if (TocMck(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
    Length = 0
    do iSym=0,TocMck(pSym)-1
      Length = Length+nTri_Elem(TocMck(pldisp+isym))
    end do
  else if (label == 'INACTIVE') then
    Length = 0
    do iS=1,TocMck(pSym)
      do jS=1,TocMck(pSym)
        ijS = Mul(iS,jS)-1
        if (btest(SymLab,ijS)) then
          jBas = TocMck(pbas+jS-1)
          iBas = TocMck(pbas+iS-1)
          if (jBas == NaN) call SysAbendMsg(TheName,'jBas == NaN at label',Label)
          if (iBas == NaN) call SysAbendMsg(TheName,'iBas == NaN at label',Label)
          Length = Length+iBas*jBas
        end if
      end do
    end do
  else if (label == 'TOTAL') then
    Length = 0
    do iS=1,TocMck(pSym)
      do jS=1,TocMck(pSym)
        ijS = Mul(iS,jS)-1
        if (btest(SymLab,ijS)) then
          jBas = TocMck(pbas+jS-1)
          iBas = TocMck(pbas+iS-1)
          if (jBas == NaN) call SysAbendMsg(TheName,'jBas == NaN at label',Label)
          if (iBas == NaN) call SysAbendMsg(TheName,'iBas == NaN at label',Label)
          Length = Length+iBas*jBas
        end if
      end do
    end do
  else
    Length = 0
    do i=1,TocMck(pSym)
      do j=1,i
        ij = Mul(i,j)-1
        if (btest(SymLab,ij)) then
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

  iData(1) = Length
  if (NoOpSiz) then
    Length = RtoI*Length
    IndDta = 1
    !IndHld = 1
    iDisk = TocMck(pOp+LenOp*(CurrOp-1)+oAddr)
    do k=0,Length-1,nBuf
      tBuf = max(0,min(nBuf,Length-k))
      call iDaFile(LuMCK,2,TmpBuf,tBuf,iDisk)
      iData(IndDta:IndDta+tBuf-1) = TmpBuf(1:tBuf)
      if (debug) then
        write(u6,'(a,z8)') ' Reading buffer to: ',IndDta
        write(u6,'(8(1x,z8))') iData(IndDta:IndDta+tBuf-1)
      end if
      IndDta = IndDta+tBuf
      !if (tBuf < iBuf) then
      !  eBuf = iBuf-tBuf
      !  !HldBuf(IndHld:IndHld+eBuf-1) = TmpBuf(IndTmp:IndTmp+eBuf-1)
      !  IndTmp = IndTmp+eBuf
      !  !IndHld = IndHld+eBuf
      !end if
    end do
    !if (.not. btest(option,sNoOri)) then
    !  iData(IndDta:IndDta+5) = HldBuf(1:6)
    !  if (debug) then
    !    write(u6,'(a,z8)') ' Reading buffer to: ',IndDta
    !    write(u6,'(8(1x,z8))') iData(IndDta:IndDta+5)
    !  end if
    !end if
    !IndDta = IndDta+6
    !if (.not. btest(option,sNoNuc)) then
    !  iData(IndDta:IndDta+1) = HldBuf(7:8)
    !  if (debug) then
    !    write(u6,'(a,z8)') ' Reading buffer to: ',IndDta
    !    write(u6,'(8(1x,z8))') iData(IndDta:IndDta+1)
    !  end if
    !end if
    !IndDta = IndDta+2
  end if
end if
!----------------------------------------------------------------------*
! Print debugging information                                          *
!----------------------------------------------------------------------*
if (Debug) then
  write(u6,*) '<<< Exiting RdMck >>>'
  write(u6,'(a,z8)') ' rc on exit:     ',rc
  write(u6,'(a,a)') ' Label on exit:  ',Label
  write(u6,'(a,z8)') ' Comp on exit:    ',Comp
  write(u6,'(a,z8)') ' SymLab on exit: ',SymLab
  write(u6,'(a,z8)') ' Option on exit: ',Option
end if

!----------------------------------------------------------------------*
! Terminate procedure                                                  *
!----------------------------------------------------------------------*
return

end subroutine RdMCK
