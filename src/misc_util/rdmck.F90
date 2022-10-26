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

subroutine RdMCK(rc,Option,InLab,iComp,data,iSymLab)
!***********************************************************************
!                                                                      *
!     Purpose: Read data from one-electron integral file               *
!                                                                      *
!     Calling parameters:                                              *
!     rc      : Return code                                            *
!     Option  : Switch to set options                                  *
!     InLab   : Identifier for the data to read                        *
!     Comp    : Composite identifier to select components              *
!     Data    : contains on output the requested data                  *
!     SymLab  : symmetry label of the requested data                   *
!                                                                      *
!     Global data declarations (Include file):                         *
!     Parm    : Table of paramaters                                    *
!     rcParm  : Table of return codes                                  *
!     Switch  : Table of options                                       *
!     Common  : Common block containing ToC                            *
!     Data    : Data definitions                                       *
!                                                                      *
!     Local data declarations:                                         *
!     Label   : character*8, used to covert incoming names             *
!     TmpBuf  : I/O buffer                                             *
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

implicit integer(A-Z)
#include "MckDat.fh"
character*(*) InLab
dimension data(*)
logical :: Debug = .false.
character(LEN=8) TmpLab, Label
dimension TmpBuf(nBuf), HldBuf(1)
character(LEN=16) :: TheName = 'RdMck'
integer :: CurrOp = 1
! Statement function
MulTab(i,j) = ieor(i-1,j-1)+1

!----------------------------------------------------------------------*
! Start procedure:                                                     *
! Pick up the file definitions                                         *
!----------------------------------------------------------------------*
icpi = itob
Len_ = rc
rc = rc0000
LuMck = AuxMck(pLu)
open = AuxMck(pOpen)
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
if (iand(option,1024) /= 0) debug = .true.
if (Debug) then
  write(6,*) '<<< Entering RdMck  >>>'
  write(6,'(a,z8)') ' rc on entry:     ',rc
  write(6,'(a,a)') ' Label on entry:  ',Label
  write(6,'(a,z8)') ' Comp on entry:   ',Comp
  write(6,'(a,z8)') ' SymLab on entry: ',SymLab
  write(6,'(a,z8)') ' Option on entry: ',Option
end if
!----------------------------------------------------------------------*
! Check reading mode                                                   *
!----------------------------------------------------------------------*
if (iand(iand(option,sRdFst),sRdNxt) /= 0) then
  call SysWarnMsg(TheName,'Invalid value',' ')
  call SysCondMsg('iAnd(iAnd(option,sRdFst),sRdNxt) /= 0',iand(iand(option,sRdFst),sRdNxt),'/=',0)
else if (iand(iand(option,sRdFst),sRdCur) /= 0) then
  call SysWarnMsg(TheName,'Invalid value',' ')
  call SysCondMsg('iAnd(iAnd(option,sRdFst),sRdCur) /= 0',iand(iand(option,sRdFst),sRdCur),'/=',0)
else if (iand(iand(option,sRdNxt),sRdCur) /= 0) then
  call SysWarnMsg(TheName,'Invalid value',' ')
  call SysCondMsg('iAnd(iAnd(option,sRdNxt),sRdCur) /= 0',iand(iand(option,sRdNxt),sRdCur),'/=',0)
end if
!----------------------------------------------------------------------*
! Read data from ToC                                                   *
!----------------------------------------------------------------------*
NoGo = sRdFst+sRdNxt+sRdCur
if ((Label == 'TITLE') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pTitle) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    call iCopy(nTitle,TocOne(pTitle+1),1,data(1),1)
    if (debug) then
      write(6,'(a,z8)') ' Reading Title:'
      write(6,'(8(1x,z8))') (data(k),k=1,nTitle)
    end if
  else
    data(1) = nTitle
    if (debug) then
      write(6,'(a,z8)') ' Reading Title:',data(1)
    end if
  end if
else if ((Label == 'CHDISP') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pCHDISP) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    Length = TocOne(pnDisp)*30/icpi+1
    call iCopy(Length,TocOne(pchdisp),1,data(1),1)
    if (debug) then
      write(6,'(a,z8)') ' Reading perturbations:'
      write(6,'(8(1x,z8))') (data(k),k=1,nTitle)
    end if
  else
    data(1) = TocOne(pnDisp)*30/icpi+1
    if (debug) then
      write(6,'(a,z8)') ' Reading perturbations:',data(1)
    end if
  end if
else if ((Label == 'NDISP') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    data(1) = TocOne(pndisp)
  else
    data(1) = 1
  end if
  if (debug) then
    write(6,'(a,z8)') ' Reading nSym: ',data(1)
  end if
else if ((Label == 'NSYM') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    data(1) = TocOne(pSym)
  else
    data(1) = 1
  end if
  if (debug) then
    write(6,'(a,z8)') ' Reading nSym: ',data(1)
  end if
else if ((Label == 'NBAS') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pBas) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    Length = TocOne(pSym)
    call iCopy(Length,TocOne(pBas),1,data,1)
    if (debug) then
      write(6,'(a,8z8)') ' Reading nBas: ',(data(k),k=1,Length)
    end if
  else
    data(1) = TocOne(pSym)
    if (debug) then
      write(6,'(a,z8)') ' Reading nBas: ',data(1)
    end if
  end if
else if ((Label == 'LDISP') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pldisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    Length = TocOne(psym)
    call iCopy(Length,TocOne(pldisp),1,data,1)
    if (debug) then
      write(6,'(a,8z8)') ' Reading ldisp: ',(data(k),k=1,Length)
    end if
  else
    data(1) = TocOne(pSym)
    if (debug) then
      write(6,'(a,z8)') ' Reading ldisp: ',data(1)
    end if
  end if
else if ((Label == 'TDISP') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(ptdisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    Length = TocOne(pndisp)
    call iCopy(Length,TocOne(ptdisp),1,data,1)
    if (debug) then
      write(6,'(a,8z8)') ' Reading nBas: ',(data(k),k=1,Length)
    end if
  else
    data(1) = TocOne(pSym)
    if (debug) then
      write(6,'(a,z8)') ' Reading nBas: ',data(1)
    end if
  end if
else if ((Label == 'NASH') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pASH) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    Length = TocOne(pSYM)
    call iCopy(Length,TocOne(pASH),1,data,1)
    if (debug) then
      write(6,'(a,8z8)') ' Reading nASH: ',(data(k),k=1,Length)
    end if
  else
    data(1) = TocOne(pSym)
    if (debug) then
      write(6,'(a,z8)') ' Reading nASH: ',data(1)
    end if
  end if
else if (label == 'PERT') then
  call icopy(16/icpi,TocOne(pPert),1,data,1)
else if (label == 'NRCTDISP') then
  if (TocOne(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocOne(pndisp)
  call iCOPY(Length,TocOne(pnrdisp),1,data,1)
else if ((Label == 'SYMOP') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    !Length = (3*TocOne(pSym)-1)/icpi+1
    Length = (3*TocOne(pSym)+ItoB-1)/ItoB
    call iCopy(Length,TocOne(pSymOp),1,data,1)
    if (debug) then
      write(6,'(a)') ' Reading symmetry operators:'
      write(6,'(8(1x,z8))') (data(k),k=1,Length)
    end if
  else
    data(1) = TocOne(pSym)
    if (debug) then
      write(6,'(a,z8)') ' Reading symmetry operators:',data(1)
    end if
  end if
else if (label == 'DEGDISP ') then
  if (TocOne(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocOne(pndisp)
  call iCOPY(Length,TocOne(pdegdisp),1,data,1)
else
!----------------------------------------------------------------------*
! Read operators from integral records                                 *
!----------------------------------------------------------------------*
  if (iand(option,sRdNxt) /= 0) then
    if (debug) then
      write(6,'(a)') ' Reading next item'
    end if
    CurrOp = CurrOp+1
    if (CurrOp > MxOp) then
      CurrOp = 0
    else if (TocOne(pOp+LenOp*(CurrOp-1)+oLabel) == Nan) then
      CurrOp = 0
    else
      i = CurrOp
      !LabTmp(1) = TocOne(pOp+LenOp*(i-1)+oLabel)
      !#ifndef _I8_
      !LabTmp(2) = TocOne(pOp+LenOp*(i-1)+oLabel+1)
      !#endif
      idx = pOp+LenOp*(i-1)+oLabel
      TmpLab = transfer(TocOne(idx:idx+iLen-1),TmpLab)
      Label = TmpLab
      InLab = Label
      SymLab = TocOne(pOp+LenOp*(i-1)+oSymLb)
      Comp = TocOne(pOp+LenOp*(i-1)+oComp)
    end if
  else if (iand(option,sRdFst) /= 0) then
    if (debug) then
      write(6,'(a)') ' Reading first item'
    end if
    CurrOp = 1
    if (TocOne(pOp+LenOp*(CurrOp-1)+oLabel) == Nan) then
      CurrOp = 0
    else
      i = CurrOp
      !LabTmp(1) = TocOne(pOp+LenOp*(i-1)+oLabel)
      !#ifndef _I8_
      !LabTmp(2) = TocOne(pOp+LenOp*(i-1)+oLabel+1)
      !#endif
      idx = pOp+LenOp*(i-1)+oLabel
      TmpLab = transfer(TocOne(idx:idx+iLen-1),TmpLab)
      Label = TmpLab
      InLab = Label
      SymLab = TocOne(pOp+LenOp*(i-1)+oSymLb)
      Comp = TocOne(pOp+LenOp*(i-1)+oComp)
    end if
  else if (iand(option,sRdCur) /= 0) then
    if (debug) then
      write(6,'(a)') ' Reading current item'
    end if
    if ((CurrOp < 1) .or. (CurrOp > MxOp)) then
      CurrOp = 0
    else if (TocOne(pOp+LenOp*(CurrOp-1)+oLabel) == Nan) then
      CurrOp = 0
    else
      i = CurrOp
      !LabTmp(1) = TocOne(pOp+LenOp*(i-1)+oLabel)
      !#ifndef _I8_
      !LabTmp(2) = TocOne(pOp+LenOp*(i-1)+oLabel+1)
      !#endif
      idx = pOp+LenOp*(i-1)+oLabel
      TmpLab = transfer(TocOne(idx:idx+iLen-1),TmpLab)
      Label = TmpLab
      InLab = Label
      SymLab = TocOne(pOp+LenOp*(i-1)+oSymLb)
      Comp = TocOne(pOp+LenOp*(i-1)+oComp)
    end if
  else
    CurrOp = 0
    do i=MxOp,1,-1
      !LabTmp(1) = TocOne(pOp+LenOp*(i-1)+oLabel)
      !#ifndef _I8_
      !LabTmp(2) = TocOne(pOp+LenOp*(i-1)+oLabel+1)
      !#endif
      idx = pOp+LenOp*(i-1)+oLabel
      TmpLab = transfer(TocOne(idx:idx+iLen-1),TmpLab)
      CmpTmp = TocOne(pOp+LenOp*(i-1)+oComp)
      TmpCmp = Comp
      if ((TmpLab == Label) .and. (CmpTmp == TmpCmp)) CurrOp = i
    end do
  end if
  if (CurrOp == 0) call SysAbendMsg(TheName,'Current Operation == 0',' ')
  SymLab = TocOne(pOp+LenOp*(CurrOp-1)+oSymLb)
  if (Label == 'MOPERT') then
    iTmp = SymLab
    do iIrr=1,TocOne(pSym)
      iTmp = iTmp/2
    end do
    na = 0
    do i=0,TocOne(psym)-1
      na = TocOne(pAsh+i)+na
    end do
    Length = na*(na+1)/2
    Length = Length*(Length+1)/2
  else if (label == 'NUCGRAD') then
    Comp = 1
    SymLab = 1
    if (TocOne(pldisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
    Length = TocOne(pldisp)
  else if (label == 'TWOGRAD') then
    Comp = 1
    SymLab = 1
    if (TocOne(pldisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
    Length = TocOne(pldisp)
  else if ((label == 'STATHESS') .or. (label == 'RESPHESS') .or. (label == 'CONNHESS') .or. (label == 'HESS')) then
    Comp = 1
    SymLab = 1
    if (TocOne(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
    Length = 0
    do iSym=0,TocOne(pSym)-1
      Length = Length+TocOne(pldisp+isym)*(TocOne(pldisp+isym)+1)/2
    end do
  else if (label == 'INACTIVE') then
    Length = 0
    do iS=1,TocOne(pSym)
      do jS=1,TocOne(pSym)
        ijS = MulTab(iS,jS)
        if (iand(2**(ijS-1),SymLab) /= 0) then
          jBas = TocOne(pbas+jS-1)
          iBas = TocOne(pbas+iS-1)
          if (jBas == NaN) call SysAbendMsg(TheName,'jBas == NaN at label',Label)
          if (iBas == NaN) call SysAbendMsg(TheName,'iBas == NaN at label',Label)
          Length = Length+iBas*jBas
        end if
      end do
    end do
  else if (label == 'TOTAL') then
    Length = 0
    do iS=1,TocOne(pSym)
      do jS=1,TocOne(pSym)
        ijS = MulTab(iS,jS)
        if (iand(2**(ijS-1),SymLab) /= 0) then
          jBas = TocOne(pbas+jS-1)
          iBas = TocOne(pbas+iS-1)
          if (jBas == NaN) call SysAbendMsg(TheName,'jBas == NaN at label',Label)
          if (iBas == NaN) call SysAbendMsg(TheName,'iBas == NaN at label',Label)
          Length = Length+iBas*jBas
        end if
      end do
    end do
  else
    Length = 0
    do i=1,TocOne(pSym)
      do j=1,i
        ij = MulTab(i,j)-1
        if (iand(2**ij,SymLab) /= 0) then
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

  data(1) = Length
  if (iand(option,sOpSiz) == 0) then
    Length = rtoi*Length
    IndDta = 1
    IndHld = 1
    iDisk = TocOne(pOp+LenOp*(CurrOp-1)+oAddr)
    do k=0,Length+nAuxDt-1,nBuf
      iBuf = max(0,min(nBuf,Length+nAuxDt-k))
      tBuf = max(0,min(nBuf,Length-k))
      eBuf = iBuf-tBuf
      call iDaFile(LuMCK,2,TmpBuf,iBuf,iDisk)
      IndTmp = 1
      call iCopy(tBuf,TmpBuf(IndTmp),1,data(IndDta),1)
      if (debug) then
        write(6,'(a,z8)') ' Reading buffer to: ',IndDta
        write(6,'(8(1x,z8))') (data(IndDta+m),m=0,tBuf-1)
      end if
      IndTmp = IndTmp+tBuf
      IndDta = IndDta+tBuf
      if (tBuf < iBuf) then
        call iCopy(eBuf,TmpBuf(IndTmp),1,HldBuf(IndHld),1)
        IndTmp = IndTmp+eBuf
        IndHld = IndHld+eBuf
      end if
    end do
    !if (iAnd(sNoOri,option) == 0) then
    !  call iCopy(6,HldBuf(1),1,Data(IndDta),1)
    !  if (debug) then
    !    write(6,'(a,z8)') ' Reading buffer to: ',IndDta
    !    write(6,'(8(1x,z8))') (Data(IndDta+m),m=0,5)
    !  end if
    !end if
    IndDta = IndDta+6
    !if (iAnd(sNoNuc,option) == 0) then
    !  call iCopy(2,HldBuf(7),1,Data(IndDta),1)
    !  if (debug) then
    !    write(6,'(a,z8)') ' Reading buffer to: ',IndDta
    !    write(6,'(8(1x,z8))') (Data(IndDta+m),m=0,1)
    !  end if
    !end if
    !IndDta = IndDta+2
  end if
end if
!----------------------------------------------------------------------*
! Print debugging information                                          *
!----------------------------------------------------------------------*
if (Debug) then
  write(6,*) '<<< Exiting RdMck >>>'
  write(6,'(a,z8)') ' rc on exit:     ',rc
  write(6,'(a,a)') ' Label on exit:  ',Label
  write(6,'(a,z8)') ' Comp on exit:    ',Comp
  write(6,'(a,z8)') ' SymLab on exit: ',SymLab
  write(6,'(a,z8)') ' Option on exit: ',Option
end if

!----------------------------------------------------------------------*
! Terminate procedure                                                  *
!----------------------------------------------------------------------*
return

end subroutine RdMCK
