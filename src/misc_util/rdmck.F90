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

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: rc, Option, iComp, iData(*), iSymLab
character(len=*) :: InLab
#include "MckDat.fh"
integer(kind=iwp) :: CmpTmp, Comp, CurrOp = 1, eBuf, i, iBas, iBuf, icpi, iDisk, idx, iIrr, ij, ijS, iLen, IndDta, IndTmp, iS, &
                     isopen, iSym, iTmp, j, jBas, jS, k, Len_, Length, LuMck, m, na, NoGo, SymLab, tBuf, TmpBuf(nBuf), TmpCmp !IFG
logical(kind=iwp) :: Debug
character(len=16), parameter :: TheName = 'RdMck'
character(len=8) TmpLab, Label

!----------------------------------------------------------------------*
! Start procedure:                                                     *
! Pick up the file definitions                                         *
!----------------------------------------------------------------------*
icpi = itob
Len_ = rc
rc = rc0000
LuMck = AuxMck(pLu)
isopen = AuxMck(pOpen)
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
Debug = btest(option,10)
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
    iData(1:nTitle) = TocOne(pTitle+1:pTitle+nTitle)
    if (debug) then
      write(u6,'(a,z8)') ' Reading Title:'
      write(u6,'(8(1x,z8))') (iData(k),k=1,nTitle)
    end if
  else
    iData(1) = nTitle
    if (debug) write(u6,'(a,z8)') ' Reading Title:',iData(1)
  end if
else if ((Label == 'CHDISP') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pCHDISP) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    Length = TocOne(pnDisp)*30/icpi+1
    iData(1:Length) = TocOne(pchdisp:pchdisp+Length-1)
    if (debug) then
      write(u6,'(a,z8)') ' Reading perturbations:'
      write(u6,'(8(1x,z8))') (iData(k),k=1,nTitle)
    end if
  else
    iData(1) = TocOne(pnDisp)*30/icpi+1
    if (debug) write(u6,'(a,z8)') ' Reading perturbations:',iData(1)
  end if
else if ((Label == 'NDISP') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    iData(1) = TocOne(pndisp)
  else
    iData(1) = 1
  end if
  if (debug) write(u6,'(a,z8)') ' Reading nSym: ',iData(1)
else if ((Label == 'NSYM') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    iData(1) = TocOne(pSym)
  else
    iData(1) = 1
  end if
  if (debug) write(u6,'(a,z8)') ' Reading nSym: ',iData(1)
else if ((Label == 'NBAS') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pBas) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    Length = TocOne(pSym)
    iData(1:Length) = TocOne(pBas:pBas+Length-1)
    if (debug) write(u6,'(a,8z8)') ' Reading nBas: ',(iData(k),k=1,Length)
  else
    iData(1) = TocOne(pSym)
    if (debug) write(u6,'(a,z8)') ' Reading nBas: ',iData(1)
  end if
else if ((Label == 'LDISP') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pldisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    Length = TocOne(psym)
    iData(1:Length) = TocOne(pldisp:pldisp+Length-1)
    if (debug) write(u6,'(a,8z8)') ' Reading ldisp: ',(iData(k),k=1,Length)
  else
    iData(1) = TocOne(pSym)
    if (debug) write(u6,'(a,z8)') ' Reading ldisp: ',iData(1)
  end if
else if ((Label == 'TDISP') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(ptdisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    Length = TocOne(pndisp)
    iData(1:Length) = TocOne(ptdisp:ptdisp+Length-1)
    if (debug) write(u6,'(a,8z8)') ' Reading nBas: ',(iData(k),k=1,Length)
  else
    iData(1) = TocOne(pSym)
    if (debug) write(u6,'(a,z8)') ' Reading nBas: ',iData(1)
  end if
else if ((Label == 'NASH') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pASH) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    Length = TocOne(pSYM)
    iData(1:Length) = TocOne(pASH:pASH+length-1)
    if (debug) write(u6,'(a,8z8)') ' Reading nASH: ',(iData(k),k=1,Length)
  else
    iData(1) = TocOne(pSym)
    if (debug) write(u6,'(a,z8)') ' Reading nASH: ',iData(1)
  end if
else if (label == 'PERT') then
  Length = 16/icpi
  iData(1:Length) = TocOne(pPert:pPert+Length-1)
else if (label == 'NRCTDISP') then
  if (TocOne(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocOne(pndisp)
  iData(1:Length) = TocOne(pnrdisp:pnrdisp+Length-1)
else if ((Label == 'SYMOP') .and. (iand(option,NoGo) == 0)) then
  if (TocOne(pSym) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  if (iand(option,sOpSiz) == 0) then
    !Length = (3*TocOne(pSym)-1)/icpi+1
    Length = (3*TocOne(pSym)+ItoB-1)/ItoB
    iData(1:Length) = TocOne(pSymOp:pSymOp+Length-1)
    if (debug) then
      write(u6,'(a)') ' Reading symmetry operators:'
      write(u6,'(8(1x,z8))') (iData(k),k=1,Length)
    end if
  else
    iData(1) = TocOne(pSym)
    if (debug) write(u6,'(a,z8)') ' Reading symmetry operators:',iData(1)
  end if
else if (label == 'DEGDISP ') then
  if (TocOne(pndisp) == NaN) call SysAbendMsg(TheName,'Undefined Label:',Label)
  Length = TocOne(pndisp)
  iData(1:Length) = TocOne(pdegdisp:pdegdisp+Length-1)
else
!----------------------------------------------------------------------*
! Read operators from integral records                                 *
!----------------------------------------------------------------------*
  if (iand(option,sRdNxt) /= 0) then
    if (debug) write(u6,'(a)') ' Reading next item'
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
    if (debug) write(u6,'(a)') ' Reading first item'
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
    if (debug) write(u6,'(a)') ' Reading current item'
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
    Length = nTri_Elem(na)
    Length = nTri_Elem(Length)
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
      Length = Length+nTri_Elem(TocOne(pldisp+isym))
    end do
  else if (label == 'INACTIVE') then
    Length = 0
    do iS=1,TocOne(pSym)
      do jS=1,TocOne(pSym)
        ijS = Mul(iS,jS)-1
        if (btest(SymLab,ijS)) then
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
        ijS = Mul(iS,jS)-1
        if (btest(SymLab,ijS)) then
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
        ij = Mul(i,j)-1
        if (btest(SymLab,ij)) then
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

  iData(1) = Length
  if (iand(option,sOpSiz) == 0) then
    Length = rtoi*Length
    IndDta = 1
    !IndHld = 1
    iDisk = TocOne(pOp+LenOp*(CurrOp-1)+oAddr)
    do k=0,Length+nAuxDt-1,nBuf
      iBuf = max(0,min(nBuf,Length+nAuxDt-k))
      tBuf = max(0,min(nBuf,Length-k))
      eBuf = iBuf-tBuf
      call iDaFile(LuMCK,2,TmpBuf,iBuf,iDisk)
      IndTmp = 1
      iData(IndDta:IndDta+tBuf-1) = TmpBuf(IndTmp:IndTmp+tBuf-1)
      if (debug) then
        write(u6,'(a,z8)') ' Reading buffer to: ',IndDta
        write(u6,'(8(1x,z8))') (iData(IndDta+m),m=0,tBuf-1)
      end if
      IndTmp = IndTmp+tBuf
      IndDta = IndDta+tBuf
      if (tBuf < iBuf) then
        !HldBuf(IndHld:IndHld+eBuf-1) = TmpBuf(IndTmp:IndTmp+eBuf-1)
        IndTmp = IndTmp+eBuf
        !IndHld = IndHld+eBuf
      end if
    end do
    !if (iAnd(sNoOri,option) == 0) then
    !  iData(IndDta:IndDta+5) = HldBuf(1:6)
    !  if (debug) then
    !    write(u6,'(a,z8)') ' Reading buffer to: ',IndDta
    !    write(u6,'(8(1x,z8))') (iData(IndDta+m),m=0,5)
    !  end if
    !end if
    IndDta = IndDta+6
    !if (iAnd(sNoNuc,option) == 0) then
    !  iData(IndDta:IndDta+1) = HldBuf(7:8)
    !  if (debug) then
    !    write(u6,'(a,z8)') ' Reading buffer to: ',IndDta
    !    write(u6,'(8(1x,z8))') (iData(IndDta+m),m=0,1)
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
