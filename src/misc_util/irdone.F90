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
!***********************************************************************

subroutine iRdOne(rc,Option,InLab,Comp,rData,SymLab)
!***********************************************************************
!                                                                      *
!     Purpose: Read data from one-electron integral file               *
!                                                                      *
!     Calling parameters:                                              *
!     rc      : Return code                                            *
!     Option  : Switch to set options                                  *
!     InLab   : Identifier for the data to read                        *
!     Comp    : Composite identifier to select components              *
!     rData   : contains on output the requested data                  *
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
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark and M.P. Fuelscher                                  *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use OneDat, only: AuxOne, LenOp, lTocOne, MxOp, NaN, nBas, nSym, oAddr, oComp, oLabel, oSymLb, pOp, rcOne, sDbg, sNoNuc, sNoOri, &
                  sOpSiz, sRdCur, sRdFst, sRdNxt, TocOne
use Definitions, only: wp, iwp, u6, RtoI, ItoB

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp), intent(in) :: Option
character(len=*), intent(inout) :: InLab
integer(kind=iwp), intent(inout) :: Comp, SymLab
integer(kind=iwp), intent(_OUT_) :: rData(*)
integer(kind=iwp), parameter :: lBuf = 1024
integer(kind=iwp) :: CmpTmp, CurrOp = 1, i, iDisk, idx, ij, iLen, IndAux, IndDta, iOpt, iRC, j, Length, LuOne, nCopy, nSave, TmpCmp
real(kind=wp) :: AuxBuf(4), TmpBuf(lBuf)
character(len=8) :: Label, TmpLab
logical(kind=iwp) :: debug, doclose
integer(kind=iwp), external :: isFreeUnit

!----------------------------------------------------------------------*
! Start procedure:                                                     *
! Pick up the file definitions                                         *
!----------------------------------------------------------------------*
rc = rcOne%good
LuOne = AuxOne%Lu
!----------------------------------------------------------------------*
! Check the file status                                                *
!----------------------------------------------------------------------*
doclose = .false.
if (.not. AuxOne%Opn) then
  !write(u6,*) ' I will open the file for you!'

  ! Well, I'll open and close it for you under the default name

  LuOne = 77
  LuOne = isFreeUnit(LuOne)
  Label = 'ONEINT  '
  !write(u6,*) 'RdOne: opening OneInt'
  iRC = -1
  iOpt = 0
  call OpnOne(iRC,iOpt,Label,LuOne)
  if (iRC /= 0) then
    write(u6,*) 'RdOne: Error opening file'
    call Abend()
  end if
  doclose = .true.
end if
!----------------------------------------------------------------------*
! Truncate the label to 8 characters and convert it to upper case      *
!----------------------------------------------------------------------*
Label = InLab
call UpCase(Label)
TmpLab = Label
iLen = len(TmpLab)/ItoB
!----------------------------------------------------------------------*
! Print debugging information                                          *
!----------------------------------------------------------------------*
debug = .false.
if (btest(option,sDbg)) debug = .true.
if (debug) then
  write(u6,*) '<<< Entering RdOne >>>'
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
  call SysWarnMsg('RdOne','Invalid value','sRdFst and sRdNxt')
else if (btest(option,sRdFst) .and. btest(option,sRdCur)) then
  call SysWarnMsg('RdOne','Invalid value','sRdFst and sRdNxt')
else if (btest(option,sRdNxt) .and. btest(option,sRdCur)) then
  call SysWarnMsg('RdOne','Invalid value','sRdNxt and sRdCur')
end if
!----------------------------------------------------------------------*
! Load back TocOne                                                     *
!----------------------------------------------------------------------*
iDisk = 0
call iDaFile(LuOne,2,TocOne,lTocOne,iDisk)
!----------------------------------------------------------------------*
! Read operators from integral records                                 *
!----------------------------------------------------------------------*
if (btest(option,sRdNxt)) then
  CurrOp = CurrOp+1
  if (CurrOp > MxOp) then
    CurrOp = 0
  else if (TocOne(pOp+LenOp*(CurrOp-1)+oLabel) == NaN) then
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
else if (btest(option,sRdFst)) then
  CurrOp = 1
  if (TocOne(pOp+LenOp*(CurrOp-1)+oLabel) == NaN) then
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
else if (btest(option,sRdCur)) then
  if ((CurrOp < 1) .or. (CurrOp > MxOp)) then
    CurrOp = 0
  else if (TocOne(pOp+LenOp*(CurrOp-1)+oLabel) == NaN) then
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
if (CurrOp == 0) then
  rc = rcOne%RD03
  !write(u6,*) 'RdOne: Information not available'
  !write(u6,*) 'Option=',Option
  !write(u6,*) 'Comp=',Comp
  !write(u6,*) 'SymLab=',SymLab
  !write(u6,*) 'Label=',Label
else
  SymLab = TocOne(pOp+LenOp*(CurrOp-1)+oSymLb)
  Length = 0
  do i=1,nSym
    do j=1,i
      ij = Mul(i,j)-1
      if (btest(SymLab,ij)) then
        if (i == j) then
          Length = Length+nTri_Elem(nBas(i))
        else
          Length = Length+nBas(i)*nBas(j)
        end if
      end if
    end do
  end do
  rData(1) = Length
  if (.not. btest(option,sOpSiz)) then
    IndAux = 0
    IndDta = 0
    iDisk = TocOne(pOp+LenOp*(CurrOp-1)+oAddr)
    do i=0,Length+3,lBuf
      nCopy = max(0,min(lBuf,Length+4-i))
      nSave = max(0,min(lBuf,Length-i))
      call dDaFile(LuOne,2,TmpBuf,nCopy,iDisk)
      call idCopy(nSave,TmpBuf,1,rData(IndDta+1),1)
      IndDta = IndDta+RtoI*nSave
      do j=nSave+1,nCopy
        IndAux = IndAux+1
        !AuxBuf(IndAux) = TmpBuf(nSave+IndAux)
        AuxBuf(IndAux) = TmpBuf(j)
      end do
    end do
    if (.not. btest(option,sNoOri)) call idCopy(3,AuxBuf,1,rData(IndDta+1),1)
    if (.not. btest(option,sNoNuc)) call idCopy(1,AuxBuf(4),1,rData(IndDta+RtoI*3+1),1)
  end if
end if

if (doclose) then
  !write(u6,*) ' I will close the file for you!'
  iRC = -1
  iOpt = 0
  call ClsOne(iRC,iOpt)
  if (iRC /= 0) then
    write(u6,*) 'RdOne: Error closing file'
    call Abend()
  end if
end if

!----------------------------------------------------------------------*
! Terminate procedure                                                  *
!----------------------------------------------------------------------*
return

! This is to allow type punning without an explicit interface
contains

subroutine idCopy(n,Src,n1,Dst,n2)

  integer(kind=iwp), intent(in) :: n, n1, n2
  real(kind=wp), intent(in) :: Src(*)
  integer(kind=iwp), intent(_OUT_), target :: Dst(*)
  real(kind=wp), pointer :: dDst(:)

  call c_f_pointer(c_loc(Dst),dDst,[n])
  call dCopy_(n,Src,n1,dDst,n2)
  nullify(dDst)

end subroutine idCopy

end subroutine iRdOne
