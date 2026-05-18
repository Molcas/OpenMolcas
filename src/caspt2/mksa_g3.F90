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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MKSA_G3(ISYM,SA,NSA,NG3,G3,idxG3)

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp, Byte
use SUPERINDEX, only: KTUV
use caspt2_module, only: NASHT, IASYM, NTUVES

implicit none
integer(kind=iwp), intent(in) :: ISYM, NSA, NG3
real(kind=wp), intent(out) :: SA(NSA)
real(kind=wp), intent(in) :: G3(NG3)
integer(kind=Byte), intent(in) :: idxG3(6,NG3)
integer(kind=iwp) iG3, iT, iU, iV, iX, iY, iZ, iST, iSU, iSV, iSX, iSY, iSZ, ituvs, ixyzs, iTU, iVX, iYZ, JSYM, ISUP, JSUP, ISADR
real(kind=wp) G3VAL

!-SVC20100831: determine indices in SA where a certain G3 value will end up
do iG3=1,NG3
  iT = idxG3(1,iG3)
  iU = idxG3(2,iG3)
  iV = idxG3(3,iG3)
  iX = idxG3(4,iG3)
  iY = idxG3(5,iG3)
  iZ = idxG3(6,iG3)
  iST = IASYM(iT)
  iSU = IASYM(iU)
  iSV = IASYM(iV)
  iSX = IASYM(iX)
  iSY = IASYM(iY)
  iSZ = IASYM(iZ)
  ituvs = Mul(IST,Mul(ISU,ISV))
  ixyzs = Mul(ISX,Mul(ISY,ISZ))
  if (ituvs /= ixyzs) cycle
  iTU = iT+NASHT*(iU-1)
  iVX = iV+NASHT*(iX-1)
  iYZ = iY+NASHT*(iZ-1)
  G3VAL = -G3(iG3)
  !-SVC20100829: 12 equivalent cases, of which the second
  !  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
  !  - G(tuvxyz) -> SA(xut,vyz)
  jSYM = Mul(IASYM(iX),Mul(IASYM(iU),IASYM(iT)))
  if (jSYM == iSYM) then
    ISUP = KTUV(iX,iU,iT)-nTUVES(jSYM)
    JSUP = KTUV(iV,iY,iZ)-nTUVES(jSYM)
    if (JSUP <= ISUP) then
      ISADR = (ISUP*(ISUP-1))/2+JSUP
      SA(ISADR) = G3VAL
    end if
  end if

  if (.not. ((iTU == iVX) .and. (iVX == iYZ))) then

    if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
      !  - G(vxtuyz) -> SA(uxv,tyz)
      jSYM = Mul(IASYM(iU),Mul(IASYM(iX),IASYM(iV)))
      if (jSYM == iSYM) then
        ISUP = KTUV(iU,iX,iV)-nTUVES(jSYM)
        JSUP = KTUV(iT,iY,iZ)-nTUVES(jSYM)
        if (JSUP <= ISUP) then
          ISADR = (ISUP*(ISUP-1))/2+JSUP
          SA(ISADR) = G3VAL
        end if
      end if
      !  - G(yzvxtu) -> SA(xzy,vtu)
      jSYM = Mul(IASYM(iX),Mul(IASYM(iZ),IASYM(iY)))
      if (jSYM == iSYM) then
        ISUP = KTUV(iX,iZ,iY)-nTUVES(jSYM)
        JSUP = KTUV(iV,iT,iU)-nTUVES(jSYM)
        if (JSUP <= ISUP) then
          ISADR = (ISUP*(ISUP-1))/2+JSUP
          SA(ISADR) = G3VAL
        end if
      end if
      !  - G(tuyzvx) -> SA(zut,yvx)
      jSYM = Mul(IASYM(iZ),Mul(IASYM(iU),IASYM(iT)))
      if (jSYM == iSYM) then
        ISUP = KTUV(iZ,iU,iT)-nTUVES(jSYM)
        JSUP = KTUV(iY,iV,iX)-nTUVES(jSYM)
        if (JSUP <= ISUP) then
          ISADR = (ISUP*(ISUP-1))/2+JSUP
          SA(ISADR) = G3VAL
        end if
      end if
    end if

    !  - G(yztuvx) -> SA(uzy,tvx)
    jSYM = Mul(IASYM(iU),Mul(IASYM(iZ),IASYM(iY)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iU,iZ,iY)-nTUVES(jSYM)
      JSUP = KTUV(iT,iV,iX)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = (ISUP*(ISUP-1))/2+JSUP
        SA(ISADR) = G3VAL
      end if
    end if
    !  - G(vxyztu) -> SA(zxv,ytu)
    jSYM = Mul(IASYM(iZ),Mul(IASYM(iX),IASYM(iV)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iZ,iX,iV)-nTUVES(jSYM)
      JSUP = KTUV(iY,iT,iU)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = (ISUP*(ISUP-1))/2+JSUP
        SA(ISADR) = G3VAL
      end if
    end if

  end if

  if ((iT == iU) .and. (iV == iX) .and. (iY == iZ)) cycle
  if ((iT == iU) .and. (iV == iZ) .and. (iX == iY)) cycle
  if ((iX == iV) .and. (iT == iZ) .and. (iU == iY)) cycle
  if ((iZ == iY) .and. (iV == iU) .and. (iX == iT)) cycle
  !  - G(utxvzy) -> SA(vtu,xzy)
  jSYM = Mul(IASYM(iV),Mul(IASYM(iT),IASYM(iU)))
  if (jSYM == iSYM) then
    ISUP = KTUV(iV,iT,iU)-nTUVES(jSYM)
    JSUP = KTUV(iX,iZ,iY)-nTUVES(jSYM)
    if (JSUP <= ISUP) then
      ISADR = (ISUP*(ISUP-1))/2+JSUP
      SA(ISADR) = G3VAL
    end if
  end if
  if ((iTU == iVX) .and. (iVX == iYZ)) cycle

  if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
    !  - G(xvutzy) -> SA(tvx,uzy)
    jSYM = Mul(IASYM(iT),Mul(IASYM(iV),IASYM(iX)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iT,iV,iX)-nTUVES(jSYM)
      JSUP = KTUV(iU,iZ,iY)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = (ISUP*(ISUP-1))/2+JSUP
        SA(ISADR) = G3VAL
      end if
    end if
    !  - G(zyxvut) -> SA(vyz,xut)
    jSYM = Mul(IASYM(iV),Mul(IASYM(iY),IASYM(iZ)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iV,iY,iZ)-nTUVES(jSYM)
      JSUP = KTUV(iX,iU,iT)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = (ISUP*(ISUP-1))/2+JSUP
        SA(ISADR) = G3VAL
      end if
    end if
    !  - G(utzyxv) -> SA(ytu,zxv)
    jSYM = Mul(IASYM(iY),Mul(IASYM(iT),IASYM(iU)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iY,iT,iU)-nTUVES(jSYM)
      JSUP = KTUV(iZ,iX,iV)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = (ISUP*(ISUP-1))/2+JSUP
        SA(ISADR) = G3VAL
      end if
    end if
  end if

  !  - G(zyutxv) -> SA(tyz,uxv)
  jSYM = Mul(IASYM(iT),Mul(IASYM(iY),IASYM(iZ)))
  if (jSYM == iSYM) then
    ISUP = KTUV(iT,iY,iZ)-nTUVES(jSYM)
    JSUP = KTUV(iU,iX,iV)-nTUVES(jSYM)
    if (JSUP <= ISUP) then
      ISADR = (ISUP*(ISUP-1))/2+JSUP
      SA(ISADR) = G3VAL
    end if
  end if
  !  - G(xvzyut) -> SA(yvx,zut)
  jSYM = Mul(IASYM(iY),Mul(IASYM(iV),IASYM(iX)))
  if (jSYM == iSYM) then
    ISUP = KTUV(iY,iV,iX)-nTUVES(jSYM)
    JSUP = KTUV(iZ,iU,iT)-nTUVES(jSYM)
    if (JSUP <= ISUP) then
      ISADR = (ISUP*(ISUP-1))/2+JSUP
      SA(ISADR) = G3VAL
    end if
  end if

end do

end subroutine MKSA_G3
