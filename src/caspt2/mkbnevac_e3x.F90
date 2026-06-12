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
! Copyright (C) 2026, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine MKBNEVAC_E3x(ISYM,nAshT,NG3,NBA,NBC,BA,BC,G1,G2,G3,Hact,Htilde,Gact,idxG3)

use caspt2_module, only: IASYM, NTUVES
use SUPERINDEX, only: KTUV
use Symmetry_Info, only: Mul
use nevpt2_mod, only: ex2, ex3
use Constants, only: Two
use Definitions, only: wp, iwp, byte

implicit none
integer(kind=iwp), intent(in) :: iSYM, nAshT, NG3, NBA, NBC
real(kind=wp), intent(inout) :: BA(NBA), BC(NBC)
real(kind=wp), intent(in) :: G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT), G3(*), Hact(nAshT,nAshT), Htilde(nAshT,nAshT), &
                             Gact(nAshT,nAshT,nAshT,nAshT)
integer(kind=byte), intent(in) :: idxG3(6,NG3)
integer(kind=iwp) :: iG3, IT, ITU, IU, IV, IVX, IX, IY, IYZ, IZ
real(kind=wp) :: G3VAL

do IG3=1,NG3
  iT = idxG3(1,iG3)
  iU = idxG3(2,iG3)
  iV = idxG3(3,iG3)
  iX = idxG3(4,iG3)
  iY = idxG3(5,iG3)
  iZ = idxG3(6,iG3)
  !iST = IASYM(iT)
  !iSU = IASYM(iU)
  !iSV = IASYM(iV)
  !iSX = IASYM(iX)
  !iSY = IASYM(iY)
  !iSZ = IASYM(iZ)
  !ituvs = MUL(IST,MUL(ISU,ISV))
  !ixyzs = MUL(ISX,MUL(ISY,ISZ))
  !if (ituvs /= ixyzs) cycle
  iTU = iT+NASHT*(iU-1)
  iVX = iV+NASHT*(iX-1)
  iYZ = iY+NASHT*(iZ-1)

  ! G3(tuvxyz) = <Psi|t+ v+ y+ z x u|Psi>
  !G3VAL = tG3(IT,IU,IV,IX,IY,IZ,G1,G2,G3(iG3))

  ! D(tuvxyz)
  !if ((MUL(iST,iSU) == iSym) .or. (MUL(iSV,ISX) == iSym) .or. (MUL(iSY,iSZ) == iSym)) then
  G3VAL = ex3(IT,IU,IV,IX,IY,IZ,G1,G2,G3(iG3))
  call Add_MKBNEVAC(iT,iU,iV,iX,iY,iZ)
  !end if

  if ((iTU /= iVX) .or. (iVX /= iYZ)) then
    if ((iTU /= iVX) .and. (iTU /= iYZ) .and. (iVX /= iYZ)) then
      ! D(vxtuyz)
      !if ((MUL(iSX,iSV) == iSym) .or. (MUL(iSX,iST) == iSym) .or. (MUL(iSU,iSY) == iSym)) then
      G3VAL = ex3(IV,IX,IT,IU,IY,IZ,G1,G2,G3(iG3))
      call Add_MKBNEVAC(iV,iX,iT,iU,iY,iZ)
      !end if
      ! D(yzvxtu)
      !if ((MUL(iSZ,iSY) == iSym) .or. (MUL(iSZ,iSV) == iSym) .or. (MUL(iSX,iST) == iSym)) then
      G3VAL = ex3(IY,IZ,IV,IX,IT,IU,G1,G2,G3(iG3))
      call Add_MKBNEVAC(iY,iZ,iV,iX,iT,iU)
      !end if
      ! D(tuyzvx)
      !if ((MUL(iSU,iST) == iSym) .or. (MUL(iSU,iSY) == iSym) .or. (MUL(iSZ,iSV) == iSym)) then
      G3VAL = ex3(IT,IU,IY,IZ,IV,IX,G1,G2,G3(iG3))
      call Add_MKBNEVAC(iT,iU,iY,iZ,iV,iX)
      !end if
    end if
    ! D(yztuvx)
    !if ((MUL(iSZ,iSY) == iSym) .or. (MUL(iSZ,iST) == iSym) .or. (MUL(iSU,iSV) == iSym)) then
    G3VAL = ex3(IY,IZ,IT,IU,IV,IX,G1,G2,G3(iG3))
    call Add_MKBNEVAC(iY,iZ,iT,iU,iV,iX)
    !end if
    ! D(vxyztu)
    !if ((MUL(iSX,iSV) == iSym) .or. (MUL(iSX,iSY) == iSym) .or. (MUL(iSZ,iST) == iSym)) then
    G3VAL = ex3(IV,IX,IY,IZ,IT,IU,G1,G2,G3(iG3))
    call Add_MKBNEVAC(iV,iX,iY,iZ,iT,iU)
    !end if
  end if

  if ((iT == iU) .and. (iV == iX) .and. (iY == iZ)) cycle
  if ((iT == iU) .and. (iV == iZ) .and. (iX == iY)) cycle
  if ((iX == iV) .and. (iT == iZ) .and. (iU == iY)) cycle
  if ((iZ == iY) .and. (iV == iU) .and. (iX == iT)) cycle

  ! D(utxvzy)
  !if ((MUL(iST,iSU) == iSym) .or. (MUL(iST,iSX) == iSym) .or. (MUL(iSV,iSZ) == iSym)) then
  G3VAL = ex3(IU,IT,IX,IV,IZ,IY,G1,G2,G3(iG3))
  call Add_MKBNEVAC(iU,iT,iX,iV,iZ,iY)
  !end if

  if ((iTU == iVX) .and. (iVX == iYZ)) cycle
  if ((iTU /= iVX) .and. (iTU /= iYZ) .and. (iVX /= iYZ)) then
    ! D(xvutzy)
    !if ((MUL(iSV,iSX) == iSym) .or. (MUL(iSV,iSU) == iSym) .or. (MUL(iST,iSZ) == iSym)) then
    G3VAL = ex3(IX,IV,IU,IT,IZ,IY,G1,G2,G3(iG3))
    call Add_MKBNEVAC(iX,iV,iU,iT,iZ,iY)
    !end if
    ! D(zyxvut)
    !if ((MUL(iSY,iSZ) == iSym) .or. (MUL(iSY,iSX) == iSym) .or. (MUL(iSV,iSU) == iSym)) then
    G3VAL = ex3(IZ,IY,IX,IV,IU,IT,G1,G2,G3(iG3))
    call Add_MKBNEVAC(iZ,iY,iX,iV,iU,iT)
    !end if
    ! D(utzyxv)
    !if ((MUL(iST,iSU) == iSym) .or. (MUL(iST,iSZ) == iSym) .or. (MUL(iSY,iSX) == iSym)) then
    G3VAL = ex3(IU,IT,IZ,IY,IX,IV,G1,G2,G3(iG3))
    call Add_MKBNEVAC(iU,iT,iZ,iY,iX,iV)
    !end if
  end if
  ! D(zyutxv)
  !if ((MUL(iSY,iSZ) == iSym) .or. (MUL(iSY,iSU) == iSym) .or. (MUL(iST,iSX) == iSym)) then
  G3VAL = ex3(IZ,IY,IU,IT,IX,IV,G1,G2,G3(iG3))
  call Add_MKBNEVAC(iZ,iY,iU,iT,iX,iV)
  !end if
  ! D(xvzyut)
  !if ((MUL(iSV,iSX) == iSym) .or. (MUL(iSV,iSZ) == iSym) .or. (MUL(iSY,iSU) == iSym)) then
  G3VAL = ex3(IX,IV,IZ,IY,IU,IT,G1,G2,G3(iG3))
  call Add_MKBNEVAC(iX,iV,iZ,iY,iU,iT)
  !end if
end do

contains

subroutine Add_MKBNEVAC(ITABS_,IUABS_,IVABS_,IXABS_,IYABS_,IZABS_)

  use Index_Functions, only: iTri

  integer(kind=iwp), intent(in) :: ITABS_, IUABS_, IVABS_, IXABS_, IYABS_, IZABS_
  integer(kind=iwp) :: IAABS, IBABS, IBADR, ISAB, ITUV_, IWABS, IXYZ_

  iTUV_ = KTUV(IVABS_,IUABS_,ITABS_)-NTUVES(ISYM)
  if (MUL(IASYM(iVABS_),MUL(IASYM(IUABS_),IASYM(ITABS_))) == iSym) then
    do IWABS=1,nAshT
      ! C: A(vut,xwz) <-- +h_{yw} Etu*Evx*Eyz
      if (MUL(IASYM(IWABS),MUL(IASYM(IXABS_),IASYM(IZABS_))) == iSym) then
        iXYZ_ = KTUV(IXABS_,IWABS,IZABS_)-NTUVES(ISYM)
        if (iTUV_ >= iXYZ_) then
          IBADR = iTri(iTUV_,iXYZ_)
          BC(IBADR) = BC(IBADR)+G3VAL*Hact(IYABS_,IWABS)
        end if
      end if
      ! C: A(vut,xyw) <-- -htilde_{wz} Etu*Evx*Eyz
      if (MUL(IASYM(IWABS),MUL(IASYM(IXABS_),IASYM(IYABS_))) == iSym) then
        iXYZ_ = KTUV(IXABS_,IYABS_,IWABS)-NTUVES(ISYM)
        if (iTUV_ >= iXYZ_) then
          IBADR = iTri(iTUV_,iXYZ_)
          BC(IBADR) = BC(IBADR)-G3VAL*Htilde(IWABS,IZABS_)
        end if
      end if
      ! C: A(vut,wyz) <-- -htilde_{wx} Etu*Evx*Eyz
      if (MUL(IASYM(IWABS),MUL(IASYM(IYABS_),IASYM(IZABS_))) == iSym) then
        iXYZ_ = KTUV(IWABS,IYABS_,IZABS_)-NTUVES(ISYM)
        if (iTUV_ >= iXYZ_) then
          IBADR = iTri(iTUV_,iXYZ_)
          BC(IBADR) = BC(IBADR)-G3VAL*Htilde(IWABS,IXABS_)
        end if
      end if
    end do

    do IAABS=1,nAshT
      do IBABS=1,nAshT
        ISAB = MUL(IASYM(IAABS),IASYM(IBABS))
        ! C: A_{vut,abx} <-- -(ab|yz) * Etu Evx Eyz
        if (MUL(IASYM(IXABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IAABS,IBABS,IXABS_)-NTUVES(ISYM)
          if (iTUV_ >= iXYZ_) then
            IBADR = iTri(iTUV_,iXYZ_)
            BC(IBADR) = BC(IBADR)-G3VAL*Gact(IAABS,IBABS,IYABS_,IZABS_)
          end if
        end if
        ! C: A_{vut,xab} <-- -(ay|bz) * Etu Evx Eyz
        if (MUL(IASYM(IXABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IXABS_,IAABS,IBABS)-NTUVES(ISYM)
          if (iTUV_ >= iXYZ_) then
            IBADR = iTri(iTUV_,iXYZ_)
            BC(IBADR) = BC(IBADR)-G3VAL*Gact(IAABS,IYABS_,IBABS,IZABS_)
          end if
        end if
        ! C: A_{vut,abz} <-- -(ax|by) * Etu Evx Eyz
        if (MUL(IASYM(IZABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IAABS,IBABS,IZABS_)-NTUVES(ISYM)
          if (iTUV_ >= iXYZ_) then
            IBADR = iTri(iTUV_,iXYZ_)
            BC(IBADR) = BC(IBADR)-G3VAL*Gact(IAABS,IXABS_,IBABS,IYABS_)
          end if
        end if
        ! C: A_{vut,ayb} <-- +(ax|bz) * Etu Evx Eyz
        if (MUL(IASYM(IYABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IAABS,IYABS_,IBABS)-NTUVES(ISYM)
          if (iTUV_ >= iXYZ_) then
            IBADR = iTri(iTUV_,iXYZ_)
            BC(IBADR) = BC(IBADR)+G3VAL*Gact(IAABS,IXABS_,IBABS,IZABS_)
          end if
        end if
      end do
    end do
  end if

  do IAABS=1,nAshT
    do IBABS=1,nAshT
      ! (1.2.11) ERI(2)
      ! A: A_{aut,byz} <-- +2(vx|ab) * Etu Evx Eyz
      if ((MUL(IASYM(IAABS),MUL(IASYM(IUABS_),IASYM(ITABS_))) == iSym) .and. &
          (MUL(IASYM(IBABS),MUL(IASYM(IYABS_),IASYM(IZABS_))) == iSym)) then
        iTUV_ = KTUV(IAABS,IUABS_,ITABS_)-NTUVES(ISYM)
        iXYZ_ = KTUV(IBABS,IYABS_,IZABS_)-NTUVES(ISYM)
        if (iTUV_ >= iXYZ_) then
          IBADR = iTri(iTUV_,iXYZ_)
          BA(IBADR) = BA(IBADR)+Two*G3VAL*Gact(IAABS,IBABS,IVABS_,IXABS_)
        end if
      end if
    end do
  end do

  iTUV_ = KTUV(IXABS_,IUABS_,ITABS_)-NTUVES(ISYM)
  if (MUL(IASYM(IXABS_),MUL(IASYM(IUABS_),IASYM(ITABS_))) == iSym) then
    do IAABS=1,nAshT
      do IBABS=1,nAshT
        ISAB = MUL(IASYM(IAABS),IASYM(IBABS))
        ! (1.2.11) ERI(3)
        ! A: A_{xut,avb} <-- +(ab|yz) * Etu Evx Eyz
        if (MUL(IASYM(IVABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IAABS,IVABS_,IBABS)-NTUVES(ISYM)
          if (iTUV_ >= iXYZ_) then
            IBADR = iTri(iTUV_,iXYZ_)
            BA(IBADR) = BA(IBADR)+G3VAL*Gact(IAABS,IBABS,IYABS_,IZABS_)
          end if
        end if
        ! (1.2.11) ERI(4)
        ! A: A_{xut,abz} <-- -(va|yb) * Etu Evx Eyz
        if (MUL(IASYM(IZABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IAABS,IBABS,IZABS_)-NTUVES(ISYM)
          if (iTUV_ >= iXYZ_) then
            IBADR = iTri(iTUV_,iXYZ_)
            BA(IBADR) = BA(IBADR)-G3VAL*Gact(IAABS,IVABS_,IYABS_,IBABS)
          end if
        end if
        ! (1.2.11) ERI(5)
        ! A: A_{xut,ayb} <-- +(va|bz) * Etu Evx Eyz
        if (MUL(IASYM(IYABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IAABS,IYABS_,IBABS)-NTUVES(ISYM)
          if (iTUV_ >= iXYZ_) then
            IBADR = iTri(iTUV_,iXYZ_)
            BA(IBADR) = BA(IBADR)+G3VAL*Gact(IAABS,IVABS_,IZABS_,IBABS)
          end if
        end if
      end do
    end do

    G3VAL = -G3VAL
    if (IVABS_ == IXABS_) G3VAL = Two*ex2(ITABS_,IUABS_,IYABS_,IZABS_,G1,G2)+G3VAL

    do IWABS=1,nAshT
      ! A: A(xut,vwz) <-- +h_{yw} EtuEvxEyz
      if (MUL(IASYM(IWABS),MUL(IASYM(IVABS_),IASYM(IZABS_))) == iSym) then
        iXYZ_ = KTUV(IVABS_,IWABS,IZABS_)-NTUVES(ISYM)
        if (iTUV_ >= iXYZ_) then
          IBADR = iTri(iTUV_,iXYZ_)
          BA(IBADR) = BA(IBADR)+G3VAL*Hact(IYABS_,IWABS)
        end if
      end if
      ! A: A(xut,vyw) <-- -h_{wz} EtuEvxEyz
      if (MUL(IASYM(IWABS),MUL(IASYM(IVABS_),IASYM(IYABS_))) == iSym) then
        iXYZ_ = KTUV(IVABS_,IYABS_,IWABS)-NTUVES(ISYM)
        if (iTUV_ >= iXYZ_) then
          IBADR = iTri(iTUV_,iXYZ_)
          BA(IBADR) = BA(IBADR)-G3VAL*Htilde(IWABS,IZABS_)
        end if
      end if
      ! A: A(xut,wyz) <-- +h_{wv} EtuEvxEyz
      if (MUL(IASYM(IWABS),MUL(IASYM(IYABS_),IASYM(IZABS_))) == iSym) then
        iXYZ_ = KTUV(IWABS,IYABS_,IZABS_)-NTUVES(ISYM)
        if (iTUV_ >= iXYZ_) then
          IBADR = iTri(iTUV_,iXYZ_)
          BA(IBADR) = BA(IBADR)+G3VAL*Hact(IWABS,IVABS_)
        end if
      end if
    end do

    do IAABS=1,nAshT
      do IBABS=1,nAshT
        ISAB = MUL(IASYM(IAABS),IASYM(IBABS))
        ! (1.2.11) ERI(1)
        ! A: A_{xut,vab} <-- +(ya|zb) * Etu Evx Eyz
        if (MUL(IASYM(IVABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IVABS_,IAABS,IBABS)-NTUVES(ISYM)
          if (iTUV_ >= iXYZ_) then
            IBADR = iTri(iTUV_,iXYZ_)
            BA(IBADR) = BA(IBADR)-G3VAL*Gact(IAABS,IYABS_,IBABS,IZABS_)
          end if
        end if
      end do
    end do
  end if

end subroutine Add_MKBNEVAC

end subroutine MKBNEVAC_E3x
