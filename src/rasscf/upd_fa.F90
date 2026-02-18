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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Upd_FA(PUVX,F,D,ExFac)
!***********************************************************************
!                                                                      *
!     compute FIA, FAA, and FAS from the integral set (pu!vx)          *
!                                                                      *
!     calling arguments:                                               *
!     PUVX    : input, array of real                                   *
!               ERIs with indices (pu!vx)                              *
!               (three active, one general)                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use general_data, only: NASH, NISH, NORB, NSYM
use Molcas, only: MxSym
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: PUVX(*), F(*), D(*), ExFac
integer(kind=iwp) :: IASH, icase, IDOFF, IFOFF, IISH, IJSYM, IORB, IP, IPU, IPUVX, IPV, IPX, ISTACK, ISYM, IU, IUV, IUX, IV, IVX, &
                     IX, JASH, JDOFF, JFOFF, JISH, JORB, JSYM, KASH, KDOFF, KL_ORB_PAIRS, KLSYM, KSYM, LASH, LSYM, &
                     off_Dmat(mxSym), off_Fmat(mxSym), off_PUVX(mxSym)
real(kind=wp) :: DUV, DUX, DVX, TEMP

! generate offsets

iStack = 0
do iSym=1,nSym
  off_Dmat(iSym) = iStack
  iAsh = nAsh(iSym)
  iStack = iStack+nTri_Elem(iAsh)
end do

iStack = 0
do iSym=1,nSym
  off_Fmat(iSym) = iStack
  iOrb = nOrb(iSym)
  iStack = iStack+nTri_Elem(iOrb)
end do

iStack = 0
do iSym=1,nSym
  off_PUVX(iSym) = iStack
  iOrb = nOrb(iSym)
  do jSym=1,nSym
    jAsh = nAsh(jSym)
    ijSym = Mul(iSym,jSym)
    do kSym=1,nSym
      kAsh = nAsh(kSym)
      do lSym=1,kSym
        lAsh = nAsh(lSym)
        klSym = Mul(kSym,lSym)
        if (ijSym == klSym) then
          kl_Orb_pairs = kAsh*lAsh
          if (kSym == lSym) kl_Orb_pairs = nTri_Elem(kAsh)
          iStack = iStack+iOrb*jAsh*kl_Orb_pairs
        end if
      end do
    end do
  end do
end do

! clear the subbocks FIA, FAA and FSA

do iSym=1,nSym
  iOrb = nOrb(iSym)
  iAsh = nAsh(iSym)
  iIsh = nIsh(iSym)
  iFoff = off_Fmat(iSym)
  do iU=iIsh+1,iIsh+iAsh
    do iP=1,iU
      F(iFoff+iTri(iU,iP)) = Zero
    end do
  end do
  do iU=iIsh+iAsh+1,iOrb
    do iP=iIsh+1,iIsh+iAsh
      F(iFoff+iTri(iP,iU)) = Zero
    end do
  end do
end do

! generate the subblocks FIA, FAA and FSA

do iSym=1,nSym
  iOrb = nOrb(iSym)
  iAsh = nAsh(iSym)
  iIsh = nIsh(iSym)
  iPUVX = off_PUVX(iSym)
  do jSym=1,iSym
    jOrb = nOrb(jSym)
    jAsh = nAsh(jSym)
    jIsh = nIsh(jSym)
    ijSym = Mul(iSym,jSym)
    do kSym=1,nSym
      kAsh = nAsh(kSym)
      do lSym=1,kSym
        lAsh = nAsh(lSym)
        klSym = Mul(kSym,lSym)

        ! find cases
        icase = 4
        if (iSym == jSym) icase = icase-2
        if (iSym == kSym) icase = icase-1

        if ((ijSym == klSym) .and. (iAsh*jAsh*kAsh*lAsh /= 0)) then

          select case (icase)

            case (1)
              ! symmetry case (II!II)
              iFoff = off_Fmat(iSym)
              iDoff = off_Dmat(iSym)
              do iV=1,kAsh
                do iX=1,iV
                  iVX = iTri(iV,iX)
                  DVX = Two*D(iDoff+iVX)
                  if (iX == iV) DVX = D(iDoff+iVX)
                  do iU=1,jAsh
                    iUV = iTri(iU,iV)
                    DUV = ExFac*Half*D(iDoff+iUV)
                    iUX = iTri(iU,iX)
                    DUX = ExFac*Half*D(iDoff+iUX)
                    if (iX == iV) then
                      DUV = ExFac*Half*DUV
                      DUX = ExFac*Half*DUX
                    end if
                    iPUVX = off_PUVX(iSym)
                    ! inactive/active block
                    do iP=1,iIsh
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPU = iTri(iIsh+iU,iP)
                      F(iFoff+iPU) = F(iFoff+iPU)+DVX*Temp
                      iPV = iTri(iIsh+iV,iP)
                      F(iFoff+iPV) = F(iFoff+iPV)-DUX*Temp
                      iPX = iTri(iIsh+iX,iP)
                      F(iFoff+iPX) = F(iFoff+iPX)-DUV*Temp
                    end do
                    ! active/active block and iP<=(iIsh+iU)
                    do iP=iIsh+1,iIsh+iU
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPU = iTri(iIsh+iU,iP)
                      F(iFoff+iPU) = F(iFoff+iPU)+DVX*Temp
                      iPV = iTri(iIsh+iV,iP)
                      F(iFoff+iPV) = F(iFoff+iPV)-ExFac*Half*DUX*Temp
                      if (iP == (iIsh+iV)) F(iFoff+iPV) = F(iFoff+iPV)-ExFac*Half*DUX*Temp
                      iPX = iTri(iIsh+iX,iP)
                      F(iFoff+iPX) = F(iFoff+iPX)-ExFac*Half*DUV*Temp
                      if (iP == (iIsh+iX)) F(iFoff+iPX) = F(iFoff+iPX)-ExFac*Half*DUV*Temp
                    end do
                    ! active/active block and iP>(iIsh+iU)
                    do iP=iIsh+iU+1,iIsh+iAsh
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPV = iTri(iIsh+iV,iP)
                      F(iFoff+iPV) = F(iFoff+iPV)-ExFac*Half*DUX*Temp
                      if (iP == (iIsh+iV)) F(iFoff+iPV) = F(iFoff+iPV)-ExFac*Half*DUX*Temp
                      iPX = iTri(iIsh+iX,iP)
                      F(iFoff+iPX) = F(iFoff+iPX)-ExFac*Half*DUV*Temp
                      if (iP == (iIsh+iX)) F(iFoff+iPX) = F(iFoff+iPX)-ExFac*Half*DUV*Temp
                    end do
                    ! active/secondary block
                    do iP=iIsh+iAsh+1,iOrb
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPU = iTri(iIsh+iU,iP)
                      F(iFoff+iPU) = F(iFoff+iPU)+DVX*Temp
                      iPV = iTri(iIsh+iV,iP)
                      F(iFoff+iPV) = F(iFoff+iPV)-DUX*Temp
                      iPX = iTri(iIsh+iX,iP)
                      F(iFoff+iPX) = F(iFoff+iPX)-DUV*Temp
                    end do
                    off_PUVX(iSym) = off_PUVX(iSym)+iOrb
                  end do
                end do
              end do

            case (2)
              ! symmetry case (II!KK)
              iFoff = off_Fmat(iSym)
              kDoff = off_Dmat(kSym)
              do iV=1,kAsh
                do iX=1,iV
                  iVX = iTri(iV,iX)
                  DVX = Two*D(kDoff+iVX)
                  if (iX == iV) DVX = D(kDoff+iVX)
                  do iU=1,jAsh
                    iPUVX = off_PUVX(iSym)
                    ! inactive/active block
                    do iP=1,iIsh
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPU = iTri(iIsh+iU,iP)
                      F(iFoff+iPU) = F(iFoff+iPU)+DVX*Temp
                    end do
                    ! active/active block and iP<=(iIsh+iU)
                    do iP=iIsh+1,iIsh+iU
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPU = iTri(iIsh+iU,iP)
                      F(iFoff+iPU) = F(iFoff+iPU)+DVX*Temp
                    end do
                    iPUVX = iPUVX+iAsh-iU
                    ! active/secondary block
                    do iP=iIsh+iAsh+1,iOrb
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPU = iTri(iIsh+iU,iP)
                      F(iFoff+iPU) = F(iFoff+iPU)+DVX*Temp
                    end do
                    off_PUVX(iSym) = off_PUVX(iSym)+iOrb
                  end do
                end do
              end do

            case (3)
              ! symmetry case (IJ!IJ)
              iFoff = off_Fmat(iSym)
              jFoff = off_Fmat(jSym)
              iDoff = off_Dmat(iSym)
              jDoff = off_Dmat(jSym)
              do iV=1,kAsh
                do iX=1,lAsh
                  do iU=1,jAsh
                    iUX = iTri(iU,iX)
                    DUX = ExFac*Half*D(jDoff+iUX)
                    iPUVX = off_PUVX(iSym)
                    ! inactive/active block
                    do iP=1,iIsh
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPV = iTri(iIsh+iV,iP)
                      F(iFoff+iPV) = F(iFoff+iPV)-DUX*Temp
                    end do
                    ! active/active block and iP<=(iIsh+iV)
                    do iP=iIsh+1,iIsh+iV
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPV = iTri(iIsh+iV,iP)
                      F(iFoff+iPV) = F(iFoff+iPV)-DUX*Temp
                    end do
                    iPUVX = iPUVX+iAsh-iV
                    ! active/secondary block
                    do iP=iIsh+iAsh+1,iOrb
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPV = iTri(iIsh+iV,iP)
                      F(iFoff+iPV) = F(iFoff+iPV)-DUX*Temp
                    end do
                    off_PUVX(iSym) = off_PUVX(iSym)+iOrb
                  end do
                  do iU=1,iAsh
                    iUV = iTri(iU,iV)
                    DUV = ExFac*Half*D(iDoff+iUV)
                    iPUVX = off_PUVX(jSym)
                    ! inactive/active block
                    do iP=1,jIsh
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPX = iTri(jIsh+iX,iP)
                      F(jFoff+iPX) = F(jFoff+iPX)-DUV*Temp
                    end do
                    ! active/active block and iP<=(jIsh+iX)
                    do iP=jIsh+1,jIsh+iX
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPX = iTri(jIsh+iX,iP)
                      F(jFoff+iPX) = F(jFoff+iPX)-DUV*Temp
                    end do
                    iPUVX = iPUVX+jAsh-iX
                    ! active/secondary block
                    do iP=jIsh+jAsh+1,jOrb
                      iPUVX = iPUVX+1
                      Temp = PUVX(iPUVX)
                      iPX = iTri(jIsh+iX,iP)
                      F(jFoff+iPX) = F(jFoff+iPX)-DUV*Temp
                    end do
                    off_PUVX(jSym) = off_PUVX(jSym)+jOrb
                  end do
                end do
              end do

            case (4)
              ! symmetry case (IJ!KL)
              do iV=1,kAsh
                do iX=1,lAsh
                  off_PUVX(iSym) = off_PUVX(iSym)+jAsh*iOrb
                  off_PUVX(jSym) = off_PUVX(jSym)+iAsh*jOrb
                end do
              end do

          end select
        end if

      end do
    end do
  end do
end do

end subroutine Upd_FA
