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

use general_data, only: NSYM, NASH, NISH, NORB
use Molcas, only: MxSym
use Constants, only: Zero, Two, Half

implicit none
real*8 PUVX(*), F(*), D(*)
integer case
integer off_PUVX(mxSym), off_Dmat(mxSym), off_Fmat(mxSym)
integer ISTACK, ISYM, IASH, IORB, JSYM, JASH, IJSYM, KSYM, KASH, LSYM, LASH, KLSYM, KL_ORB_PAIRS, IISH, IFOFF, IU, IP, IPU, IPUVX, &
        JORB, JISH, IDOFF, IV, IX, IVX, IUV, IUX, IPV, IPX, KDOFF, JFOFF, JDOFF
real*8 DVX, DUV, EXFAC, DUX, TEMP
! Statement function
integer i, iTri
iTri(i) = (i*i-i)/2

! generate offsets

iStack = 0
do iSym=1,nSym
  off_Dmat(iSym) = iStack
  iAsh = nAsh(iSym)
  iStack = iStack+(iAsh*iAsh+iAsh)/2
end do

iStack = 0
do iSym=1,nSym
  off_Fmat(iSym) = iStack
  iOrb = nOrb(iSym)
  iStack = iStack+(iOrb*iOrb+iOrb)/2
end do

iStack = 0
do iSym=1,nSym
  off_PUVX(iSym) = iStack
  iOrb = nOrb(iSym)
  do jSym=1,nSym
    jAsh = nAsh(jSym)
    ijSym = 1+ieor(iSym-1,jSym-1)
    do kSym=1,nSym
      kAsh = nAsh(kSym)
      do lSym=1,kSym
        lAsh = nAsh(lSym)
        klSym = 1+ieor(kSym-1,lSym-1)
        if (ijSym == klSym) then
          kl_Orb_pairs = kAsh*lAsh
          if (kSym == lSym) kl_Orb_pairs = (kAsh*kAsh+kAsh)/2
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
      iPU = iP+iTri(iU)
      F(iFoff+iPU) = Zero
    end do
  end do
  do iU=iIsh+iAsh+1,iOrb
    do iP=iIsh+1,iIsh+iAsh
      iPU = iP+iTri(iU)
      F(iFoff+iPU) = Zero
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
    ijSym = 1+ieor(iSym-1,jSym-1)
    do kSym=1,nSym
      kAsh = nAsh(kSym)
      do lSym=1,kSym
        lAsh = nAsh(lSym)
        klSym = 1+ieor(kSym-1,lSym-1)

        ! find cases
        case = 4
        if (iSym == jSym) case = case-2
        if (iSym == kSym) case = case-1

        if ((ijSym == klSym) .and. (iAsh*jAsh*kAsh*lAsh /= 0)) then

          goto(100,200,300,400) case

          ! symmetry case (II!II)
100       continue
          iFoff = off_Fmat(iSym)
          iDoff = off_Dmat(iSym)
          do iV=1,kAsh
            do iX=1,iV
              iVX = iTri(iV)+iX
              DVX = Two*D(iDoff+iVX)
              if (iX == iV) DVX = D(iDoff+iVX)
              do iU=1,jAsh
                iUV = iTri(iU)+iV
                if (iV > iU) iUV = iTri(iV)+iU
                DUV = ExFac*Half*D(iDoff+iUV)
                iUX = iTri(iU)+iX
                if (iX > iU) iUX = iTri(iX)+iU
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
                  iPU = iTri(iIsh+iU)+iP
                  F(iFoff+iPU) = F(iFoff+iPU)+DVX*Temp
                  iPV = iTri(iIsh+iV)+iP
                  F(iFoff+iPV) = F(iFoff+iPV)-DUX*Temp
                  iPX = iTri(iIsh+iX)+iP
                  F(iFoff+iPX) = F(iFoff+iPX)-DUV*Temp
                end do
                ! active/active block and iP<=(iIsh+iU)
                do iP=iIsh+1,iIsh+iU
                  iPUVX = iPUVX+1
                  Temp = PUVX(iPUVX)
                  iPU = iTri(iIsh+iU)+iP
                  F(iFoff+iPU) = F(iFoff+iPU)+DVX*Temp
                  iPV = iTri(iIsh+iV)+iP
                  if (iP > (iIsh+iV)) iPV = iTri(iP)+iIsh+iV
                  F(iFoff+iPV) = F(iFoff+iPV)-ExFac*Half*DUX*Temp
                  if (iP == (iIsh+iV)) F(iFoff+iPV) = F(iFoff+iPV)-ExFac*Half*DUX*Temp
                  iPX = iTri(iIsh+iX)+iP
                  if (iP > (iIsh+iX)) iPX = iTri(iP)+iIsh+iX
                  F(iFoff+iPX) = F(iFoff+iPX)-ExFac*Half*DUV*Temp
                  if (iP == (iIsh+iX)) F(iFoff+iPX) = F(iFoff+iPX)-ExFac*Half*DUV*Temp
                end do
                ! active/active block and iP>(iIsh+iU)
                do iP=iIsh+iU+1,iIsh+iAsh
                  iPUVX = iPUVX+1
                  Temp = PUVX(iPUVX)
                  iPV = iTri(iIsh+iV)+iP
                  if (iP > (iIsh+iV)) iPV = iTri(iP)+iIsh+iV
                  F(iFoff+iPV) = F(iFoff+iPV)-ExFac*Half*DUX*Temp
                  if (iP == (iIsh+iV)) F(iFoff+iPV) = F(iFoff+iPV)-ExFac*Half*DUX*Temp
                  iPX = iTri(iIsh+iX)+iP
                  if (iP > (iIsh+iX)) iPX = iTri(iP)+iIsh+iX
                  F(iFoff+iPX) = F(iFoff+iPX)-ExFac*Half*DUV*Temp
                  if (iP == (iIsh+iX)) F(iFoff+iPX) = F(iFoff+iPX)-ExFac*Half*DUV*Temp
                end do
                ! active/secondary block
                do iP=iIsh+iAsh+1,iOrb
                  iPUVX = iPUVX+1
                  Temp = PUVX(iPUVX)
                  iPU = iTri(iP)+iIsh+iU
                  F(iFoff+iPU) = F(iFoff+iPU)+DVX*Temp
                  iPV = iTri(iP)+iIsh+iV
                  F(iFoff+iPV) = F(iFoff+iPV)-DUX*Temp
                  iPX = iTri(iP)+iIsh+iX
                  F(iFoff+iPX) = F(iFoff+iPX)-DUV*Temp
                end do
                off_PUVX(iSym) = off_PUVX(iSym)+iOrb
              end do
            end do
          end do
          goto 500

          ! symmetry case (II!KK)
200       continue
          iFoff = off_Fmat(iSym)
          kDoff = off_Dmat(kSym)
          do iV=1,kAsh
            do iX=1,iV
              iVX = iTri(iV)+iX
              DVX = Two*D(kDoff+iVX)
              if (iX == iV) DVX = D(kDoff+iVX)
              do iU=1,jAsh
                iPUVX = off_PUVX(iSym)
                ! inactive/active block
                do iP=1,iIsh
                  iPUVX = iPUVX+1
                  Temp = PUVX(iPUVX)
                  iPU = iTri(iIsh+iU)+iP
                  F(iFoff+iPU) = F(iFoff+iPU)+DVX*Temp
                end do
                ! active/active block and iP<=(iIsh+iU)
                do iP=iIsh+1,iIsh+iU
                  iPUVX = iPUVX+1
                  Temp = PUVX(iPUVX)
                  iPU = iTri(iIsh+iU)+iP
                  F(iFoff+iPU) = F(iFoff+iPU)+DVX*Temp
                end do
                iPUVX = iPUVX+iAsh-iU
                ! active/secondary block
                do iP=iIsh+iAsh+1,iOrb
                  iPUVX = iPUVX+1
                  Temp = PUVX(iPUVX)
                  iPU = iTri(iP)+iIsh+iU
                  F(iFoff+iPU) = F(iFoff+iPU)+DVX*Temp
                end do
                off_PUVX(iSym) = off_PUVX(iSym)+iOrb
              end do
            end do
          end do
          goto 500

          ! symmetry case (IJ!IJ)
300       continue
          iFoff = off_Fmat(iSym)
          jFoff = off_Fmat(jSym)
          iDoff = off_Dmat(iSym)
          jDoff = off_Dmat(jSym)
          do iV=1,kAsh
            do iX=1,lAsh
              do iU=1,jAsh
                iUX = iTri(iU)+iX
                if (iX > iU) iUX = iTri(iX)+iU
                DUX = ExFac*Half*D(jDoff+iUX)
                iPUVX = off_PUVX(iSym)
                ! inactive/active block
                do iP=1,iIsh
                  iPUVX = iPUVX+1
                  Temp = PUVX(iPUVX)
                  iPV = iTri(iIsh+iV)+iP
                  F(iFoff+iPV) = F(iFoff+iPV)-DUX*Temp
                end do
                ! active/active block and iP<=(iIsh+iV)
                do iP=iIsh+1,iIsh+iV
                  iPUVX = iPUVX+1
                  Temp = PUVX(iPUVX)
                  iPV = iTri(iIsh+iV)+iP
                  F(iFoff+iPV) = F(iFoff+iPV)-DUX*Temp
                end do
                iPUVX = iPUVX+iAsh-iV
                ! active/secondary block
                do iP=iIsh+iAsh+1,iOrb
                  iPUVX = iPUVX+1
                  Temp = PUVX(iPUVX)
                  iPV = iTri(iP)+iIsh+iV
                  F(iFoff+iPV) = F(iFoff+iPV)-DUX*Temp
                end do
                off_PUVX(iSym) = off_PUVX(iSym)+iOrb
              end do
              do iU=1,iAsh
                iUV = iTri(iU)+iV
                if (iV > iU) iUV = iTri(iV)+iU
                DUV = ExFac*Half*D(iDoff+iUV)
                iPUVX = off_PUVX(jSym)
                ! inactive/active block
                do iP=1,jIsh
                  iPUVX = iPUVX+1
                  Temp = PUVX(iPUVX)
                  iPX = iTri(jIsh+iX)+iP
                  F(jFoff+iPX) = F(jFoff+iPX)-DUV*Temp
                end do
                ! active/active block and iP<=(jIsh+iX)
                do iP=jIsh+1,jIsh+iX
                  iPUVX = iPUVX+1
                  Temp = PUVX(iPUVX)
                  iPX = iTri(jIsh+iX)+iP
                  F(jFoff+iPX) = F(jFoff+iPX)-DUV*Temp
                end do
                iPUVX = iPUVX+jAsh-iX
                ! active/secondary block
                do iP=jIsh+jAsh+1,jOrb
                  iPUVX = iPUVX+1
                  Temp = PUVX(iPUVX)
                  iPX = iTri(iP)+jIsh+iX
                  F(jFoff+iPX) = F(jFoff+iPX)-DUV*Temp
                end do
                off_PUVX(jSym) = off_PUVX(jSym)+jOrb
              end do
            end do
          end do
          goto 500

          ! symmetry case (IJ!KL)
400       continue
          do iV=1,kAsh
            do iX=1,lAsh
              off_PUVX(iSym) = off_PUVX(iSym)+jAsh*iOrb
              off_PUVX(jSym) = off_PUVX(jSym)+iAsh*jOrb
            end do
          end do

500       continue
        end if

      end do
    end do
  end do
end do

end subroutine Upd_FA
