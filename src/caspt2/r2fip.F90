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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************
      Subroutine R2FIP(CHSPC,NCHSPC,WRK,ipWRK,NUMV,nBasT,iSym,          &
     &                iSkip,irc,JREDC)

      Use Cholesky, only: INFVEC, nDimRS, nnBstR
      use definitions, only: wp, iwp
      use Constants, only: Zero

      implicit none

      integer(kind=iwp), intent(in) :: NCHSPC, ipWRK(8), NUMV, nBasT,   &
     &  iSym, iSkip(8)
      real(kind=wp), intent(inout) :: CHSPC(NCHSPC), WRK(nBasT*nBasT)
      integer(kind=iwp), intent(inout) :: irc, JREDC

      integer(kind=iwp) :: kloc, iVec, lscr, JREDL, ipVecL, jloc,       &
     &  l_NDIMRS
!
!     Transform the reduced form to the full form in place
!
      l_NDIMRS = size(nDIMRS)
      kloc = 0
      Do iVec = 1, NUMV
        If (l_NDIMRS < 1) Then
          lscr  = NNBSTR(iSym,3)
        Else
          JREDL = INFVEC(iVec,2,iSym)
          lscr  = nDimRS(iSym,JREDL) !! JRED?
        End If
        kloc = kloc + lscr
      End Do

      ipVecL = 1 + kloc !! lscr*(JNUM-1)
      jloc = (NUMV-1)*nBasT**2+1
      Do iVec = NUMV, 1, -1
        If (l_NDIMRS < 1) Then
          lscr  = NNBSTR(iSym,3)
        Else
          JREDL = INFVEC(iVec,2,iSym)
          lscr  = nDimRS(iSym,JREDL) !! JRED?
        End If
        ipVecL = ipVecL - lscr
        WRK(1:nBasT**2) = Zero
        Call Cho_ReOrdr(irc,CHSPC(ipVecL),lscr,1,                       &
     &                  1,1,1,iSym,JREDC,2,ipWRK,WRK,                   &
     &                  iSkip)
        CHSPC(jloc:jloc+nBasT**2-1) = WRK(1:nBasT**2)
        jloc = jloc-nBasT**2
      End Do

      Return

      End Subroutine R2FIP
