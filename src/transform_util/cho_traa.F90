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
! Copyright (C) 2004,2005, Giovanni Ghigo                              *
!               2021, Roland Lindh                                     *
!***********************************************************************
!  Cho_TraA
!
!> @brief
!>   Routine for the transformation of the Cholesky vectors in MO-based TCVx for case \c Sym(i) &ne; \c Sym(j)
!> @author Giovanni Ghigo
!>
!> @details
!> In the inner batch the Cholesky Full Vectors are transformed and
!> stored in memory. Adresses (``1``) and length (``2``) are stored in the
!> matrix <tt>TCVX(iType, Sym(i), Sym(j)) of allocatable 2D arrays</tt>.
!>
!> - \c iType = ``1``: TCVA
!> - \c iType = ``2``: TCVB
!> - \c iType = ``3``: TCVC
!> - \c iType = ``4``: TCVD
!> - \c iType = ``5``: TCVE
!> - \c iType = ``6``: TCVF
!> - \c iType = ``7``: TCVG
!>
!> Types ``1``, ``2`` and ``4``--``7`` are generated only if \c DoTCVA = ``.True.``
!> TCVC is always generated.
!>
!> In the first half-transformation the vectors are contracted
!> only with the occupied (inactive and active) MO coefficients
!> for \c Sym(j). In the second half-transformation the vectors are
!> contracted with all MO coefficients.
!>
!> @note
!> The logical matrix \c TCVXist must be defined.
!>
!> @param[in] iSym        Symmetry(``i``) of the Cholesky Full Vector
!> @param[in] jSym        Symmetry(``j``) of the Cholesky Full Vector
!> @param[in] NumV        Number of Cholesky vectors to transform in the current batch
!> @param[in] CMO         MO coefficients
!> @param[in] NCMO        Total number of MO coefficients
!> @param[in] lUCHFV      Unit number of the Cholesky Full Vector to transform (``CHFV``)
!> @param[in] iStrtVec_AB Current initial disk pointer of the Cholesky Full Vector to transform (``CHFV``)
!> @param[in] nFVec       Number of Cholesky vectors to transform in the inner batch procedure
!***********************************************************************

subroutine Cho_TraA(iSym,jSym,NumV,CMO,NCMO,lUCHFV,iStrtVec_AB,nFVec)
!***********************************************************************
!  This is the routine for the transformation from AO basis to MO      *
!  basis of the Cholesky Full Vectors when  iSym /= jSym.              *
!  The new Transformed Cholesky Full Vectors L are :                   *
!                                                                      *
!  IF (.not. DoFull) (CASPT2, MCLR,MBPT2)                              *
!   TCVA : L_ij, L_ji=T(L_ij) if DoTCVA=.True.                         *
!   TCVB : L_tj, L_ui         if DoTCVA=.True.                         *
!   TCVC : L_aj, L_bi                                                  *
!   TCVD : L_tu, L_ut=T(L_tu) if DoTCVA=.True.                         *
!   TCVE : L_au, L_bt         if DoTCVA=.True.                         *
!   TCVF : L_ab               if DoTCVA=.True.                         *
!   TCVG : L_jt, L_iu         if DoTCVA=.True.                         *
!  For generation of <pk|ql>  p,q: All MO, k,l: Occupied (i & t)       *
!  MO Indices  i,j: Inactive;   t,u: Active;   a,b: Secondary          *
!                                                                      *
!----------------------------------------------------------------------*
!  Author  :  Giovanni Ghigo                                           *
!             Lund University, Sweden & University di Torino, Italia   *
!  Written :  October 2004                                             *
!  Modified:  July 2005                                                *
!***********************************************************************

use Cho_Tra, only: nAsh, nBas, nFro, nIsh, nSsh, TCVX, TCVXist
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSym, jSym, NumV, NCMO, lUCHFV, iStrtVec_AB, nFVec
real(kind=wp), intent(in) :: CMO(NCMO)
integer(kind=iwp) :: i, iFBatch, iiVec, iStrt, iStrt0MO, iStrtVec_FAB, iVec, j, jStrt, jStrt0MO, jVec, Len_XAb, Len_XAj, Len_XAu, &
                     Len_XBi, Len_XBt, Nab, Naj, Nau, Nbi, Nbt, NFAB, Nij, Niu, Nji, Njt, Ntj, Ntu, Nui, NumFV, Nut
logical(kind=iwp) :: TCVA, TCVB, TCVBt, TCVC, TCVCt, TCVD, TCVE, TCVEt, TCVF
real(kind=wp), allocatable :: FAB(:,:), XAb(:), XAj(:), XAu(:), XBi(:), XBt(:)

! Memory to allocate & Nr. of Cholesky vectors transformable
! A=Alpha(AO);  B=Beta(AO)
TCVA = .false.
TCVB = .false.
TCVBt = .false.
TCVC = .false.
TCVCt = .false.
TCVD = .false.
TCVE = .false.
TCVEt = .false.
TCVF = .false.

NFAB = 0
Nij = 0  ! A
Nji = 0  ! A"
Ntj = 0  ! B
Njt = 0  ! G
Nui = 0  ! B"
Niu = 0  ! G"
Naj = 0  ! C
Nbi = 0  ! C"
Ntu = 0  ! D
Nut = 0  ! D"
Nau = 0  ! E
Nbt = 0  ! E"
Nab = 0  ! F

Len_XAj = 0   ! A/A", B, C
Len_XAu = 0   ! D/D", E
Len_XAb = 0   ! F
Len_XBi = 0   ! B",C"
Len_XBt = 0   ! E"

NFAB = nBas(iSym)*nBas(jSym)

! Allocate memory for Transformed Cholesky Vectors - TCVx
! TCV-A :
if (TCVXist(1,iSym,jSym)) then
  TCVA = .true.
  Len_XAj = nBas(iSym)*nIsh(jSym)
  Nij = nIsh(iSym)*nIsh(jSym)
  Nji = nIsh(jSym)*nIsh(iSym)
  call mma_allocate(TCVX(1,iSym,jSym)%A,Nij,NumV,Label='TCVA')
  call mma_allocate(TCVX(1,jSym,iSym)%A,Nji,NumV,Label='TCVA')
end if

! TCV-B :
if (TCVXist(2,iSym,jSym)) then
  TCVB = .true.
  Len_XAj = nBas(iSym)*nIsh(jSym)
  Ntj = nAsh(iSym)*nIsh(jSym)
  Njt = nIsh(jSym)*nAsh(iSym)
  call mma_allocate(TCVX(2,iSym,jSym)%A,Ntj,NumV,Label='TCVB')
  call mma_allocate(TCVX(7,jSym,iSym)%A,Njt,NumV,Label='TCVB')
end if
if (TCVXist(2,jSym,iSym)) then
  TCVBt = .true.
  Len_XBi = nBas(jSym)*nIsh(iSym)
  Nui = nAsh(jSym)*nIsh(iSym)
  Niu = nIsh(iSym)*nAsh(jSym)
  call mma_allocate(TCVX(2,jSym,iSym)%A,Nui,NumV,Label='TCVBt')
  call mma_allocate(TCVX(7,iSym,jSym)%A,Niu,NumV,Label='TCVBt')
end if

! TCV-C :
if (TCVXist(3,iSym,jSym)) then
  TCVC = .true.
  Len_XAj = nBas(iSym)*nIsh(jSym)
  Naj = nSsh(iSym)*nIsh(jSym)
  call mma_allocate(TCVX(3,iSym,jSym)%A,Naj,NumV,Label='TCVC')
end if
if (TCVXist(3,jSym,iSym)) then
  TCVCt = .true.
  Len_XBi = nBas(jSym)*nIsh(iSym)
  Nbi = nSsh(jSym)*nIsh(iSym)
  call mma_allocate(TCVX(3,jSym,iSym)%A,Nbi,NumV,Label='TCVCt')
end if

! TCV-D :
if (TCVXist(4,iSym,jSym)) then
  TCVD = .true.
  Len_XAu = nBas(iSym)*nAsh(jSym)
  Ntu = nAsh(iSym)*nAsh(jSym)
  Nut = nAsh(jSym)*nAsh(iSym)
  call mma_allocate(TCVX(4,iSym,jSym)%A,Ntu,NumV,Label='TCVD')
  call mma_allocate(TCVX(4,jSym,iSym)%A,Nut,NumV,Label='TCVD')
end if

! TCV-E :
if (TCVXist(5,iSym,jSym)) then
  TCVE = .true.
  Len_XAu = nBas(iSym)*nAsh(jSym)
  Nau = nSsh(iSym)*nAsh(jSym)
  call mma_allocate(TCVX(5,iSym,jSym)%A,Nau,NumV,Label='TCVE')
end if
if (TCVXist(5,jSym,iSym)) then
  TCVEt = .true.
  Len_XBt = nBas(jSym)*nAsh(iSym)
  Nbt = nSsh(jSym)*nAsh(iSym)
  call mma_allocate(TCVX(5,jSym,iSym)%A,Nbt,NumV,Label='TCVEt')
end if

! TCV-F :
if (TCVXist(6,iSym,jSym)) then
  TCVF = .true.
  Len_XAb = nBas(iSym)*nSsh(jSym)
  Nab = nSsh(iSym)*nSsh(jSym)
  call mma_allocate(TCVX(6,jSym,iSym)%A,Nab,NumV,Label='TCVF')
end if

iStrt = 1
do i=1,iSym-1
  iStrt = iStrt+nBas(i)*nBas(i)
end do

jStrt = 1
do j=1,jSym-1
  jStrt = jStrt+nBas(j)*nBas(j)
end do

! START LOOP iiVec
do iiVec=1,NumV,nFVec
  NumFV = max(nFVec,NumV-iiVec+1)
  iFBatch = (iiVec+nFVec-1)/nFVec

  ! Allocate memory & Load Full Cholesky Vectors - CHFV

  iStrtVec_FAB = iStrtVec_AB+nFVec*(iFBatch-1)

  call mma_allocate(FAB,NFAB,NumFV,Label='FAB')
  call RdChoVec(FAB,NFAB,NumFV,iStrtVec_FAB,lUCHFV)

  ! Start Loop jVec
  do jVec=iiVec,iiVec+NumFV-1   ! Loop  jVec
    iVec = jVec-iiVec+1

    ! 1st Half-Transformation  iBeta(AO) -> q(MO) only occupied
    jStrt0MO = jStrt+nFro(jSym)*nBas(jSym)

    ! From CHFV A(Alpha,Beta) to XAj(Alpha,jMO)
    if (TCVA .or. TCVB .or. TCVC) then
      call mma_allocate(XAj,Len_XAj,Label='XAj')
      call ProdsA_2(FAB(:,iVec),nBas(iSym),nBas(jSym),CMO(jStrt0MO),nIsh(jSym),XAj)
    end if
    jStrt0MO = jStrt0MO+nIsh(jSym)*nBas(jSym)

    ! From CHFV A(Alpha,Beta) to XAu(Alpha,uMO)
    if (TCVD .or. TCVE) then
      call mma_allocate(XAu,Len_XAu,Label='XAu')
      call ProdsA_2(FAB(:,iVec),nBas(iSym),nBas(jSym),CMO(jStrt0MO),nAsh(jSym),XAu)
    end if
    jStrt0MO = jStrt0MO+nAsh(jSym)*nBas(jSym)

    ! From CHFV A(Alpha,Beta) to XAb(Alpha,bMO)
    if (TCVF) then
      call mma_allocate(XAb,Len_XAb,Label='XAb')
      call ProdsA_2(FAB(:,iVec),nBas(iSym),nBas(jSym),CMO(jStrt0MO),nSsh(jSym),XAb)
    end if

    iStrt0MO = iStrt+nFro(iSym)*nBas(iSym)

    ! From CHFV A(Alpha,Beta) to XBi(Beta,iMO)
    if (TCVBt .or. TCVCt) then
      call mma_allocate(XBi,Len_XBi,Label='XBi')
      call ProdsA_2t(FAB(:,iVec),nBas(iSym),nBas(jSym),CMO(iStrt0MO),nIsh(iSym),XBi)
    end if
    iStrt0MO = iStrt0MO+nIsh(iSym)*nBas(iSym)

    ! From CHFV A(Alpha,Beta) to XBt(Beta,tMO)
    if (TCVEt) then
      call mma_allocate(XBt,Len_XBt,Label='XBt')
      call ProdsA_2t(FAB(:,iVec),nBas(iSym),nBas(jSym),CMO(iStrt0MO),nAsh(iSym),XBt)
    end if

    ! 2nd Half-Transformation  iAlpha(AO) -> p(MO)
    iStrt0MO = iStrt+nFro(iSym)*nBas(iSym)

    ! From XAj(Alpha,jMO) to ij(i,j)
    if (TCVA) then
      call ProdsA_1(XAj,nBas(iSym),nIsh(jSym),CMO(iStrt0MO),nIsh(iSym),TCVX(1,iSym,jSym)%A(:,jVec))
      call Trnsps(nIsh(iSym),nIsh(jSym),TCVX(1,iSym,jSym)%A(:,jVec),TCVX(1,jSym,iSym)%A(:,jVec))
    end if
    iStrt0MO = iStrt0MO+nIsh(iSym)*nBas(iSym)

    ! From XAj(Alpha,jMO) to tj(t,j)
    if (TCVB) then
      call ProdsA_1(XAj,nBas(iSym),nIsh(jSym),CMO(iStrt0MO),nAsh(iSym),TCVX(2,iSym,jSym)%A(:,jVec))
      call Trnsps(nAsh(iSym),nIsh(jSym),TCVX(2,iSym,jSym)%A(:,jVec),TCVX(7,jSym,iSym)%A(:,jVec))
    end if

    ! From XAu(Alpha,jMO) to tu(t,u)
    if (TCVD) then
      call ProdsA_1(XAu,nBas(iSym),nAsh(jSym),CMO(iStrt0MO),nAsh(iSym),TCVX(4,iSym,jSym)%A(:,jVec))
      call Trnsps(nAsh(iSym),nAsh(jSym),TCVX(4,iSym,jSym)%A(:,jVec),TCVX(4,jSym,iSym)%A(:,jVec))
    end if
    iStrt0MO = iStrt0MO+nAsh(iSym)*nBas(iSym)

    ! From XAj(Alpha,jMO) to aj(a,j)
    if (TCVC) then
      call ProdsA_1(XAj,nBas(iSym),nIsh(jSym),CMO(iStrt0MO),nSsh(iSym),TCVX(3,iSym,jSym)%A(:,jVec))
    end if

    ! From XAu(Alpha,jMO) to au(a,u)
    if (TCVE) then
      call ProdsA_1(XAu,nBas(iSym),nAsh(jSym),CMO(iStrt0MO),nSsh(iSym),TCVX(5,iSym,jSym)%A(:,jVec))
    end if

    ! From XAb(Alpha,jMO) to ab(a,b)
    if (TCVF) then
      call ProdsA_1(XAb,nBas(iSym),nSsh(jSym),CMO(iStrt0MO),nSsh(iSym),TCVX(6,jSym,iSym)%A(:,jVec))
    end if

    jStrt0MO = jStrt+(nFro(jSym)+nIsh(jSym))*nBas(jSym)

    ! From XBi(Beta,jMO) to ui(u,i)
    if (TCVBt) then
      call ProdsA_1(XBi,nBas(jSym),nIsh(iSym),CMO(jStrt0MO),nAsh(jSym),TCVX(2,jSym,iSym)%A(:,jVec))
      call Trnsps(nAsh(jSym),nIsh(iSym),TCVX(2,jSym,iSym)%A(:,jVec),TCVX(7,iSym,jSym)%A(:,jVec))
    end if
    jStrt0MO = jStrt0MO+nAsh(jSym)*nBas(jSym)

    ! From XBi(Beta,jMO) to bi(b,i)
    if (TCVCt) then
      call ProdsA_1(XBi,nBas(jSym),nIsh(iSym),CMO(jStrt0MO),nSsh(jSym),TCVX(3,jSym,iSym)%A(:,jVec))
    end if

    ! From XBt(Beta,jMO) to bt(b,t)
    if (TCVEt) then
      call ProdsA_1(XBt,nBas(jSym),nAsh(iSym),CMO(jStrt0MO),nSsh(jSym),TCVX(5,jSym,iSym)%A(:,jVec))
    end if

    ! End of Transformations

    if (allocated(XAj)) call mma_deallocate(XAj)
    if (allocated(XAu)) call mma_deallocate(XAu)
    if (allocated(XAb)) call mma_deallocate(XAb)
    if (allocated(XBi)) call mma_deallocate(XBi)
    if (allocated(XBt)) call mma_deallocate(XBt)

  end do
  ! End Loop  jVec

  call mma_deallocate(FAB)
end do
! END LOOP iiVec

return

end subroutine Cho_TraA
