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

subroutine MKBNEVE(nAshT,Hact,Gact,G1,G2)

use Index_Functions, only: iTri, nTri_Elem
use caspt2_global, only: LUSBT
use general_data, only: NASH
use caspt2_module, only: NAES, NINDEP, NSYM
use EQSOLV, only: IDBMAT
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAshT
real(kind=wp), intent(in) :: Hact(nAshT,nAshT), Gact(nAshT,nAshT,nAshT,nAshT), G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT)
integer(kind=iwp) :: IBE, IDISK, ISYM, iSymU, iSymUV, iSymV, iSymY, IT, ITABS, IU, IUABS, IV, IVABS, IX, IXABS, IY, IYABS, NAS, &
                     NBE, NINM
real(kind=wp) :: Val
real(kind=wp), allocatable :: BE(:)

do ISYM=1,NSYM
  if (NINDEP(ISYM,6) == 0) cycle
  NINM = NINDEP(ISYM,7)
  NAS = NASH(ISYM)
  NBE = nTri_Elem(NAS)
  if (NBE > 0) call mma_Allocate(BE,NBE,LABEL='BE')

  ! Eq.(A3): K(t,x) = h(u,x)*tR(u,t) + (ce|da)*(tE(td)E(ce))
  !                 = h(u,x)*(2del(u,t) - R(t,u)) + (uv|yx)*(2*del(ty)*D(uv)-D(yu,tv)-del(tu)D(yv))
  ! Or, considering using Eq.(A4)
  do IT=1,NAS
    ITABS = IT+NAES(ISYM)
    do IX=1,IT
      IXABS = IX+NAES(ISYM)
      Val = Zero
      !do IU=1,NAS
      !  IUABS = IU+NAES(ISYM)
      !  Val = Val-Hact(IUABS,IXABS)*G1(ITABS,IUABS)
      !  if (ITABS == IUABS) Val = Val+Two*Hact(IUABS,IXABS)
      !  do IV=1,NAS
      !    IVABS = IV+NAES(ISYM)
      !    do IY=1,NAS
      !      IYABS = IY+NAES(ISYM)
      !      Val = Val-Gact(IUABS,IVABS,IYABS,IXABS)*G2(IYABS,ITABS,IUABS,IVABS)
      !      if (ITABS == IYABS) Val = Val+Two*Gact(IUABS,IVABS,IYABS,IXABS)*G1(IUABS,IVABS)
      !      if (ITABS == IUABS) Val = Val-Gact(IUABS,IVABS,IYABS,IXABS)*G1(IYABS,IVABS)
      !    end do
      !  end do
      !end do
      do IU=1,NAS
        IUABS = IU+NAES(iSym)
        Val = Val-Hact(IUABS,IXABS)*G1(ITABS,IUABS)
        if (ITABS == IUABS) Val = Val+Two*Hact(IUABS,IXABS)
      end do
      do iSymU=1,nSym
        do iSymV=1,nSym
          iSymUV = MUL(iSymU,iSymV)
          do iSymY=1,nSym
            if (iSym == MUL(iSymY,iSymUV)) then
              do IU=1,NASH(iSymU)
                IUABS = IU+NAES(iSymU)
                do IV=1,NASH(iSymV)
                  IVABS = IV+NAES(iSymV)
                  do IY=1,NASH(iSymY)
                    IYABS = IY+NAES(iSymY)
                    Val = Val-Gact(IUABS,IVABS,IYABS,IXABS)*G2(IYABS,ITABS,IUABS,IVABS)
                    if (ITABS == IYABS) Val = Val+Two*Gact(IUABS,IVABS,IYABS,IXABS)*G1(IUABS,IVABS)
                    if (ITABS == IUABS) Val = Val-Gact(IUABS,IVABS,IYABS,IXABS)*G1(IYABS,IVABS)
                  end do
                end do
              end do
            end if
          end do
        end do
      end do
      IBE = iTri(IT,IX)
      BE(IBE) = Val
    end do
  end do

  if ((NBE > 0) .and. (NINDEP(ISYM,6) > 0)) then
    IDISK = IDBMAT(ISYM,6)
    call DDAFILE(LUSBT,1,BE,NBE,IDISK)
    if ((NINM > 0) .and. (NINDEP(ISYM,7) > 0)) then
      IDISK = IDBMAT(ISYM,7)
      call DDAFILE(LUSBT,1,BE,NBE,IDISK)
    end if
    call mma_deallocate(BE)
  end if
end do

end subroutine MKBNEVE
