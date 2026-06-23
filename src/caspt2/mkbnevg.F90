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

subroutine MKBNEVG(nAshT,Hact,Gact,G1,G2)

use Index_Functions, only: iTri, nTri_Elem
use caspt2_global, only: LUSBT
use caspt2_module, only: NAES, NASH, NINDEP, NSYM
use EQSOLV, only: IDBMAT
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAshT
real(kind=wp), intent(in) :: Hact(nAshT,nAshT), Gact(nAshT,nAshT,nAshT,nAshT), G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT)
integer(kind=iwp) :: IBG, IDISK, ISYM, iSymU, iSymUV, iSymV, iSymY, IT, ITABS, IU, IUABS, IV, IVABS, IX, IXABS, IY, IYABS, NAS, &
                     NBG, NINM
real(kind=wp) :: Val
real(kind=wp), allocatable :: BG(:)

do ISYM=1,NSYM
  if (NINDEP(ISYM,10) == 0) cycle
  NINM = NINDEP(ISYM,11)
  NAS = NASH(ISYM)
  NBG = nTri_Elem(NAS)
  if (NBG > 0) call mma_Allocate(BG,NBG,LABEL='BG')

  ! Eq.(27): K(t,x) = - h(x,u)*R(t,u) - (uy|xv)*D(tu,vy)
  do IT=1,NAS
    ITABS = IT+NAES(ISYM)
    do IX=1,IT
      IXABS = IX+NAES(ISYM)
      Val = Zero
      !do IU=1,NAS
      !  IUABS = IU+NAES(ISYM)
      !  Val = Val+Hact(IXABS,IUABS)*G1(ITABS,IUABS)
      !  do IV=1,NAS
      !    IVABS = IV+NAES(ISYM)
      !    do IY=1,NAS
      !      IYABS = IY+NAES(ISYM)
      !      Val = Val+Gact(IUABS,IYABS,IXABS,IVABS)*G2(ITABS,IVABS,IUABS,IYABS)
      !    end do
      !  end do
      !end do
      do IU=1,NAS
        IUABS = IU+NAES(iSym)
        Val = Val+Hact(IXABS,IUABS)*G1(ITABS,IUABS)
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
                    Val = Val+Gact(IUABS,IYABS,IXABS,IVABS)*G2(ITABS,IVABS,IUABS,IYABS)
                  end do
                end do
              end do
            end if
          end do
        end do
      end do
      IBG = iTri(IT,IX)
      BG(IBG) = -Val
      !write(u6,'("BG",2i3,f20.10)') it,ix,bg(ibg)
    end do
  end do

  if (NBG > 0) then
    if (NINDEP(ISYM,10) > 0) then
      IDISK = IDBMAT(ISYM,10)
      call DDAFILE(LUSBT,1,BG,NBG,IDISK)
    end if
    if ((NINM > 0) .and. (NINDEP(ISYM,11) > 0)) then
      IDISK = IDBMAT(ISYM,11)
      call DDAFILE(LUSBT,1,BG,NBG,IDISK)
    end if
    call mma_deallocate(BG)
  end if
end do

end subroutine MKBNEVG
