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
  use caspt2_global, only: LUSBT
  use caspt2_module, only: NASH, NAES, NINDEP, NSYM
  use EQSOLV, only: IDBMAT
  use stdalloc, only: mma_allocate, mma_deallocate
  use definitions, only: iwp,wp
  use Constants, only: Zero, Two
  use Symmetry_Info, only: Mul

  implicit none

  integer(kind=iwp), intent(in) :: nAshT
  real(kind=wp), intent(in) :: Hact(nAshT,nAshT), Gact(nAshT,nAshT,nAshT,nAshT), &
                               G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT)

  integer(kind=iwp) :: IBE, IDISK, ISYM, IT, ITABS, IU, IUABS, IV, IVABS, IX, IXABS, IY, IYABS, &
                       NAS, NBE, NINM
  integer(kind=iwp) :: iSymU, iSymUV, iSymV, iSymY
  real(kind=wp) :: VALUE

  real(kind=wp), allocatable :: BE(:)

  DO ISYM=1,NSYM
    IF(NINDEP(ISYM,6) == 0) cycle
    NINM=NINDEP(ISYM,7)
    NAS=NASH(ISYM)
    NBE=(NAS*(NAS+1))/2
    if (NBE > 0) then
      CALL mma_Allocate(BE,NBE,LABEL='BE')
    end if

    ! Eq.(A3): K(t,x) = h(u,x)*tR(u,t) + (ce|da)*(tE(td)E(ce))
    !                 = h(u,x)*(2del(u,t) - R(t,u)) + (uv|yx)*(2*del(ty)*D(uv)-D(yu,tv)-del(tu)D(yv))
    ! Or, considering using Eq.(A4)
    DO IT=1,NAS
      ITABS=IT+NAES(ISYM)
      DO IX=1,IT
        IXABS=IX+NAES(ISYM)
        VALUE=Zero
!       DO IU=1,NAS
!         IUABS=IU+NAES(ISYM)
!         VALUE = VALUE - Hact(IUABS,IXABS)*G1(ITABS,IUABS)
!         if (ITABS == IUABS) VALUE = VALUE + Two*Hact(IUABS,IXABS)
!         DO IV=1,NAS
!           IVABS=IV+NAES(ISYM)
!           DO IY=1,NAS
!             IYABS=IY+NAES(ISYM)
!             VALUE = VALUE - Gact(IUABS,IVABS,IYABS,IXABS)*G2(IYABS,ITABS,IUABS,IVABS)
!             if (ITABS == IYABS) VALUE = VALUE + Two*Gact(IUABS,IVABS,IYABS,IXABS)*G1(IUABS,IVABS)
!             if (ITABS == IUABS) VALUE = VALUE - Gact(IUABS,IVABS,IYABS,IXABS)*G1(IYABS,IVABS)
!           END DO
!         END DO
!       END DO
        do IU = 1, NAS
          IUABS = IU + NAES(iSym)
          VALUE = VALUE - Hact(IUABS,IXABS)*G1(ITABS,IUABS)
          if (ITABS == IUABS) VALUE = VALUE + Two*Hact(IUABS,IXABS)
        end do
        do iSymU = 1, nSym
          do iSymV = 1, nSym
            iSymUV = MUL(iSymU,iSymV)
            do iSymY = 1, nSym
              if (iSym == MUL(iSymY,iSymUV)) then
                do IU = 1, NASH(iSymU)
                  IUABS = IU + NAES(iSymU)
                  do IV = 1, NASH(iSymV)
                    IVABS = IV + NAES(iSymV)
                    do IY = 1, NASH(iSymY)
                      IYABS = IY + NAES(iSymY)
                      VALUE = VALUE - Gact(IUABS,IVABS,IYABS,IXABS)*G2(IYABS,ITABS,IUABS,IVABS)
                      if (ITABS == IYABS) VALUE = VALUE + Two*Gact(IUABS,IVABS,IYABS,IXABS)*G1(IUABS,IVABS)
                      if (ITABS == IUABS) VALUE = VALUE - Gact(IUABS,IVABS,IYABS,IXABS)*G1(IYABS,IVABS)
                    end do
                  end do
                end do
              end if
            end do
          end do
        end do
        IBE=(IT*(IT-1))/2+IX
        BE(IBE)=VALUE
      END DO
    END DO

    if (NBE > 0 .and. NINDEP(ISYM,6) > 0) then
      IDISK=IDBMAT(ISYM,6)
      CALL DDAFILE(LUSBT,1,BE,NBE,IDISK)
      IF(NINM > 0 .and. NINDEP(ISYM,7) > 0) THEN
        IDISK=IDBMAT(ISYM,7)
        CALL DDAFILE(LUSBT,1,BE,NBE,IDISK)
      END IF
      CALL mma_deallocate(BE)
    end if
  end do

end subroutine MKBNEVE
