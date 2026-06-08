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
  use caspt2_global, only: LUSBT
  use caspt2_module, only: NAES, NASH, NINDEP, NSYM
  use EQSOLV, only: IDBMAT
  use stdalloc, only: mma_allocate, mma_deallocate
  use definitions, only: iwp,wp
  use Constants, only: Zero
  use Symmetry_Info, only: Mul

  implicit none

  integer(kind=iwp), intent(in) :: nAshT
  real(kind=wp), intent(in) :: Hact(nAshT,nAshT), Gact(nAshT,nAshT,nAshT,nAshT), &
                               G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT)

  integer(kind=iwp) :: IBG, IDISK, ISYM, IT, ITABS, IU, IUABS, IV, IVABS, IX, IXABS, IY, IYABS, &
                       NAS, NBG, NINM
  integer(kind=iwp) :: iSymU, iSymUV, iSymV, iSymY
  real(kind=wp) :: VALUE

  real(kind=wp), allocatable :: BG(:)

  DO ISYM=1,NSYM
    IF(NINDEP(ISYM,10) == 0) cycle
    NINM=NINDEP(ISYM,11)
    NAS=NASH(ISYM)
    NBG=(NAS*(NAS+1))/2
    if (NBG > 0) then
      CALL mma_Allocate(BG,NBG,LABEL='BG')
    end if

    ! Eq.(27): K(t,x) = - h(x,u)*R(t,u) - (uy|xv)*D(tu,vy)
    DO IT=1,NAS
      ITABS=IT+NAES(ISYM)
      DO IX=1,IT
        IXABS=IX+NAES(ISYM)
        VALUE=Zero
!       DO IU=1,NAS
!         IUABS=IU+NAES(ISYM)
!         VALUE = VALUE + Hact(IXABS,IUABS)*G1(ITABS,IUABS)
!         DO IV=1,NAS
!           IVABS=IV+NAES(ISYM)
!           DO IY=1,NAS
!             IYABS=IY+NAES(ISYM)
!             VALUE = VALUE + Gact(IUABS,IYABS,IXABS,IVABS)*G2(ITABS,IVABS,IUABS,IYABS)
!           END DO
!         END DO
!       END DO
        do IU = 1, NAS
          IUABS = IU + NAES(iSym)
          VALUE = VALUE + Hact(IXABS,IUABS)*G1(ITABS,IUABS)
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
                      VALUE = VALUE + Gact(IUABS,IYABS,IXABS,IVABS)*G2(ITABS,IVABS,IUABS,IYABS)
                    end do
                  end do
                end do
              end if
            end do
          end do
        end do
        IBG=(IT*(IT-1))/2+IX
        BG(IBG)=-VALUE
!       write (6,'("BG",2i3,f20.10)') it,ix,bg(ibg)
      END DO
    END DO

    if (NBG > 0) then
      IF(NINDEP(ISYM,10) > 0) THEN
       IDISK=IDBMAT(ISYM,10)
       CALL DDAFILE(LUSBT,1,BG,NBG,IDISK)
      END IF
      IF(NINM > 0 .and. NINDEP(ISYM,11) > 0) THEN
        IDISK=IDBMAT(ISYM,11)
        CALL DDAFILE(LUSBT,1,BG,NBG,IDISK)
      END IF
      CALL mma_deallocate(BG)
    end if
  end do

end subroutine MKBNEVG
