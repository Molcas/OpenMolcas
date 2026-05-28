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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine NADIAG()
! Set up non-active diagonal elements of H0.

use Symmetry_Info, only: Mul
use SUPERINDEX, only: MAGEB, MAGTB, MIGEJ, MIGTJ
use EQSOLV, only: IDBMAT
use caspt2_global, only: LUSBT
use caspt2_module, only: EPSE, EPSI, NAGEB, NAGEBES, NAGEBES, NAGTB, NAGTBES, NASUP, NIES, NIGEJ, NIGEJES, NIGTJ, NIGTJES, NINDEP, &
                         NISH, NISUP, NSES, NSSH, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: I2, I2ABS, IA, IAABS, IAB, IABQ, IBABS, ICASE, IDID, II, IIABS, IIJ, IIJQ, IIQ, IIS, IJABS, ISYM, ISYMA, &
                     ISYMAB, ISYMI, ISYMIJ, NAS, NIN, NIS
real(kind=wp) :: Dummy(1)
real(kind=wp), allocatable :: BD(:), ID(:)

do ICASE=1,13
  do ISYM=1,NSYM

    NIN = NINDEP(ISYM,ICASE)
    if (NIN == 0) cycle
    NIS = NISUP(ISYM,ICASE)
    NAS = NASUP(ISYM,ICASE)
    if (ICASE > 11) call mma_allocate(BD,NAS,LABEL='BD')
    call mma_allocate(ID,NIS,LABEL='ID')

    select case (ICASE)
      case (1)
        ! VJTU CASE:
        ID(1:NIS) = -EPSI(NIES(ISYM)+1:NIES(ISYM)+NIS)

      case (2)
        ! VJTIP CASE:
        do IIS=1,NIS
          IIQ = IIS+NIGEJES(ISYM)
          IIABS = MIGEJ(1,IIQ)
          IJABS = MIGEJ(2,IIQ)
          ID(IIS) = -EPSI(IIABS)-EPSI(IJABS)
        end do

      case (3)
        ! VJTIM CASE:
        do IIS=1,NIS
          IIQ = IIS+NIGTJES(ISYM)
          IIABS = MIGTJ(1,IIQ)
          IJABS = MIGTJ(2,IIQ)
          ID(IIS) = -EPSI(IIABS)-EPSI(IJABS)
        end do

      case (4)
        ! ATVX  CASE:
        ID(1:NIS) = EPSE(NSES(ISYM)+1:NSES(ISYM)+NIS)

      case (5)
        ! AIVX  CASE:
        IIS = 0
        do ISYMA=1,NSYM
          ISYMI = Mul(ISYMA,ISYM)
          do IA=1,NSSH(ISYMA)
            IAABS = IA+NSES(ISYMA)
            do II=1,NISH(ISYMI)
              IIABS = II+NIES(ISYMI)
              IIS = IIS+1
              ID(IIS) = -EPSI(IIABS)+EPSE(IAABS)
            end do
          end do
        end do

      case (6)
        ! VJAIP CASE:
        IIS = 0
        do ISYMA=1,NSYM
          ISYMIJ = Mul(ISYMA,ISYM)
          do I2=1,NIGEJ(ISYMIJ)
            I2ABS = I2+NIGEJES(ISYMIJ)
            IIABS = MIGEJ(1,I2ABS)
            IJABS = MIGEJ(2,I2ABS)
            do IA=1,NSSH(ISYMA)
              IAABS = IA+NSES(ISYMA)
              IIS = IIS+1
              ID(IIS) = -EPSI(IIABS)-EPSI(IJABS)+EPSE(IAABS)
            end do
          end do
        end do

      case (7)
        ! VJAIM CASE:
        IIS = 0
        do ISYMA=1,NSYM
          ISYMIJ = Mul(ISYMA,ISYM)
          do I2=1,NIGTJ(ISYMIJ)
            I2ABS = I2+NIGTJES(ISYMIJ)
            IIABS = MIGTJ(1,I2ABS)
            IJABS = MIGTJ(2,I2ABS)
            do IA=1,NSSH(ISYMA)
              IAABS = IA+NSES(ISYMA)
              IIS = IIS+1
              ID(IIS) = -EPSI(IIABS)-EPSI(IJABS)+EPSE(IAABS)
            end do
          end do
        end do

      case (8)
        ! BVATP CASE:
        do IIS=1,NIS
          IIQ = IIS+NAGEBES(ISYM)
          IAABS = MAGEB(1,IIQ)
          IBABS = MAGEB(2,IIQ)
          ID(IIS) = EPSE(IAABS)+EPSE(IBABS)
        end do

      case (9)
        ! BVATM CASE:
        do IIS=1,NIS
          IIQ = IIS+NAGTBES(ISYM)
          IAABS = MAGTB(1,IIQ)
          IBABS = MAGTB(2,IIQ)
          ID(IIS) = EPSE(IAABS)+EPSE(IBABS)
        end do

      case (10)
        ! BJATP CASE:
        IIS = 0
        do ISYMI=1,NSYM
          ISYMAB = Mul(ISYMI,ISYM)
          do I2=1,NAGEB(ISYMAB)
            I2ABS = I2+NAGEBES(ISYMAB)
            IAABS = MAGEB(1,I2ABS)
            IBABS = MAGEB(2,I2ABS)
            do II=1,NISH(ISYMI)
              IIABS = II+NIES(ISYMI)
              IIS = IIS+1
              ID(IIS) = -EPSI(IIABS)+EPSE(IAABS)+EPSE(IBABS)
            end do
          end do
        end do

      case (11)
        ! BJATM CASE:
        IIS = 0
        do ISYMI=1,NSYM
          ISYMAB = Mul(ISYMI,ISYM)
          do I2=1,NAGTB(ISYMAB)
            I2ABS = I2+NAGTBES(ISYMAB)
            IAABS = MAGTB(1,I2ABS)
            IBABS = MAGTB(2,I2ABS)
            do II=1,NISH(ISYMI)
              IIABS = II+NIES(ISYMI)
              IIS = IIS+1
              ID(IIS) = -EPSI(IIABS)+EPSE(IAABS)+EPSE(IBABS)
            end do
          end do
        end do

      case (12)
        ! BJAIP CASE:
        do IAB=1,NAGEB(ISYM)
          IABQ = IAB+NAGEBES(ISYM)
          IAABS = MAGEB(1,IABQ)
          IBABS = MAGEB(2,IABQ)
          IIS = IIS+1
          BD(IAB) = EPSE(IAABS)+EPSE(IBABS)
        end do
        do IIJ=1,NIGEJ(ISYM)
          IIJQ = IIJ+NIGEJES(ISYM)
          IIABS = MIGEJ(1,IIJQ)
          IJABS = MIGEJ(2,IIJQ)
          ID(IIJ) = -EPSI(IIABS)-EPSI(IJABS)
        end do

      case (13)
        ! BJAIM CASE:
        do IAB=1,NAGTB(ISYM)
          IABQ = IAB+NAGTBES(ISYM)
          IAABS = MAGTB(1,IABQ)
          IBABS = MAGTB(2,IABQ)
          IIS = IIS+1
          BD(IAB) = EPSE(IAABS)+EPSE(IBABS)
        end do
        do IIJ=1,NIGTJ(ISYM)
          IIJQ = IIJ+NIGTJES(ISYM)
          IIABS = MIGTJ(1,IIJQ)
          IJABS = MIGTJ(2,IIJQ)
          ID(IIJ) = -EPSI(IIABS)-EPSI(IJABS)
        end do

      case default
        write(u6,*) 'NADIAG: illegal case number'
        call ABEND()
    end select

    ! NOTE: BDIAG elements used in cases 12 & 13.
    IDID = IDBMAT(ISYM,ICASE)
    if (ICASE > 11) then
      call DDAFILE(LUSBT,1,BD,NAS,IDID)
      call mma_deallocate(BD)
    else
      ! Dummy read the BDIAG elements. NOTE: NAS, not NIN.
      call DDAFILE(LUSBT,0,Dummy,NAS,IDID)
    end if
    call DDAFILE(LUSBT,1,ID,NIS,IDID)
    call mma_deallocate(ID)

  end do
end do

end subroutine NADIAG
