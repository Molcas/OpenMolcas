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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine COMMWEW(IVEC,JVEC,DCOM)

use Symmetry_Info, only: Mul
use constants, only: Zero, One, Two
use SUPERINDEX, only: KTUV, KTGEU, KTGTU, KTU
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_global, only: LUSBT
use EQSOLV, only: IDSMAT
use caspt2_module, only: NSYM, NASUP, NISUP, NASHT, NTUVES, NASH, NAES, IASYM, NTGEUES, NTGTUES, NTUES
use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC
real(kind=wp), intent(inout) :: DCOM(NASHT,NASHT)
real(kind=wp), allocatable :: CBLK(:), TBLK(:), SMAT(:)
integer(kind=iwp) ICASE, ISYM, NAS, NIS, NCBLK, NS, IDS
integer(kind=iwp) K000, IIS, ISYMT, ISYMTU, ISYMU, ISYMX, ITABS, ITUX, ITUY, ITXU, ITYU, IU, IUABS, IX, IXABS, IXTU, IY, IYABS, &
                  IYTU, NAX
real(kind=wp) SUM
integer(kind=iwp) IT, IXT, IYT
real(kind=wp) PARTSUM
real(kind=wp) SGN
integer(kind=iwp) ITX1, ITX2, ITY1, ITY2, IXT1, IXT2, IYT1, IYT2, NAS1

! This subroutine is one of the components needed to compute the active/active
! transition density matrix elements for the two first-order vectors IVEC and
! JVEC.
! Adds, into the matrix DCOM, a correction obtained by commutation relations.
! Present assumption: The two vectors nr. IVEC and JVEC, stored on LUSOLV,
! are both in contravariant representation. Possibly, IVEC equals JVEC.

do ICASE=1,11
  do ISYM=1,NSYM
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    NCBLK = NAS*NIS
    if (NCBLK == 0) cycle
    ! Allocate CBLK, TBLK
    call MMA_ALLOCATE(CBLK,NCBLK,Label='CBLK')
    call MMA_ALLOCATE(TBLK,NCBLK,Label='TBLK')
    ! First, read in the ICASE, ISYM block of coefficients from JVEC into CBLK:
    ! Note carefully, this is not a mistake: vector JVEC into CBLK it is!
    call RDBLKC(ISYM,ICASE,JVEC,CBLK,NCBLK)
    ! Allocate overlap matrix:
    NS = (NAS*(NAS+1))/2
    call MMA_ALLOCATE(SMAT,NS)
    IDS = IDSMAT(ISYM,ICASE)
    call DDAFILE(LUSBT,2,SMAT,NS,IDS)
    ! Compute TBLK as the covariant representation of vector JVEC, by multiplying
    ! with the overlap matrix. Then get rid of the overlap matrix.
    call DCOPY_(NCBLK,[Zero],0,TBLK,1)
    call TRIMUL(NAS,NIS,One,SMAT,CBLK,NAS,TBLK,NAS)
    call MMA_DEALLOCATE(SMAT)
    ! Finally, if IVEC not equals JVEC, read in the contravariant block of vector
    ! IVEC into CBLK:
    if (IVEC /= JVEC) call RDBLKC(ISYM,ICASE,IVEC,CBLK,NCBLK)
    ! Finally, branch to the appropriate code section:

    select case (ICASE)
      case (1)
        ! Case 1 code section:
        K000 = NTUVES(ISYM)

        do ISYMX=1,NSYM
          NAX = NASH(ISYMX)
          do IX=1,NAX
            IXABS = NAES(ISYMX)+IX
            do IY=1,NAX
              IYABS = NAES(ISYMX)+IY

              SUM = Zero
              ISYMTU = Mul(ISYMX,ISYM)
              do ITABS=1,NASHT
                ISYMT = IASYM(ITABS)
                ISYMU = Mul(ISYMT,ISYMTU)
                do IU=1,NASH(ISYMU)
                  IUABS = NAES(ISYMU)+IU

                  IYTU = KTUV(IYABS,ITABS,IUABS)-K000
                  IXTU = KTUV(IXABS,ITABS,IUABS)-K000
                  ITYU = KTUV(ITABS,IYABS,IUABS)-K000
                  ITXU = KTUV(ITABS,IXABS,IUABS)-K000
                  ITUY = KTUV(ITABS,IUABS,IYABS)-K000
                  ITUX = KTUV(ITABS,IUABS,IXABS)-K000

                  do IIS=1,NIS
                    SUM = SUM+CBLK(IXTU+NAS*(IIS-1))*TBLK(IYTU+NAS*(IIS-1))
                    SUM = SUM+CBLK(ITXU+NAS*(IIS-1))*TBLK(ITYU+NAS*(IIS-1))
                    SUM = SUM-CBLK(ITUY+NAS*(IIS-1))*TBLK(ITUX+NAS*(IIS-1))
                  end do

                end do
              end do
              DCOM(IXABS,IYABS) = DCOM(IXABS,IYABS)+SUM

            end do
          end do
        end do

      case (2)
        ! Case 2 code section:
        do ISYMX=1,NSYM
          NAX = NASH(ISYMX)
          do IX=1,NAX
            IXABS = NAES(ISYMX)+IX
            do IY=1,NAX
              IYABS = NAES(ISYMX)+IY

              SUM = Zero
              ISYMT = Mul(ISYMX,ISYM)
              do IT=1,NASH(ISYMT)
                ITABS = NAES(ISYMT)+IT
                if (ITABS >= IXABS) then
                  IXT = KTGEU(ITABS,IXABS)-NTGEUES(ISYM)
                else
                  IXT = KTGEU(IXABS,ITABS)-NTGEUES(ISYM)
                end if
                if (ITABS >= IYABS) then
                  IYT = KTGEU(ITABS,IYABS)-NTGEUES(ISYM)
                else
                  IYT = KTGEU(IYABS,ITABS)-NTGEUES(ISYM)
                end if
                PARTSUM = Zero
                do IIS=1,NIS
                  PARTSUM = PARTSUM+CBLK(IXT+NAS*(IIS-1))*TBLK(IYT+NAS*(IIS-1))
                end do
                if (ITABS == IXABS) PARTSUM = Two*PARTSUM
                SUM = SUM+PARTSUM
              end do
              DCOM(IXABS,IYABS) = DCOM(IXABS,IYABS)+SUM

            end do
          end do
        end do

      case (3)
        ! Case 3 code section:
        do ISYMX=1,NSYM
          NAX = NASH(ISYMX)
          do IX=1,NAX
            IXABS = NAES(ISYMX)+IX
            do IY=1,NAX
              IYABS = NAES(ISYMX)+IY

              SUM = Zero
              ISYMT = Mul(ISYMX,ISYM)
              do IT=1,NASH(ISYMT)
                ITABS = NAES(ISYMT)+IT
                if (ITABS == IXABS) cycle
                if (ITABS == IYABS) cycle
                if (ITABS > IXABS) then
                  IXT = KTGTU(ITABS,IXABS)-NTGTUES(ISYM)
                  SGN = One
                else
                  IXT = KTGTU(IXABS,ITABS)-NTGTUES(ISYM)
                  SGN = -One
                end if
                if (ITABS > IYABS) then
                  IYT = KTGTU(ITABS,IYABS)-NTGTUES(ISYM)
                else
                  IYT = KTGTU(IYABS,ITABS)-NTGTUES(ISYM)
                  SGN = -SGN
                end if
                PARTSUM = Zero
                do IIS=1,NIS
                  PARTSUM = PARTSUM+CBLK(IXT+NAS*(IIS-1))*TBLK(IYT+NAS*(IIS-1))
                end do
                SUM = SUM+SGN*PARTSUM

              end do
              DCOM(IXABS,IYABS) = DCOM(IXABS,IYABS)+SUM

            end do
          end do
        end do

      case (4)
        ! Case 4 code section:
        K000 = NTUVES(ISYM)

        do ISYMX=1,NSYM
          NAX = NASH(ISYMX)
          do IX=1,NAX
            IXABS = NAES(ISYMX)+IX
            do IY=1,NAX
              IYABS = NAES(ISYMX)+IY

              SUM = Zero
              ISYMTU = Mul(ISYMX,ISYM)
              do ITABS=1,NASHT
                ISYMT = IASYM(ITABS)
                ISYMU = Mul(ISYMT,ISYMTU)
                do IU=1,NASH(ISYMU)
                  IUABS = NAES(ISYMU)+IU

                  IYTU = KTUV(IYABS,ITABS,IUABS)-K000
                  IXTU = KTUV(IXABS,ITABS,IUABS)-K000
                  ITYU = KTUV(ITABS,IYABS,IUABS)-K000
                  ITXU = KTUV(ITABS,IXABS,IUABS)-K000
                  ITUY = KTUV(ITABS,IUABS,IYABS)-K000
                  ITUX = KTUV(ITABS,IUABS,IXABS)-K000

                  do IIS=1,NIS
                    SUM = SUM-CBLK(IYTU+NAS*(IIS-1))*TBLK(IXTU+NAS*(IIS-1))
                    SUM = SUM+CBLK(ITXU+NAS*(IIS-1))*TBLK(ITYU+NAS*(IIS-1))
                    SUM = SUM-CBLK(ITUY+NAS*(IIS-1))*TBLK(ITUX+NAS*(IIS-1))
                  end do

                end do
              end do
              DCOM(IXABS,IYABS) = DCOM(IXABS,IYABS)+SUM

            end do
          end do
        end do

      case (5)
        ! Case 5 code section:
        NAS1 = NAS/2

        do ISYMX=1,NSYM
          NAX = NASH(ISYMX)
          do IX=1,NAX
            IXABS = NAES(ISYMX)+IX
            do IY=1,NAX
              IYABS = NAES(ISYMX)+IY

              SUM = Zero
              ISYMT = Mul(ISYMX,ISYM)
              do IT=1,NASH(ISYMT)
                ITABS = NAES(ISYMT)+IT
                IXT1 = KTU(IXABS,ITABS)-NTUES(ISYM)
                IYT1 = KTU(IYABS,ITABS)-NTUES(ISYM)
                ITX1 = KTU(ITABS,IXABS)-NTUES(ISYM)
                ITY1 = KTU(ITABS,IYABS)-NTUES(ISYM)
                IXT2 = IXT1+NAS1
                IYT2 = IYT1+NAS1
                ITX2 = ITX1+NAS1
                ITY2 = ITY1+NAS1
                do IIS=1,NIS
                  SUM = SUM+CBLK(IXT1+NAS*(IIS-1))*TBLK(IYT1+NAS*(IIS-1))
                  SUM = SUM-CBLK(ITY1+NAS*(IIS-1))*TBLK(ITX1+NAS*(IIS-1))
                  SUM = SUM+CBLK(IXT2+NAS*(IIS-1))*TBLK(IYT2+NAS*(IIS-1))
                  SUM = SUM-CBLK(ITY2+NAS*(IIS-1))*TBLK(ITX2+NAS*(IIS-1))
                end do
              end do
              DCOM(IXABS,IYABS) = DCOM(IXABS,IYABS)+SUM

            end do
          end do
        end do

      case (6)
        ! Case 6 code section:
        NAX = NASH(ISYM)
        do IX=1,NAX
          IXABS = NAES(ISYM)+IX
          do IY=1,NAX
            IYABS = NAES(ISYM)+IY

            SUM = Zero
            do IIS=1,NIS
              SUM = SUM+CBLK(IX+NAS*(IIS-1))*TBLK(IY+NAS*(IIS-1))
            end do
            DCOM(IXABS,IYABS) = DCOM(IXABS,IYABS)+SUM

          end do
        end do

      case (7)
        ! Case 7 code section:
        NAX = NASH(ISYM)
        do IX=1,NAX
          IXABS = NAES(ISYM)+IX
          do IY=1,NAX
            IYABS = NAES(ISYM)+IY

            SUM = Zero
            do IIS=1,NIS
              SUM = SUM+CBLK(IX+NAS*(IIS-1))*TBLK(IY+NAS*(IIS-1))
            end do
            DCOM(IXABS,IYABS) = DCOM(IXABS,IYABS)+SUM

          end do
        end do

      case (8)
        ! Case 8 code section:
        do ISYMX=1,NSYM
          NAX = NASH(ISYMX)
          do IX=1,NAX
            IXABS = NAES(ISYMX)+IX
            do IY=1,NAX
              IYABS = NAES(ISYMX)+IY

              SUM = Zero
              ISYMT = Mul(ISYMX,ISYM)
              do IT=1,NASH(ISYMT)
                ITABS = NAES(ISYMT)+IT
                if (ITABS >= IXABS) then
                  IXT = KTGEU(ITABS,IXABS)-NTGEUES(ISYM)
                else
                  IXT = KTGEU(IXABS,ITABS)-NTGEUES(ISYM)
                end if
                if (ITABS >= IYABS) then
                  IYT = KTGEU(ITABS,IYABS)-NTGEUES(ISYM)
                else
                  IYT = KTGEU(IYABS,ITABS)-NTGEUES(ISYM)
                end if
                PARTSUM = Zero
                do IIS=1,NIS
                  PARTSUM = PARTSUM-CBLK(IYT+NAS*(IIS-1))*TBLK(IXT+NAS*(IIS-1))
                end do
                if (ITABS == IYABS) PARTSUM = Two*PARTSUM
                SUM = SUM+PARTSUM

              end do
              DCOM(IXABS,IYABS) = DCOM(IXABS,IYABS)+SUM

            end do
          end do
        end do

      case (9)
        ! Case 9 code section:
        do ISYMX=1,NSYM
          NAX = NASH(ISYMX)
          do IX=1,NAX
            IXABS = NAES(ISYMX)+IX
            do IY=1,NAX
              IYABS = NAES(ISYMX)+IY

              SUM = Zero
              ISYMT = Mul(ISYMX,ISYM)
              do IT=1,NASH(ISYMT)
                ITABS = NAES(ISYMT)+IT
                if (ITABS == IXABS) cycle
                if (ITABS == IYABS) cycle
                if (ITABS > IXABS) then
                  IXT = KTGTU(ITABS,IXABS)-NTGTUES(ISYM)
                  SGN = One
                else
                  IXT = KTGTU(IXABS,ITABS)-NTGTUES(ISYM)
                  SGN = -One
                end if
                if (ITABS > IYABS) then
                  IYT = KTGTU(ITABS,IYABS)-NTGTUES(ISYM)
                else
                  IYT = KTGTU(IYABS,ITABS)-NTGTUES(ISYM)
                  SGN = -SGN
                end if
                PARTSUM = Zero
                do IIS=1,NIS
                  PARTSUM = PARTSUM-CBLK(IYT+NAS*(IIS-1))*TBLK(IXT+NAS*(IIS-1))
                end do
                SUM = SUM+SGN*PARTSUM

              end do
              DCOM(IXABS,IYABS) = DCOM(IXABS,IYABS)+SUM

            end do
          end do
        end do

      case (10)
        ! Case 10 code section:
        NAX = NASH(ISYM)
        do IX=1,NAX
          IXABS = NAES(ISYM)+IX
          do IY=1,NAX
            IYABS = NAES(ISYM)+IY

            SUM = Zero
            do IIS=1,NIS
              SUM = SUM-CBLK(IY+NAS*(IIS-1))*TBLK(IX+NAS*(IIS-1))
            end do
            DCOM(IXABS,IYABS) = DCOM(IXABS,IYABS)+SUM

          end do
        end do

      case (11)
        ! Case 11 code section:
        NAX = NASH(ISYM)
        do IX=1,NAX
          IXABS = NAES(ISYM)+IX
          do IY=1,NAX
            IYABS = NAES(ISYM)+IY

            SUM = Zero
            do IIS=1,NIS
              SUM = SUM-CBLK(IY+NAS*(IIS-1))*TBLK(IX+NAS*(IIS-1))
            end do
            DCOM(IXABS,IYABS) = DCOM(IXABS,IYABS)+SUM

          end do
        end do

      case default
        call ABEND()
    end select

    call MMA_DEALLOCATE(CBLK)
    call MMA_DEALLOCATE(TBLK)

    ! Here ends the loops over ISYM and ICASE.
  end do
end do

end subroutine COMMWEW
