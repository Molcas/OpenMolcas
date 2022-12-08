!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module guga_global

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: MXVERT = 1000

integer(kind=iwp) :: IA(MXVERT), IADD10, IADD11, IAF(MXVERT), IB(MXVERT), IBF(MXVERT), ICH(55), ICOUP(55), ICOUP1(55), IFIRST, &
                     IJ(55), IJF(55), ILIM, IOUT, IPO(0:MXVERT), IPRINT, IRC(4), ISPIN, IV0, IVF0, IWAY(55), IX(4*MXVERT), &
                     IY(4*MXVERT,3), J1(55), J2(55), JM(55), JM1(55), JRC(4), K0(4*MXVERT), K0F(0:MXVERT), K1(4*MXVERT), &
                     K1F(0:MXVERT), K2(4*MXVERT), K2F(0:MXVERT), K3(4*MXVERT), K3F(0:MXVERT), LN, LNP, Lu_10, Lu_11, N, NBUF, &
                     NIORB, NMAT, NSM(55), NSYM
real(kind=wp) :: COUP(55), COUP1(55), S
integer(kind=iwp), allocatable :: ICASE(:), JNDX(:)
real(kind=wp), allocatable :: BL1(:), BL2(:), BS1(:), BS2(:), BS3(:), BS4(:)

public :: BL1, BL2, BS1, BS2, BS3, BS4, COUP, COUP1, free_all, IA, IADD10, IADD11, IAF, IB, IBF, ICASE, ICH, ICOUP, ICOUP1, &
          IFIRST, IJ, IJF, ILIM, IOUT, IPO, IPRINT, IRC, ISPIN, IV0, IVF0, IWAY, IX, IY, J1, J2, JM, JM1, JNDX, JRC, K0, K0F, K1, &
          K1F, K2, K2F, K3, K3F, LN, LNP, Lu_10, Lu_11, MXVERT, N, NBUF, NIORB, NMAT, NSM, NSYM, S

contains

subroutine free_all()
  use stdalloc, only: mma_deallocate

  if (allocated(BL1)) call mma_deallocate(BL1)
  if (allocated(BL2)) call mma_deallocate(BL2)
  if (allocated(BS1)) call mma_deallocate(BS1)
  if (allocated(BS2)) call mma_deallocate(BS2)
  if (allocated(BS3)) call mma_deallocate(BS3)
  if (allocated(BS4)) call mma_deallocate(BS4)
  if (allocated(ICASE)) call mma_deallocate(ICASE)
  if (allocated(JNDX)) call mma_deallocate(JNDX)
end subroutine free_all

end module guga_global
