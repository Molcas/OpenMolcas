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

integer(kind=iwp), parameter :: LIX = 500000, MAXB = 8, MXCASE = 14000, MXVERT = 1000

integer(kind=iwp) :: IA(MXVERT), IADD11, IAF(MXVERT), IB(MXVERT), IBF(MXVERT), ICASE(MXCASE), ICH(55), ICOUP(55), ICOUP1(55), &
                     IFIRST, IJ(55), IJF(55), ILIM, IOUT, IPO(0:MXVERT), IPRINT, IRC(4), ISPA, ISPIN, IV0, IVF0, IWAY(55), &
                     IX(4*MXVERT), IY(4*MXVERT,3), J1(55), J2(55), JM(55), JM1(55), JNDX(LIX), JRC(4), K0(4*MXVERT), &
                     K0F(0:MXVERT), K1(4*MXVERT), K1F(0:MXVERT), K2(4*MXVERT), K2F(0:MXVERT), K3(4*MXVERT), K3F(0:MXVERT), LN, &
                     LNP, Lu_10, Lu_11, N, NBUF, NIORB, NMAT, NSM(55), NSYM
real(kind=wp) :: BL1(MAXB+3), BL2(MAXB+3), BS1(MAXB+3), BS2(MAXB+3), BS3(MAXB+3), BS4(MAXB+3), COUP(55), COUP1(55), S

public :: BL1, BL2, BS1, BS2, BS3, BS4, COUP, COUP1, IA, IADD11, IAF, IB, IBF, ICASE, ICH, ICOUP, ICOUP1, IFIRST, IJ, IJF, ILIM, &
          IOUT, IPO, IPRINT, IRC, ISPA, ISPIN, IV0, IVF0, IWAY, IX, IY, J1, J2, JM, JM1, JNDX, JRC, K0, K0F, K1, K1F, K2, K2F, K3, &
          K3F, LN, LNP, Lu_10, Lu_11, MAXB, MXVERT, N, NBUF, NIORB, NMAT, NSM, NSYM, S

end module guga_global
