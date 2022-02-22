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

module cpf_global

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
private

! Lu_CIGuga - SYMBOLIC FORMULAS
! Lu_TraInt - TRANSFORMED MO 2-EL INTEGRALS
! Lu_TiABIJ - SORTED AIBJ, ABIJ AND AIJK INTEGRALS
! Lu_TiABCI - SORTED IJKL AND ABCI INTEGRALS
! Lu_TiABCD - SORTED ABCD INTEGRALS
! Lu_TraOne - ONE ELECTRON INTEGRALS
! Lu_CPFORB - (Formatted sequential!) CPF-ORBITALS OUT
! Lu_25     - FOCK MATRIX AND DIAGONAL CSF MATRIX ELEMENTS
! Lu_CI     - CI VECTOR
! Lu_27     - SCRATCH IN IIJJ
! Lu_30

#include "Molcas.fh"
#include "tratoc.fh"

integer(kind=iwp), parameter :: MADR = 20000, NTIBUF = nTraBuf
real(kind=wp), parameter :: SQ2 = sqrt(Two)
integer(kind=iwp) :: IPS(200), Lu_CIGuga, Lu_TraInt, Lu_TraOne, Lu_CI, Lu_TiABIJ, Lu_TiABCI, Lu_TiABCD, Lu_25, Lu_27, Lu_30, &
                     Lu_CPFORB, N, LN, NDIAG(MXORB), LSYM, LIC, KBUF, JJS(18), LW(99), NNS(8), IDENS, IREST, NCONF, ICH(MXORB), &
                     JBUF, IROW(MXORB+1), NSYM, MUL(8,8), IPASS, NSM(MXORB), IPRINT, IFIRST, IRC(4), ISC(4), JSC(4), LBUF, ITER, &
                     IV0, IV1, NSYS(9), MAXIT, ICPF, ISDCI, ITPUL, ICONV, IREF0, IDIIS, MAXITP, KBUFF1, INCPF, NREF, NPFRO(8), &
                     NFRO(8), NISH(8), NASH(8), NVIR(8), NORB(8), NBAS(8), NVIRT, NORBT, LASTAD(MADR), IAD25S, IADDP(79), &
                     IADC(79), IADABCI, ITOC17(64)
real(kind=wp) :: POTNUC, ETHRE, CTRSH, ETOT, WLEV, DETOT
logical(kind=iwp) :: LWSP
character(len=LenIn8) :: BNAME(MXBAS)

public :: IPS, LWSP, Lu_CIGuga, Lu_TraInt, Lu_TraOne, Lu_CI, Lu_TiABIJ, Lu_TiABCI, Lu_TiABCD, Lu_25, Lu_27, Lu_30, Lu_CPFORB, &
          POTNUC, ETHRE, CTRSH, ETOT, WLEV, DETOT, N, LN, NDIAG, LSYM, LIC, KBUF, JJS, LW, NNS, IDENS, IREST, NCONF, ICH, JBUF, &
          IROW, NSYM, MUL, IPASS, NSM, IPRINT, IFIRST, IRC, ISC, JSC, LBUF, ITER, IV0, IV1, NSYS, MAXIT, ICPF, ISDCI, ITPUL, &
          ICONV, IREF0, IDIIS, MAXITP, KBUFF1, INCPF, NREF, NPFRO, NFRO, NISH, NASH, NVIR, NORB, NBAS, NVIRT, NORBT, BNAME, SQ2, &
          MADR, LASTAD, NTIBUF, IAD25S, IADDP, IADC, IADABCI, ITOC17

end module CPF_global
