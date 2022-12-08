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

integer(kind=iwp), parameter :: KBUFF1 = 2*9600, MADR = 20000, NTIBUF = nTraBuf
real(kind=wp), parameter :: SQ2 = sqrt(Two)
integer(kind=iwp) :: IAD25S, IADABCI, IADDP(79), ICH(MXORB), ICONV, ICPF, IDENS, IDIIS, IFIRST, ILIM, INCPF, IPASS, IPRINT, &
                     IPS(200), IR1, IRC(4), IREF0, IREST, IROW(MXORB+1), ISC(4), ISDCI, ISMAX, ITER, ITOC17(64), ITPUL, IV0, IV1, &
                     JBUF, JJS(18), JMAX, JSC(4), KBUF, LASTAD(MADR), LBUF, LIC, LN, LSYM, Lu_25, Lu_27, Lu_30, Lu_CI, Lu_CIGuga, &
                     Lu_CPFORB, Lu_TiABCD, Lu_TiABCI, Lu_TiABIJ, Lu_TraInt, Lu_TraOne, MAX11, MAXIT, MAXITP, MX1, MX2, N, NASH(8), &
                     NBAS(8), NCONF, NDIAG(MXORB), NFRO(8), NISH(8), NNS(8), NORB(8), NORBT, NOV, NOV1, NPFRO(8), NREF, &
                     NSM(MXORB), NSYM, NSYS(9), NTMAX, NVIR(8), NVIRT, NVMAX, NVT5
real(kind=wp) :: CTRSH, DETOT, ETHRE, ETOT, POTNUC, WLEV
logical(kind=iwp) :: LWSP
character(len=LenIn8) :: BNAME(MXBAS)
integer(kind=iwp), allocatable :: INDX(:), ISAB(:), JSY(:), ICASE(:)

public :: BNAME, CTRSH, DETOT, ETHRE, ETOT, IAD25S, IADABCI, IADDP, ICASE, ICH, ICONV, ICPF, IDENS, IDIIS, IFIRST, ILIM, INCPF, &
          INDX, IPASS, IPRINT, IPS, IR1, IRC, IREF0, IREST, IROW, ISAB, ISC, ISDCI, ISMAX, ITER, ITOC17, ITPUL, IV0, IV1, JBUF, &
          JJS, JMAX, JSC, JSY, KBUF, KBUFF1, LASTAD, LBUF, LIC, LN, LSYM, Lu_25, Lu_27, Lu_30, Lu_CI, Lu_CIGuga, Lu_CPFORB, &
          Lu_TiABCD, Lu_TiABCI, Lu_TiABIJ, Lu_TraInt, Lu_TraOne, LWSP, MADR, MAX11, MAXIT, MAXITP, MX1, MX2, N, NASH, NBAS, NCONF, &
          NDIAG, NFRO, NISH, NNS, NORB, NORBT, NOV, NOV1, NPFRO, NREF, NSM, NSYM, NSYS, NTIBUF, NTMAX, NVIR, NVIRT, NVMAX, NVT5, &
          POTNUC, SQ2, WLEV

end module CPF_global
