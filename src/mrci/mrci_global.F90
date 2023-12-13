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

module mrci_global

use Constants, only: Two, Half
use Definitions, only: wp, iwp

implicit none
private

#include "Molcas.fh"
#include "tratoc.fh"

integer(kind=iwp), parameter :: IDVER = 1, IVVER = 0, MCHAIN = 40000, MXPROP = 30, MXREF = 1000, MXVEC = 50, MXZ = 4*MXROOT, &
                                NSECT = 256, NSRTMX = nTraBuf, NTIBUF = nTraBuf
real(kind=wp), parameter :: SQ2 = sqrt(Two), SQ2INV = sqrt(Half)

integer(kind=iwp) :: IAD25S, IADABCI, ICH(MXORB), ICPF, IDFREE, IDISKC(MXVEC), IDISKD, IDISKS(MXVEC), IFIRST, INDSRT(NSRTMX+2), &
                     IORB(MXORB), IPASS, IPCOMP(MXPROP), IPRINT, IRC(4), IREFCI, IREFX(MXREF), IREST, IROOT(MXROOT), &
                     IROW(MXORB+1), ISMAX, ITER, ITOC17(64), ITRANS, JJS(18), JSC(4), KBUFF1, LASTAD(MCHAIN), LN, LSYM, Lu_25, &
                     Lu_27, Lu_60, Lu_70, Lu_80, LUEIG, LUONE, LUREST, LUSYMB, LUTRA, LUVEC, MAXIT, MBUF, MEMPRM, MEMTOT, MEMWRK, &
                     NASH(8), NBAS(8), NBAST, NBITM1, NBITM2, NBITM3, NBMAX, NBMN, NBTRI, NCHN1, NCHN2, NCHN3, NCMO, &
                     NCOMP(MXROOT), NCONF, NCSHT, NCSPCK, NCVAL, NDEL(8), NDIAG(MXORB), NDMO(8), NELEC, NFMO(8), NFRO(8), NISH(8), &
                     NNEW, NORB(8), NORBT, NPROP, NREF, NRROOT, NSEL, NSM(MXORB), NSTOT, NSYM, NVEC, NVIR(8), NVIRP(8), NVIRT, &
                     NVMAX, NVPAIR(8), NVSQ, NVTOT
real(kind=wp) :: CSEL(50), CTRSH, ENGY(MXROOT,3), ENP, ESHIFT, ESMALL(MXVEC), ETHRE, ETRSH, GFAC, PNUC(MXPROP), PORIG(3,MXPROP), &
                 POTNUC, SQNLIM, THRORB, TIBUF(NTIBUF), VALSRT(NSRTMX), VSMALL(MXVEC,MXVEC)
character(len=LenIn8) :: BNAME(MXBAS)
character(len=20) :: SSEL(50)
character(len=8) :: PNAME(MXPROP)
character(len=4) :: PTYPE(MXPROP)
integer(kind=iwp), allocatable :: CSPCK(:), INDX(:), INTSY(:), ISAB(:,:), JREFX(:)
real(kind=wp), allocatable :: CISEL(:,:), DMO(:), FIJKL(:), FOCK(:), HZERO(:,:), SZERO(:,:), TDMO(:,:), VZERO(:,:)

public :: BNAME, CISEL, CSEL, CSPCK, CTRSH, DMO, ENGY, ENP, ESHIFT, ESMALL, ETHRE, ETRSH, FIJKL, FOCK, GFAC, HZERO, IAD25S, &
          IADABCI, ICH, ICPF, IDFREE, IDISKC, IDISKD, IDISKS, IDVER, IFIRST, INDSRT, INDX, INTSY, IORB, IPASS, IPCOMP, IPRINT, &
          IRC, IREFCI, IREFX, IREST, IROOT, IROW, ISAB, ISMAX, ITER, ITOC17, ITRANS, IVVER, JJS, JREFX, JSC, KBUFF1, LASTAD, LN, &
          LSYM, Lu_25, Lu_27, Lu_60, Lu_70, Lu_80, LUEIG, LUONE, LUREST, LUSYMB, LUTRA, LUVEC, MAXIT, MBUF, MCHAIN, MEMPRM, &
          MEMTOT, MEMWRK, MXREF, MXVEC, MXZ, NASH, NBAS, NBAST, NBITM1, NBITM2, NBITM3, NBMAX, NBMN, NBTRI, NCHN1, NCHN2, NCHN3, &
          NCMO, NCOMP, NCONF, NCSHT, NCSPCK, NCVAL, NDEL, NDIAG, NDMO, NELEC, NFMO, NFRO, NISH, NNEW, NORB, NORBT, NPROP, NREF, &
          NRROOT, NSECT, NSEL, NSM, NSRTMX, NSTOT, NSYM, NTIBUF, NVEC, NVIR, NVIRP, NVIRT, NVMAX, NVPAIR, NVSQ, NVTOT, PNAME, &
          PNUC, PORIG, POTNUC, PTYPE, SQ2, SQ2INV, SQNLIM, SSEL, SZERO, TDMO, THRORB, TIBUF, VALSRT, VSMALL, VZERO

end module mrci_global
