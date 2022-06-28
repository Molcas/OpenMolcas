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

module ccsort_global

use Definitions, only: wp, iwp

implicit none
private

!1.   parameters for expansion of orbitals
! nsize  nmbas
!
!3.   names of TEMP files and status matrix for TEMP files
! tmpnam stattemp lrectemp nrectemp
!
!4.   arrays for expanding of orbitals
! valh jh kh lh nshow
!
!5.   arrays for expanding of orbitals
!     reflecting permutation
! np nq nr ns typ idis
!
!6.   four mapd,mapi matrices and corresponding initial positions variables
!     mapd and mapi for R_i matrix, required for making T3 integrals
! pos10 pos20 pos30 mapd1 mapi1 mapd2 mapi2 mapd3 mapi3 posri0 mapdri mapiri
!
!7.   lun for INTA1 <ma||ef>aaaa, <ma||ef>baab
!     lun for INTA2 <ma||ef>bbbb, <ma||ef>abab
!     lun for INTA3 <ma||ej>aaaa, <ma||ej>baab, <ma||ej>baba
!     lun for INTA4 <ma||ej>bbbb, <ma||ej>abab, <ma||ej>abba
!     lun for INTAB _a_b(p,q)
!     lunt3 - Lun for t3nam file
! luna1 luna2 luna3 luna4 lunab lunt3
!
!8.   parameters for direct access file
!     lun and reclen (in R8 words) for direct access TEMPDA1,TEMPDA2
! lunda1 lunda2 reclen
!
!10.  input keys
!     cckey - key for doing CCSD integrals
!     t3key - key for doing T3 integrals
!     clopkey - closed/open key
!     nfror - forzen orbitals per symmetry in Reorg
!     ndelr - deleted orbitals per symmetry in Reorg
!     maxspace - maximal allowed allocatable area
!     fullprint - output printing control key
!     noop - no operation key
!     iokey - disk handling control key
!     zrkey - key for reading I values and indices simultanously
! cckey t3key clopkey nfror ndelr maxspace fullprint noop iokey zrkey
!
!     disk addresses for MOLCAS DA file handling
! daddr

integer(kind=iwp), parameter :: mbas = 1024, nsize = 8192, reclen = 100
integer(kind=iwp) :: cckey, clopkey, daddr(128), fullprint, IADR15(64), idis(8,8,8), iokey, IPT2, ISCF, ISPIN, jh(nsize), JOBIPH, &
                     kh(nsize), lh(nsize), lrectemp(mbas), LROOT, LSYM, LUINTM, luna1, luna2, luna3, luna4, lunab, lunda1, lunda2, &
                     lunpublic, lunt3, mapd1(0:512,6), mapd2(0:512,6), mapd3(0:512,6), mapdri(0:512,6), mapi1(8,8,8), &
                     mapi2(8,8,8), mapi3(8,8,8), mapiri(8,8,8), maxspace, NACTEL, NASH(8), NASHT, NBAS(8), nBasX(8), NCONF, &
                     NDEL(8), ndelr(8), nDelX(8), NELE3, NFRO(8), nfror(8), nFroX(8), NHOLE1, NISH(8), NISHT, noa(8), nob(8), &
                     noop, NORB(8), np(8,8,8), nq(8,8,8), nr(8,8,8), nrectemp(mbas), NROOTS, ns(8,8,8), nshow(mbas), NSSH(8), &
                     NSSHT, NSYM, nSymX, nva(8), nvb(8), pos10, pos20, pos30, posri0, stattemp(mbas), t3key, typ(8,8,8), zrkey
real(kind=wp) :: EScf, valh(nsize)
character(len=7) :: tmpnam(mbas)

public :: cckey, clopkey, daddr, Escf, fullprint, IADR15, idis, iokey, IPT2, ISCF, ISPIN, jh, JOBIPH, kh, lh, lrectemp, LROOT, &
          LSYM, LUINTM, luna1, luna2, luna3, luna4, lunab, lunda1, lunda2, lunpublic, lunt3, mapd1, mapd2, mapd3, mapdri, mapi1, &
          mapi2, mapi3, mapiri, maxspace, mbas, NACTEL, NASH, NASHT, NBAS, nBasX, NCONF, NDEL, ndelr, nDelX, NELE3, NFRO, nfror, &
          nFroX, NHOLE1, NISH, NISHT, noa, nob, noop, NORB, np, nq, nr, nrectemp, NROOTS, ns, nshow, nsize, NSSH, NSSHT, NSYM, &
          nSymX, nva, nvb, pos10, pos20, pos30, posri0, reclen, stattemp, t3key, tmpnam, typ, valh, zrkey

end module ccsort_global
