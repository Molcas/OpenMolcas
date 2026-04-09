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

module CASPT2_Module

use Molcas, only: LenIn, MxAct, MxBas, MxIna, MxOrb, MxRoot
use RASDim, only: MxTit
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: MXCASE = 13, MxExt = MxBas-MxAct, MxState = MxRoot

integer(kind=iwp) :: DWTYPE, IAD1M(64), IASYM(MxAct), IEOF1M, IEXTIS(MXEXT), IFQCAN, IINAIS(MxIna), IISYM(MxIna), &
                     IOFFRHS(8,MXCASE), iRlxRoot, IROOT(MxRoot), ISCF, ISPIN, JSTATE, LROOTS, MAXIT, MSTATE(MXSTATE), NACTEL, &
                     NAES(8), NAGEB(8), NAGEBES(8), NAGTB(8), NAGTBES(8), NAMX, NASH(8), NASHT, NASUP(8,MXCASE), NBAS(8), NBAST, &
                     NBMX, NBSQT, NBTCH(8), NBTCHES(8), NBTRI, NCASES, NCONF, NDEL(8), NDET, NELE3, NEXCES(8,MXCASE), NFRO(8), &
                     NFROT, NGROUP, NGROUPSTATE(MXSTATE), NHOLE1, NIAES(8), NIES(8), NIGEJ(8), NIGEJES(8), NIGTJ(8), NIGTJES(8), &
                     NIMX, NINABX, NINDEP(8,MXCASE), NISH(8), NISHT, NISUP(8,MXCASE), NLYGROUP, NLYROOT, NOMX, NORB(8), NORBT, &
                     NOSH(8), NOSQT, NOTRI, NRAS1(8), NRAS1T, NRAS2(8), NRAS2T, NRAS3(8), NRAS3T, NROOTS, NSECBX, NSES(8), NSMX, &
                     NSSH(8), NSSHT, NSTATE, NSYM, NTGEU(8), NTGEUES(8), NTGTU(8), NTGTUES(8), NTU(8), NTUES(8), NTUV(8), &
                     NTUVES(8), ROOT2STATE(MxRoot), STSYM
real(kind=wp) :: CPUEIG, CPUFG3, CPUFMB, CPUGIN, CPUGRD, CPUINT, CPULCS, CPUNAD, CPUOVL, CPUPCG, CPUPRP, CPUPT2, CPURHS, CPUSBM, &
                 CPUSCA, CPUSER, CPUSGM, CPUSIN, CPUVEC, DENORM, E2CORR, E2TOT, EASUM, ENERGY(MXSTATE), EPS(MXORB), EPSA(MxAct), &
                 EPSE(MXEXT), EPSI(MxIna), EREF, ERFSELF, POTNUC, REFENE(MXSTATE), RNORM, THRCONV, THRENE, THROCC, THRSHN, THRSHS, &
                 TIOEIG, TIOFG3, TIOFMB, TIOGIN, TIOGRD, TIOINT, TIOLCS, TIONAD, TIOOVL, TIOPCG, TIOPRP, TIOPT2, TIORHS, TIOSBM, &
                 TIOSCA, TIOSER, TIOSGM, TIOSIN, TIOVEC, ZETA
logical(kind=iwp) :: DMRG, DoCumulant, G1SECIN, IfChol, IFDENS, IFDORTHO, IFDW, IFMIX, IFMSCOUP, IFPROP, IFRMS, IFSADREF, IFXMS, &
                     JMS, PRORB, PRSD, RFPERT, RHSDIRECT
character(len=LenIn+8) :: BNAME(2*MXORB)
character(len=32) :: HZERO
character(len=8) :: BMATRIX, BSPECT, BTRANS, CASES(MXCASE), FOCKTYPE, ISNAM(MXEXT), ORBIN, ORBNAM(MXORB), OUTFMT, SDECOM, SMATRIX
character(len=4) :: TITLE(18,MxTit)
character(len=2) :: HEADER(72)

public :: BMATRIX, BNAME, BSPECT, BTRANS, CASES, CPUEIG, CPUFG3, CPUFMB, CPUGIN, CPUGRD, CPUINT, CPULCS, CPUNAD, CPUOVL, CPUPCG, &
          CPUPRP, CPUPT2, CPURHS, CPUSBM, CPUSCA, CPUSER, CPUSGM, CPUSIN, CPUVEC, DENORM, DMRG, DoCumulant, DWTYPE, E2CORR, E2TOT, &
          EASUM, ENERGY, EPS, EPSA, EPSE, EPSI, EREF, ERFSELF, FOCKTYPE, G1SECIN, HEADER, HZERO, IAD1M, IASYM, IEOF1M, IEXTIS, &
          IfChol, IFDENS, IFDORTHO, IFDW, IFMIX, IFMSCOUP, IFPROP, IFQCAN, IFRMS, IFSADREF, IFXMS, IINAIS, IISYM, IOFFRHS, &
          iRlxRoot, IROOT, ISCF, ISNAM, ISPIN, JMS, JSTATE, LROOTS, MAXIT, MSTATE, MXCASE, MXEXT, NACTEL, NAES, NAGEB, NAGEBES, &
          NAGTB, NAGTBES, NAMX, NASH, NASHT, NASUP, NBAS, NBAST, NBMX, NBSQT, NBTCH, NBTCHES, NBTRI, NCASES, NCONF, NDEL, NDET, &
          NELE3, NEXCES, NFRO, NFROT, NGROUP, NGROUPSTATE, NHOLE1, NIAES, NIES, NIGEJ, NIGEJES, NIGTJ, NIGTJES, NIMX, NINABX, &
          NINDEP, NISH, NISHT, NISUP, NLYGROUP, NLYROOT, NOMX, NORB, NORBT, NOSH, NOSQT, NOTRI, NRAS1, NRAS1T, NRAS2, NRAS2T, &
          NRAS3, NRAS3T, NROOTS, NSECBX, NSES, NSMX, NSSH, NSSHT, NSTATE, NSYM, NTGEU, NTGEUES, NTGTU, NTGTUES, NTU, NTUES, NTUV, &
          NTUVES, ORBIN, ORBNAM, OUTFMT, POTNUC, PRORB, PRSD, REFENE, RFPERT, RHSDIRECT, RNORM, ROOT2STATE, SDECOM, SMATRIX, &
          STSYM, THRCONV, THRENE, THROCC, THRSHN, THRSHS, TIOEIG, TIOFG3, TIOFMB, TIOGIN, TIOGRD, TIOINT, TIOLCS, TIONAD, TIOOVL, &
          TIOPCG, TIOPRP, TIOPT2, TIORHS, TIOSBM, TIOSCA, TIOSER, TIOSGM, TIOSIN, TIOVEC, TITLE, ZETA

end module CASPT2_Module
