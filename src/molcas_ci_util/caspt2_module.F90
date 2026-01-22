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
Module CASPT2_Module
use definitions, only: wp, iwp
use Molcas
use rasdim
      implicit none
      integer(kind=iwp), Parameter:: MxState=MxRoot, MxExt=MxBas-MxAct
      real(kind=wp)ENERGY(MXSTATE), REFENE(MXSTATE)
      logical(kind=iwp)  lOPTO,RHSDIRECT,                               &
          RFPERT,G1SECIN,PRORB,PRSD,Found2,IfDirect,IfChol,             &
          MODE_DGA_COMPAT,DoCumulant,JMS
      real(kind=wp)THRCONV,THRSHN,THRSHS,THRSHF,THRENE,THROCC,          &
          POTNUC,ECORE,EREF,ERFSELF

      logical(kind=iwp) IFDENS,IFMIX,IFPROP,IFMSCOUP,                   &
          IFXMS,IFRMS,IFDW,IFSILPRROT,IFSADREF,IFDORTHO,                &
          IFINVAR,DMRG

! Note: MUL(8,8) is just a copy of Mul in the Symmetry_Info module
      integer(kind=iwp)MAXIT,NTIT,NACTEL,ISPIN,NSYM,STSYM,NCONF,NUNIQAT,&
          NFRO(8),NFROT,NISH(8),NIES(8),NISHT,NRAS1(8),NRAS1T,NDET,     &
          NRAS2(8),NRAS2T,NRAS3(8),NRAS3T,NASH(8),NAES(8),NASHT,        &
          NOSH(8),NOSHT,NSSH(8),NSES(8),NSSHT,NORB(8),NORBT,NOTRI,      &
          NDEL(8),NDELT,NBAS(8),NBAST,NBTRI,NIMX,NAMX,NSMX,NOMX,NBMX,   &
          NINABX,NSECBX,                                                &
          NOSQT,NBSQT,MUL(8,8),ISCF,                                    &
          IISYM(MXINA),IASYM(MXACT),IESYM(MXEXT),IINAIS(MXINA),         &
          IACTIS(MXACT),IEXTIS(MXEXT),                                  &
          NSTATE,MSTATE(MXSTATE),JSTATE,LROOTS,NROOTS,IROOT(MXROOT),    &
          ROOT2STATE(MXROOT),                                           &
          iRlxRoot,NGROUP,NGROUPSTATE(MXSTATE),                         &
          IAD1M(64),IEOF1M,NELE3,NHOLE1,IFQCAN,                         &
          NLYROOT,NLYGROUP,                                             &
          DWTYPE
      real(kind=wp) EPS(MXORB),EPSI(MXINA),EPSA(MXACT),EPSE(MXEXT),     &
          EASUM,DENORM,ZETA,                                            &
          ECOMM,E2CORR,E2TOT,RNORM,REFWGT

      integer(kind=iwp) NSTUV(8),NSTU(8),NSTGEU(8),NSTGTU(8),           &
          NTUV(8),NTU(8),NTGEU(8),NTGTU(8),                             &
          NIGEJ(8),NIGTJ(8),NAGEB(8),NAGTB(8),NTUVES(8),NTUES(8),       &
          NTGEUES(8),NTGTUES(8),NIGEJES(8),NIGTJES(8),NAGEBES(8),       &
          NAGTBES(8),NIAES(8)
! Excitation operators, sizes and offsets
      integer(kind=iwp), Parameter:: MXCASE=13
      integer(kind=iwp) NCASES,NASUP(8,MXCASE),NISUP(8,MXCASE),         &
          NINDEP(8,MXCASE),NEXC(8,MXCASE),NEXCES(8,MXCASE),             &
          NBTCH(8),NBTCHES(8),                                          &
          NSCP,NSCQ,NBATCH_TOT,NJSCT,NJSCT_TOT,                         &
          IOFFRHS(8,MXCASE)
      CHARACTER(Len=32) HZERO
      CHARACTER(LEN=4) TITLE(18,MxTIT)
      CHARACTER(LEN=LenIn+8) NAME (2*MXORB)
      CHARACTER(LEN=2) HEADER (72)

      CHARACTER(Len=8) CASES(MXCASE),FOCKTYPE,SMATRIX,SDECOM,           &
          BMATRIX,BTRANS,BSPECT,ORBIT,ORBIN,OUTFMT,                     &
          ORBNAM(MXORB),IINAM(MXINA),IANAM(MXACT),ISNAM(MXEXT)

      real(kind=wp) CPUGIN,CPUINT,CPUFMB,                               &
          CPUSIN,CPUFG3,                                                &
          CPUPT2,CPUSBM,CPUEIG,CPUNAD,CPURHS,CPUSER,                    &
          CPUPCG,CPUSCA,CPULCS,CPUOVL,CPUVEC,CPUSGM,                    &
          CPUPRP,CPUGRD,                                                &
          TIOGIN,TIOINT,TIOFMB,                                         &
          TIOSIN,TIOFG3,                                                &
          TIOPT2,TIOSBM,TIOEIG,TIONAD,TIORHS,TIOSER,                    &
          TIOPCG,TIOSCA,TIOLCS,TIOOVL,TIOVEC,TIOSGM,                    &
          TIOPRP,TIOGRD
End Module CASPT2_Module
