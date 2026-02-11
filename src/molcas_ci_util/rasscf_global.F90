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

module rasscf_global

use Molcas, only: LenIn, MxAct, MxOrb, MxRoot
use RASDim, only: MxIter, MxRef
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: ITRIM = mxAct*(mxAct+1)/2

integer(kind=iwp) :: hfocc(mxact), hRoots, IADR15(30), iAlphaBeta, IBLB(8), IBLBM, ICI(mxRoot,mxRef), ICICH, ICICP, iCIonly, &
                     ICIRFROOT, ICIRST, ICMSIterMax, ICMSIterMin, ICMSP, IEXPAND, IFORDE, INOCALC, iOrbOnly, iOrbTyp, IORDEM, &
                     iOverwr, IPCMROOT, IPR, IPT2, iRlxRoot, IROOT(mxRoot), IRotPsi, ISAVE_EXP, ISCF, ISPDEN, ISTORD(9), &
                     ISTORP(9), ISUPSM, ISYMBB, ITCORE, ITER, ITERCI, ITERSX, ITMAX, iToc(64), ITRI(ITRIM), IXMSP, IXSYM(mxOrb), &
                     IZROT(ITRIM), JBLB(8), JBLBM, JCJ(mxRoot,mxRef), KTIGHT, LOWMS, LROOTS, MAXIT, MAXJT, MAXORBOUT, MxDMRG, &
                     n_Det, n_keep, NAC, NACPAR, NACPR2, NDIMSX, NewFock, NFINT, NFR, NIN, NO2M, NORBT, NQUNE, NROOT, NROOTS, &
                     NSEC, NSM(mxOrb), NSXS, NTIT, NTOT3, NTOT4
real(kind=wp) :: CBLB(8), CBLBM = Zero, CCI(mxRoot,mxRef), CMAX, CMSThreshold, CONV(6,mxIter+2), CoreShift, DE, E2act, &
                 ECAS = Zero, EMY, ENER(mxRoot,mxIter+2), ESX, ExFac, FDIAG(mxOrb) = Zero, HALFQ = Zero, HALFQ1 = Zero, LVSHFT, &
                 POTNUC, PRETHR, PROTHR, PRWTHR, RLXGRD, ROTMAX, S, SXSHFT = Zero, THFACT, THRE, THREN, THRSX, THRTE, TMIN, &
                 Tot_Charge, Tot_El_Charge, Tot_Nuc_Charge, VIA = Zero, VIA_DFT = Zero, WEIGHT(mxRoot)
logical(kind=iwp) :: DoBlockDMRG, doDMRG, DoFaro, DOFCIDUMP, IfCRPR, kIvo, l_casdft, lSquare, NonEq, RFpert, Start_Vectors
character(len=LenIn+8) :: BName(mxOrb)
character(len=256) :: CMSStartMat
character(len=80) :: KSDFT, KSDFT_TEMP, TITLE(18)
character(len=8) :: IPHNAME, OutFmt1, OutFmt2, PURIFY, SXSEL
character(len=4) :: DFTFOCK
character(len=3) :: QNUPDT
character(len=2) :: HEADER(72), QNSTEP

public :: BName, CBLB, CBLBM, CCI, CMAX, CMSStartMat, CMSThreshold, CONV, CoreShift, DE, DFTFOCK, DoBlockDMRG, doDMRG, DoFaro, &
          DOFCIDUMP, E2act, ECAS, EMY, ENER, ESX, ExFac, FDIAG, HALFQ, HALFQ1, HEADER, hfocc, hRoots, IADR15, iAlphaBeta, IBLB, &
          IBLBM, ICI, ICICH, ICICP, iCIonly, ICIRFROOT, ICIRST, ICMSIterMax, ICMSIterMin, ICMSP, IEXPAND, IfCRPR, IFORDE, INOCALC, &
          iOrbOnly, iOrbTyp, IORDEM, iOverwr, IPCMROOT, IPHNAME, IPR, IPT2, iRlxRoot, IROOT, IRotPsi, ISAVE_EXP, ISCF, ISPDEN, &
          ISTORD, ISTORP, ISUPSM, ISYMBB, ITCORE, ITER, ITERCI, ITERSX, ITMAX, iToc, ITRI, ITRIM, IXMSP, IXSYM, IZROT, JBLB, &
          JBLBM, JCJ, kIvo, KSDFT, KSDFT_TEMP, KTIGHT, l_casdft, LOWMS, LROOTS, lSquare, LVSHFT, MAXIT, MAXJT, MAXORBOUT, MxDMRG, &
          n_Det, n_keep, NAC, NACPAR, NACPR2, NDIMSX, NewFock, NFINT, NFR, NIN, NO2M, NonEq, NORBT, NQUNE, NROOT, NROOTS, NSEC, &
          NSM, NSXS, NTIT, NTOT3, NTOT4, OutFmt1, OutFmt2, POTNUC, PRETHR, PROTHR, PRWTHR, PURIFY, QNSTEP, QNUPDT, RFpert, RLXGRD, &
          ROTMAX, S, Start_Vectors, SXSEL, SXSHFT, THFACT, THRE, THREN, THRSX, THRTE, TITLE, TMIN, Tot_Charge, Tot_El_Charge, &
          Tot_Nuc_Charge, VIA, VIA_DFT, WEIGHT

#ifdef _DMRG_
integer(kind=iwp) :: MPSCompressM
logical(kind=iwp) :: DoDelChk, domcpdftDMRG, DoNEVPT2Prep, twordm_qcm

public :: DoDelChk, domcpdftDMRG, DoNEVPT2Prep, MPSCompressM, twordm_qcm
#endif

#if defined (_ENABLE_CHEMPS2_DMRG_)
integer(kind=iwp) :: chemps2_lrestart, max_canonical, max_sweep
real(kind=wp) :: chemps2_blb, chemps2_noise, davidson_tol
logical(kind=iwp) :: chemps2_restart

public :: chemps2_blb, chemps2_lrestart, chemps2_noise, chemps2_restart, davidson_tol, max_canonical, max_sweep
#endif

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_)
logical(kind=iwp) :: Do3RDM

public :: Do3RDM
#endif

#ifdef _ENABLE_DICE_SHCI_
integer(kind=iwp) :: dice_iter, dice_sampleN, nref_dice
real(kind=wp) :: dice_eps1, dice_eps2
logical(kind=iwp) :: dice_restart, dice_stoc
character(LEN=500) :: diceocc(20)

public :: dice_eps1, dice_eps2, dice_iter, dice_restart, dice_sampleN, dice_stoc, diceocc, nref_dice
#endif

end module rasscf_global
