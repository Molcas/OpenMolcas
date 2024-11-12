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
!
! (TODO) the rasscf_data module should be merge into this module
Module rasscf_global
  use definitions, only: wp, iwp
  Private
#include "rasdim.fh"
  REAL(kind=wp), Public::D0,D1,D2,D4
!
  INTEGER(kind=iwp), PARAMETER, Public :: ITRIM=mxAct*(mxAct+1)/2
  INTEGER(kind=iwp), Public::  IW,IPR,IPRINT,IPRDIA,ITERCI,ITERSX,ITER,MAXIT,NROOTS,IROOT(mxRoot),NAC,       &
                               ITRI(ITRIM),IADR15(30), IBLBM,JBLBM,ISYMBB,ISCF,NQUNE
  CHARACTER(LEN=LENIN8), Public:: NAME(mxOrb)
  CHARACTER(LEN=2), Public:: HEADER(72)
  CHARACTER(LEN=2), Public:: QNSTEP
  CHARACTER(LEN=3), Public:: QNUPDT
  Character(LEN=80), Public:: KSDFT, KSDFT_TEMP
  Character(LEN=4), Public:: DFTFOCK
!
  REAL(kind=wp), Public:: ENER(mxRoot,mxIter+2),CONV(6,mxIter+2),FDIAG(mxOrb)=0.0D0,THRE,THRTE,THRSX,ROTMAX,       &
                          ECAS=0.0D0,CMAX,WEIGHT(mxRoot),DE,CBLBM=0.0D0,THREN,THFACT,TMIN,PRETHR,PROTHR,                 &
                          Tot_Charge,Tot_Nuc_Charge,Tot_El_Charge,ExFac,E2act,      &
                          CoreShift
  REAL(kind=wp), Public:: VIA_DFT=0.0D0, HALFQ=0.0D0, VIA=0.0D0, HALFQ1=0.0D0
!
  INTEGER(kind=iwp), Public:: LROOTS,ICICH,IDIAG,ICIRST,KAVER,KSYM(4),MAXJT,    &
                              IXSYM(mxOrb),IPRSEC(7),MAXORBOUT,                 &
                              ICI(mxRoot,mxRef),JCJ(mxRoot,mxRef),              &
                              NFR,NIN,NSEC,NTIT,NO2M,                           &
                              NACPAR,NACPR2,ISTORD(9),NFINT,                    &
                              NORBT,NTOT3,ISTORP(9),NTOT4,ICICP,                &
                              ITMAX,IPT2,ISPDEN,LOWMS,                          &
                              ISUPSM,IORDEM,IFORDE,IPCMROOT,ICIRFROOT,          &
                              iRlxRoot,n_Det,iAlphaBeta,                        &
                              MxDMRG,INOCALC,ISAVE_EXP,IEXPAND,ITCORE,          &
                              hRoots,n_keep,hfocc(mxact)
!
  Logical(kind=iwp), Public:: RFpert,lSquare, Start_Vectors, NonEq, DoFaro,     &
                              DOFCIDUMP,DoBlockDMRG,DoGradPDFT,DoNOGRAD,        &
                              DoGSOR, kIvo, IfCRPR, doDMRG, l_casdft

!
  REAL(kind=wp), Public:: PRWTHR,POTNUC,CCI(mxRoot,mxRef),ECAS1=0.0D0,RLXGRD,EVAC=0.0D0
  CHARACTER(LEN=80), Public:: TITLE(18)
!
  Character(LEN=8), Public:: IPHNAME,OutFmt1,OutFmt2,SXSEL,PURIFY

!
  INTEGER(kind=iwp), Public:: iCIonly,NSM(mxOrb),KTIGHT
!
  INTEGER(kind=iwp), Public:: IRotPsi,IXMSP,ICMSP,ICMSIterMax,ICMSIterMin
  REAL(kind=wp), Public:: CMSThreshold
  CHARACTER(LEN=256),Public::  CMSStartMat
!
  REAL(kind=wp), Public::          EMY,S
!
  INTEGER(kind=iwp), Public:: IBLB(8),JBLB(8),NSXS,NROOT,NDIMSX,IZROT(ITRIM),NewFock
!
  REAL(kind=wp), Public:: CBLB(8),LVSHFT,ESX,SXSHFT=0.0D0
!
  INTEGER(kind=iwp), Public:: iOrbTyp,iOrbOnly,iOrbRoot
!
  Character(LEN=8), Public:: FnJob,FnOrb
  Integer(kind=iwp), Public:: LuJob,LuOrb
!
  INTEGER(kind=iwp), Public:: iToc(64)
!
  INTEGER(kind=iwp), Public:: iOverwr
!
  REAL(kind=wp), Public:: Acc,Bcc,Aoo,Boo,Avv,Bvv
!

#ifdef _DMRG_
  Logical(kind=iwp), Public:: domcpdftDMRG, twordm_qcm
#endif

#ifdef _ENABLE_CHEMPS2_DMRG_
  REAL(kind=wp), Public:: davidson_tol,chemps2_blb,chemps2_noise
  INTEGER(kind=iwp), Public:: max_sweep,chemps2_lrestart,max_canonical
  LOGICAL(kind=iwp), Public:: chemps2_restart
#endif
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_)
  LOGICAL(kind=iwp), Public:: Do3RDM
#endif
#ifdef _ENABLE_DICE_SHCI_
  REAL(kind=wp), Public:: dice_eps1,dice_eps2
  INTEGER(kind=iwp), Public:: nref_dice,dice_sampleN,dice_iter
  LOGICAL(kind=iwp), Public:: dice_stoc,dice_restart
  CHARACTER(LEN=500), Public:: diceocc(20)
#endif
#ifdef _DMRG_
!DMRG-NEVPT2 variables: MPS compression, 4-RDM evaluation
Integer(kind=iwp), Public :: MPSCompressM
Logical(kind=iwp), Public :: DoNEVPT2Prep, DoDelChk
#endif
End Module rasscf_global
