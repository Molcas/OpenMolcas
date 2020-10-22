************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996, Anders Bernhardsson                              *
************************************************************************
      SubRoutine CIDIA_TD(iSym)
      use Exp, only: nexp, nexp_max
      Implicit Real*8 (a-h,o-z)
#include "detdim.fh"
#include "crun_mclr.fh"
#include "cicisp_mclr.fh"
#include "spinfo_mclr.fh"
#include "csfbas_mclr.fh"
#include "incdia.fh"
#include "WrkSpc.fh"
#include "Input.fh"
#include "Pointers.fh"
#include "negpre.fh"
      Integer iSM(1),LSPC(1),iSPC(1),IDUM(1)
*
*     This is just a interface to hide Jeppe from the rest of the world
*     we dont want to let world see the work of the danish
*     (I hope he never reads that)
*     Anyway concerning the CSF/SD stuff.
*     If we work with spin dependent perturbations
*     we never use CSF's (to complicated), instead we use
*     SD in all parts of the program,
*     otherwise we will switch to SD representation in this routine
*
      NSPC=1
      ISPC(1)=1
      iSM(1)=iSym
      IAMCMP=0
      ICISTR=1
      i=2
      If (isym.eq.state_sym) i=1
      If (NOCSF.eq.0) Then
        ncsf1=NCSASM(ISYM)
        nsd=max(ncsf(isym),nint(XISPSM(ISYM,1)))
*       Call GetMem('DIAGCSF','ALLO','REAL',ipDCSF,nSD)
        ipdcsfi=ipget(nsd)
        ipdcsf=ipin(ipdcsfi)
        Call GetMem('DIAGSD','ALLO','REAL',ipDSD,nSD)
      Else
        nsd=max(ncsf(isym),nint(XISPSM(ISYM,1)))
*       Call GetMem('DIAGSD','ALLO','REAL',ipDSD,nSD)
        ipDSDi=ipGet(nsd)
        ipdsd=ipin(ipdsdi)
      End If
      If (nocsf.eq.0) Then
         nD=NCSASM(ISYM)
         ipDia=ipDCSF
         ipdiai=ipdcsfi
      Else
        nD=idint(XISPSM(ISYM,ISPC(1)))
        ipDIA=ipDSD
        ipdiai=ipdsdi
      End If
      LSPC(1)=nSD
      Call IntDia(Work(ipDSD),NSPC,ISPC,ISM,LSPC,
     &           IAMCMP,rin_ene+potnuc)
      If (Nocsf.ne.1)
     &Call CSDIAG(Work(ipDCSF),Work(ipDSD),
     &              NCNATS(1,ISYM),NTYP,
     &              iWORK(KICTS(i)),NDPCNT,NCPCNT,0,
     &              0,IDUM,IPRNT)

      If (nocsf.eq.0)
     & Call GetMem('DIAGSD','FREE','REAL',ipDSD,nSD)
*     Calculate explicit part of hamiltonian
*
      np2=Min(nd,nexp_max)
      np1=0
      nq=0
      If (np2.ne.0) Then
       irc=ipnout(ipdiai)
       call h0(Work(ipin(ipdiai)),np1,nexp_max,nq,isym,nexp,TimeDep)
      Else
        nexp=0
      End if
*     Call Getmem('KICONF1','FREE','INTEGER',kiconf(i),lConf)

      ip=ipin(ipdiai)

      ECAS=ERASSCF(1)
      Do iC=1,nD
        Work(ip+iC-1)=(Work(ip+ic-1)-ECAS)
      End do

      ipdia=ipdiai
      RETURN
      END
