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
      SubRoutine CIDIA_sa(iSym,ralp,S)
      use Arrays, only: CNSM
      use ipPage, only: W
      use negpre
      Implicit Real*8 (a-h,o-z)
#include "detdim.fh"
#include "crun_mclr.fh"
#include "cicisp_mclr.fh"
#include "spinfo_mclr.fh"
#include "incdia.fh"

#include "Input.fh"
#include "Pointers.fh"
#include "sa.fh"
      Integer iSM(1),LSPC(1),iSPC(1),IDUM(1)
      Real*8 ralp(*),S(*)
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
         ipdcsfi=ipget(nsd)
         ipdcsf=ipin(ipdcsfi)
         ipDSDi=ipGet(nSD)
      Else
         nsd=max(ncsf(isym),nint(XISPSM(ISYM,1)))
         ipDSDi=ipGet(nsd)
         ipdsd=ipin(ipdsdi)
      End If

      If (nocsf.eq.0) Then
         nD=NCSASM(ISYM)
         ipdiai=ipdcsfi
      Else
         nD=idint(XISPSM(ISYM,ISPC(1)))
         ipdiai=ipdsdi
      End If
      LSPC(1)=nSD

      irc=ipin(ipDSDi)
      Call IntDia(W(ipDSDi)%Vec,NSPC,ISPC,ISM,LSPC,IAMCMP,
     &            rin_ene+potnuc)

      If (Nocsf.ne.1) Call CSDIAG(W(ipDCSFi)%Vec,W(ipDSDi)%Vec,
     &                            NCNATS(1,ISYM),NTYP,
     &                            CNSM(i)%ICTS,NDPCNT,NCPCNT,0,
     &                            0,IDUM,IPRNT)

      If (nocsf.eq.0) irc=ipClose(ipDSDi)
*     Calculate explicit part of hamiltonian
*
      ipdia=ipdiai

      If (FANCY_PRECONDITIONER) Then
         irc=ipin(ipdia)
         Call SA_PREC(S,W(ipdia)%Vec)
      Else
         irc=ipin(ipdiai)
         irc=ipin(ipCI)
         ip2=1
         Do j=1,nroots
            ECAS=ERASSCF(j)
            We=Weight(j)
            ralp(j)=0.0d0
            Do i=1,ncsf(State_SYM)
               ralp(j)=ralp(j)+1.0d0/(W(ipdiai)%Vec(i)-ECAS)*We*
     &                   W(ipCI)%Vec(ip2)**2
               ip2=ip2+1
            End Do
         End Do
      End If

      RETURN
      END

