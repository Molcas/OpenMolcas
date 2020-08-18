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
* Copyright (C) 2010, Francesco Aquilante                              *
************************************************************************
      SUBROUTINE Tw_corr_drv(EOrb,nEO,CMO,nCMO,Ecorr)

      Implicit Real*8 (a-h,o-z)
#include "mxdm.fh"
#include "infscf.fh"
#include "WrkSpc.fh"
      Integer nEO, nCMO
      Real*8 EOrb(nEO), CMO(nCMO), Ecorr

      Call GetMem('Eov','Allo','Real',ipEOkk,nEO)

      ipEVir=ipEOkk+nnOc
      iOff=0
      jOff=0
      kOff=0
      Do iSym=1,nSym
         nOkk=nFro(iSym)+nOcc(iSym,1)
         nExt=nBas(iSym)-nDel(iSym)-nOkk
         jOrb=1+jOff
         jOkk=ipEOkk+iOff
         Do i=0,nOkk-1
            Work(jOkk+i)=EOrb(jOrb+i)
         End Do
         jOrb=jOrb+nOkk
         jVir=ipEVir+kOff
         Do i=0,nExt-1
            Work(jVir+i)=EOrb(jOrb+i)
         End Do
         iOff=iOff+nOkk
         jOff=jOff+nOrb(iSym)
         kOff=kOff+nExt
      End Do

      Call Tw_corr(irc,Ecorr,CMO,Work(ipEOkk),Work(ipEVir))

      Call GetMem('Eov','Free','Real',ipEOkk,nEO)
      Return
      End
*****************************************************************************
*                                                                           *
*****************************************************************************
      SUBROUTINE Tw_corr(irc,DeTW,CMOI,EOcc,EVir)

#include "implicit.fh"
#include "WrkSpc.fh"
      Real*8 DeTW, CMOI(*), EOcc(*), EVir(*)
C
      Logical DoDens_
      Integer ChoAlg_
      Integer nExt(8)
#include "mxdm.fh"
#include "infscf.fh"
#include "chomp2_cfg.fh"
      Dimension Grad(1)

      DoDens_= DoDens
      DoDens = .false.
      ChoAlg_= ChoAlg
      ChoAlg = 2
*
      CALL GETMEM('DMAT','ALLO','REAL',ip_DM0,2*nBT)
      ip_DM=ip_DM0+nBT

      nElk=0
      Do i=1,nSym
         nExt(i)=nBas(i)-nDel(i)-nOcc(i,1)-nFro(i)
         nElk=nElk+2*(nFro(i)+nOcc(i,1))
      End Do

      CALL DM_FNO_RHF(irc,nSym,nBas,nFro,nOcc(1,1),nExt,nDel,
     &                    CMOI,EOcc,EVir,Work(ip_DM0),Work(ip_DM))
      If (irc .ne. 0) Then
         Write(6,*) 'DM_FNO_RHF returned ',irc
         Call SysAbendMsg('DM_FNO_RHF',
     &                    'Non-zero return code from DM_FNO_RHF',
     &                    ' ')
      EndIf

      CALL GETMEM('FMAT','ALLO','REAL',ipF_DFT,nBT)
*
      Call Fold_tMat(nSym,nBas,Work(ip_DM),Work(ip_DM))
      call dscal_(nBT,0.5d0,Work(ip_DM),1)
      Call Fold_tMat(nSym,nBas,Work(ip_DM0),Work(ip_DM0))
      call dscal_(nBT,0.5d0,Work(ip_DM0),1)
      Grad=0.0d0

      Call wrap_DrvNQ('HUNTER',Work(ipF_DFT),1,TW,
     &                      Work(ip_DM),nBT,1,
     &                      .false.,
     &                      Grad,1,'SCF ')

      Call wrap_DrvNQ('HUNTER',Work(ipF_DFT),1,TW0,
     &                      Work(ip_DM0),nBT,1,
     &                      .false.,
     &                      Grad,1,'SCF ')
      DeTW=(TW-TW0)/dble(nElk)
*
      CALL GETMEM('FMAT','FREE','REAL',ipF_DFT,nBT)
      CALL GETMEM('DMAT','FREE','REAL',ip_DM0,2*nBT)
*
      Return
      End
*****************************************************************************
*                                                                           *
*****************************************************************************

      SUBROUTINE DM_FNO_RHF(irc,nSym,nBas,nFro,nIsh,nSsh,nDel,
     &                          CMOI,EOcc,EVir,DM0,DM)
*****************************************************************************
*                                                                           *
*     Purpose:  setup of FNO density matrix calculation (RHF-based)         *
*                                                                           *
*     Author:   F. Aquilante  (Geneva, Sep 2010)                            *
*                                                                           *
*****************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
*
      Integer nBas(nSym),nFro(nSym),nIsh(nSym),nSsh(nSym),
     &        nDel(nSym)
      Real*8  CMOI(*), EOcc(*), EVir(*), DM0(*), DM(*)
      Real*8  DeMP2
      Logical MP2_small
      Common / ChFNOPT/ DeMP2, MP2_small
*
      Integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
*
*
      irc=0
      MP2_small=.false.
*
*----------------------------------------------------------------------*
*     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
*----------------------------------------------------------------------*
*
      nBasT=0
      ntri=0
      nSQ=0
      nBmx=0
      nOrb=0
      nVV=0
      Do i=1,nSym
        nBasT=nBasT+nBas(i)
        nOrb=nOrb+nFro(i)+nIsh(i)+nSsh(i)+nDel(i)
        ntri=ntri+nBas(i)*(nBas(i)+1)/2
        nSQ=nSQ+nBas(i)**2
        nVV=nVV+nSsh(i)**2
        nBmx=Max(nBmx,nBas(i))
      End Do
      IF(nBasT.GT.mxBas) then
       Write(6,'(/6X,A)')
     & 'The number of basis functions exceeds the present limit'
       Call Abend
      Endif
*
      NCMO=nSQ
      CALL GETMEM('LCMO','ALLO','REAL',LCMO,2*NCMO)
      iCMO=LCMO+NCMO
      CALL DCOPY_(NCMO,CMOI,1,WORK(LCMO),1)
*
      nOA=0
      Do iSym=1,nSym  ! setup info
         lnFro(iSym)=nFro(iSym)
         lnOcc(iSym)=nIsh(iSym)
         nOA=nOA+lnOcc(iSym)
         lnVir(iSym)=nSsh(iSym)
         lnOrb(iSym)=lnOcc(iSym)+lnVir(iSym)
         lnDel(iSym)=nDel(iSym)
      End Do
*
      Call GetMem('Eorb','Allo','Real',ipOrbE,4*nOrb)
      jOff=0
      kOff=0
      lOff=0
      Do iSym=1,nSym
         jp=ipOrbE+lOff+nFro(iSym)
         jOcc=jOff+1
         call dcopy_(nIsh(iSym),EOcc(jOcc),1,Work(jp),1)
         jVir=kOff+1
         jp=jp+nIsh(iSym)
         call dcopy_(nSsh(iSym),EVir(jVir),1,Work(jp),1)
         jOff=jOff+nIsh(iSym)
         kOff=kOff+nSsh(iSym)
         lOff=lOff+nBas(iSym)
      End Do
      ip_ZZ=ipOrbE
      ipEorb=ipOrbE+nOrb
      ip_Z=ipEorb
      kEOcc=ipEorb+nOrb
      kEVir=kEOcc+nOrb
      ioff=0
      joff=0
      koff=0
      Do iSym=1,nSym
         ifr=ipOrbE+ioff+nFro(iSym)
         ito=kEOcc+joff
         call dcopy_(nIsh(iSym),Work(ifr),1,Work(ito),1)
         ifr=ipOrbE+ioff+nFro(iSym)+nIsh(iSym)
         ito=kEVir+koff
         call dcopy_(nSsh(iSym),Work(ifr),1,Work(ito),1)
         ioff=ioff+nBas(iSym)
         joff=joff+nIsh(iSym)
         koff=koff+nSsh(iSym)
      End Do
      Call GetMem('Dmat','Allo','Real',ip_X,nVV+nOA)
      ip_Y=ip_X+nVV
      Call FZero(Work(ip_X),nVV+nOA)
*
      Call FnoSCF_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,ip_X,ip_Y)
      Call FZero(Work(iCMO),NCMO)
      iOff=0
      Do iSym=1,nSym
         kfr=LCMO+iOff+nBas(iSym)*nFro(iSym)
         kto=iCMO+iOff+nBas(iSym)*lnFro(iSym)
         call dcopy_(nBas(iSym)*lnOcc(iSym),Work(kfr),1,Work(kto),1)
         kfr=LCMO+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
         kto=kto+nBas(iSym)*lnOcc(iSym)
         call dcopy_(nBas(iSym)*lnVir(iSym),Work(kfr),1,Work(kto),1)
         iOff=iOff+nBas(iSym)**2
      End Do
      Call Check_Amp_SCF(nSym,lnOcc,lnVir,iSkip)
      If (iSkip.gt.0) Then
         Call ChoMP2_Drv(irc,Dummy,Work(iCMO),Work(kEOcc),Work(kEVir))
         If(irc.ne.0) then
           write(6,*) 'MP2 pseudodensity calculation failed !'
           Call Abend
         Endif
      Else
         write(6,*)
         write(6,*)'There are ZERO amplitudes T(ai,bj) with the given '
         write(6,*)'combinations of inactive and virtual orbitals !! '
         write(6,*)'Check your input and rerun the calculation! Bye!!'
         Call Abend
      Endif
*
*
*     Compute the correlated density in AO basis
*     -------------------------------------------------------------
      jOcc=ip_X+nVV
c           write(6,*) ' Occ    : ',(Work(jOcc+j),j=0,nOA-1)
c           write(6,*) ' Sum    : ',ddot_(nOA,1.0d0,0,Work(jOcc),1)
      call dscal_(nOA,2.0d0,Work(jOcc),1)
      Call daxpy_(nOA,2.0d0,[1.0d0],0,Work(jOcc),1)
*
      iOff=0
      jOff=0
      kDM=1
      Do iSym=1,nSym
*
         kto=LCMO+jOff
         nOkk=nFro(iSym)+nIsh(iSym)
         Call DGEMM_Tri('N','T',nBas(iSym),nBas(iSym),nOkk,
     &                      2.0d0,Work(kto),nBas(iSym),
     &                            Work(kto),nBas(iSym),
     &                      0.0d0,DM0(kDM),nBas(iSym))
*
         sqocc=sqrt(2.0d0)
         call dscal_(nBas(iSym)*nFro(iSym),sqocc,Work(kto),1)
         Do j=0,nIsh(iSym)-1
             sqocc=sqrt(Work(jOcc+j))
             ito=kto+nBas(iSym)*j
             call dscal_(nBas(iSym),sqocc,Work(ito),1)
         End Do
         Call DGEMM_Tri('N','T',nBas(iSym),nBas(iSym),nOkk,
     &                      1.0d0,Work(kto),nBas(iSym),
     &                            Work(kto),nBas(iSym),
     &                      0.0d0,DM(kDM),nBas(iSym))
*
         if (nSsh(iSym).gt.0) then
           jD=ip_X+iOff
*     Eigenvectors will be in increasing order of eigenvalues
           Call Eigen_Molcas(nSsh(iSym),Work(jD),Work(ip_Z),Work(ip_ZZ))
*     Reorder to get relevant eigenpairs first
           Do j=1,nSsh(iSym)/2
              Do i=1,nSsh(iSym)
                 lij=jD-1+nSsh(iSym)*(j-1)+i
                 kij=jD-1+nSsh(iSym)**2-(nSsh(iSym)*j-i)
                 tmp=Work(lij)
                 Work(lij)=Work(kij)
                 Work(kij)=tmp
              End Do
              tmp=Work(ip_Z-1+j)
              Work(ip_Z-1+j)=Work(ip_Z+nSsh(iSym)-j)
              Work(ip_Z+nSsh(iSym)-j)=tmp
           End Do
*
*     Compute new MO coeff. : X=C*U
           kfr=iCMO+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
           kto=LCMO+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
           Call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),
     &                        1.0d0,Work(kfr),nBas(iSym),
     &                              Work(jD),nSsh(iSym),
     &                        0.0d0,Work(kto),nBas(iSym))

c           write(6,*) ' Occ_vir: ',(Work(ip_Z+j),j=0,nSsh(iSym)-1)
c           write(6,*) ' Sum_vir: ',ddot_(nSsh(iSym),1.0d0,0,Work(ip_Z),1)
           Do j=0,nSsh(iSym)-1
              sqocc=sqrt(2.0d0*Work(ip_Z+j))
              jto=kto+nBas(iSym)*j
              call dscal_(nBas(iSym),sqocc,Work(jto),1)
           End Do
           Call DGEMM_Tri('N','T',nBas(iSym),nBas(iSym),nSsh(iSym),
     &                        1.0d0,Work(kto),nBas(iSym),
     &                              Work(kto),nBas(iSym),
     &                        1.0d0,DM(kDM),nBas(iSym))

           iOff=iOff+nSsh(iSym)**2
         endif
         jOff=jOff+nBas(iSym)**2
         kDM =kDM +nBas(iSym)*(nBas(iSym)+1)/2
         jOcc=jOcc+nIsh(iSym)
      End Do
*
      Call GetMem('Eorb','Free','Real',ipOrbE,4*nOrb)
      Call GetMem('Dmat','Free','Real',ip_X,nVV+nOA)
      CALL GETMEM('LCMO','FREE','REAL',LCMO,2*NCMO)
*
      Return
      End
************************************************************************
*                                                                      *
************************************************************************
      Subroutine Check_Amp_SCF(nSym,nOcc,nVir,iSkip)

      Implicit Real*8 (a-h,o-z)
      Integer nSym, nOcc(nSym), nVir(nSym), iSkip
      Integer nT1amTot, nT1am(8)

      MulD2h(i,j)=iEor(i-1,j-1) + 1

      iSkip=0
      nT1amTot=0
      Do iSym = 1,nSym
         nT1am(iSym) = 0
         Do iSymi = 1,nSym
            iSyma = MulD2h(iSymi,iSym)
            nT1am(iSym) = nT1am(iSym)
     &                  + nVir(iSyma)*nOcc(iSymi)
         End Do
         nT1amTot = nT1amTot + nT1am(iSym)
      End Do

      If (nT1amTot .gt. 0) iSkip=1
      Return
      End
************************************************************************
*                                                                      *
************************************************************************
      SubRoutine FnoSCF_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,
     &                         ip_X,ip_Y)
C
C     Purpose: put info in MP2 common blocks.
C
#include "implicit.fh"
      Integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
      Integer ip_X, ip_Y
C
#include "corbinf.fh"
#include "chomp2_cfg.fh"
C
C
      nSym = mSym
C
      Do iSym = 1,nSym
         nOrb(iSym) = lnOrb(iSym)
         nOcc(iSym) = lnOcc(iSym)
         nFro(iSym) = lnFro(iSym)
         nDel(iSym) = lnDel(iSym)
         nExt(iSym) = lnVir(iSym)
      End Do
C
      DoFNO=.true.
      ip_Dab=ip_X
      ip_Dii=ip_Y
      l_Dab=nExt(1)
      l_Dii=nOcc(1)
      Do iSym=2,nSym
         l_Dab=l_Dab+nExt(iSym)**2
         l_Dii=l_Dii+nOcc(iSym)
      End Do
C
      Return
      End
