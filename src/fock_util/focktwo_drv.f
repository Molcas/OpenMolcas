************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine FockTwo_Drv(nSym,nBas,nAux,Keep,
     &                       DLT,DSQ,FLT,nFLT,
     &                       ExFac,nBSQT,nBMX)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 DLT(*),DSQ(*),FLT(nFLT)
      Integer nBas(8), nAux(8), Keep(8)
      Logical DoCholesky,GenInt
      Integer ALGO
      Logical REORD,DECO
      Real*8 CMO_DUMMY(1)

      Common /CHORAS / REORD,DECO,ALGO
*
* nAux is the number of occupied orbitals
      GenInt=.false.
      DoCholesky=.false.
      if(ALGO.eq.0) GenInt=.true. !use GenInt to regenerate integrals

      Call DecideOnCholesky(DoCholesky)

      Call GetMem('LWFSQ','Allo','Real',LWFSQ,NBSQT)
      call dcopy_(NBSQT,Zero,0,Work(LWFSQ),1)

      if((.not.DoCholesky).or.(GenInt)) then
      Call GetMem('LW2','Allo','Real',LW2,NBMX*NBMX)
      end if
*
      Call Allocate_Work(ipTemp,nFlt)
      Call FZero(Work(ipTemp),nFlt)
*
      Call GetMem('LW1','MAX','Real',LW1,LBUF)
*
* Standard building of the Fock matrix from Two-el integrals
*

      IF (.not.DoCholesky) THEN
         Call GetMem('LW1','Allo','Real',LW1,LBUF)

      If (LBUF.LT.1+NBMX**2) Then
         WRITE(6,*)' FockTwo_Drv Error: Too little memory remains for'
     &     //' the call to FOCKTWO.'
         WRITE(6,*)' Largest allocatable array size LBUF=',LBUF
         WRITE(6,*)' Max nr of bf in any symmetry,  NBMX=',NBMX
         WRITE(6,*)' Required minimum size     1+NBMX**2=',1+NBMX**2
         WRITE(6,*)'    (All in Real*8-size words)'
         Call QTRACE()
         Call  ABEND()
      End If
*
      Call FOCKTWO(nSym,nBas,nAux,Keep,
     &             DLT,DSQ,Work(ipTemp),nFlt,
     &             Work(LWFSQ),LBUF,Work(LW1),Work(LW2),ExFac)

      ENDIF
*
*
* Building of the Fock matrix regenerating the integrals on the fly
*
      IF (DoCholesky.and.GenInt) THEN ! save some space for GenInt
           LBUF = MAX(LBUF-LBUF/10,0)
           Call GetMem('LW1','Allo','Real',LW1,LBUF)

      If (LBUF.LT.1+NBMX**2) Then
         WRITE(6,*)' FockTwo_Drv Error: Too little memory remains for'
     &     //' the call to FOCKTWO.'
         WRITE(6,*)' Largest allocatable array size LBUF=',LBUF
         WRITE(6,*)' Max nr of bf in any symmetry,  NBMX=',NBMX
         WRITE(6,*)' Required minimum size     1+NBMX**2=',1+NBMX**2
         WRITE(6,*)'    (All in Real*8-size words)'
         Call QTRACE()
         Call  ABEND()
      End If
*
      Call FOCKTWO(nSym,nBas,nAux,Keep,
     &             DLT,DSQ,Work(ipTemp),nFlt,
     &             Work(LWFSQ),LBUF,Work(LW1),Work(LW2),ExFac)

      ENDIF

*
* Building of the Fock matrix directly from Cholesky vectors
*
      IF (DoCholesky .and. .not.GenInt) THEN

*
* CMO_DUMMY is required call argument of choras_drv:
* (Not used, see logical flags in choras_drv)
*      SUBROUTINE CHORAS_DRV(nSym,nBas,nOcc,DSQ,DLT,FLT,
*     &                      ExFac,LWFSQ,CMO)
       CALL CHOras_drv(nSym,nBas,nAux,DSQ,DLT,
     &                 Work(ipTemp),ExFac,LWFSQ,CMO_DUMMY)
*
      ENDIF
*


      Call DaXpY_(nFlt,One,Work(ipTemp),1,FLT,1)
*
      Call Free_Work(ipTemp)

      if(.not.DoCholesky)then
      Call GetMem('LW1','Free','Real',LW1,LBUF)
      Call GetMem('LW2','Free','Real',LW2,NBMX*NBMX)
      endif

      Call GetMem('LWFSQ','Free','Real',LWFSQ,NBSQT)
*
      Return
      End
