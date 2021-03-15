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
      Subroutine FockTwo_Drv(nSym,nBas,nAux,Keep,DLT,DSQ,FLT,nFLT,
     &                       ExFac,nBMX)
      Use Data_Structures, only: DSBA_Type, Allocate_DSBA,
     &                           Deallocate_DSBA
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8 DLT(*),DSQ(*),FLT(nFLT)
      Integer nBas(8), nAux(8), Keep(8)
      Logical DoCholesky,GenInt
      Real*8 CMO_DUMMY(1)

      Type (DSBA_Type) WFSQ

      Real*8, Allocatable:: W1(:), W2(:), Temp(:)

#include "choras.fh"
*
* nAux is the number of occupied orbitals
      GenInt=.false.
      DoCholesky=.false.
      if(ALGO.eq.0) GenInt=.true. !use GenInt to regenerate integrals

      Call DecideOnCholesky(DoCholesky)

      Call Allocate_DSBA(WFSQ,nBas,nBas,nSym)
      WFSQ%A0(:)=Zero

      if((.not.DoCholesky).or.(GenInt)) then
        Call mma_allocate(W2,NBMX**2,Label='W2')
      end if
*
      Call mma_allocate(Temp,nFlt,Label='Temp')
      Temp(:)=Zero
*
      Call mma_maxDBLE(LBUF)
*                                                                      *
************************************************************************
*                                                                      *
*     Standard building of the Fock matrix from Two-el integrals
*
      IF (.not.DoCholesky) THEN
*                                                                      *
************************************************************************
*                                                                      *
         Call mma_allocate(W1,LBUF,Label='W1')

         If (LBUF.LT.1+NBMX**2) Then
            WRITE(6,*)' FockTwo_Drv Error: Too little memory remains'
     &     //'for the call to FOCKTWO.'
            WRITE(6,*)' Largest allocatable array size LBUF=',LBUF
            WRITE(6,*)' Max nr of bf in any symmetry,  NBMX=',NBMX
            WRITE(6,*)' Required minimum size     1+NBMX**2=',1+NBMX**2
            WRITE(6,*)'    (All in Real*8-size words)'
            Call  ABEND()
         End If
*
         Call FOCKTWO(nSym,nBas,nAux,Keep,DLT,DSQ,Temp,nFlt,
     &                WFSQ%A0,LBUF,W1,W2,ExFac)

*                                                                      *
************************************************************************
*                                                                      *
*     Building of the Fock matrix regenerating the integrals on the fly
*
      Else IF (DoCholesky.and.GenInt) THEN ! save some space for GenInt
*                                                                      *
************************************************************************
*                                                                      *
         LBUF = MAX(LBUF-LBUF/10,0)
         Call mma_allocate(W1,LBUF,Label='W1')

         If (LBUF.LT.1+NBMX**2) Then
            WRITE(6,*)' FockTwo_Drv Error: Too little memory remains'
     &     //'for the call to FOCKTWO.'
            WRITE(6,*)' Largest allocatable array size LBUF=',LBUF
            WRITE(6,*)' Max nr of bf in any symmetry,  NBMX=',NBMX
            WRITE(6,*)' Required minimum size     1+NBMX**2=',1+NBMX**2
            WRITE(6,*)'    (All in Real*8-size words)'
            Call  ABEND()
         End If
*
         Call FOCKTWO(nSym,nBas,nAux,Keep,DLT,DSQ,Temp,nFlt,
     &               WFSQ%A0,LBUF,W1,W2,ExFac)

*                                                                      *
************************************************************************
*                                                                      *
*     Building of the Fock matrix directly from Cholesky vectors
*
      Else IF (DoCholesky .and. .not.GenInt) THEN
*                                                                      *
************************************************************************
*                                                                      *
*        CMO_DUMMY is required call argument of choras_drv:
*        (Not used, see logical flags in choras_drv)
*        SUBROUTINE CHORAS_DRV(nSym,nBas,nOcc,DSQ,DLT,FLT,ExFac,WFSQ,
*    &                         CMO)
          CALL CHOras_drv(nSym,nBas,nAux,DSQ,DLT,Temp,ExFac,WFSQ,
     &                    CMO_DUMMY)
*                                                                      *
************************************************************************
*                                                                      *
      ENDIF
*                                                                      *
************************************************************************
*                                                                      *
      Call DaXpY_(nFlt,One,Temp,1,FLT,1)
*
      Call mma_deallocate(Temp)
      If (Allocated(W1))   Call mma_deallocate(W1)
      If (Allocated(W2))   Call mma_deallocate(W2)
      Call Deallocate_DSBA(WFSQ)
*
      Return
      End
