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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      SubRoutine RInt_SP(rkappa,rmos,rmoa,Focki,Sigma)
      use Arrays, only: FAMO_SpinP, FAMO_SpinM, SFock,
     &                  G2mm, G2mp, G2pp, Fp, Fm, G1p, G1m
*
*                          ^   ~
*     Constructs  F  = <0|[Q  ,H]|0>
*                  pq       pq
*
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "spin.fh"
      Real*8 rkappa(nDensC),Sigma(nDensC), Focki(ndens2),rMOs(*),rmoa(*)
      Real*8, Allocatable:: MT1(:), MT2(:), MT3(:), Scr(:)
*
*     D,FA used in oit of FA
*     p11 -> (~~|  )
*     p12 -> (  |~~)
*     sign1 k**t=sign K
*     sign2=Sign K * Lambda
*
*
      Call mma_allocate(MT1,nmba,Label='MT1')
      Call mma_allocate(MT2,nmba,Label='MT2')
      Call mma_allocate(MT3,nmba,Label='MT3')
      Call mma_allocate(Scr,ndensc,Label='Scr')
*
      Call Oit_sp(rkappa,Scr,-1,1,
     &            -One,G2mp,One,Fm,
     &            G1m,FAMO_Spinm,
     &            MT1,MT2,Focki)
*     Call DYAX(ndensc,rbetaA/Two,SCR,1,sigma,1)
      Call DYAX(ndensc,One,SCR,1,sigma,1)
      Call Recprt(' ',' ',SCR,ndensc,1)
*
*     kappa_S
*
      Call Oit_sp(rkappa,Scr,1,1,
     &            -One,G2mp,One,Fm,
     &            G1m,FAMO_Spinm,
     &            MT1,MT2,Focki)
*     call daxpy_(ndensc,-rbetaA/Two,SCR,1,sigma,1)
      call daxpy_(ndensc,-One,SCR,1,sigma,1)
      Call Recprt(' ',' ',SCR,ndensc,1)
*
*     alpha_S
*
      If (rbetas.ne.Zero) Then
        Call Oit_sp(rkappa,Scr,-1,-1,
     &            One,G2pp,One,Fp,
     &            G1p,FAMO_Spinp,
     &            MT1,MT2,Focki)
*     call daxpy_(ndensc,rbetaS,SCR,1,sigma,1)
      call daxpy_(ndensc,One,SCR,1,sigma,1)
      Call Recprt(' ',' ',SCR,ndensc,1)
      End if
      Call Oit_sp(rkappa,Scr,-1,-1,
     &            One,G2pp,One,G2mm,
     &            G1p,Famo_spinp,
     &            MT1,MT2,Focki)
*     call daxpy_(ndensc,ralphas,SCR,1,sigma,1)
      Call Recprt(' ',' ',SCR,ndensc,1)

      Call AddGrad_sp(rKappa,Scr,SFock,1,One,One)
      Call Recprt(' ',' ',SCR,ndensc,1)
      call daxpy_(ndensc,rbetaA/Two,SCR,1,sigma,1)
*
      Call DZAXPY(nmba,One,MT1,1,MT2,1,MT3,1)
      Call PickMO_MCLR(MT3,rmos,1)
*
      Call DZAXPY(nmba,-One,MT2,1,MT1,1,MT3,1)
      Call PickMO_MCLR(MT3,rmoa,1)

      Call mma_deallocate(Scr)
      Call mma_deallocate(MT3)
      Call mma_deallocate(MT2)
      Call mma_deallocate(MT1)
*
      return
      End
