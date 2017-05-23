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
*
*                          ^   ~
*     Constructs  F  = <0|[Q  ,H]|0>
*                  pq       pq
*
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "glbbas_mclr.fh"
#include "spin.fh"
      Real*8 rkappa(nDensC),Sigma(nDensC),
     &       Focki(ndens2),rMOs(*),rmoa(*)
*
*     D,FA used in oit of FA
*     p11 -> (~~|  )
*     p12 -> (  |~~)
*     sign1 k**t=sign K
*     sign2=Sign K * Lambda
*
*
      Call GetMem('MOTemp2','ALLO','REAL',ipMT1,nmba)
      Call GetMem('MOTemp1','ALLO','REAL',ipMT2,nmba)
      Call GetMem('MOTemp3','ALLO','REAL',ipMT3,nmba)
      Call GetMem('Sigma','ALLO','REAL',ipScr,ndensc)
*
      Call Oit_sp(rkappa,Work(ipScr),-1,1,
     &            -1.0d0,Work(ipG2mp),1.0d0,Work(ipfm),
     &            Work(ipg1m),Work(ipfamo_spinm),
     &            Work(ipMT1),Work(ipMT2),Focki)
*     Call DYAX(ndensc,rbetaA/2.0d0,Work(ipSCR),1,sigma,1)
      Call DYAX(ndensc,1.0d0,Work(ipSCR),1,sigma,1)
      Call Recprt(' ',' ',Work(ipSCR),ndensc,1)
*
*     kappa_S
*
      Call Oit_sp(rkappa,Work(ipScr),1,1,
     &            -1.0d0,Work(ipG2mp),1.0d0,Work(ipfm),
     &            Work(ipg1m),Work(ipfamo_spinm),
     &            Work(ipMT1),Work(ipMT2),Focki)
*     call daxpy_(ndensc,-rbetaA/2.0d0,Work(ipSCR),1,sigma,1)
      call daxpy_(ndensc,-1.0d0,Work(ipSCR),1,sigma,1)
      Call Recprt(' ',' ',Work(ipSCR),ndensc,1)
*
*     alpha_S
*
      If (rbetas.ne.0.0d0) Then
        Call Oit_sp(rkappa,Work(ipScr),-1,-1,
     &            1.0d0,Work(ipG2pp),1.0d0,Work(ipfp),
     &            Work(ipg1p),Work(ipfamo_spinp),
     &            Work(ipMT1),Work(ipMT2),Focki)
*     call daxpy_(ndensc,rbetaS,Work(ipSCR),1,sigma,1)
      call daxpy_(ndensc,1.0d0,Work(ipSCR),1,sigma,1)
      Call Recprt(' ',' ',Work(ipSCR),ndensc,1)
      End if
      Call Oit_sp(rkappa,Work(ipScr),-1,-1,
     &            1.0d0,Work(ipG2pp),1.0d0,Work(ipG2mm),
     &            Work(ipg1p),Work(ipFamo_spinp),
     &            Work(ipMT1),Work(ipMT2),Focki)
*     call daxpy_(ndensc,ralphas,Work(ipSCR),1,sigma,1)
      Call Recprt(' ',' ',Work(ipSCR),ndensc,1)

      Call AddGrad_sp(rKappa,Work(ipScr),Work(ipFS),1,1.0d0,1.0d0)
      Call Recprt(' ',' ',Work(ipSCR),ndensc,1)
      call daxpy_(ndensc,rbetaA/2.0d0,Work(ipSCR),1,sigma,1)
*
      Call DZAXPY(nmba,1.0d0,Work(ipMT1),1,
     &           Work(ipMT2),1,Work(ipMT3),1)
      Call PickMO_MCLR(Work(ipMT3),rmos,1)
*
      Call DZAXPY(nmba,-1.0d0,Work(ipMT2),1,
     &           Work(ipMT1),1,Work(ipMT3),1)
      Call PickMO_MCLR(Work(ipMT3),rmoa,1)

      Call GetMem('Sigma','FREE','REAL',ipScr,ndensc)
      Call GetMem('MOTemp3','FREE','REAL',ipMT3,nmba)
      Call GetMem('MOTemp2','FREE','REAL',ipMT2,nmba)
      Call GetMem('MOTemp1','FREE','REAL',ipMT1,nmba)
*
      return
      end
