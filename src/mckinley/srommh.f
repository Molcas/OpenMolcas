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
* Copyright (C) 1993, Roland Lindh                                     *
************************************************************************
      Subroutine sroMmH(nHer,MmSROH,la,lb,lr)
************************************************************************
*                                                                      *
*  Object: to compute the number of real*8 the kernal routine will     *
*          need for the computation of a matrix element between two    *
*          cartesian Gaussin functions with the total angular momentum *
*          of la and lb (la=0 s-function, la=1 p-function, etc.)       *
*          lr is the order of the operator (this is only used when the *
*          integrals are computed with the Hermite-Gauss quadrature).  *
*                                                                      *
*  Called from: OneEl                                                  *
*                                                                      *
************************************************************************
*
      use Basis_Info, only: dbsc, nCnttp, Shells
#include "itmax.fh"
#include "info.fh"
*
      nElem(i) = (i+1)*(i+2)/2
*
      nOrder = 0
      nordop=lr
      ld=2
      MmSROH = 0
      Do 1960 iCnttp = 1, nCnttp
         If (.Not.dbsc(iCnttp)%ECP) Go To 1960
         Do 1966 iAng = 0, nSRO_Shells(iCnttp)-1
            iShll = ipSRO(iCnttp) + iAng
            nExpi=Shells(iShll)%nExp
            If (nExpi.eq.0) Go To 1966
*
            ip = 0
            nac = nElem(la)*nElem(iAng)
            ncb = nElem(iAng)*nElem(lb)
            ip = ip + nElem(la)*nElem(lb)*21 ! Final

            ip = ip + nExpi*nExpi ! tmp

            ip=ip+10*nac*nExpi ! FA1 & FA2
            ip=ip+10*ncb*nExpi ! FB1 & FB2

            nHer = (la+1+iAng+1+ld)/2
            nOrder = Max(nHer,nOrder)
            iacore=6+3*nHer*(la+1+ld)+3*nHer*(iAng+1)+
     &           3*nHer*(nOrdOp+1)+3*(la+1+ld)*(iAng+1)*(nOrdOp+1)+1

            nHer = (lb+1+iAng+1+ld)/2
            nOrder = Max(nHer,nOrder)
            icoreb=6+3*nHer*(lb+1+ld)+3*nHer*(iAng+1)+
     &           3*nHer*(nOrdOp+1)+3*(lb+1+ld)*(iAng+1)*(nOrdOp+1)+1

            icores = MAX(icoreb,iacore)*nExpi
            MmSROH = Max(MmSROH,ip+icores)
*
 1966    Continue
 1960 Continue
      nHer = nOrder
*
      Return
      End
