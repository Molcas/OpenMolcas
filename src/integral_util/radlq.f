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
* Copyright (C) 1993, Kurt Pfingst                                     *
************************************************************************
      SubRoutine Radlq(Zeta,nZeta,lsum,Rnr,icop)
************************************************************************
*                                                                      *
* Object: to compute the radial part of the continuum  Coulomb         *
*         integrals outside the  R-matrix sphere                       *
*                                                                      *
* Called from: KneInt                                                  *
*                                                                      *
* Author: K.Pfingst 21/5/93                                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "rmat.fh"
#include "real.fh"
#include "nrmf.fh"
      external fradf
      Character*80 Label
      Real*8 Zeta(nZeta), Rnr(nZeta,0:lsum)
      Parameter(limit=200,lenw=4*limit)
      Integer iScrt(limit)
      Real*8 Scrt(lenw)
*
      iRout = 122
      iPrint = nPrint(iRout)
ccccccccccccccccccccccccccccccccccccccc
c     iPrint = 99
ccccccccccccccccccccccccccccccccccccccc
*                                                                      *
************************************************************************
*                                                                      *
      Do ir=0,lsum
         Do iZeta=1,nZeta
            expsum=Zeta(iZeta)
            ier=0
            l=ir-icop
            Call dqagi(fradf,Rmatr,1,Epsabs,Epsrel,result,abser,neval,
     &                 ier,
     &                 limit,lenw,last,iScrt,Scrt)
            If (ier.gt.0) Then
               Call WarningMessage(1,
     &         ' WARNING in Radlq; Consult output for details!')
               write(6,*) ' ier=',ier,
     &                    ' Error in Dqagi called from Radlq.'
               write(6,*) ' result=',result
               write(6,*) ' abser =',abser
               write(6,*) ' neval =',neval
               write(6,*) ' WARNING in Radlq'
            End If
            Rnr(iZeta,ir)=result
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.99) Then
         Write (6,*) ' Result in Radlq'
         Write (Label,'(A)') ' Rnr'
         Call RecPrt(Label,' ',Rnr(1,0),nZeta,lsum+1)
      End If
*
      Return
      End
