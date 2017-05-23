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
* Copyright (C) 2000, Jonna Stalring                                   *
************************************************************************
*
      Subroutine calcerr(kappa,iestate)
*
*-------------------------------------------------
* Jonna 000411
*
* Calculates the derivative d E_j
*                           ------ = 2\sum_pq Fpq kappa_pq
*                            d w
*
* which estimates the error in the SA.
*---------------------------------------------------
*
      Implicit Real*8 (a-h,o-z)
#include "Input.fh"
#include "SysDef.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "glbbas_mclr.fh"
#include "Files_mclr.fh"
#include "sa.fh"
      Real*8 kappa(*)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
      ng1=itri(ntash,ntash)
      ng2=itri(ng1,ng1)
*
      Call Getmem('G1Q ','ALLO','REAL',ipG1Q,ng1)
      Call Getmem('G1R ','ALLO','REAL',ipG1R,ntash**2)
      Call Getmem('G2R ','ALLO','REAL',ipG2R,ntash**4)
      Call Getmem('T   ','ALLO','REAL',ipT,ndens2)
      Call Getmem('pq  ','ALLO','REAL',ipq,ndens2)
c
c Get densities of iestate in triangular and rectangular form.
c
      Call rddj(Work(ipG1R),Work(ipG1Q),Work(ipG2R),iestate)
c
c Get Fock matrix ipT
c
      Call FockGen(1.0d0,Work(ipG1R),Work(ipG2R),
     &             Work(ipT),Work(ipq),1)
c
c Print energy
c
c     Call prnte(Work(ipG1q),Work(ipT))
c
c     Call Recprt('KAPP',' ',kappa,nDens2,1)
c     Call Recprt('FOCK',' ',Work(ipT),nDens2,1)
c
      dejdw = 0.0d0
      Do iS=1,nsym
         dejdw = dejdw + ddot_(nBas(is)**2,Work(ipT+ipmat(is,is)-1),
     &           1,kappa(ipmat(is,is)),1)
      End do
      dejdw = 2.0d0*dejdw
c     is=1
c     write(*,*)'Work(ipT+ipmat(is,is)-1)',Work(ipT+ipmat(is,is)-1)
c     write(*,*)'kappa(ipmat(is,is)+1)',kappa(ipmat(is,is)+1)
c
c     d E_j / d w_i, without constraint
c     Write(6,200) iestate,istate,dejdw
c
      If (iestate.eq.istate) Then
c        d E_i / d w_i, with the constraint sum(w_i)=1,
c        multiplied by 1-w_i... which happens to be equal to
c        d E_i / d w_i, without constraint
c        change sign, because this is the *error*
         Write(6,100) iestate,-dejdw
      End If
100   Format(' **********'/'                 ',
     &       ' Estimated error in the energy of state ',I5,': ',ES12.5/
     &       ' **********')
*200   Format(' Derivative of the energy of state ',I5,
*     &       ' with the weight ',I5,': ',ES12.5)
*
      Call Getmem('pq  ','FREE','REAL',ipq,ndens2)
      Call Getmem('T   ','FREE','REAL',ipT,ndens2)
      Call Getmem('G2R ','FREE','REAL',ipG2R,ntash**4)
      Call Getmem('G1R ','FREE','REAL',ipG1R,ntash**2)
      Call Getmem('G1Q ','FREE','REAL',ipG1Q,ng1)
*
      Return
      End
