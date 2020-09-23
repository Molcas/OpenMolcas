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
      Subroutine NEMO_Opt1()
      use Basis_Info
      Implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "itmax.fh"
#include "info.fh"
#include "warnings.fh"
#include "rinfo.fh"
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
      Integer nBas_Prim(0:7), nBas_cont(0:7), lOper(3)
      Parameter(MxMltPl=10)
      Integer ipMP((MxMltPl+1)*(MxMltPl+2)*(MxMltPl+3)/6),
     &        iSm((MxMltPl+1)*(MxMltPl+2)*(MxMltPl+3)/6)
      Integer ip(3), iSml(3)
      Character*8 Label
      Real*8, Dimension(:), Allocatable :: P_Matrix, MP_Matrix
      Dimension Length(1),nInt(1)
*
      iRout=77
      iPrint=nPrint(iRout)
      Call QEnter('NEMO_Opt1')
*                                                                      *
************************************************************************
*                                                                      *
*     Save basis set info from contracted run
*
      if(iprint.ge.10) write(6,*) ' In NEMO_Opt1', ncnttp
      kCof=0
      kAng=0
      kExp=0
      kC=0
*
*     Normalize coefficients
*
      do iCnttp=1,nCnttp
*
*-- Make a check that no cartesian d or higher have been used.
*   The reason for this restriction is found in tr_prm_cnt.
*   If that routine is generalized, then remove this check.
*
        lSh=0
        kShStr = dbsc(iCnttp)%iVal
        kShEnd = dbsc(iCnttp)%iVal+dbsc(iCnttp)%nVal-1
        Do kSh = kShStr, kShEnd
          If(.not.Shells(kSh)%Transf.and.lSh.ge.2) then
            Call WarningMessage(2,'   NEMO Error')
            Write(6,*)
            Write(6,*)
            Write(6,*)'Error! The NEMO keyword does not work with'
     &//' cartesian d-functions or higher.'
            Write(6,*)'Request spherical functions to proceed.'
            Write(6,*)
            Write(6,*)
            Call Quit(_RC_INPUT_ERROR_)
          Endif
          lSh=lSh+1
        Enddo
*
*-- End check.
*
        Do icnt = 1, dbsc(iCnttp)%nCntr
        kC=kC+1
c           Do iAngr=0,nAngr(icnt)
           Do iAngr=0,nAngr(kC)
c              rI=iAngr+1.0d0+half
               rI=DBLE(iAngr)+One+Half
              kAng=kAng+1
              Do iBas=1,nBasisr(kAng)
                 Sum=Zero
                 kExpi=kExp
                 kCofi=kCof
                 Do iExp=1,nPrimr(kAng)
                    kExpi=kExpi+1
                    kCofi=kCofi+1
                    rExpi=rExp(kExpi)
c                    write(6,'(a11,f20.8)') ' Exponents',rExpi
                    rCofi=rCof(kCofi)
                    kExpj=kExp
                    kCofj=kCof
                    Do jExp=1,nPrimr(kAng)
                       kExpj=kExpj+1
                       kCofj=kCofj+1
                       rExpj=rExp(kExpj)
                       rCofj=rCof(kCofj)
                       Sum=Sum+rCofi*rCofj*
     &                 (Two*sqrt(rExpi*rExpj)/(rExpi+rExpj))**rI
                    End Do
                 End Do
                 rNorm=One/sqrt(Sum)
                 if(iprint.ge.10) write(6,*) ' rNorm', kAng,rNorm
                   Do iExp=1,nPrimr(kAng)
                      rCof(kCof+iExp)=rCof(kCof+iExp)*rNorm
                     if(iprint.ge.10) then
                         write(6,'(a24,f20.6)')
     &                   ' normalized coefficients',
     &                   rCof(kCof+iExp)
                     endif
                   End Do
                 kCof=kCof+nPrimr(kAng)
              End Do
              kExp=kExp+nPrimr(kAng)
           End Do
        End Do
      End Do
*
      If (iPrint.ge.10) Then
      i=0
      Do L=1,nrSym
         write(6,*) ' Irreducible representation', L
         Do ibas=1,nrBas(L)
            i=i+1
            Write (6,'(20i4)') i, icent(i),lnang(i),lmag(i)
         End Do
      End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Close ONEINT and re-open ONEREL
*
      Call iCopy(8,nBas,1,nBas_Cont,1)
*
      nSym=nIrrep
      iOpt=0
      Call ClsOne(iRC,iOpt)
      iOpt = 0
      iRC = -1
      Lu_One=2
      Call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
      If (iRC.ne.0) Go To 9999
*
      Call OneBas('PRIM')
      Call Get_iArray('nBas_Prim',nBas_Prim,nIrrep)
*
      If(iPrint.ge.10) then
         write(6,'(a,8i5)') ' Symmetries          ', nSym
         write(6,'(a,8i5)') ' Primitive basis fcns',
     &                          (nBas_Prim(i),i=0,nSym-1)
      Endif
*                                                                      *
************************************************************************
*                                                                      *
*     Read P_Matrix from ONEREL
*
      nComp=3
      nLength_Tot=0
      Do iComp = 1, nComp
         iOpt=1
         iRC=-1
         Call iRdOne(iRC,iOpt,'P_matrix',iComp,Length,iSmLbl)
         If (iRC.ne.0) Then
            Call WarningMessage(2,'Error reading length of P-Matrix')
            Write (6,*) 'iComp=',iComp
            Call Abend()
         End If
         iSml(iComp)=iSmLbl
         lOper(iComp)=1
         ip(iComp)=1 + nLength_Tot
         nLength_Tot = nLength_Tot + Length(1) + 4
      End Do
*
      Call mma_allocate(P_Matrix,nLength_Tot,label='P_Matrix')
      Call FZero(P_Matrix,nLength_Tot)
*
      Do iComp = 1, nComp
         iOpt=0
         iRC=-1
         iSmLbl=iSml(iComp)
         ip(iComp) = ip(iComp)
         Call RdOne(iRC,iOpt,'P_matrix',iComp,P_Matrix(ip(iComp)),
     &              iSmLbl)
         If (iRC.ne.0) Then
            Call WarningMessage(2,'Error reading P-Matrix')
            Write (6,*) 'iComp=',iComp
            Call Abend()
         End If
      End Do
      If (iPrint.ge.10) Call PrMtrx('P_matrix',lOper,nComp,ip,P_Matrix)
*                                                                      *
************************************************************************
*                                                                      *
*     Read multipole integrals from ONEREL
*
      nip = 0
      nInt_Tot = 0
      Do iMltPl = 0, MxMltPl
         Write (Label,'(a,i2)') 'MLTPL ',iMltPl
         nComp = (iMltPl+1)*(iMltPl+2)/2
         Do iComp = 1, nComp
            iRC = -1
            iOpt = 1
            nInt=0
            Call iRdOne(iRC,iOpt,Label,iComp,nInt,iSmLbl)
            If (iRC.ne.0) Then
               If (iComp.ne.1) Then
                  Call WarningMessage(2,' Error reading length!')
                  Write (6,*) ' Label=', Label,' Comp=',iComp
                  Call Abend()
               Else
                  Go To 100
               End If
            End If
            nip = nip + 1
            iSm(nip)=iSmLbl
            ipMP(nip)= 1 + nInt_Tot
            nInt_Tot = nInt_Tot + nInt(1)+4
         End Do
      End Do
 100  Continue
*
      Call mma_allocate(MP_Matrix,nInt_Tot,label='MP_Matrix')
*
      iip=0
      Do iMltPl = 0, MxMltPl
         Write (Label,'(a,i2)') 'MLTPL ',iMltPl
         nComp = (iMltPl+1)*(iMltPl+2)/2
         Do iComp = 1, nComp
            iip = iip + 1
            if (iip.gt.nip) Go To 200
            iRC=-1
            iOpt=0
            ipMP(iip) = ipMP(iip)
            iSmLbl=iSm(iip)
            Call RdOne(iRC,iOpt,Label,iComp,MP_Matrix(ipMP(iip)),iSmLbl)
            If (iRC.ne.0) Then
               Call WarningMessage(2,' Error reading integrals!')
               Write (6,*) ' Label=', Label,' Comp=',iComp
               Call Abend()
            End If
         End Do
      End Do
 200  Continue
*                                                                      *
************************************************************************
*                                                                      *
*     Close ONEREL and re-open ONEINT
*
      iOpt = 0
      iRC = -1
      Call ClsOne(iRC,iOpt)
      If (iRC.ne.0) Go To 9999
      iOpt = 0
      iRC = -1
      Lu_One=2
      Call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
      If (iRC.ne.0) Go To 9999
*
      Call OneBas('PRIM')
*                                                                      *
************************************************************************
*                                                                      *
*     Put the transformation matrix on RUNFILE
*
      idbg=0
      Call tr_prm_cnt(idbg,nBas_Cont,nBas_Prim)
*                                                                      *
************************************************************************
*                                                                      *
*     Process the P-matrix
*
      nComp=3
      iOpt=0
      Do iComp = 1, nComp
         iRC=-1
         Call WrOne(iRC,iOpt,'P_matrix',iComp,P_Matrix(ip(iComp)),
     &              iSml(iComp))
         If (iRC.ne.0) Then
            Call WarningMessage(2,'Error reading P-Matrix')
            Write (6,*) 'iComp=',iComp
            Call Abend()
         End If
      End Do
      Call mma_deallocate(P_Matrix)
*                                                                      *
************************************************************************
*                                                                      *
*     Process multipole integrals
*
      iip = 0
      Do iMltPl = 0, MxMltPl
         Write (Label,'(a,i2)') 'PLTPL ',iMltpl
         nComp = (iMltpl+1)*(iMltpl+2)/2
         Do iComp = 1, nComp
            iip = iip + 1
            If (iip.gt.nip) Go To 300
            iRC = -1
            iOpt = 0
            Call WrOne(iRC,iOpt,Label,iComp,MP_Matrix(ipMP(iip)),
     &                 iSm(iip))
            If (iRC.ne.0) Then
               Call WarningMessage(2,' Error writing integrals!')
               Write (6,*) ' Label=', Label,' Comp=',iComp
               Call Abend()
            End If
         End Do
      End Do
 300  Continue
      Call mma_deallocate(MP_Matrix)
*                                                                      *
************************************************************************
*                                                                      *
*     And now change it back!
*
      Call OneBas('CONT')
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('NEMO_Opt1')
      Return
*
 9999 Continue
      Call WarningMessage(2,
     &            ' *** Error in subroutine NEMO_Opt1 ***;'
     &          //'     Abend in subroutine OpnOne or ClsOne')
      Call Abend()
      End
