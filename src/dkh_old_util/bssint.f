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
      Subroutine BSSint
      use Basis_Info
      use Symmetry_Info, only: nIrrep
      Implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "rinfo.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "wldata.fh"
#include "oneswi.fh"
#include "RelLight.fh"
      integer ipaddr(3)
      Character*8 Label
      Logical IfTest
      Data IfTest/.False./
*
      iRout=77
      iPrint=nPrint(iRout)
      Call QEnter('BSSInt')
#ifdef _DEBUGPRINT_
      IfTest=.True.
#endif
*
*     Save basis set info from contracted run
*
      if(iprint.ge.10) write(6,*) ' In dkint', ncnttp
      kCof=0
      kAng=0
      kExp=0
      kC=0
*
*     Normalize coefficients
*
      do iCnttp=1,nCnttp
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
                    Enddo
                 Enddo
                 rNorm=One/sqrt(Sum)
                 if(iprint.ge.10) write(6,*) ' rNorm', kAng,rNorm
                   Do iExp=1,nPrimr(kAng)
                      rCof(kCof+iExp)=rCof(kCof+iExp)*rNorm
                     if(iprint.ge.10) then
                         write(6,'(a24,f20.6)')
     &                   ' normalized coefficients',
     &                   rCof(kCof+iExp)
                     endif
                   Enddo
                 kCof=kCof+nPrimr(kAng)
              Enddo
              kExp=kExp+nPrimr(kAng)
           Enddo
        Enddo
      Enddo
*
      i=0
      If(iPrint.ge.10) then
      Do L=1,nrSym
         write(6,*) ' Irreducible representation', L
         Do ibas=1,nrBas(L)
         i=i+1
         write(6,'(20i4)') i, icent(i),lnang(i),lmag(i)
         enddo
      enddo
      endif
*
*     Close ONEINT and re-open ONEREL
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
      Call Get_iArray('nBas_Prim',nbas,nSym)
*
      If(iPrint.ge.10) then
         write(6,'(a11,10i5)') ' Symmetries', nSym
         write(6,'(a11,10i5)') ' Primitive basis fcns',
     &                          (nBas(i),i=0,nSym-1)
      Endif
*
*
*     Allocate memory for relativistic part
*
*
      LenIntf=0
      LenIntf1=0
      Do L=0,nSym-1
      n=nBas(L)
      LenIntf=LenIntf+n*n+4
      LenIntf1=LenIntf1+n+4
      Enddo
*
*
      CALL GetMem('Eigf    ','ALLO','REAL',iEigf,LenIntf)
      CALL GetMem('Sinvf   ','ALLO','REAL',iSinvf,LenIntf)
      CALL GetMem('Revtf   ','ALLO','REAL',iRevtf,LenIntf)
      CALL GetMem('Aaf     ','ALLO','REAL',iAaf,LenIntf1)
      CALL GetMem('Rrf     ','ALLO','REAL',iRrf,LenIntf1)
      call dcopy_(LenIntf,[Zero],0,Work(iEigf),1)
      call dcopy_(LenIntf,[Zero],0,Work(iSinvf),1)
      call dcopy_(LenIntf,[Zero],0,Work(iRevtf),1)
      call dcopy_(LenIntf1,[Zero],0,Work(iAaf),1)
      call dcopy_(LenIntf1,[Zero],0,Work(iRrf),1)
*
      VELIT=CLightAU
      iSizep=0
      Do L=0,nSym-1
         iSizep=iSizep + nBas(L)*(nBas(L)+1)/2
      Enddo
        If(iPrint.ge.10) write(6,*) ' iSizep', iSizep
*
      CALL GetMem('Kin     ','ALLO','REAL',iK,iSizep+4)
      CALL GetMem('SS      ','ALLO','REAL',iSS,iSizep+4)
      CALL GetMem('V       ','ALLO','REAL',iV,iSizep+4)
      CALL GetMem('pVp     ','ALLO','REAL',ipVp,iSizep+4)
*
      If (iprint.ge.20) write(6,*) '  indices', iss, ik, iv, ipvp
      Label='Mltpl  0'
      iComp=1
      iOpt=0
      iRC = -1
      Call RdOne(iRC,iOpt,Label,1,Work(iSS),lOper)
      If (iRC.ne.0) then
         Write (6,*) 'BSSInt: Error reading from ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend
      End If
c         write(6,'(8f9.4)') (Work(iSS+k),k=0,iSizep-1)
      nComp=1
      ipaddr(1)=iSS
      If (iPrint.ge.20) Call PrMtrx(Label,[lOper],nComp,ipaddr,Work)
      Label='Attract '
      iRC = -1
      Call RdOne(iRC,iOpt,Label,1,Work(iV),lOper)
      If (iRC.ne.0) then
         Write (6,*) 'BSSInt: Error reading from ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend
      End If
      Label='Kinetic '
      iRC = -1
      Call RdOne(iRC,iOpt,Label,1,Work(iK),lOper)
      If (iRC.ne.0) then
         Write (6,*) 'BSSInt: Error reading from ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend
      End If
      Label='pVp     '
      iRC = -1
      Call RdOne(iRC,iOpt,Label,1,Work(ipVp),lOper)
      If (iRC.ne.0) then
         Write (6,*) 'BSSInt: Error reading from ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend
      End If
*
      nComp=3
      Call GetMem ('lOper1  ','ALLO','INTE',ip2,nComp)
      Call GetMem ('ipl     ','ALLO','INTE',ipV,nComp)
         ixyz=1
         iWork(ip2  ) = 2**IrrFnc(ixyz)
         ixyz=2
         iWork(ip2+1) = 2**IrrFnc(ixyz)
         ixyz=4
         iWork(ip2+2) = 2**IrrFnc(ixyz)
*
      Call GetMem ('iiVpf   ','ALLO','INTE',iVpf,nComp)
      Call GetMem ('iipVf   ','ALLO','INTE',ipVf,nComp)
*
*
C     Read pV matrix , elements <iSyma|pV|iSymb>,
C     iSymb.le.iSyma !
*
      DO iComp=1,nComp
*
      Label='pV      '
      iRC = -1
      iSmlbl=iWork(ip2+iComp-1)
      LenInt=n2Tri(iSmlbl)
      Call GetMem ('pV      ','ALLO','REAL',iWork(ipV+iComp-1),LenInt+4)
      ip4=iWork(ipV+iComp-1)
      Call RdOne(iRC,iOpt,Label,iComp,Work(ip4),iSmlbl)
*
*
      ENDDO
*
*     Call PrMtrx('pV      ',iWork(ip2),nComp,iWork(ipV),Work)
*
*
*
C     Read Vp matrix , elements <iSyma|Vp|iSymb>,
C     iSymb.le.iSyma !
*
*
C     Read Vp matrix
*
      Label='Vp      '
*
      Call GetMem ('ipk     ','ALLO','INTE',iVp,nComp)
*
      DO iComp=1,nComp
*
      iRC = -1
      iSmlbl=iWork(ip2+iComp-1)
      LenInt=n2Tri(iSmlbl)
      Call GetMem ('Vp      ','ALLO','REAL',iWork(iVp+iComp-1),LenInt+4)
      ip5=iWork(iVp+iComp-1)
      Call RdOne(iRC,iOpt,Label,iComp,Work(ip5),iSmlbl)
      ENDDO


*     Call PrMtrx(Label,iWork(ip2),nComp,iWork(iVp),Work)
*
*
C     Build the whole matrix <iSyma|pV|iSymb> lower triangel
C     and put into a Work(ip6+ip1+k) and the upper triange
C     and put into a Work(ip6+LenInt+iip1+k).
C     For 'Vp' matrix the same , p7 instead of p6
*
*
      Do iComp =1,nComp
*
        iSmlbl=iWork(ip2+iComp-1)
        LenInt=n2Tri(iSmlbl)
*
        Call GetMem ('pVf     ','ALLO','REAL',
     &              iWork(ipVf+iComp-1),2*LenInt+4)
        Call GetMem ('Vpf     ','ALLO','REAL',
     &              iWork(iVpf+iComp-1),2*LenInt+4)
*
        ip4=iWork(ipV+iComp-1)
        ip5=iWork(iVp+iComp-1)
        ip6=iWork(ipVf+iComp-1)
        ip7=iWork(iVpf+iComp-1)
        ip1=0
        iip1=0
*
      Do iSyma=0,nSym-1
        na=nBas(iSyma)
        If (na.eq.0) go to 81
*
      Do iSymb=0,iSyma
        nb=nBas(iSymb)
        If (nb.eq.0) go to 82
        If (iAnd(iSmLbl,2**iEor(iSyma,iSymb)).eq.0)
     &  Go to 82
      If (iSyma.eq.iSymb) Then
        iSize=na*(na+1)/2
C
C     Elements <iSyma|pV|iSyma> and <iSyma|Vp|iSyma>
C
      Do k=0,iSize-1
        Work(ip6+ip1+k)=Work(ip4+ip1+k)
        Work(ip7+ip1+k)=Work(ip5+ip1+k)
      Enddo
        ip1=ip1+na*(na+1)/2
      Else
        iSize=na*nb
C
C     Elements <iSyma|pV|iSymb> and <iSymb|pV|iSyma> = -<iSyma|Vp|iSymb>
C     AND
C     Elements <iSyma|Vp|iSymb> and <iSymb|Vp|iSyma> = -<iSyma|pV|iSymb>
C
      Do k=0,iSize-1
        Work(ip6+ip1+k)=Work(ip4+ip1+k)
        Work(ip7+ip1+k)=Work(ip5+ip1+k)
        Work(ip6+LenInt+4+iip1+k)=-Work(ip5+ip1+k)
        Work(ip7+LenInt+4+iip1+k)=-Work(ip4+ip1+k)
      Enddo
        ip1=ip1+na*nb
        iip1=iip1+na*nb
      End If
*
 82   Continue
      Enddo
 81   Continue
      Enddo
*
*
C
C
*     WRITE(*,*) 'Created  whole matrix (NAxNB) PV '
C
*     WRITE(*,*) 'pV matrix <iSyma|pV|iSymb> '
*     WRITE(*,*)
*     WRITE(*,*) (Work(iWork(ipVf+iComp-1)+i),i=0,LenInt-1)
*     WRITE(*,*)
*     WRITE(*,*)
*     WRITE(*,*) 'pV matrix <iSymb|pV|iSyma> if iSyma.ne.iSymb
*    &            otherwise zero matrix '
*     WRITE(*,*)
*     WRITE(*,*) (Work(iWork(ipVf+iComp-1)+LenInt+4+i),i=0,LenInt-1)
*
*
C     Create the whole matrix (NAxNB) VP
C
*     WRITE(*,*) 'Create the whole matrix (NAxNB) VP '
*
*     WRITE(*,*)
*     WRITE(*,*) 'Vp matrix <iSyma|Vp|iSymb> '
*     WRITE(*,*)
*     WRITE(*,*) (Work(iWork(iVpf+iComp-1)+i),i=0,LenInt-1)
*     WRITE(*,*)
*     WRITE(*,*)
*     WRITE(*,*)
*     WRITE(*,*) 'Vp matrix <iSymb|Vp|iSyma> if iSyma.ne.iSymb
*    &            otherwise zero matrix '
*     WRITE(*,*)
*     WRITE(*,*)
*     WRITE(*,*) (Work(iWork(iVpf+iComp-1)+LenInt+4+i),i=0,LenInt-1)
*
*
C     End of iComp loop
*
      Enddo
*
*
*
*     Main loop  1
      epsilon=1.d-10
      k=0
      L=0
      k1=0
      k2=0
      Do L = 0, nSym-1
*
*
         n=nBas(L)
         iSize=n*(n+1)/2
CAJS protection against zero dimension representation
         if (iSize.eq.0) goto 9
CAJS
*   put zeroes
         Len=4*(iSize+4)+5*(n*n+4)+5*(n+4)
         CALL GetMem('Scratch ','ALLO','REAL',iCr,Len+4)
         call dcopy_(Len,[Zero],0,Work(iCr),1)
         CALL GetMem('Scratch ','FREE','REAL',iCr,Len+4)
*
*        Allocate
*
         CALL GetMem('Bu      ','ALLO','REAL',iBu,isize+4)
         CALL GetMem('P       ','ALLO','REAL',iP,isize+4)
         CALL GetMem('G       ','ALLO','REAL',iG,isize+4)
         CALL GetMem('Ev2     ','ALLO','REAL',iEv2,isize+4)
         CALL GetMem('Eig     ','ALLO','REAL',iEig,n*n+4)
         CALL GetMem('Sinv    ','ALLO','REAL',iSinv,n*n+4)
         CALL GetMem('Revt    ','ALLO','REAL',iRevt,n*n+4)
         CALL GetMem('Aux     ','ALLO','REAL',iAux,n*n+4)
         CALL GetMem('Ove     ','ALLO','REAL',iOve,n*n+4)
         CALL GetMem('Ew      ','ALLO','REAL',iEw,n+4)
         CALL GetMem('E       ','ALLO','REAL',iE,n+4)
         CALL GetMem('Aa      ','ALLO','REAL',iAa,n+4)
         CALL GetMem('Rr      ','ALLO','REAL',iRr,n+4)
         CALL GetMem('Tt      ','ALLO','REAL',iTt,n+4)
*
*        Debug output on unit idbg
         If(IfTest)Then
           idbg=41
         Else
           idbg=-1
         Endif
*
*
*   call to package relsewb
         CALL SCFCLI2(idbg,epsilon,Work(iSS+k),Work(iK+k),Work(iV+k),
     *   Work(ipVp+k),n,iSize,VELIT,Work(iBu),Work(iP),
     *   Work(iG),Work(iEv2),Work(iEig),Work(iSinv),Work(iRevt),
     *        Work(iAux),Work(iOve),Work(iEw),Work(iE),Work(iAa),
     *        Work(iRr),Work(iTt))
*
*
*
*
*
         call dcopy_(n*n+4,Work(iEig),1,Work(iEigf+k1),1)
         call dcopy_(n*n+4,Work(iSinv),1,Work(iSinvf+k1),1)
         call dcopy_(n*n+4,Work(iRevt),1,Work(iRevtf+k1),1)
         call dcopy_(n+4,Work(iAa),1,Work(iAaf+k2),1)
         call dcopy_(n+4,Work(iRr),1,Work(iRrf+k2),1)
*
*
*
*
*
         CALL GetMem('Bu      ','FREE','REAL',iBu,isize+4)
         CALL GetMem('P       ','FREE','REAL',iP,isize+4)
         CALL GetMem('G       ','FREE','REAL',iG,isize+4)
         CALL GetMem('Ev2     ','FREE','REAL',iEv2,isize+4)
         CALL GetMem('Eig     ','FREE','REAL',iEig,n*n+4)
         CALL GetMem('Sinv    ','FREE','REAL',iSinv,n*n+4)
         CALL GetMem('Revt    ','FREE','REAL',iRevt,n*n+4)
         CALL GetMem('Aux     ','FREE','REAL',iAux,n*n+4)
         CALL GetMem('Ove     ','FREE','REAL',iOve,n*n+4)
         CALL GetMem('Ew      ','FREE','REAL',iEw,n+4)
         CALL GetMem('E       ','FREE','REAL',iE,n+4)
         CALL GetMem('Aa      ','FREE','REAL',iAa,n+4)
         CALL GetMem('Rr      ','FREE','REAL',iRr,n+4)
         CALL GetMem('Tt      ','FREE','REAL',iTt,n+4)
 9       k=k+isize
         k1=k1+n*n+4
         k2=k2+n+4
      Enddo
*
*
*
*     Main loop  2
      nComp=3
*
      Do iComp=1,nComp
*
      ip6=iWork(ipVf+iComp-1)
      ip7=iWork(iVpf+iComp-1)
      ip1=0
      iip1=0
      epsilon=1.d-10
      kh=0
      k1a=0
      k2a=0
      iSmlbl=iWork(ip2+iComp-1)
      LenInt=n2Tri(iSmlbl)
*
      Do iSyma = 0, nSym-1
*
         na=nBas(iSyma)
         iSizea=na*(na+1)/2
CAJS protection against zero dimension representation
         if (iSizea.le.0) goto 19
CAJS
*
       k1b=0
       k2b=0
*
       Do iSymb=0,nSym-1
*
       nb=nBas(iSymb)
       iSizeb=nb*(nb+1)/2
       If (iSizeb.le.0) goto 29
       iSizeab=na*nb
       If (iAnd(iSmlbl,2**iEor(iSyma,iSymb)).eq.0)
     *    Go to 29
*
*
       Call GetMem('ipVa    ','ALLO','REAL',ifa,iSizea)
       Call GetMem('iVpa    ','ALLO','REAL',if2a,iSizea)
*
       If (iSyma.eq.iSymb) Then
*
       call dcopy_(iSizea,[Zero],0,Work(ifa),1)
       call dcopy_(iSizea,[Zero],0,Work(if2a),1)
       Do k=0,iSizea-1
       Work(ifa+k)=Work(ip6+ip1+k)
       Work(if2a+k)=Work(ip7+ip1+k)
       Enddo
       ip1=ip1+iSizea
       End If
*
*
       Call GetMem('ifpV    ','ALLO','REAL',if,iSizeab)
       Call GetMem('ifVp    ','ALLO','REAL',if2,iSizeab)
       Call GetMem('ScpV    ','ALLO','REAL',iScpV,iSizeab)
       Call GetMem('ScVp    ','ALLO','REAL',iScVp,iSizeab)
       call dcopy_(iSizeab,[Zero],0,Work(if),1)
       call dcopy_(iSizeab,[Zero],0,Work(if2),1)
       call dcopy_(iSizeab,[Zero],0,Work(iScpV),1)
       call dcopy_(iSizeab,[Zero],0,Work(iScVp),1)
*
       If (iSyma.gt.iSymb) Then

       Do k=0,iSizeab-1
       Work(if+k) = Work(ip6+ip1+k)
       Work(if2+k) = Work(ip7+ip1+k)
       Enddo
       ip1=ip1+na*nb
*
*
       End If
*
       If (iSyma.lt.iSymb) Then
*
       Do k=0,iSizeab-1
       Work(if+k)=Work(ip6+iip1+LenInt+4+k)
       Work(if2+k)=Work(ip7+iip1+LenInt+4+k)
       Enddo
       CALL DCOPY_(iSizeab,Work(if),1,Work(iScpV),1)
       CALL DCOPY_(iSizeab,Work(if2),1,Work(iScVp),1)
*
       iip1=iip1+na*nb
       End If
*
*
*
*
*
*
*
*
*        Allocate
*
         CALL GetMem('Bu2     ','ALLO','REAL',iBu2,na*nb+4)
         CALL GetMem('Bu4     ','ALLO','REAL',iBu4,na*nb+4)
         CALL GetMem('Bu6     ','ALLO','REAL',iBu6,na*na+4)
         CALL GetMem('Bu      ','ALLO','REAL',iBu,iSizea+4)
         CALL GetMem('P       ','ALLO','REAL',iP,iSizea+4)
         CALL GetMem('G2      ','ALLO','REAL',iG2,na*nb+4)
         CALL GetMem('Ev4     ','ALLO','REAL',iEv4,isizea+4)
         CALL GetMem('Eig4    ','ALLO','REAL',iEig4,na*na+4)
*        CALL GetMem('Sinv    ','ALLO','REAL',iSinv,n*n+4)
*        CALL GetMem('Revt    ','ALLO','REAL',iRevt,n*n+4)
         CALL GetMem('Aux2    ','ALLO','REAL',iAux2,na*nb+4)
         CALL GetMem('Cmm1    ','ALLO','REAL',iCmm1,na*nb+4)
         CALL GetMem('Cmm2    ','ALLO','REAL',iCmm2,na*nb+4)
         CALL GetMem('Ew4     ','ALLO','REAL',iEw4,na+4)
*        CALL GetMem('E       ','ALLO','REAL',iE,n+4)
*        CALL GetMem('Tt      ','ALLO','REAL',iTt,n+4)

*
*
*
      call dcopy_(na*nb+4,[Zero],0,Work(iBu2),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iBu4),1)
      call dcopy_(na*na+4,[Zero],0,Work(iBu6),1)
      call dcopy_(na*na+4,[Zero],0,Work(iEig4),1)
      call dcopy_(na+4,[Zero],0,Work(iEw4),1)
      call dcopy_(iSizea+4,[Zero],0,Work(iEv4),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iG2),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iAux2),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iCmm1),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iCmm2),1)
*
*
*
*
*
*
*
*     Calculate the matrix elements of the
*     following commutator :
*     Revta*AAA*<a|[bp,V]|b>store in CMM1(iSyma,iSymb) matrix,
*
*
*
         CALL VPBMBPV(idbg,epsilon,Work(iEigf+k1a),
     *   Work(iEigf+k1b),
     *   Work(iRevtf+k1a),Work(iSinvf+k1a),
     *   Work(iSinvf+k1b),na,nb,iSizea,iSizeb,VELIT,
     *   Work(iAaf+k2a),Work(iAaf+k2b),
     *   Work(iRrf+k2a),Work(iRrf+k2b),
     *   Work(if ),Work(if2),iSyma,iSymb,Work(iBu2),
     *   Work(iG2),Work(iAux2),Work(iCmm1),Work(iBu4),
     *   Work(iCmm2),Work(ifa),Work(if2a),Work(iScpV),
     *   Work(iScVp))
*
*
*
      call dcopy_(na*nb+4,[Zero],0,Work(iBu2),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iBu4),1)
      call dcopy_(na*na+4,[Zero],0,Work(iBu6),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iG2),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iAux2),1)
      call dcopy_(na*na+4,[Zero],0,Work(iEig4),1)
*
*
*
*        call to package relsewc, BSS up to the fourth order in alpha
C
         CALL SCFCLI4(idbg,epsilon,
     *   Work(iSS+kh),Work(iK+kh),
     *   Work(iRevtf+k1a),Work(iSinvf+k1a),
     *   na,nb,iSizea,iSizeb,VELIT,
     *   Work(iAaf+k2a),Work(iAaf+k2b),
     *   iSyma,iSymb,Work(iCmm1),
     *   Work(iCmm2),Work(iEv4),Work(iBu2),Work(iBu6),
     *   Work(iEig4),Work(iEw4),Work(iP))
*
*
*
*
*
*        Free a space
*
*
*
C
         CALL GetMem('Bu2     ','FREE','REAL',iBu2,na*nb+4)
         CALL GetMem('Bu4     ','FREE','REAL',iBu4,na*nb+4)
         CALL GetMem('Bu6     ','FREE','REAL',iBu6,na*na+4)
         CALL GetMem('Bu      ','FREE','REAL',iBu,iSizea+4)
         CALL GetMem('P       ','FREE','REAL',iP,iSizea+4)
         CALL GetMem('G2      ','FREE','REAL',iG2,na*nb+4)
         CALL GetMem('Ev4     ','FREE','REAL',iEv4,isizea+4)
         CALL GetMem('Eig4    ','FREE','REAL',iEig4,na*na+4)
*        CALL GetMem('Sinv    ','FREE','REAL',iSinv,n*n+4)
*        CALL GetMem('Revt    ','FREE','REAL',iRevt,n*n+4)
         CALL GetMem('Aux2    ','FREE','REAL',iAux2,na*nb+4)
         CALL GetMem('Cmm1    ','FREE','REAL',iCmm1,na*nb+4)
         CALL GetMem('Cmm2    ','FREE','REAL',iCmm2,na*nb+4)
         CALL GetMem('Ew4     ','FREE','REAL',iEw4,na+4)
*        CALL GetMem('E       ','FREE','REAL',iE,n+4)
*        CALL GetMem('Tt      ','FREE','REAL',iTt,n+4)
*
*
*
         Call GetMem('ifpV    ','FREE','REAL',if,iSizeab)
         Call GetMem('ifVp    ','FREE','REAL',if2,iSizeab)
         Call GetMem('ScpV    ','FREE','REAL',iScpV,iSizeab)
         Call GetMem('ScVp    ','FREE','REAL',iScVp,iSizeab)
*
*
*
         Call GetMem('ipVa    ','FREE','REAL',ifa,iSizea)
         Call GetMem('iVpa    ','FREE','REAL',if2a,iSizea)
*
 29   Continue
      k1b=k1b+nb*nb+4
      k2b=k2b+nb+4
      Enddo
 19   kh=kh+iSizea
      k1a=k1a+na*na+4
      k2a=k2a+na+4
      Enddo
      Enddo
*
*
      Do iComp=1,nComp
      iSmlbl=iWork(ip2+iComp-1)
      LenInt=n2Tri(iSmlbl)
      Call GetMem ('pVf     ','FREE','REAL',
     &              iWork(ipVf+iComp-1),2*LenInt+4)
      Call GetMem ('Vpf     ','FREE','REAL',
     &              iWork(iVpf+iComp-1),2*LenInt+4)
      Enddo
*
      do iComp=1,nComp
      iSmlbl=iWork(ip2+iComp-1)
      LenInt=n2Tri(iSmlbl)
      Call GetMem ('Vp      ','FREE','REAL',iWork(iVp+iComp-1),LenInt+4)
      enddo
      Call GetMem ('ipk     ','FREE','INTE',iVp,nComp)
*
      do iComp=1,nComp
      iSmlbl=iWork(ip2+iComp-1)
      LenInt=n2Tri(iSmlbl)
      Call GetMem ('pV      ','FREE','REAL',iWork(ipV+iComp-1),LenInt+4)
      enddo
*
*
      Call GetMem ('iiVpf   ','FREE','INTE',iVpf,nComp)
      Call GetMem ('iipVf   ','FREE','INTE',ipVf,nComp)
*
*
*
      Call GetMem ('ipl     ','FREE','INTE',ipV,nComp)
      Call GetMem ('lOper1  ','FREE','INTE',ip2,nComp)

      CALL GetMem('pVp     ','FREE','REAL',ipVp,iSizep+4)
*
*     open arrays in contracted basis
*
      iSizec=0
      do L=1,nSym
         iSizec=iSizec+nrBas(L)*(nrBas(L)+1)/2
      enddo
      CALL GetMem('H       ','ALLO','REAL',iH,iSizec+4)
      Call FZero(Work(iH),iSizec+4)
      CALL GetMem('H_nr    ','ALLO','REAL',iH_nr,iSizec+4)
      Call FZero(Work(iH_nr),iSizec+4)
      CALL GetMem('H_temp  ','ALLO','REAL',iH_temp,iSizec+4)
      Call FZero(Work(iH_temp),iSizec+4)
*
*     compute stripped non-relativistic H
*
      Label='Kinetic '
      Call RdOne(iRC,iOpt,Label,1,Work(iSS),lOper)
      If (iRC.ne.0) then
         Write (6,*) 'BSSInt: Error reading from ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend
      End If
C
C
      Label='Attract '
      Call RdOne(iRC,iOpt,Label,1,Work(iV),lOper)
      If (iRC.ne.0) then
         Write (6,*) 'BSSInt: Error reading from ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend
      End If
C
C
      Call DaXpY_(iSizep+4,One,Work(iSS),1,Work(iV),1)
C
C
C
C
C
      call dcopy_(4,[Zero],0,Work(iH_temp+iSizec),1)
      Call repmat(idbg,Work(iV),Work(iH_temp))
*
C
C
C
      CALL GetMem('V       ','FREE','REAL',iV,iSizep+4)
      CALL GetMem('SS      ','FREE','REAL',iSS,iSizep+4)
*
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
*
*     Transform to contracted basis
*
*     The Hamiltonian is now in Kin.
*


C
      Call repmat(idbg,Work(iK),Work(iH))
*
      CALL GetMem('Kin     ','FREE','REAL',iK,iSizep+4)
*
      iOpt = 0
      iRC = -1
      Label='OneHam 0'
      Call RdOne(iRC,iOpt,Label,1,Work(iH_nr),lOper)
      If (iRC.ne.0) then
         Write (6,*) 'BSSInt: Error reading from ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend
      End If
      iOpt = 0
      iRC = -1
*
*     final Hamiltonian computed as H(nrel) + ( Hrel(s) - Hnrel(s))
*     where (s) is stripped and with full charge
*
      Call DaXpY_(iSizec+4,-One,Work(iH_temp),1,Work(iH),1)
      Call DaXpY_(iSizec+4,One,Work(iH_nr),1,Work(iH),1)
*
      Call Get_iArray('nBas',nBas,nSym)
      If(iPrint.ge.10) then
         write(6,'(a11,10i5)') ' Symmetries', nSym
         write(6,'(a11,10i5)') ' Contracted',(nBas(i),i=0,nSym-1)
      Endif
      Label='OneHam 0'
      lOper=1
      nComp=1
      ipaddr(1)=iH
      If (iPrint.ge.20) Call PrMtrx(Label,[lOper],nComp,ipaddr,Work)
*
*     Replace 1-el Hamiltonian on ONEINT
*
      iRC = -1
      Call WrOne(iRC,iOpt,Label,1,Work(iH),lOper)
      Label='OneHam  '
      Call WrOne(iRC,iOpt,Label,1,Work(iH),lOper)
      If (iRC.ne.0) then
         Write (6,*) 'BSSInt: Error writing to ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend
      End If
*
      CALL GetMem('OneHam  ','FREE','REAL',iH,iSizec+4)
      CALL GetMem('H_nr    ','FREE','REAL',iH_nr,iSizec+4)
      CALL GetMem('H_temp  ','FREE','REAL',iH_temp,iSizec+4)
*
*
      CALL GetMem('Eigf    ','FREE','REAL',iEigf,LenIntf)
      CALL GetMem('Sinvf   ','FREE','REAL',iSinvf,LenIntf)
      CALL GetMem('Revtf   ','FREE','REAL',iRevtf,LenIntf)
      CALL GetMem('Aaf     ','FREE','REAL',iAaf,LenIntf1)
      CALL GetMem('Rrf     ','FREE','REAL',iRrf,LenIntf1)
      Call QExit('BSSInt')
      Return
*
 9999 Continue
      Call qTrace
      Write (6,*) ' *** Error in subroutine BSSint ***'
      Write (6,*) '     Abend in subroutine OpnOne'
      Call Abend
      End
