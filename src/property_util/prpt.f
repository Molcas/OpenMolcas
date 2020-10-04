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
* Copyright (C) 1991, Roland Lindh                                     *
*               2018, Sijia S. Dong                                    *
************************************************************************
      Subroutine Prpt()
************************************************************************
c
c     Purpose: To set up all calling arguments for the subroutine
c              Prpt_ . For RASSCF to work MO coefficents and Occupation
c              numbers must be available on TMPORB.
c
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
      Integer nBas(8)
      Character*8 Method
      Logical var, Short, ifallorb
      Character*81 note, lbl*2, PrpLst*4
      Dimension Dummy(1),iDummy(1)
*
      Call GetEnvf("MOLCAS_PROPERTIES",PrpLst)
      Call UpCase(PrpLst)
      If (PrPlst(1:3).eq.'LON') Then
         Short=.False.
*         ifallorb=.True.
      Else
         Short=.True.
         ifallorb=.False.
      End If
*
*     This variable is used so we know if the density we search for is labeled
*     variational or not.
      var =.false.
*
      Call Get_cArray('Relax Method',Method,8)
*
      Call Get_iScalar('nSym', nIrrep)
*
      Call Get_iArray('nBas',nBas,nIrrep)
*
      nDim = 0
      n2Dim = 0
      nTriDim = 0
      n2Tot= 0
      Do iSym = 1,nIrrep
         nDim = nDim + nBas(iSym)
         nTriDim = nTriDim + nBas(iSym)*(nBas(iSym)+1)/2
         n2Tot=n2Tot+nBas(iSym)**2
      End Do
*
      ipOcc    = ip_Dummy ! dummy initialization
      ipOcc_ab = ip_Dummy ! dummy initialization
      ipVec    = ip_Dummy ! dummy initialization
      ipVec_ab = ip_Dummy ! dummy initialization
*
      If((Method.eq.'RHF-SCF ').or.(Method.eq.'IVO-SCF ').or.
     &   (Method.eq.'KS-DFT  ').or.Method.eq.'UHF-SCF ') Then
         Call Get_iScalar('SCF mode',iUHF)
      Else
         iUHF=0
      End If
*
      If (iUHF.eq.1.or.Method.eq.'RASSCFSA') Then
         Call GetMem('Occ','Allo','Real',ipOcc,2*nDim)
         ipOcc_ab=ipOcc + nDim
      Else
         Call GetMem('Occ','Allo','Real',ipOcc,nDim)
      End If
      If (Short) Then
         ipVec = ip_Dummy
         lbl='O '
         n2Tot=1
      Else
         If (iUHF.ne.1.or.Method.eq.'RASSCFSA') Then
            Call GetMem('Vec','Allo','Real',ipVec,n2Tot)
         Else
            Call GetMem('Vec','Allo','Real',ipVec,2*n2Tot)
            ipVec_ab=ipVec + n2Tot
         End If
         lbl='CO'
      End If
*
      Lu = 10
      Lu=IsFreeUnit(Lu)
      If((Method.eq.'RHF-SCF ').or.(Method.eq.'IVO-SCF ').or.
     &   (Method.eq.'KS-DFT  ').or.Method.eq.'UHF-SCF ') Then
         If(iUHF.ne.1) then
            Call RdVec('SCFORB',Lu,Lbl,nIrrep,nBas,
     &                 nBas,Work(ipVec),Work(ipOcc),Dummy,
     &                 iDummy,'',0,iError)
         Else
            Call RdVec_('UHFORB',Lu,Lbl,iUHF,nIrrep,nBas,
     &                  nBas,Work(ipVec),Work(ipVec_ab),Work(ipOcc),
     &                  Work(ipOcc_ab),Dummy,Dummy,iDummy,'',1,iError,
     &                  iWFtype)
            If (Short) Then
               Do i = 0,nDim-1
                  Work(ipOcc+i) = Work(ipOcc+i) + Work(ipOcc_ab+i)
               End Do
            End If
         End If
      Else If((Method.eq.'RASSCF  ').or.(Method.eq.'CASSCF  ').or.
     &        (Method.eq.'CASDFT  ').or.(Method.eq.'CASSCFSA').or.
     &        (Method.eq.'CASPT2  ').or.(Method.eq.'RASSCFSA')) Then
         If (Method.eq.'RASSCFSA') Then
            Call RdVec_('TMPORB',Lu,Lbl,iUHF,nIrrep,nBas,
     &                  nBas,Work(ipVec),Work(ipVec_ab),Work(ipOcc),
     &                  Work(ipOcc_ab),Dummy,Dummy,iDummy,'',1,iError,
     &                  iWFtype)
            If (Short) Then
               Do i = 0,nDim-1
                  Work(ipOcc+i) = Work(ipOcc+i) + Work(ipOcc_ab+i)
               End Do
            End If
            var = .False.
         Else
            Call RdVec('TMPORB',Lu,Lbl,nIrrep,nBas,
     &                 nBas,Work(ipVec),Work(ipOcc),Dummy,
     &                 iDummy,note,0,iError)
            If(Note(2:4).eq.'var') var = .true.
         End If
      Else If(Method.eq.'MBPT2   ') then
*     MBPT2 has no occupation-numbers at the moment.
         Call FZero(Work(ipOcc),nDim)
         var = .true.
      Else
         Write (6,*) 'Properties not supported for ',Method
      End If
*
      MaxScr = nTriDim + nDim*(nDim + 1)/2 + 10 + 480 + 4*10
      Call GetMem('Scr','Allo','Real',ipScr,MaxScr)
      Call FZero(Work(ipScr),MaxScr)
*
      Call Prpt_(nIrrep,nBas,n2Dim,
     &            nDim,Work(ipOcc),n2Tot,Work(ipVec),MaxScr,
     &            Work(ipScr),var,Short,iUHF,ifallorb)
*
      Call GetMem('Scr','Free','Real',ipScr,MaxScr)
      Call GetMem('Occ','Free','Real',ipOcc,nDim)
      If (.Not.Short) Call GetMem('Vec','Free','Real',ipVec,n2Tot)
      Return
      End
*
      Subroutine Prpt_(nIrrep,nBas,n2Dim,nDim,Occ,n2Tot,Vec,
     &                 MaxScr,Scr,var,Short,iUHF,ifallorb)
************************************************************************
c
c     Purpose: calculation of expectation values of different
c              operators as available on the 'ONEINT' file
c
c     Caution: before calling this subroutine one needs to
c              open the ONEINT file
c
c     Calling parameters:
c
c       nIrRep            number of irreducible representations
c
c       nBas(0:nIrRep)    number of basis functions in each repre-
c                         sentation
c
c       ndim              total number of basis functions
c
c       occ(1:ndim)       occupation number for all eigenvectors
c
c       maxscr            the maximum available size of the
c                         scratch area
c
c       scr(1:maxscr)     a scratch area whose size, nscr, can be
c                         calculated in the following way:
c
c                         nscr =sum(i,i=0,nirrep-1)(nbas(i)*
c                                  (nbas(i)+1)/2)
c                               +3
c                               +3
c                               +nComp
c                               +nComp
c                               +2*nComp
c                               +ndim*(ndim+1)/2+4
c                               +2*ntComp*(ntComp+1)
c
c                         and should not exceed maxscr.
c                         nComp is the number of cartresian compo-
c                         nets for the given operator
c                         ntComp is the number of components of
c                         the l-th moment opartors which are trans-
c                         formed into l-pole moment. Currently,
c                         ntComp=15 (hexadecapole moments)
c
c     ifallorb        logical option for whether the property of
c                     all orbitals are printed (and not weighted by
c                     occupation number)in property calculation
c                     (S.S.Dong, 2018)
c
c
* 1991 R. Lindh, Dept. of Theor. Chem. Univ. of Lund, Sweden.          *
* Modified by S.S.Dong, 2018, Univ. of Minnesota                       *
* - Enable properties to be printed for all orbitals                   *
* (including virtuals) and not weighted by occupation numbers          *
************************************************************************
      Implicit real*8 (a-h,o-z)
*
#include "real.fh"
#include "WrkSpc.fh"
c
      Character*8 label
      Logical short, NxtOpr, var, Reduce_Prt, ifallorb
      External Reduce_Prt
      Integer nBas(0:nirrep-1), mBas(0:7)
      Real*8  occ(1:ndim), scr(1:maxscr), Vec(n2Tot)
      Dimension idum(1)
c
      Call Prpt_Internal(Scr)
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(n2Dim)
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine Prpt_Internal(Scr)
      Use Iso_C_Binding
      Real*8, Target :: Scr(*)
      Character, Pointer :: cScr(:)
#ifdef _DEBUG_
#endif
*                                                                      *
************************************************************************
*                                                                      *
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0

      If (iPL.ge.2) Then
         Write (6,*)
         Call CollapseOutput(1,'   Molecular properties:')
         Write (6,'(3X,A)')    '   ---------------------'
         Write (6,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Thrs=1.0D-6
      Call ICopy(nIrrep,nBas,1,mBas,1)
      If (Short) Then
         mDim = 1
      Else
         mDim = nDim
         If (iUHF.eq.1) Then
            mDim = 2*mDim
            Do iIrrep = 0, nIrrep-1
               mBas(iIrrep)=2*mBas(iIrrep)
            End Do
         End If
      End If
      iadopr=-1   ! dummy initialize
      iadtmt=-1   ! dummy initialize
      iadtmp=-1   ! dummy initialize
      iadDen_ab=1 ! dummy initialize
      iOcc_ab=1   ! dummy initialize
*                                                                      *
************************************************************************
*                                                                      *
      nblock=0
      do 1 i=0,nirrep-1
        nblock=nblock+nbas(i)*(nbas(i)+1)/2
1     continue
      nfblock=ndim*(ndim+1)/2
c
c     estimate the minimum scratch area needed for calculatig
c     the average values of a 3 component operator
c
      nscr=nblock+nfblock
      iscr=nscr+76
c
      if (iscr.gt.maxscr) Go To 999
*                                                                      *
************************************************************************
*                                                                      *
      If (Short) Then
         iadDen=1
c
c        calculate the density matrix with all off-diagonal elements
c        multipled by 2
c
*        If (iUHF.eq.0) then
            call dcopy_(nblock,[Zero],0,Scr(iadDen),1)
            If (var) Then
               Call GetMem('D1ao','Allo','Real',ipD1ao,nBlock)
               Call Get_D1ao_Var(Work(ipD1ao),nBlock)
            Else
               Call Get_D1ao(ipD1ao,nBlock)
            End If
            Do i = 1,nBlock
               SCR(i) = Work(ipD1ao+i-1)
            End Do
            Call GetMem('Dens','Free','Real',ipD1ao,nBlock)
*        End if
         iadC1 =iadDen+nblock
*
      Else
*
*        Make iadDen to point at the MO vectors
         ipVec=ip_of_Work(Vec(1))
         ipScr=ip_of_Work(Scr(1))
         iadDen=ipVec-ipScr+1
         iadC1=1
         If (iUHF.eq.1) Then
            iadDen_ab=iadDen+nDim**2
            iOcc_ab=1+nDim
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
c     Scan the ONEINT file for multipole moment operators
c
*
      iadC2 =iadC1+3
      iadNuc=iadC2+3
c
C     Write (*,*) ' Starting scan of ONEINT for multipole moments'
      do 100 i=0,99
        NxtOpr = .False.
        nComp=(i+1)*(i+2)/2
c
        iadEl =iadNuc+nComp
        iadLab=iadEl +nComp
        If (.Not.Short) Then
           Call GetMem('iadEl1','Allo','Real',iadEl_Work,nComp*mDim)
           ip_Scr=ip_of_Work(Scr(1))
           iadEl=iadEl_Work-(ip_Scr-1)
        End If
        iscr=nscr+10+4*nComp
        If (iscr.gt.maxscr) Then
           If (.Not.Short) Call Free_Work(iadEl_Work)
           Go To 999
        End If
c
        call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
        call dcopy_(nComp*mDim,[Zero],0,Scr(iadEl ),1)
        call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
        write (label,'(a,i2)') 'MLTPL ',i
        do 101 iComp=1,nComp
          irc=-1
          iopt=1
          Call iRdOne (irc,iopt,label,iComp,idum,iSmLbl)
          if (irc.eq.0) mInt=idum(1)
          if (irc.ne.0) go to 101
          NxtOpr = .True.
          irc=-1
          iopt=0
          iadOpr=iadLab+2*nComp
          Call RdOne (irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
          if (irc.ne.0) go to 101
          If (mInt.ne.0)
     &       Call CmpInt(scr(iadOpr),mInt,nBas,nIrrep,iSmLbl)
          scr(iadNuc+icomp-1)=scr(iadOpr+mInt+3)
          if (iComp.eq.1) then
            do 102 k=0,2
              scr(iadC1+k)=scr(iadOpr+mInt+k)
              scr(iadC2+k)=scr(iadOpr+mInt+k)
102         continue
          endif
          If (mInt.eq.0) Go To 101
          Call Xprop(short,ifallorb,
     &               nIrrep,nBas,
     &               nBlock,Scr(iadDen),nDim,Occ,Thrs,
     &               nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim))
          If (.Not.Short.and.iUHF.eq.1) Call Xprop(short,ifallorb,
     &               nIrrep,nBas,
     &               nBlock,Scr(iadDen_ab),nDim,Occ(iOcc_ab),Thrs,
     &               nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim+nDim))
101     continue
        If (.Not.NxtOpr) Then
           If (.Not.Short) Call Free_Work(iadEl_Work)
           Go To 199
        End If
        iadTmt=iadOpr+nblock
        if (i.le.4) then
          iadTmp=iadTmt+nComp**2
          iscr=iadTmp+nComp
          If (iscr.gt.maxscr) Then
             If (.Not.Short) Call Free_Work(iadEl_Work)
             Go To 999
          End If
        else
          iadTmp=iadTmt
        endif
c
        Call C_F_Pointer(C_Loc(scr(iadLab)),cScr,[1])
        Call prop (short,label,scr(iadC1),scr(iadC2),
     &             nirrep,mBas,mDim,occ,Thrs,
     &             scr(iadEl),scr(iadNuc),i,cScr,
     &             scr(iadTmt),scr(iadTmp),ifallorb)
        Nullify(cScr)
        If (.Not.Short) Call Free_Work(iadEl_Work)
100   continue
c
c     Scan 'ONEINT' for electric field integrals
c
199   continue
*                                                                      *
************************************************************************
*                                                                      *
*
C     Write (*,*) ' Starting scan of ONEINT for various elec. field integrals'
*
      do 210 iEF=0,2
         nComp=(iEF+1)*(iEF+2)/2
*
         iadEl =iadNuc+nComp
         iadLab=iadEl +nComp
         If (.Not.Short) Then
            Call GetMem('iadEl2','Allo','Real',iadEl_Work,nComp*mDim)
            ip_Scr=ip_of_Work(Scr(1))
            iadEl=iadEl_Work-(ip_Scr-1)
         End If
*        create vectors to store the sums of electronic and nuclear components over all centers
         Call GetMem('ElSum','Allo','Real',iadElSum,nComp)
         Call GetMem('NucSum','Allo','Real',iadNucSum,nComp)
         Call DZero(Work(iadElSum),nComp)
         Call DZero(Work(iadNucSum),nComp)
*
*        loop over different operator origins (max.99999)
*
         maxCen=99999
         nCen=0
         Do 200 i=1,maxCen
            call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
            call dcopy_(nComp*mDim,[Zero],0,Scr(iadEl ),1)
            call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
            Write (label,'(a,i1,i5)') 'EF',iEF,i
            NxtOpr=.False.
            Do 201 iComp=1,nComp
               irc=-1
               iopt=1
               Call iRdOne (irc,iopt,label,iComp,idum,iSmLbl)
               If (irc.eq.0) mInt=idum(1)
               If (irc.ne.0) go to 201
               NxtOpr = .True.
               irc=-1
               iopt=0
               iadOpr=iadLab+2*nComp
               Call RdOne (irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
               If (irc.ne.0) Go To 201
               If (mInt.ne.0)
     &            Call CmpInt(scr(iadOpr),mInt,nBas,nIrrep,iSmLbl)
               scr(iadNuc+icomp-1)=scr(iadOpr+mInt+3)
               If (iComp.eq.1) Then
                  Do k=0,2
                    scr(iadC1+k)=scr(iadOpr+mInt+k)
                    scr(iadC2+k)=scr(iadOpr+mInt+k)
                  End Do
               Endif
               If (mInt.eq.0) Go To 201
               Call Xprop(short,ifallorb,
     &                    nIrrep,nBas,
     &                    nBlock,Scr(iadDen),nDim,Occ,Thrs,
     &                    nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim))
               If (.Not.Short.and.iUHF.eq.1) Call Xprop(short,ifallorb,
     &                    nIrrep,nBas,
     &                    nBlock,Scr(iadDen_ab),nDim,Occ(iOcc_ab),Thrs,
     &                    nblock,Scr(iadOpr),
     &                    Scr(iadEl+(iComp-1)*mDim+nDim))
201         Continue
            If (.Not.NxtOpr) Then
               If (.Not.Short) Call Free_Work(iadEl_Work)
               Go to 299
            End If
*
            Call C_F_Pointer(C_Loc(scr(iadLab)),cScr,[1])
            Call Prop (short,label,scr(iadC1),scr(iadC2),
     &                 nirrep,mBas,mDim,occ,Thrs,
     &                 scr(iadEl),scr(iadNuc),iEF,cScr,
     &                 scr(iadTmt),scr(iadTmp),ifallorb)
            Nullify(cScr)
*           add the components to the sums, and update the total number of centers
            Do iComp=0,nComp-1
              iInd1=iadElSum+iComp
              Do iOcc=0,mDim-1
                Work(iInd1)=Work(iInd1)+Scr(iadEl+iComp*mDim+iOcc)
              End Do
            End Do
            Call DaXpY_(nComp,One,Scr(iadNuc),1,Work(iadNucSum),1)
            nCen=i
200      Continue
         If (.Not.Short) Call Free_Work(iadEl_Work)
*
299      Continue
         If (nCen.gt.0) Then
*          set the tolerance according to the total number of centers
*          (assuming error scales with sqrt(ncen))
           iTol=5
           iTol=iTol-NInt(Half*Log10(Dble(nCen)))
           Write (label,'(a,i1,a)') 'EF',iEF,'   el'
           Call Add_Info(label,Work(iadElSum),nComp,iTol)
           Write (label,'(a,i1,a)') 'EF',iEF,'  nuc'
           Call Add_Info(label,Work(iadNucSum),nComp,iTol)
         End If
         Call GetMem('ElSum','Free','Real',iadElSum,nComp)
         Call GetMem('NucSum','Free','Real',iadNucSum,nComp)
210   Continue
*
*                                                                      *
************************************************************************
*                                                                      *
*
C     Write (*,*) ' Starting scan of ONEINT for various contact term integrals'
*
      nComp=1
*
      iadEl =iadNuc+nComp
      iadLab=iadEl +nComp
      If (.Not.Short) Then
         Call GetMem('iadEl2','Allo','Real',iadEl_Work,nComp*mDim)
         ip_Scr=ip_of_Work(Scr(1))
         iadEl=iadEl_Work-(ip_Scr-1)
      End If
*     create vectors to store the sums of electronic and nuclear components over all centers
      Call GetMem('ElSum','Allo','Real',iadElSum,nComp)
      Call GetMem('NucSum','Allo','Real',iadNucSum,nComp)
      Call DZero(Work(iadElSum),nComp)
      Call DZero(Work(iadNucSum),nComp)
*
*     loop over different operator origins (max.99999)
*
      maxCen=99999
      nCen=0
      Do 300 i=1,maxCen
         call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
         call dcopy_(nComp*mDim,[Zero],0,Scr(iadEl ),1)
         call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
         Write (label,'(a,i5)') 'CNT',i
         NxtOpr=.False.
*
         iComp=1
         irc=-1
         iopt=1
         Call iRdOne (irc,iopt,label,iComp,idum,iSmLbl)
         If (irc.eq.0) mInt=idum(1)
         If (irc.ne.0) go to 301
         NxtOpr = .True.
         irc=-1
         iopt=0
         iadOpr=iadLab+2*nComp
         Call RdOne (irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
         If (irc.ne.0) Go To 301
         If (mInt.ne.0) Call CmpInt(scr(iadOpr),mInt,nBas,nIrrep,iSmLbl)
         scr(iadNuc+icomp-1)=scr(iadOpr+mInt+3)
         Do k=0,2
            scr(iadC1+k)=scr(iadOpr+mInt+k)
            scr(iadC2+k)=scr(iadOpr+mInt+k)
         End Do
         If (mInt.eq.0) Go To 301
         Call Xprop(short,ifallorb,
     &              nIrrep,nBas,
     &              nBlock,Scr(iadDen),nDim,Occ,Thrs,
     &              nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim))
         If (.Not.Short.and.iUHF.eq.1) Call Xprop(short,ifallorb,
     &              nIrrep,nBas,
     &              nBlock,Scr(iadDen_ab),nDim,Occ(iOcc_ab),Thrs,
     &              nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim+nDim))
 301     Continue
         If (.Not.NxtOpr) Then
            If (.Not.Short) Call Free_Work(iadEl_Work)
            Go To 399
         End If
*
         Call C_F_Pointer(C_Loc(scr(iadLab)),cScr,[1])
         Call Prop (short,label,scr(iadC1),scr(iadC2),
     &              nirrep,mBas,mDim,occ,Thrs,
     &              scr(iadEl),scr(iadNuc),iEF,cScr,
     &              scr(iadTmt),scr(iadTmp),ifallorb)
         Nullify(cScr)
*        add the components to the sums, and update the total number of centers
         Do iComp=0,nComp-1
           iInd1=iadElSum+iComp
           iInd2=iadEl+iComp*mDim
           Do iOcc=0,mDim-1
             Work(iInd1)=Work(iInd1)+Scr(iInd2+iOcc)
           End Do
         End Do
         Call DaXpY_(nComp,One,Scr(iadNuc),1,Work(iadNucSum),1)
         nCen=i
300      Continue
         If (.Not.Short) Call Free_Work(iadEl_Work)
*
399      Continue
         If (nCen.gt.0) Then
*          set the tolerance according to the total number of centers
*          (assuming error scales with sqrt(ncen))
           iTol=5
           iTol=iTol-NInt(Half*Log10(Dble(nCen)))
           Write (label,'(a,a)') 'CNT','   el'
           Call Add_Info(label,Work(iadElSum),nComp,iTol)
           Write (label,'(a,a)') 'CNT','  nuc'
           Call Add_Info(label,Work(iadNucSum),nComp,iTol)
         End If
         Call GetMem('ElSum','Free','Real',iadElSum,nComp)
         Call GetMem('NucSum','Free','Real',iadNucSum,nComp)
*                                                                      *
************************************************************************
*                                                                      *
C     Write (*,*) ' Starting scan of ONEINT diamagnetic shielding'
c
      nComp=9
      lpole=2
c
      iadEl =iadNuc+nComp
      iadLab=iadEl +nComp
      If (.Not.Short) Then
         Call GetMem('iadEl3','Allo','Real',iadEl_Work,nComp*mDim)
         ip_Scr=ip_of_Work(Scr(1))
         iadEl=iadEl_Work-(ip_Scr-1)
      End If
c
      maxGG =99
      maxCen=99
c     loop over different gauge origins (max.99)
      do 400 j=1,maxGG
         call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
         call dcopy_(nComp*mDim,[Zero],0,Scr(iadEl ),1)
         call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
        jRC = 0
c       loop over different operator origins (max.99)
        do 401 i=1,maxCen
          write (label,'(a,i2,i2)') 'DMS ',j,i
          NxtOpr = .False.
          do 402 iComp=1,nComp
            irc=-1
            iopt=1
            Call iRdOne (irc,iopt,label,iComp,idum,iSmLbl)
            if (irc.eq.0) mInt=idum(1)
            if (irc.ne.0) go to 402
            NxtOpr = .True.
            irc=-1
            iopt=0
            iadOpr=iadLab+2*nComp
            Call RdOne (irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
            if (irc.ne.0) go to 402
            If (mInt.ne.0)
     &         Call CmpInt(scr(iadOpr),mInt,nBas,nIrrep,iSmLbl)
            scr(iadNuc+icomp-1)=scr(iadOpr+mInt+3)
            if (iComp.eq.1) then
              do 403 k=0,2
                scr(iadC1+k)=scr(iadOpr+mInt+k)
403           continue
            endif
            if (iComp.eq.2) then
              do 404 k=0,2
                scr(iadC2+k)=scr(iadOpr+mInt+k)
404           continue
            endif
            If (mInt.eq.0) Go To 402
            Call Xprop(short,ifallorb,
     &                 nIrrep,nBas,
     &                 nBlock,Scr(iadDen),nDim,Occ,Thrs,
     &                 nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim))
            If (.Not.Short.and.iUHF.eq.1) Call Xprop(short,ifallorb,
     &                 nIrrep,nBas,
     &                 nBlock,Scr(iadDen_ab),nDim,Occ(iOcc_ab),Thrs,
     &                 nblock,Scr(iadOpr),
     &                 Scr(iadEl+(iComp-1)*mDim+nDim))
402       continue
          If (.Not.NxtOpr) Go To 4000
c
          Call C_F_Pointer(C_Loc(scr(iadLab)),cScr,[1])
          call prop (short,label,scr(iadC1),scr(iadC2),
     &               nirrep,mBas,mDim,occ,Thrs,
     &               scr(iadEl),scr(iadNuc),lpole,cScr,
     &               scr(iadTmt),scr(iadTmp),ifallorb)
          Nullify(cScr)
           jRC = 1
401     continue
4000    If (jRC.eq.0) Then
           If (.Not.Short) Call Free_Work(iadEl_Work)
           Go To 499
        End If
400   continue
      If (.Not.Short) Call Free_Work(iadEl_Work)
*                                                                      *
************************************************************************
*                                                                      *
499   continue
#ifdef _DEBUG_
#endif
      If (iPL.ge.2) Then
         Call CollapseOutput(0,'   Molecular properties:')
         Write(6,*)
      End If
      Return
*                                                                      *
************************************************************************
*                                                                      *
  999 continue
      Write (6,'(//1x,a/1x,a/1x,a,i8,a,i8)')
     &  ' Warning:',
     &  ' Not enough scratch area to perform calculations',
     &  ' Needed at least:',iscr,'   available:',maxscr
*                                                                      *
************************************************************************
*                                                                      *
      If (iPL.ge.2) Then
         Call CollapseOutput(0,'   Molecular properties:')
         Write(6,*)
      End IF
      Return
      End Subroutine Prpt_Internal
*
      End

