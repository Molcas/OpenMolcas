!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************
      Subroutine Prpt_old(nirrep,nbas,ndim,n2dim,vec,occ,maxscr,scr)
!***********************************************************************
!
!     Purpose: calculation of expectation values of different
!              operators as available on the 'ONEINT' file
!
!     Caution: before calling this subroutine one needs to
!              open the ONEINT file
!
!     Calling parameters:
!
!       nIrRep            number of irreducible representations
!
!       nBas(0:nIrRep)    number of basis functions in each repre-
!                         sentation
!
!       ndim              total number of basis functions
!
!       n2dim             sum(i,i=0,nirrep-1)(nbas(i)**2): elements
!                         of all vectors in all representations
!
!       vec(1:n2dim)      eigenvectors, all for each representation
!
!       occ(1:ndim)       occupation number for all eigenvectors
!
!       maxscr            the maximum available size of the
!                         scratch area
!
!       scr(1:maxscr)     a scratch area whose size, nscr, can be
!                         calculated in the following way:
!
!                         nscr =sum(i,i=0,nirrep-1)(nbas(i)*
!                                  (nbas(i)+1)/2)
!                               +3
!                               +3
!                               +nComp
!                               +nComp
!                               +2*nComp
!                               +ndim*(ndim+1)/2+4
!                               +2*ntComp*(ntComp+1)
!
!                         and should not exceed maxscr.
!                         nComp is the number of cartresian compo-
!                         nets for the given operator
!                         ntComp is the number of components of
!                         the l-th moment opartors which are trans-
!                         formed into l-pole moment. Currently,
!                         ntComp=15 (hexadecpole moments)
!
!
!
! 1991 R. Lindh, Dept. of Theor. Chem. Univ. of Lund, Sweden.          *
!***********************************************************************
      Implicit real*8 (a-h,o-z)
!
#include "real.fh"
!
      Character*8 label
      Logical short, NxtOpr, ifallorb
      Integer nbas(0:nirrep)
      Real*8 vec(1:n2dim), occ(1:ndim), scr(1:maxscr)
      Dimension idum(1)
!
      Call Prpt_old_Internal(Scr)
!
!     This is to allow type punning without an explicit interface
      Contains
      Subroutine Prpt_old_Internal(Scr)
      Use Iso_C_Binding
      Real*8, Target :: Scr(*)
      Character, Pointer :: cScr(:)
      Write (6,*)
      Call CollapseOutput(1,'   Molecular properties:')
      Write (6,'(3X,A)')    '   ---------------------'
      Write (6,*)
      short=.true.
      ifallorb=.false.
      mDim = 1
      iadopr=-1 ! dummy initialize
      iadtmt=-1 ! dummy initialize
      iadtmp=-1 ! dummy initialize
!
      nblock=0
      do 1 i=0,nirrep-1
        nblock=nblock+nbas(i)*(nbas(i)+1)/2
1     continue
      nfblock=ndim*(ndim+1)/2
!
!     estimate the minimum scratch area needed for calculatig
!     the average values of a 3 component operator
!
      nscr=nblock+nfblock
      iscr=nscr+76
!
      if (iscr.gt.maxscr) Go To 999
!
      iadDen=1
!
!     calculate the density matrix with all off-diagonal elements
!     multipled by 2
!
      call dcopy_(nblock,[Zero],0,Scr(iadDen),1)
      iCount=iadDen
      iVec  =0
      iOcc  =0
      do i=0,nIrrep-1
         do ii=1,nBas(i)
            iOcc=iOcc+1
            jCount=iCount
            do 4 il=1,nBas(i)
               Do ir=1,il-1
                  Scr(jCount)=Scr(jCount)                               &
     &                       +Two*Vec(iVec+il)*Vec(iVec+ir)*Occ(iOcc)
                   jCount=jCount+1
                End Do
                Scr(jCount)=Scr(jCount)                                 &
     &                     +Vec(iVec+il)*Vec(iVec+il)*Occ(iOcc)
                jCount=jCount+1
    4        continue
            iVec=iVec+nBas(i)
         End Do
         iCount=iCount+nbas(i)*(nBas(i)+1)/2
      End Do
!                                                                      *
!***********************************************************************
!                                                                      *
!     Scan the ONEINT file for multipole moment operators
!
      iadC1 =iadDen+nblock
      iadC2 =iadC1+3
      iadNuc=iadC2+3
!
!     Write (*,*) ' Starting scan of ONEINT for multipole moments'
      do 100 i=1,99
        NxtOpr = .False.
        nComp=(i+1)*(i+2)/2
!
        iadEl =iadNuc+nComp
        iadLab=iadEl +nComp
        iscr=nscr+10+4*nComp
        if (iscr.gt.maxscr) Go To 999
!
        call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
        call dcopy_(nComp,[Zero],0,Scr(iadEl ),1)
        call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
        write (label,'(a,i2)') 'MLTPL ',i
        do 101 iComp=1,nComp
          irc=-1
          iopt=1
          Call iRdOne (irc,iopt,label,iComp,idum,iSmLbl)
          if (irc.eq.0) nInt=idum(1)
          if (irc.ne.0) go to 101
          NxtOpr = .True.
          irc=-1
          iopt=0
          iadOpr=iadLab+2*nComp
          Call RdOne (irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
          if (irc.ne.0) go to 101
          If (nInt.ne.0)                                                &
     &       Call CmpInt(scr(iadOpr),nInt,nBas,nIrrep,iSmLbl)
          scr(iadNuc+icomp-1)=scr(iadOpr+nInt+3)
          if (iComp.eq.1) then
            do 102 k=0,2
              scr(iadC1+k)=scr(iadOpr+nInt+k)
              scr(iadC2+k)=scr(iadOpr+nInt+k)
102         continue
          endif
          If (nInt.eq.0) Go To 101
          Call Xprop(short, ifallorb,                                   &
     &               nIrrep,nBas,                                       &
     &               nBlock,Scr(iadDen),nDim,Occ,dummy,                 &
     &               nblock,Scr(iadOpr),Scr(iadEl+iComp-1))
101     continue
        if (.Not.NxtOpr) Go To 199
        iadTmt=iadOpr+nblock
        if (i.le.4) then
          iadTmp=iadTmt+nComp**2
          iscr=iadTmp+nComp
          if (iscr.gt.maxscr) Go To 999
        else
          iadTmp=iadTmt
        endif
!
        Call C_F_Pointer(C_Loc(Scr(iadLab)),cScr,[1])
        call prop (short,label,scr(iadC1),scr(iadC2),                   &
     &             nirrep,nBas,mDim,occ,dummy,                          &
     &             scr(iadEl),scr(iadNuc),i,cScr,                       &
     &             scr(iadTmt),scr(iadTmp),ifallorb)
        Nullify(cScr)
100   continue
!                                                                      *
!***********************************************************************
!                                                                      *
!     Scan 'ONEINT' for electric field integrals
!
199   continue
!
!
!     Write (*,*) ' Starting scan of ONEINT for various elec. field integrals'
!
      do 210 iEF=0,2
        nComp=(iEF+1)*(iEF+2)/2
!
        iadEl =iadNuc+nComp
        iadLab=iadEl +nComp
!
!
!       loop over differnt operator origins (max.9999)
!
        maxCen=9999
        do 200 i=1,maxCen
          call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
          call dcopy_(nComp,[Zero],0,Scr(iadEl ),1)
          call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
          write (label,'(a,i1,i5)') 'EF',iEF,i
          NxtOpr=.False.
          do 201 iComp=1,nComp
            irc=-1
            iopt=1
            Call iRdOne (irc,iopt,label,iComp,idum,iSmLbl)
            if (irc.eq.0) nInt=idum(1)
            if (irc.ne.0) go to 201
            NxtOpr = .True.
            irc=-1
            iopt=0
            iadOpr=iadLab+2*nComp
            Call RdOne (irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
            if (irc.ne.0) go to 201
            If (nInt.ne.0)                                              &
     &         Call CmpInt(scr(iadOpr),nInt,nBas,nIrrep,iSmLbl)
            scr(iadNuc+icomp-1)=scr(iadOpr+nInt+3)
            if (iComp.eq.1) then
              do 202 k=0,2
                scr(iadC1+k)=scr(iadOpr+nInt+k)
                scr(iadC2+k)=scr(iadOpr+nInt+k)
202           continue
            endif
            If (nInt.eq.0) Go To 201
            Call Xprop(short, ifallorb,                                 &
     &                 nIrrep,nBas,                                     &
     &                 nBlock,Scr(iadDen),nDim,Occ,dummy,               &
     &                 nblock,Scr(iadOpr),Scr(iadEl+iComp-1))
201       continue
          if (.Not.NxtOpr) go to 299
!
          Call C_F_Pointer(C_Loc(Scr(iadLab)),cScr,[1])
          call prop (short,label,scr(iadC1),scr(iadC2),                 &
     &               nirrep,nBas,mDim,occ,dummy,                        &
     &               scr(iadEl),scr(iadNuc),i,cScr,                     &
     &               scr(iadTmt),scr(iadTmp),ifallorb)
          Nullify(cScr)
200     continue
!
299      continue
210   continue
!                                                                      *
!***********************************************************************
!                                                                      *
!     Write (*,*) ' Starting scan of ONEINT diamagnetic shielding'
!
      nComp=9
!
      iadEl =iadNuc+nComp
      iadLab=iadEl +nComp
!
!
!
      maxGG =99
      maxCen=99
!     loop over differnt gauge origins (max.99)
      do 400 j=1,maxGG
         call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
         call dcopy_(nComp,[Zero],0,Scr(iadEl ),1)
         call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
        jRC = 0
!       loop over differnt operator origins (max.99)
        do 401 i=1,maxCen
          write (label,'(a,i2,i2)') 'DMS ',j,i
          NxtOpr = .False.
          do 402 iComp=1,nComp
            irc=-1
            iopt=1
            Call iRdOne (irc,iopt,label,iComp,idum,iSmLbl)
            if (irc.eq.0) nInt=idum(1)
            if (irc.ne.0) go to 402
            NxtOpr = .True.
            irc=-1
            iopt=0
            iadOpr=iadLab+2*nComp
            Call RdOne (irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
            if (irc.ne.0) go to 402
            If (nInt.ne.0)                                              &
     &         Call CmpInt(scr(iadOpr),nInt,nBas,nIrrep,iSmLbl)
            scr(iadNuc+icomp-1)=scr(iadOpr+nInt+3)
            if (iComp.eq.1) then
              do 403 k=0,2
                scr(iadC1+k)=scr(iadOpr+nInt+k)
403           continue
            endif
            if (iComp.eq.2) then
              do 404 k=0,2
                scr(iadC2+k)=scr(iadOpr+nInt+k)
404           continue
            endif
            If (nInt.eq.0) Go To 402
            Call Xprop(short,ifallorb,                                  &
     &                 nIrrep,nBas,                                     &
     &                 nBlock,Scr(iadDen),nDim,Occ,dummy,               &
     &                 nblock,Scr(iadOpr),Scr(iadEl+iComp-1))
402       continue
          If (.Not.NxtOpr) Go To 4000
!
          Call C_F_Pointer(C_Loc(Scr(iadLab)),cScr,[1])
          call prop (short,label,scr(iadC1),scr(iadC2),                 &
     &               nirrep,nBas,mDim,occ,dummy,                        &
     &               scr(iadEl),scr(iadNuc),i,cScr,                     &
     &               scr(iadTmt),scr(iadTmp),ifallorb)
          Nullify(cScr)
           jRC = 1
401     continue
4000    If (jRC.eq.0) Go To 499
400   continue
!
!                                                                      *
!***********************************************************************
!                                                                      *
499   continue
!

      Call CollapseOutput(0,'   Molecular properties:')
      Write(6,*)
      Return
!
  999 continue
      write (6,'(//1x,a/1x,a/1x,a,i8,a,i8)')                            &
     &  ' Warrning:',                                                   &
     &  ' Not enough scratch area to perform calculations',             &
     &  ' Needed at least:',iscr,'   available:',maxscr
!
      Call CollapseOutput(0,'   Molecular properties:')
      Write(6,*)
      Return
      End Subroutine Prpt_old_Internal
!
      end
