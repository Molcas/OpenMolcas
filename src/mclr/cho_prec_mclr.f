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
* Copyright (C) Mickael G. Delcey                                      *
************************************************************************
      SUBROUTINE CHO_Prec_MCLR(CMO,nIsh,nAsh,LuAChoVec,LuChoInt)
************************************************************************
*                                                                      *
*  Author : M. G. Delcey                                               *
*                                                                      *
*  Purpose: form 2-electron integrals needed for the preconditioner    *
*           those are (ii|pq) and (ip|iq) with i inactive              *
*           and p, q active+virtual                                    *
*           as well as (tu|pq) and (tp|uq) with t, u active            *
*           and p, q all molecular orbitals                            *
*                                                                      *
*           For small active spaces, this is much less than            *
*           the full list of integrals!                                *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Real*8 CMO(*)
#include "warnings.fh"
      Character*13 SECNAM
      Parameter (SECNAM = 'CHO_PREC_MCLR')
      Integer   ISTLT(8),ISTSQ(8),ISSQ(8,8),ipLpq(8)
      Integer   LuAChoVec(8),ipCMOt(8),LuChoInt(2)
      Integer  nAsh(8),nIsh(8),nIshb(8),nIshe(8),iSkip(8),
     &        nAshb(8),nAshe(8)
      Real*8    tread(2),ttran(2),tform(2) ,tform2(2) ,
     &                            tforma(2),tforma2(2),tMO(2)
      Logical timings
#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
      Character*50 CFmt
      parameter ( N2 = InfVec_N2 )
      parameter (zero = 0.0D0, one = 1.0D0, xone=-1.0D0)
      Character*6 mode
      Logical DoRead,taskleft
      parameter (DoRead = .false. )
      Integer   Cho_LK_MaxVecPerBatch
      External  Cho_LK_MaxVecPerBatch
************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      InfVec(i,j,k) = iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
******
      nDimRS(i,j) = iWork(ip_nDimRS-1+nSym*(j-1)+i)
************************************************************************
      timings=.false.
      CALL CWTIME(TCstart1,TWstart1)
      do i=1,2            ! 1 --> CPU   2 --> Wall
         tread(i) = zero  !time read vectors
         ttran(i) = zero  !time transform vectors
         tform(i) = zero  !time form integrals
         tform2(i) = zero  !time form integrals
         tforma(i) = zero  !time form integrals
         tforma2(i) = zero  !time form integrals
         tMO(i)     = zero  !time for final MO transform
      end do
      MaxVecPerBatch=Cho_LK_MaxVecPerBatch()
      iLoc = 3
*
      ISTLT(1)=0
      ISTSQ(1)=0
      DO ISYM=2,NSYM
        NBB=NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NBB
        ISTSQ(iSYM)=ISTSQ(iSYM-1)+NBAS(ISYM-1)**2
      END DO
*
**    Start with big loop over symmetries
*
      Do jsym=1,nsym
        iAdr=0
*
**      Compute some sizes
*
        nip=0
        ntp=0
        npq=0
        maxpq=0
        maxtpq=0
        ntoti=0
        ntota=0
        Do isymb=1,nsym
          iSyma=MulD2h(iSymb,jsym)

          npq  = npq + nBas(iSymb)**2
          maxpq  =max(npq,nBas(iSymb)**2)
          maxtpq  =max(maxtpq,nAsh(iSyma)*nBas(iSymb)**2)

!         For inactive half-transformed Cho vector + Lii^J
          nip  =nip + nIsh(iSyma)*(nBas(isymb)+1)
!         For active half-transformed Cho vector + Lij^J
          ntp  =ntp + nAsh(iSyma)*(nBas(isymb)+nAsh(iSyma))

          ntoti=ntoti+nIsh(isymb)
          ntota=ntota+nAsh(isymb)
        End Do
*
        NumCV=NumCho(jSym)
        Call GAIGOP_SCAL(NumCV,'max')
        If (NumCV .lt. 1) GOTO 999
*
*
        maxRS=0
        JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
        JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec
        Call GAIGOP_SCAL(JRED1,'min')
        Call GAIGOP_SCAL(JRED2,'max')
        Do Jred=JRED1,JRED2
          maxRS = max(maxRS,nDimRS(JSYM,JRED))
        End Do
*
**      Check for maxmem to see if integrals would fit in memory
**      or if batching is required
*
        Call GetMem('MaxM','Max','Real',KDUM,LWORK)
        Call GAIGOP_SCAL(LWORK,'min')
*
        Do i=1,nsym
          nIshb(i)=0
          nAshb(i)=0
        EndDo
*
**      Loop over i and t batches
**      each batch will increase the i/o, therefore we preferably want to
**      have the whole iiab and tupq in-core!
*
*
**      First, do we have enough memory at all!
*
        memneeded=maxRS+max(nip,ntp)                  ! for 1 Jbatch
        memneeded=max(memneeded,maxpq)                ! for MO transform
        If (ntota.gt.0) Then
           memneeded=memneeded+ maxtpq                ! for 1 (ta|ub)
!          for 1 (tu|ab) and 1 reduced set
           If (jsym.eq.1)
     &        memneeded=memneeded+ ntota*(maxpq+maxRS)
        ElseIf (ntoti.gt.0) Then
           memneeded=memneeded+ maxpq                 ! for 1 (ia|ib)
!          for 1 (ii|ab) and 1 reduced set
           If (jsym.eq.1) memneeded=memneeded+npq+maxRS
        EndIf
*
        If (memneeded.gt.lWork) Then
          WRITE(6,*) SECNAM//': Insufficient memory for I/T batch'
          WRITE(6,*) 'LWORK= ',LWORK
          WRITE(6,*) 'min. mem. need= ',memneeded
          WRITE(6,*) 'maxRS+nip/ntp= ',maxRS+max(nip,ntp)
          WRITE(6,*) 'maxpq        = ',maxpq
          If (ntota.gt.0) Then
            WRITE(6,*) 'maxtpq       = ',maxtpq
            If (jsym.eq.1)
     &       WRITE(6,*) 'ntota*(maxpq,maxRS) = ',ntota*maxpq,ntota*maxRS
          Else
            If (jsym.eq.1)
     &        WRITE(6,*) 'npq,maxRS = ',npq,maxRS
          EndIf
          CALL Quit(_RC_MEMORY_ERROR_)
        EndIf
        lWork=lWork-max(maxpq,maxRS+max(nip,ntp))
*
**      How many nab + nrs can be stored in memory while stil having place for
**      at least 1 nRS, 1 maxpq and 1 nIP
**      (actually nIP will be smaller but that does not matter very much)
**
**      nIshe is the number of inactive orbitals dealt in the batch
**      nIshb is inactive dealt with after the previous batches
*
 50     Continue
        lWorke=lWork
        Do i=1,nsym
          nIshe(i)=0
          nAshe(i)=0
        End Do
        taskleft=.false.
        Libatch=0
        Do i=1,nsym
          k=MulD2h(i,jsym)
*
          nab=0
          nRS=0
          If (jsym.eq.1) Then
            nab=npq
            nRS=maxRS
          EndIf
          nab2=nBas(k)**2
          nileft=nIsh(i)-nIshb(i)
          If ((nab+nab2)*nileft.gt.0) Then
            nIshe(i)=min(lWorke/(nrs+nab+nab2),nileft)
            lWorke=lWorke-(nrs+nab+nab2)*nIshe(i)
            libatch=libatch+(nrs+nab+nab2)*nIshe(i)
            If (nIshe(i).ne.nileft)  Then
              taskleft=.true.
              Go to 10
            EndIf
          EndIf
*
        End Do
 10     Continue
*
**      Update nip and compute sum(nIshe)
*
        nip=0
        ntotie=0
        Do i=1,nsym
          k=MulD2h(i,jsym)
          nip=nip+nIshe(i)*nBas(k)
          ntotie=ntotie+nIshe(i)
        End Do
        nip=nip+ntotie ! for Lii^J
*
        Call GetMem('iiab','Allo','Real',ipiiab,libatch)
        nab=0
        If ((jsym.eq.1).and.(ntotie.gt.0)) Then
          nab=npq
          Call GetMem('iirs','Allo','Real',ipiirs,ntotie*maxRS)
        EndIf
        ipiaib=ipiiab+nab*ntotie
        Call dcopy_(libatch,[0.0d0],0,Work(ipiiab),1)

*init for compilers
        iptupq=ip_Dummy
        iptpuq=iptupq
*
        If (taskleft) Then
           ntotae=0
           ntue=0
           ntp=0      ! do not allocate those
           labatch=0
           write(6,*) 'Batching loop i'
        Else
*
**        Batching T loop
*
           labatch=0
           Do i=1,nsym
             k=MulD2h(i,jsym)
             nab=0
             nRS=0
             If (jsym.eq.1) Then
               nab=ntota*npq
               nRS=ntota*maxRS
             EndIf
             nab2=nAsh(i)*nBas(k)**2
             naleft=nAsh(i)-nAshb(i)
             If ((nab+nab2)*naleft.gt.0) Then
               nAshe(i)=min(lWorke/(nrs+nab+nab2),naleft)
               lWorke=lWorke-(nrs+nab+nab2)*nAshe(i)
               labatch=labatch+(nrs+nab+nab2)*nAshe(i)
               If (nAshe(i).ne.naleft)  Then
                 taskleft=.true.
                 Go to 11
               EndIf
             EndIf
           End Do
 11        Continue
*
           ntotae=0
           ntue=0
           ntp=0
           Do i=1,nsym
             k=MulD2h(i,jsym)
             ntotae=ntotae+nAshe(i)
             ntue=ntue+nAshe(i)*nAsh(i)
             ntp  =ntp + nAsh(i)*(nBas(k)+nAshe(i))
           End Do

           If (ntotae.gt.0)
     &        Call GetMem('tupq','Allo','Real',iptupq,labatch)
           nab=0
           If (jsym.eq.1) Then
             nab=npq
             If (ntue.gt.0)
     &          Call GetMem('turs','Allo','Real',ipturs,ntue*maxRS)
           EndIf
           iptpuq=iptupq+nab*ntue
           call dcopy_(labatch,[0.0d0],0,Work(iptupq),1)
           If (taskleft) Then
              write(6,*) 'Batching loop a'
           EndIf
        EndIf
*FIXME: can iptpuq be uninitialized?
*
**    Transpose CMO
*
        nCMO=0
        Do i=1,nsym
          nCMO=nCMO+nIshe(i)*nBas(i)
        End Do
        Call GetMem('CMOt','Allo','Real',ipCMOt(1),nCMO)
        ioff =0
        ioff2=0
        Do iSym=1,nsym
          ipCMOt(iSym)=ipCMOt(1)+ioff2
          Do j=1,nIshe(iSym)
            ioff3=ioff+nBas(iSym)*(nIshb(iSym)+j-1)
            call dcopy_(nBas(iSym),CMO(1+ioff3),1,
     &                 Work(ipCMOt(isym)+j-1),nIshe(iSym))
          End Do
          ioff =ioff +nBas(iSym)**2
          ioff2=ioff2+nBas(iSym)*nIshe(iSym)
        End Do
*
**      Loop over reduced sets
*
        Do JRED=JRED1,JRED2
          CALL Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

          If (nVrs.eq.0) GOTO 998  ! no vectors in that (jred,isym)

          if (nVrs.lt.0) then
            Write(6,*)SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
            call Abend
          endif

          Call Cho_X_SetRed(irc,iLoc,JRED)
c         !set index arrays at iLoc
          if(irc.ne.0)then
            Write(6,*)SECNAM//' cho_X_setred non-zero return code.'//
     &                        ' rc= ',irc
            call Abend
          endif

          IREDC=JRED

          nRS = nDimRS(JSYM,JRED)

          If (jSym.eq.1) Then
            If (ntotie.gt.0)
     &         call dcopy_(ntotie*nRS,[0.0d0],0,Work(ipiirs),1)
            If (ntue.gt.0)
     &         call dcopy_(ntue*nRS,[0.0d0],0,Work(ipturs),1)
          EndIf

          Call GetMem('MaxM','Max','Real',KDUM,LWORKe)
          nVec= min(LWORKE/(nRS+max(nip,ntp)),min(nVrs,MaxVecPerBatch))
          If (nVec.lt.1) Then
            WRITE(6,*) SECNAM//': Insufficient memory for J batch'
            WRITE(6,*) 'That should not happen here'
            WRITE(6,*) 'Contact the developers'
            CALL Quit(_RC_MEMORY_ERROR_)
            nBatch = -9999  ! dummy assignment
          EndIf
          LREAD = nRS*nVec
*
          Call GetMem('rsL','Allo','Real',ipLrs,LREAD)
          Call GetMem('ChoT','Allo','Real',ipChoT,max(nIP,ntp)*nVec)
*
          nBatch = (nVrs-1)/nVec + 1
*
**        Loop over batches of Cholesky vectors
*
          DO jBatch=1,nBatch

            If (jBatch.eq.nBatch) Then
               JNUM = nVrs - nVec*(nBatch-1)
            else
               JNUM = nVec
            endif
            JVEC = nVec*(jBatch-1) + iVrs
            IVEC2 = JVEC - 1 + JNUM
************************************************************************
*                                                                      *
*                Let's start the real work                             *
*                                                                      *
************************************************************************
*
**          Read Cholesky vectors
*
            CALL CWTIME(TCR1,TWR1)

            CALL CHO_VECRD(Work(ipLrs),LREAD,JVEC,IVEC2,JSYM,
     &                     NUMV,IREDC,MUSED)

            If (NUMV.le.0 .or.NUMV.ne.JNUM ) then
               rc=77
               RETURN
            End If

            CALL CWTIME(TCR2,TWR2)
            tread(1) = tread(1) + (TCR2 - TCR1)
            tread(2) = tread(2) + (TWR2 - TWR1)
*
************************************************************************
**          MO Half-transformation
**          Liq^J= sum_p Lpq^J Xip
*
            lChoI=0
            Do i=1,nSym

               k = Muld2h(i,JSYM)
               iSkip(k) = Min(1,
     &              nBas(k)*nIshe(i))

               ipLpq(k) = ipChoT + lChoI       ! Lvb,J

               lChoI= lChoI + nIshe(i)*nBas(k)*JNUM

            End Do

            iSwap = 1 ! Lqi,J are returned
            kMOs = 1  !
            nMOs = 1  ! Active MOs (1st set)
*
            CALL CHO_X_getVtra(irc,Work(ipLrs),LREAD,jVEC,JNUM,
     &                        JSYM,iSwap,IREDC,nMOs,kMOs,ipCMOt,
     &                        nIshe,ipLpq,iSkip,DoRead)

            if (irc.ne.0) then
               rc = irc
               RETURN
            endif

            CALL CWTIME(TCR1,TWR1)
            ttran(1) = ttran(1) + (TCR1 - TCR2)
            ttran(2) = ttran(2) + (TWR1 - TWR2)
*
************************************************************************
**          Integral formation
**          (i p | i q ) = sum_J Lip^J  Liq^J
*
            ip1=ipiaib
            Do isym=1,nsym
              ksym=MulD2h(iSym,jsym)
              Do ii=1,nIshe(isym)
*                 ipLip=ipLpq(isym)+nBas(kSym)*(ii-1)
                 ipLip=ipLpq(ksym)+nBas(kSym)*(ii-1)
*
                 Call DGEMM_('N','T',nBas(kSym),nBas(kSym),JNUM,
     &                        1.0d0,Work(ipLip),nBas(kSym)*nIshe(iSym),
     &                              Work(ipLip),nBas(kSym)*nIshe(iSym),
     &                        1.0d0,Work(ip1)  ,nBas(kSym))
                 ip1=ip1+nBas(kSym)**2
              End Do
            End Do
            CALL CWTIME(TCR2,TWR2)
            tform(1) = tform(1) + (TCR2 - TCR1)
            tform(2) = tform(2) + (TWR2 - TWR1)
*
************************************************************************
**          Second MO transformation
**          Lii^J = sum_q Liq^J Xiq
*
            If (jSym.eq.1) Then
*             The end of ipChoT is allocated for Lii^J
              ipLii=ipChoT+(nIP-ntotie)*nVec
              ipMO=1
              Do isym=1,nsym
                Do ii=1,nIshe(iSym)
                  ipMO=ipMO+ISTSQ(iSym)
                  ipLip=ipLpq(iSym)+nBas(iSym)*(ii-1)
                  ipMOi=ipMO+(nIshb(isym)+ii-1)*nBas(iSym)
                  ipLiii=ipLii+(ii-1)*JNUM
                  Call dGeMV_('T',nBas(iSym),JNUM,
     &                       1.0d0,Work(ipLip),nBas(iSym)*nIshe(iSym),
     &                             CMO(ipMOi),1,
     &                       0.0d0,Work(ipLiii),1)
                End Do
*                ipMO=ipMO+nIsh(iSym)*nBas(iSym)
                ipLii=ipLii+JNUM*nIshe(iSym)
              End Do
            CALL CWTIME(TCR1,TWR1)
            ttran(1) = ttran(1) + (TCR1 - TCR2)
            ttran(2) = ttran(2) + (TWR1 - TWR2)
*
************************************************************************
**            Integral formation
**            (i i | p q ) = sum_J Lii^J  Lpq^J
*
              ipLii=ipChoT+(nIP-ntotie)*nVec
              Call DGEMM_('N','N',nRS,ntotie,JNUM,
     &                    1.0d0,Work(ipLrs),nRS,
     &                          Work(ipLii),JNUM,
     &                    1.0d0,Work(ipiirs),nRS)
            CALL CWTIME(TCR2,TWR2)
            tform2(1) = tform2(1) + (TCR2 - TCR1)
            tform2(2) = tform2(2) + (TWR2 - TWR1)

            EndIf
*
************************************************************************
**          Read half-transformed active vectors
*
            lChoA=0
            Do i=1,nSym
               k = Muld2h(i,JSYM)
               ipLpq(k) = ipChoT + lChoA       ! Lvb,J
               lChoA= lChoA + nAsh(k)*nBas(i)*JNUM
            End Do
*
            ioff=0
            Do i=1,nSym
              k = Muld2h(i,JSYM)
              lvec=nAsh(k)*nBas(i)*JNUM
              iAdr2=(JVEC-1)*nAsh(k)*nBas(i)+ioff
              call DDAFILE(LuAChoVec(Jsym),2,Work(ipLpq(k)),
     &                     lvec,iAdr2)
            End Do
*
************************************************************************
**          Form (tp|uq) integrals
*
            Do j=1,JNUM
              ioff=0
              Do i=1,nsym
                k = Muld2h(i,JSYM)
                Do it=0,nAshe(k)-1
                  ipLtp=ipLpq(k)+nAshb(k)+it+nAsh(k)*nBas(i)*(j-1)
                  Do iu=0,nAshb(k)+it
                    itu=it*(2*nAshb(k)+it+1)/2+iu
                    ipLuq=ipLpq(k)+iu+nAsh(k)*nBas(i)*(j-1)
                    ipInt=iptpuq+ioff+itu*nBas(i)**2
                    Call DGER(nBas(i),nBas(i),
     &                        1.0d0,Work(ipLtp),nAsh(k),
     &                              Work(ipLuq),nAsh(k),
     &                              Work(ipInt),nBas(i))
                  End Do
                End Do
                ioff=ioff+nAshe(k)*(2*nAshb(k)+nAshe(k)+1)/2*
     &               nBas(i)**2
              End Do
            End Do
            CALL CWTIME(TCR1,TWR1)
            tforma(1) = tforma(1) + (TCR1 - TCR2)
            tforma(2) = tforma(2) + (TWR1 - TWR2)
*
************************************************************************
**          Second MO transformation
*
           If (jsym.eq.1) Then
             ipLtu=ipChoT+(ntp-ntue)*nVec
             ioff=0
             Do i=1,nsym
               Do j=1,JNUM
                 ipLtp=ipLpq(i)+(j-1)*(nBas(i)*nAsh(i))
                 ipMO=1+ioff+nBas(i)*(nIsh(i)+nAshb(i))
                 Do k=0,nAshe(i)-1
                   Call dGeMV_('N',nAshb(i)+k+1,nBas(i),
     &                         1.0d0,Work(ipLtp),nAsh(i),
     &                               CMO(ipMO+k*nBas(i)),1,
     &                         0.0d0,Work(ipLtu),1)
                   ipLtu=ipLtu+(nAshb(i)+k+1)
                 EndDo
               End Do
               ioff=ioff+nBas(i)**2
             End Do
            CALL CWTIME(TCR2,TWR2)
*
************************************************************************
**           Formation of the (tu|rs) integral
*
             ipInt=ipturs
             ipLtu=ipChoT+(ntp-ntue)*nVec
             Do i=1,nsym
               na2=nAshe(i)*nAshb(i)+nAshe(i)*(nAshe(i)+1)/2
               If (na2.gt.0) Then
                 Call DGEMM_('N','T',nRS,na2,JNUM,
     &                       1.0d0,Work(ipLrs),nRS,
     &                             Work(ipLtu),na2,
     &                       1.0d0,Work(ipInt),nRS)
                 ipLtu=ipLtu+na2*JNUM
                 ipInt=ipInt+nRS*na2
               End If
             End Do
            CALL CWTIME(TCR1,TWR1)
            tforma2(1) = tforma2(1) + (TCR1 - TCR2)
            tforma2(2) = tforma2(2) + (TWR1 - TWR2)
           EndIf
************************************************************************
*                                                                      *
*                Cholesky loop is over!                                *
*                                                                      *
************************************************************************

          End Do ! J batch
*
**        Transform to full storage, use Lrs as temp storage
*
          If (jsym.eq.1) Then
            Do i=1,ntotie
              ip1=ipiiab+nab*(i-1)
              ipRS1=ipiirs+nRS*(i-1)
              mode = 'tofull'
              Call play_rassi_sto(irc,iLoc,JSYM,ISTSQ,ISSQ,
     &                                   ip1,ipRS1,mode)
            End Do
            Do i=1,ntue
              ip1=iptupq+npq*(i-1)
              ipRS1=ipturs+nrs*(i-1)
              mode = 'tofull'
*              Call play_rassi_sto(irc,iLoc,JSYM,ISTLT,ISSQ,
              Call play_rassi_sto(irc,iLoc,JSYM,ISTSQ,ISSQ,
     &                                   ip1,ipRS1,mode)
            End Do
          EndIf
*
          Call GetMem('rsL','Free','Real',ipLrs,LREAD)
          Call GetMem('ChoT','Free','Real',ipChoT,max(nIP,ntp)*nVec)
 998      Continue
        End Do ! reduced set JRED
*
        If (jsym.eq.1) Then
          If (ntotie.gt.0)
     &      Call GetMem('iirs','Free','Real',ipiirs,ntotie*maxRS)
          If (ntue.gt.0)
     &      Call GetMem('turs','Free','Real',ipturs,ntue*maxRS)
        EndIf
*
** MGD  Gather integrals from parallel runs
*

*
        CALL CWTIME(TCR1,TWR1)
        Call GetMem('MOt-int','Allo','Real',ipInt,maxpq)
        ip1=ipiaib
        ip3=iptpuq
        isum=0
        isum2=0
        Do isym=1,nSym
          ksym=MulD2h(iSym,jsym)
          nvirt=nBas(kSym)-nIsh(kSym)
          ipMO=1+ISTSQ(kSym)+nBas(kSym)*nIsh(kSym)

          Do ii=1,nIshe(isym)
*
**          MO transform ( i i | p q)
*
            If (jsym.eq.1) Then
              isum=isum+1
              ip2=ipiiab+nab*(isum-1)
              ioff2=0
              Do ksym2=1,nsym
                ipMO2=1+ioff2+nBas(kSym2)*nIsh(kSym2)
                nvirt2=nBas(kSym2)-nIsh(kSym2)
                Do j=1,nvirt2
                  ipMOj=ipMO2+(j-1)*nBas(kSym2)
                  ipIntj=ipInt+(j-1)*nBas(kSym2)
                  Call DSPMV_('U',nBas(kSym2),1.0d0,Work(ip2),
     &                       CMO(ipMOj),1,0.0d0,Work(ipIntj),1)
                End Do
                Call DGEMM_('T','N',nvirt2,nvirt2,nBas(kSym2),
     &                      1.0d0,Work(ipInt),nBas(kSym2),
     &                            CMO(ipMO2),nBas(kSym2),
     &                      0.0d0,Work(ip2)  ,nvirt2)
*
                call DDAFILE(LuChoInt(1),1,Work(ip2),nvirt2**2,iAdr)
*
                ip2=ip2+nBas(kSym2)**2
                ioff2 =ioff2 +nBas(kSym2)**2
              End Do
            Else
              Do i=1,nsym
                iAdr=iAdr+(nBas(i)-nIsh(i))**2
              End Do
            EndIf
*
**      MO transform ( i p | i q)
*
            Call DGEMM_('N','N',nBas(kSym),nvirt,nBas(kSym),
     &                  1.0d0,Work(ip1)  ,nBas(kSym),
     &                        CMO(ipMO),nBas(kSym),
     &                  0.0d0,Work(ipInt),nBas(kSym))
            Call DGEMM_('T','N',nvirt,nvirt,nBas(kSym),
     &                  1.0d0,Work(ipInt),nBas(kSym),
     &                        CMO(ipMO),nBas(kSym),
     &                  0.0d0,Work(ip1)  ,nvirt)
*
            Do i=1,ksym-1
                iAdr=iAdr+(nBas(i)-nIsh(i))**2
            End Do
            call DDAFILE(LuChoInt(1),1,Work(ip1),nvirt**2,iAdr)
            Do i=ksym+1,nsym
                iAdr=iAdr+(nBas(i)-nIsh(i))**2
            End Do
*
            ip1=ip1+nBas(kSym)**2
          End Do
*
**        MO transform ( t u | p q)
*
*compute address
          iAdrtu=0
          Do i=1,ksym-1
            na2=nAshe(i)*nAshb(i)+nAshe(i)*(nAshe(i)+1)/2
            Do j=1,na2
              iAdrtu=iAdrtu+npq
              Do k=1,nsym
                iAdrtu=iAdrtu+nBas(k)**2
              End Do
            End Do
          End Do
*
          ipMO=1+ISTSQ(iSym)
          na2=nAshe(ksym)*nAshb(ksym)+nAshe(ksym)*(nAshe(ksym)+1)/2
          Do itu=1,na2
            If (jsym.eq.1) Then
              isum2=isum2+1
              ip2=iptupq+npq*(isum2-1)
              ioff2=0
              Do ksym2=1,nsym
                ipMO2=1+ioff2
                Do j=1,nBas(ksym2)
                  ipMOj=ipMO2+(j-1)*nBas(kSym2)
                  ipIntj=ipInt+(j-1)*nBas(kSym2)
                  Call DSPMV_('U',nBas(kSym2),1.0d0,Work(ip2),
     &                       CMO(ipMOj),1,0.0d0,Work(ipIntj),1)
                End Do
                If (nBas(ksym2).gt.0) Then
                  Call DGEMM_('T','N',nBas(ksym2),nBas(kSym2),
     &                              nBas(kSym2),1.0d0,
     &                              Work(ipInt),nBas(kSym2),
     &                              CMO(ipMO2),nBas(kSym2),
     &                        0.0d0,Work(ip2)  ,nBas(kSym2))
                  call DDAFILE(LuChoInt(2),1,Work(ip2),nBas(kSym2)**2,
     &            iAdrtu)
*
                  ip2=ip2+nBas(kSym2)**2
                  ioff2 =ioff2 +nBas(kSym2)**2
                EndIf
              EndDo
            Else
              iAdrtu=iAdrtu+npq
            EndIf
*
**        MO transform ( t p | u q)
*
            Call DGEMM_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),
     &                  1.0d0,Work(ip3)  ,nBas(iSym),
     &                        CMO(ipMO),nBas(iSym),
     &                  0.0d0,Work(ipInt),nBas(iSym))
            Call DGEMM_('T','N',nBas(iSym),nBas(iSym),nBas(iSym),
     &                  1.0d0,Work(ipInt),nBas(iSym),
     &                        CMO(ipMO),nBas(iSym),
     &                  0.0d0,Work(ip3) ,nBas(iSym))
            Do i=1,isym-1
                iAdrtu=iAdrtu+nBas(i)**2
            End Do

            call DDAFILE(LuChoInt(2),1,Work(ip3),nBas(iSym)**2,iAdrtu)

            Do i=isym+1,nsym
*               j=MulD2h(i,jsym)
                iAdrtu=iAdrtu+nBas(i)**2
            End Do
            ip3=ip3+nBas(iSym)**2
*
*            Do i=1,nsym
*              iAdrtu=iAdrtu+nBas(i)**2
*            enddo
          End Do

        End Do
        CALL CWTIME(TCR2,TWR2)
        tMO(1) = tMO(1) + (TCR2 - TCR1)
        tMO(2) = tMO(2) + (TWR2 - TWR1)
*
        Call GetMem('MOt-int','Free','Real',ipInt,maxpq)
*
        Do i=1,nsym
          nIshb(i)=nIshb(i)+nIshe(i)  ! now those are done!
          nAshb(i)=nAshb(i)+nAshe(i)  ! now those are done!
        EndDo
        Call GetMem('CMOt','Free','Real',ipCMOt(1),nCMO)
        Call GetMem('iiab','Free','Real',ipiiab,libatch)
        If (ntotae.gt.0)
     &     Call GetMem('tupq','Free','Real',iptupq,labatch)
        If (taskleft) Go to 50  ! loop over i/t batches
*
 999    Continue
*
      End Do ! jsym
*
      CALL CWTIME(TCstart2,TWstart2)
      TOTCPU=TCstart2-TCstart1
      TOTWALL=TWstart2-TWstart1
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky MCLR timing from '//SECNAM
      Write(6,CFmt)'----------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'Integral construction           CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ VECTORS                     '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'TRANSFORMATION                   '
     &                           //'         ',ttran(1),ttran(2)
         Write(6,'(2x,A26,2f10.2)')'(IA|IB) FORMATION                '
     &                           //'         ',tform(1),tform(2)
         Write(6,'(2x,A26,2f10.2)')'(II|AB) FORMATION                '
     &                           //'         ',tform2(1),tform2(2)
         Write(6,'(2x,A26,2f10.2)')'(TP|UQ) FORMATION                '
     &                           //'         ',tforma(1),tforma(2)
         Write(6,'(2x,A26,2f10.2)')'(TU|PQ) FORMATION                '
     &                           //'         ',tforma2(1),tforma2(2)
         Write(6,'(2x,A26,2f10.2)')'MO TRANSFORMATION                '
     &                           //'         ',tMO(1),tMO(2)
         Write(6,*)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif

      end
