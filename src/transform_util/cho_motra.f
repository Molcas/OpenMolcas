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
        Subroutine Cho_MOtra(CMO,nCMOs,Do_int,ihdf5)
        Implicit None
        Integer nCMOs, ihdf5
        Real*8  CMO(nCMOs)
        Logical Do_int

        Character*6 BName

        Integer iSym
        Integer nSym
        Integer nBas(8), nOrb(8)
        Integer nFro(8), nIsh(8), nAsh(8)
        Integer nSsh(8), nDel(8)

        Logical InitChoEnv

        Call Get_iScalar('nSym',nSym)
        Call Get_iArray('nBas',nBas,nSym)
        Call Get_iArray('nOrb',nOrb,nSym)
        Call Get_iArray('nFro',nFro,nSym)
        Call Get_iArray('nIsh',nIsh,nSym)
        Call Get_iArray('nAsh',nAsh,nSym)
        Call Get_iArray('nDel',nDel,nSym)
        Do iSym=1,nSym
           nSsh(iSym)=nBas(iSym)-nDel(iSym)-nAsh(iSym)
     &               -nIsh(iSym)-nFro(iSym)
        End Do
        BName='_CHMOT'
        InitChoEnv=.true.
        Call Cho_MOTra_(CMO,nCMOs,nSym,nBas,nOrb,
     &                  nFro,nIsh,nAsh,nSsh,nDel,
     &                  BName,Do_Int,ihdf5,InitChoEnv)

        End
************************************************************************
        Subroutine Cho_MOTra_(CMO,nCMOs,nSym,nBas,nOrb,
     &                        nFro,nIsh,nAsh,nSsh,nDel,
     &                        BName,Do_int,ihdf5,Do_ChoInit)
C
C       Note: frozen and deleted orbitals are not included in the
C             transformation.
C
        Implicit Real*8 (a-h,o-z)
        Integer nCMOs, ihdf5
        real*8  CMO(nCMOs)
        Integer nSym
        Integer nBas(nSym), nOrb(nSym)
        Integer nFro(nSym), nIsh(nSym), nAsh(nSym)
        Integer nSsh(nSym), nDel(nSym)
        Character*6 BName
        Logical Do_int
        Logical Do_ChoInit

        Logical timings
        COMMON /CHOTIME /timings
#include "WrkSpc.fh"
**************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
**************************************************

        n=nBas(1)**2
        Do iSym=2,nSym
           n=n+nBas(iSym)**2
        End Do
        If (n.ne.nCMOs) Then
           ! The dimension of CMO is assumed to be nBas**2 in Transp_MOs
           ! with each symmetry block starting at
           !  sym=1: 1
           !  sym=2: 1+nBas(1)**2
           !  sym=3: 1+nBas(1)**2+nBas(2)**2
           ! etc.
           ! This differs from, e.g., subroutine PriMO where each
           ! symmetry block starts at
           !  sym=1: 1
           !  sym=2: 1+nBas(1)*nOrb(1)
           !  sym=3: 1+nBas(1)*nOrb(1)+nBas(2)*nOrb(2)
           ! etc.
           ! This is only a potential problem if orbitals were removed
           ! due to linear dependence (not deleted virtual orbitals)
           Call WarningMessage(2,'Cho_MOTra_: n != nCMOs')
           Write(6,*) 'n,nCMOs=',n,nCMOs
           Call Abend()
        End If
        Call GetMem('CHMOs','Allo','Real',ipCMO,nCMOs)
        Call Transp_MOs(CMO,Work(ipCMO),nSym,nFro,nIsh,nAsh,nSsh,nBas)
c
        timings=.True.
c
        If (Do_int) Then
          Lu_Xint = 80
          Lu_Xint = isfreeunit(Lu_Xint)
          call DaName_mf_wa(Lu_Xint,'DIAGINT')
          lXint=0
          Do jSym=1,nSym
             Do iSymq=1,nSym
                nOrbq=nIsh(iSymq)+nAsh(iSymq)+nSsh(iSymq)
                iSymp=MulD2h(iSymq,jSym)
                If (iSymp.eq.iSymq) Then
                   lXint=lXint+nOrbq*(nOrbq+1)/2
                ElseIf (iSymp.lt.iSymq) Then
                   nOrbp=nIsh(iSymp)+nAsh(iSymp)+nSsh(iSymp)
                   lXint=lXint+nOrbp*nOrbq
                EndIf
             End Do
          End Do
          Call GetMem('DIAGON','Allo','Real',ipXint,lXint)
        EndIf
c
        If (Do_ChoInit) Then
           FracMem=0.0d0 ! in a parallel run set it to a sensible value
           irc=0
           Call Cho_X_Init(irc,FracMem) ! initialize cholesky info
           If (irc.ne.0) Then
              Call WarningMessage(2,
     &                        'Cho_MOTra_: non-zero rc from Cho_X_Init')
              Write(6,*) 'rc=',irc
              Call Abend()
           End If
        End If
        call CHO_TR_drv(irc,nIsh,nAsh,nSsh,ipCMO,BName,
     &                      Do_int,ihdf5,Work(ipXint),lXint)
        If (Do_ChoInit) Then
           Call Cho_X_final(irc)
           If (irc.ne.0) Then
              Call WarningMessage(2,
     &                       'Cho_MOTra_: non-zero rc from Cho_X_Final')
              Write(6,*) 'rc=',irc
              Call Abend()
           End If
        End If
c
        If (Do_int) Then
           Call GADSum(Work(ipXint),lXint)
           kdisk=0
           Call ddafile(Lu_Xint,1,Work(ipXint),lXint,kdisk)
           Call daclos(Lu_Xint)
           Call GetMem('DIAGON','Free','Real',ipXint,lXint)
        EndIf
        Call GetMem('CHMOs','Free','Real',ipCMO,nCMOs)
c
        return
c Avoid unused argument warnings
        If (.False.) Then
           Call Unused_integer_array(nOrb)
           Call Unused_integer_array(nDel)
        End If
        end
************************************************************************
      SUBROUTINE CHO_TR_drv(rc,nIsh,nAsh,nSsh,ipPorb,BName,Do_int,ihdf5,
     &                      Xint,lXint)
**********************************************************************
C
C      a,b,g,d:  AO-index
C      p,q,r,s:  MO-indeces belonging to all (fro and del excluded)
C
**********************************************************************
#ifdef _HDF5_QCM_
      use hdf5_utils
#endif
      Implicit Real*8 (a-h,o-z)

      Integer   rc,nIsh(*),nAsh(*),nSsh(*),lXint, ihdf5
      Real*8    Xint(0:lXint-1)
      Character*6 BName

      Real*8    tread(2),tmotr1(2),tmotr2(2)
      Logical   Debug,timings,DoRead,Do_int
      Integer   nPorb(8),ipOrb(8)
      Integer   ipLpb(8),iSkip(8)
      Integer   LunChVF(8),kOff(8),iOffB(8),nOB(8)

      Character*7  Fnam
      Character*50 CFmt
      Character*10 SECNAM
      Parameter (SECNAM = 'CHO_TR_drv')

      COMMON    /CHOTIME /timings
      Character*3  tv2disk
      COMMON    /CHOTRAW /tv2disk

      parameter (zero = 0.0D0, one = 1.0D0)

#ifdef _HDF5_QCM_
      integer(HID_T)  :: choset_id
      integer(HID_T)  :: space_id
#endif

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      parameter ( N2 = InfVec_N2 )

      integer isfreeunit

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      InfVec(i,j,k) = iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
******
      nDimRS(i,j) = iWork(ip_nDimRS-1+nSym*(j-1)+i)
************************************************************************


#ifdef _DEBUG_
      Debug=.true.
#else
      Debug=.false.
#endif

      Call QEnter(SECNAM)

#ifdef _HDF5_QCM_
      ! Leon 13.6.2017: Avoid opening a regular file if HDF5 is used
      if (ihdf5/=1) then
#endif
        Do i=1,nSym
          LunChVF(i) = 80
          LunChVF(i) = isfreeunit(LunChVF(i))
          write(Fnam,'(A6,I1)') BName,i
          call DaName_mf_wa(LunChVF(i),Fnam)
        End Do
#ifdef _HDF5_QCM_
      End If
#endif

      DoRead  = .false.
      IREDC = -1  ! unknown reduced set in core

      iSwap = 0  ! Lpb,J are returned by cho_x_getVtra
      kMOs = 1
      nMOs = 1


      CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

      do i=1,2            ! 1 --> CPU   2 --> Wall
         tread(i) = zero   !time read/write vectors
         tmotr1(i) = zero  !time 1st MO half-transf.
         tmotr2(i) = zero  !time 2nd MO half-transf.
      end do

      If (Do_int) Call Fzero(Xint(0),lXint)

c --- Define MO space used
c -----------------------------------
      do i=1,nSym
         nPorb(i) = nIsh(i) + nAsh(i) + nSsh(i)
      end do


C ==================================================================

c --- Various offsets & pointers
c ------------------------------
      ipOrb(1)=ipPorb
      DO ISYM=2,NSYM
        NB=NBAS(ISYM-1)
        NP=NPORB(ISYM-1)
        NP2=NP*NB
        ipOrb(ISYM)=ipOrb(ISYM-1)+NP2 !  MO coeff. symm pointers
      END DO

      iLoc = 3 ! use scratch location in reduced index arrays

C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
c
c
      Mpq=0

      DO jSym=1,nSym

! Init the HDF5 section
#ifdef _HDF5_QCM_
         if(ihdf5 == 1)then
           ! max size of the cholesky vector for now is the
           ! max(nPorb)*(max(nPorb)+1)/2
           ! probably this can be chosen more efficiently,
           ! but would matter only if we use symmetry
           call hdf5_init_wr_cholesky(file_id(1), JSym,
     &       maxval(nPorb(1:nSym))*(maxval(nPorb(1:nSym))+1)/2,
     &       NumCho(JSym), choset_id, space_id)
         end if
#endif
         If (NumCho(jSym).lt.1) GOTO 1000

C --- Set up the skipping flags + some initializations --------
C -------------------------------------------------------------
         Do i=1,nSym
            ipLpb(i) = -6666
            k=Muld2h(i,JSYM)
            iSkip(i) = Min(1,nPorb(i)*nPorb(k)) ! skip Lik vector
            If (i.lt.k) Then
               kOff(i)=Mpq
               nOB(i)=nPorb(i)*nPorb(k)*NumCho(jSym)
               Mpq=Mpq+nPorb(i)*nPorb(k)
            ElseIf (k.eq.i) Then
               kOff(i)=Mpq
               nOB(i)=nPorb(i)*(nPorb(i)+1)/2*NumCho(jSym)
               Mpq=Mpq+nPorb(i)*(nPorb(i)+1)/2
            Else
               nOB(i)=0
            EndIf
            iOffB(i)=0
         End Do
*
         Do i=2,nSym
            iOffB(i)=iOffB(i-1)+nOB(i-1)
         End Do
         Do i=1,nSym
            If (nOB(i).eq.0) Then
               k=Muld2h(i,JSYM)
               iOffB(i)=iOffB(k)
               kOff(i)=kOff(k)
            EndIf
         End Do
C -------------------------------------------------------------


C ****************     MEMORY MANAGEMENT SECTION    *****************
C ------------------------------------------------------------------
C --- compute memory needed to store at least 1 vector of JSYM
C --- and do all the subsequent calculations
C ------------------------------------------------------------------
         mTvec  = 0  ! mem for storing half-transformed vec Laq,J
         mTTvec = 0  ! mem for storing transformed vec Lpq,J

         do l=1,nSym
            k=Muld2h(l,JSYM)
            mTvec = mTvec + nPorb(k)*nBas(l)
            mTTvec = Max(mTTvec,nPorb(k)*nPorb(l))
         end do

         mvec = mTvec + mTTvec

C ------------------------------------------------------------------
C ------------------------------------------------------------------

         JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
         JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec

         Do JRED=JRED1,JRED2

            CALL Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

            If (nVrs.eq.0) GOTO 999  ! no vectors in that (jred,jsym)

            if (nVrs.lt.0) then
               Write(6,*)SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
               call qtrace()
               call abend()
            endif

            Call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
            if(irc.ne.0)then
              Write(6,*)SECNAM//'cho_X_setred non-zero return code.',
     &                         ' rc= ',irc
              call qtrace()
              call abend()
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            Call GetMem('MaxM','Max','Real',KDUM,LWORK)

            nVec  = Min(LWORK/(nRS+mvec+1),nVrs)

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) 'LWORK= ',LWORK
               WRITE(6,*) 'Min. mem. need= ',nRS+mvec+1
               WRITE(6,*) 'Reading ',nRS,' and then MO-transform.'
               WRITE(6,*) 'In jsym= ',jsym,' and JRED= ',JRED
               rc = 33
               CALL QTrace()
               CALL Abend()
               nBatch = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call GetMem('rsL','Allo','Real',ipLrs,LREAD)
            Call GetMem('ChoT','Allo','Real',ipChoT,(mvec+1)*nVec)

C --- BATCH over the vectors ----------------------------

            nBatch = (nVrs-1)/nVec + 1

            DO iBatch=1,nBatch

               If (iBatch.eq.nBatch) Then
                  JNUM = nVrs - nVec*(nBatch-1)
               else
                  JNUM = nVec
               endif

               JVEC = nVec*(iBatch-1) + iVrs
               IVEC2 = JVEC - 1 + JNUM

               CALL CWTIME(TCR1,TWR1)

               CALL CHO_VECRD(Work(ipLrs),LREAD,JVEC,IVEC2,JSYM,
     &                        NUMV,IREDC,MUSED)

               If (NUMV.le.0 .or.NUMV.ne.JNUM) then
                  rc=77
                  RETURN
               End If

               CALL CWTIME(TCR2,TWR2)
               tread(1) = tread(1) + (TCR2 - TCR1)
               tread(2) = tread(2) + (TWR2 - TWR1)

               lChoT=JNUM
               Do iSymp=1,nSym

                  iSymb = MulD2h(jSym,iSymp)

                  ipLpb(iSymp) = ipChoT + lChoT
                  lChoT = lChoT + nPorb(iSymp)*nBas(iSymb)*JNUM

               End Do

               ipLpq = ipChoT + lChot

C --------------------------------------------------------------------
C --- First half MO transformation  Lpb,J = sum_a  C(p,a) * Lab,J
C --------------------------------------------------------------------

               CALL CWTIME(TCM1,TWM1)

               CALL CHO_X_getVtra(irc,Work(ipLrs),LREAD,jVEC,JNUM,
     &                           JSYM,iSwap,IREDC,nMOs,kMOs,ipOrb,nPorb,
     &                           ipLpb,iSkip,DoRead)

               if (irc.ne.0) then
                  rc = irc
                  RETURN
               endif

               CALL CWTIME(TCM2,TWM2)
               tmotr1(1) = tmotr1(1) + (TCM2 - TCM1)
               tmotr1(2) = tmotr1(2) + (TWM2 - TWM1)

C --------------------------------------------------------------------
C --- 2nd half of MO transformation  Lpq,J = sum_b  Lpb,J * C(q,b)
C --------------------------------------------------------------------

               IF (JSYM.eq.1) THEN   !  Lpq,J in LT-storage

                  Do iSymb=1,nSym

                     NAp = nPorb(iSymb)
                     NApq=NAp*(NAp+1)/2

                     CALL CWTIME(TCM3,TWM3)

                     If (NApq.ne.0) Then

                      Do JVC=1,JNUM

                       ipLJpb = ipLpb(iSymb)
     &                        + nPorb(iSymb)*nBas(iSymb)*(JVC-1)
                       ipLJpq = ipLpq
     &                        + NApq*(JVC-1)

                       CALL DGEMM_Tri('N','T',NAp,NAp,nBas(iSymb),
     &                            One,Work(ipLJpb),NAp,
     &                                Work(ipOrb(iSymb)),NAp,
     &                           Zero,Work(ipLJpq),NAp)

                      End Do

                     EndIf

                     CALL CWTIME(TCM4,TWM4)
                     tmotr2(1) = tmotr2(1) + (TCM4 - TCM3)
                     tmotr2(2) = tmotr2(2) + (TWM4 - TWM3)

                     CALL CWTIME(TCR3,TWR3)

                     If (NApq.ne.0) Then

                        If (tv2disk.eq.'PQK') Then
#ifdef _HDF5_QCM_
                         if (ihdf5/=1) then
#endif
                           Call ddafile(LunChVF(jSym),1,Work(ipLpq),
     &                                                NApq*JNUM,
     &                                                iOffB(iSymb))
#ifdef _HDF5_QCM_
                         else
                           ! this should never happen, this case should be caught in motra.f
                           Write(6,*)' Writing of Cholesky vectors'//
     &                       'in HDF5 format as (pq,k) is not'//
     &                       'supported.'
                           call Abend()
                         end if
#endif
                           If (Do_int) Then
                              Do ipq=0,NApq-1
                                 kt=kOff(iSymb)+ipq
                                 kLpq=ipLpq+ipq
                                Xint(kt)=Xint(kt)+ddot_(JNUM,Work(kLpq),
     &                                             NApq,Work(kLpq),NApq)
                              End Do
                           EndIf
                        Else
                           Do ipq=0,NApq-1
                              call dcopy_(JNUM,Work(ipLpq+ipq),NApq,
     &                                        Work(ipChoT),1)
                              If (Do_int) Then
                                 kt=kOff(iSymb)+ipq
                                 Xint(kt)=Xint(kt)+ddot_(JNUM,
     &                                                  Work(ipChoT),1,
     &                                                  Work(ipChoT),1)
                              EndIf
                              idisk=iOffB(iSymb)+NumCho(jSym)*ipq
#ifdef _HDF5_QCM_
                             ! Leon 13.6.2017: Do not write Cholesky vectors to the regular file
                             ! if the hdf5 file is written. It becomes counterproductive to write
                             ! the same content twice for large basis sets
                             if (ihdf5/=1) then
#endif
                              Call ddafile(LunChVF(jSym),1,Work(ipChoT),
     &                                                   JNUM,
     &                                                   idisk)

#ifdef _HDF5_QCM_
                             ! Write the transformed Cholesky batch to the hdf5 dataset
                             ! The ordering in HDF5 is in column-major order, corresponding to
                             ! the 'Kpq' storage
                             ! This way all the elements needed to compute one integral can be
                             ! read with one read operation.

                             ! TODO: eventually row-major order storage + chunked dataset might
                             ! improve the performance -- but probably it's irrelevant.
                             ! Leon 22.4.2016 -- modified the write_cholesky call below to
                             ! account for multiple reduced sets
                             else
                               call hdf5_write_cholesky(choset_id,
     &                           space_id,ipq,nVec*(iBatch-1)+iVrs-1,
     &                                             JNUM,Work(ipChoT))
                             end if
#endif

                           End Do
                           iOffB(iSymb)=iOffB(iSymb)+JNUM
                        EndIf

                     EndIf

                     CALL CWTIME(TCR4,TWR4)
                     tread(1) = tread(1) + (TCR4 - TCR3)
                     tread(2) = tread(2) + (TWR4 - TWR3)

                  End Do

               ELSE

                  Do iSymb=1,nSym

                     iSymp = MulD2h(JSYM,iSymb)
                     NAp = nPorb(iSymp)
                     NAq = nPorb(iSymb) ! iSymb=iSymq
                     NApq=NAp*NAq

                     CALL CWTIME(TCM3,TWM3)

                     If(NApq.ne.0.and.iSymp.lt.iSymb)Then

                      Do JVC=1,JNUM

                       ipLJpb = ipLpb(iSymp)
     &                        + nPorb(iSymp)*nBas(iSymb)*(JVC-1)
                       ipLJpq = ipLpq
     &                        + nPorb(iSymp)*nPorb(iSymb)*(JVC-1)

                       CALL DGEMM_('N','T',NAp,NAq,nBas(iSymb),
     &                            One,Work(ipLJpb),NAp,
     &                                Work(ipOrb(iSymb)),NAq,
     &                           Zero,Work(ipLJpq),NAp)

                      End Do

                     EndIf

                     CALL CWTIME(TCM4,TWM4)
                     tmotr2(1) = tmotr2(1) + (TCM4 - TCM3)
                     tmotr2(2) = tmotr2(2) + (TWM4 - TWM3)

                     CALL CWTIME(TCR3,TWR3)

                     If (NApq.ne.0) Then

                        If (tv2disk.eq.'PQK') Then
                           Call ddafile(LunChVF(jSym),1,Work(ipLpq),
     &                                                NApq*JNUM,
     &                                                iOffB(iSymp))
                           If (Do_int) Then
                              Do ipq=0,NApq-1
                                 kt=kOff(iSymp)+ipq
                                 kLpq=ipLpq+ipq
                                Xint(kt)=Xint(kt)+ddot_(JNUM,Work(kLpq),
     &                                             NApq,Work(kLpq),NApq)
                              End Do
                           EndIf

                        Else
                           Do ipq=0,NApq-1
                              call dcopy_(JNUM,Work(ipLpq+ipq),NApq,
     &                                        Work(ipChoT),1)
                              If (Do_int) Then
                                 kt=kOff(iSymp)+ipq
                                 Xint(kt)=Xint(kt)+ddot_(JNUM,
     &                                                  Work(ipChoT),1,
     &                                                  Work(ipChoT),1)
                              EndIf
                              idisk=iOffB(iSymp)+NumCho(jSym)*ipq
#ifdef _HDF5_QCM_
                             ! Write the transformed Cholesky batch to the hdf5 dataset
                             ! The ordering in HDF5 is in column-major order, corresponding to
                             ! the 'Kpq' storage
                             ! See above for more explanation
                             if (ihdf5==1) then
                               call hdf5_write_cholesky(choset_id,
     &                           space_id,ipq,nVec*(iBatch-1),JNUM,
     &                                                    Work(ipChoT))
                             end if
#endif

                              Call ddafile(LunChVF(jSym),1,Work(ipChoT),
     &                                                   JNUM,
     &                                                   idisk)
                           End Do
                           iOffB(iSymp)=iOffB(iSymp)+JNUM
                        EndIf

                     EndIf

                     CALL CWTIME(TCR4,TWR4)
                     tread(1) = tread(1) + (TCR4 - TCR3)
                     tread(2) = tread(2) + (TWR4 - TWR3)

                  End Do

               EndIf

C --------------------------------------------------------------------
C --------------------------------------------------------------------

            END DO  ! end batch loop

C --- free memory
            Call GetMem('ChoT','Free','Real',ipChoT,(mvec+1)*nVec)
            Call GetMem('rsL','Free','Real',ipLrs,LREAD)

999         CONTINUE

         END DO   ! loop over red sets

#ifdef _HDF5_QCM_
         ! close Cholesky HDF5 stuff
         if (ihdf5==1) then
           call hdf5_close_cholesky(choset_id, space_id)
         else
#endif
           call daclos(LunChVF(jSym))
#ifdef _HDF5_QCM_
         endif
#endif
1000     CONTINUE

      END DO   !loop over JSYM


      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1
*
*---- Write out timing information
      if(timings)then

      CFmt='(6x,A)'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'Cholesky-MOTRA timings            CPU       WALL '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

      If (Do_int) Then
         Write(6,'(6x,A28,2f10.2)')'I/O vectors + diag ERIs step     '
     &                           //'         ',tread(1),tread(2)
      Else
         Write(6,'(6x,A28,2f10.2)')'I/O vectors                      '
     &                           //'         ',tread(1),tread(2)
      EndIf
         Write(6,'(6x,A28,2f10.2)')'1st half-transf.                 '
     &                           //'         ',tmotr1(1),tmotr1(2)
         Write(6,'(6x,A28,2f10.2)')'2nd half-transf.                 '
     &                           //'         ',tmotr2(1),tmotr2(2)
         Write(6,*)
         Write(6,'(6x,A28,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif

      rc  = 0

      write(6,*)
      If (tv2disk.eq.'PQK') Then
         tv2disk(1:2)='pq'
         write(6,*)'     Transformed Cholesky vectors stored as L(',
     &                                tv2disk(1:2),',',tv2disk(3:3),')'
      Else
         tv2disk(2:3)='pq'
         write(6,*)'     Transformed Cholesky vectors stored as L(',
     &                                tv2disk(1:1),',',tv2disk(2:3),')'
      EndIf
      write(6,*)

      CAll QExit(SECNAM)

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(ihdf5)
      END
