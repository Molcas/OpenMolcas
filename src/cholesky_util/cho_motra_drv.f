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

      SUBROUTINE CHO_motra_drv(rc,nIsh,nAsh,nSsh,nFr,nFrVir,
     &                         ipPorb)

**********************************************************************
C
C      a,b,g,d:  AO-index
C      p,q,r,s:  MO-indeces (probably frozen/deleted are excluded)
C
**********************************************************************

      Implicit Real*8 (a-h,o-z)

      Integer   rc,nIsh(*),nAsh(*),nSsh(*), nFr(*), nFrVir(*)

      Real*8    tread(2),tmotr1(2),tmotr2(2)
      Logical   timings,DoRead
      Integer   nPorb(8),ipOrb(8),nPvir(8),nPocc(8)
      Integer   ipLpb(8),iSkip(8),LuLTra(4)
      Integer   kOff1(8),kOff1ij(8),kOff1ia(8),kOff1ai(8),kOff1ab(8)
      Character*6  Fname

      Character*50 CFmt
      Character*13 SECNAM
      Parameter (SECNAM = 'CHO_motra_drv')

      COMMON    /CHOTIME /timings

      parameter (zero = 0.0D0, one = 1.0D0)

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      parameter ( N2 = InfVec_N2 )

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      InfVec(i,j,k) = iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
******
      nDimRS(i,j) = iWork(ip_nDimRS-1+nSym*(j-1)+i)

************************************************************************

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

c --- Define MOs to be used by a given code
c -----------------------------------------
        do i=1,nSym
           nPorb(i) = nIsh(i) + nAsh(i) + nSsh(i) + nFr(i) + nFrVir(i)
           nPocc(i) = nFr(i) + nIsh(i) + nAsh(i)
           nPvir(i) = nSsh(i) + nFrVir(i)
        end do

C ==================================================================

c --- Various offsets & pointers
c ------------------------------
      ipOrb(1)=ipPorb
      DO ISYM=2,NSYM
        NB=NBAS(ISYM-1)
        NP=NPORB(ISYM-1)
        NP2=NB*NP
        ipOrb(ISYM)=ipOrb(ISYM-1)+NP2 !  MO coeff. symm pointers
      END DO

      iLoc = 3 ! use scratch location in reduced index arrays
      nType = 4

C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
      DO jSym=1,nSym
         If (NumCho(jSym).lt.1) GOTO 1000
*        Type refers to occ-occ 2 x vir-occ (with vir or occ as outer
*        index) and vir-vir transformed cholesky vectors. All integrals
*        we need can be constructed from these and the only redundancy
*        we get is frozen orbitals in an inner index.

         Do iType=1,nType
            iSeed = 7 + jSym-1
            LuLTra(iType) = IsFreeUnit(iSeed)
            Write(Fname,'(A4,I1,I1)') 'LTra',jSym,iType
            Call DANAME_MF_WA(LuLTra(iType),Fname)
         End Do
         Call ChoMP2_OpenF(1,1,jSym)
C --- Set up the skipping flags + some initializations --------
C -------------------------------------------------------------
         Do i=1,nSym
            k = Muld2h(i,JSYM)
            iSkip(i) = Min(1,nBas(i)*nBas(k)) ! skip Lik vector
            kOff1(i) = 0
            kOff1ij(i) = 0
            kOff1ia(i) = 0
            kOff1ai(i) = 0
            kOff1ab(i) = 0
            ipLpb(i) = -6666
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
               call abend()
            endif

            Call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
            if(irc.ne.0)then
              Write(6,*)SECNAM//'cho_X_setred non-zero return code.'
              Write(6,*)        'rc= ',irc
              call abend()
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            Call GetMem('MaxM','Max','Real',KDUM,LWORK)

            nVec  = Min(LWORK/(nRS+mvec),nVrs)

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) 'LWORK= ',LWORK
               WRITE(6,*) 'min. mem. need= ',nRS+mTvec
               WRITE(6,*) 'reading ',nRS,' and transforming to ',mvec
               WRITE(6,*) 'of jsym= ',jsym,' and JRED= ',JRED
               rc = 33
               CALL Abend()
               NumBat = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call GetMem('rsL','Allo','Real',ipLrs,LREAD)
            Call GetMem('ChoT','Allo','Real',ipChoT,mvec*nVec)

C --- BATCH over the vectors ----------------------------

            NumBat = (nVrs-1)/nVec + 1
            DO iBatch=1,NumBat
               iAdr = 0
               If (iBatch.eq.NumBat) Then
                  JNUM = nVrs - nVec*(NumBat-1)
               else
                  JNUM = nVec
               endif

               JVEC = nVec*(iBatch-1) + iVrs
               IVEC2 = JVEC - 1 + JNUM
*
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

*              The total number of elements in all cholesky-vectors
*              for this batch

               lChoT=0
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
               kOff2 = 0
               kOff2ij = 0
               kOff2ia = 0
               kOff2ai = 0
               kOff2ab = 0
               Do iSymb=1,nSym

                  iSymp = MulD2h(JSYM,iSymb)
                  NAp = nPorb(iSymp)
                  NAq = nPorb(iSymb) ! iSymb=iSymq

                  CALL CWTIME(TCM3,TWM3)

                  If(NAp*NAq.ne.0)Then
                     Do JVC=1,JNUM
                        ipLJpb = ipLpb(iSymp)
     &                         + nPorb(iSymp)*nBas(iSymb)*(JVC-1)
                        ipLJpq = ipLpq
     &                         + nPorb(iSymp)*nPorb(iSymb)*(JVC-1)
*


*                        ipLJai = ipLJia + nPocc(iSymP)*nPvir(iSymB)
*                        ipLJab = ipLJai + nPvir(iSymP)*nPocc(iSymB)
*
                        CALL DGEMM_('N','T',NAp,NAq,nBas(iSymb),
     &                             One,Work(ipLJpb),NAp,
     &                             Work(ipOrb(iSymb)),NAq,
     &                             Zero,Work(ipLJpq),NAp)
                        lChoMO = NAp*NAq
                        iAdr = 1 + kOff2 + kOff1(iSymB)
                        Call ddaFile(lUnit_F(jSym,1),1,
     &                               Work(ipLJpq),lChoMO,iAdr)
                        kOff1(iSymB) = kOff1(iSymB) + nPQ_prod(jSym)

*                        Write(6,*) 'Long ChoVec', JVC
*                        Do i = 0, lChoMO-1
*                           Write(6,*) Work(ipLJpq+i)
*                        End Do
*
                        Do iI = 1, nPocc(iSymP)
                           Do iJ = 1,nPocc(iSymB)
                              Work(ip_TmpL+(iJ-1)+(iI-1)*nPocc(iSymB)) =
     &                         Work(ipLJpq+(iJ-1)+(iI-1)*nPorb(iSymB))
                           End Do
                        End Do
                        lChoMOij = nPocc(iSymP)*nPocc(iSymB)
                        iType = 1
                        iAdrij = 1 + kOff2ij + kOff1ij(iSymB)
                        Call dDaFile(LuLtra(iType),1,
     &                               Work(ip_TmpL),
     &                               lChoMOij, iAdrij)

*                        Write(6,*) 'ChoVec nr', JVC
*                        Do i = 0, lChoMOij-1
*                          Write(6,*) Work(ip_TmpL+i)
*                        End Do

*                       Put L-vectors ia on disk
                        Do iI = 1, nPocc(iSymP)
                           Do iA = 1,nPvir(iSymB)
                              Work(ip_TmpL+(iA-1)+(iI-1)*nPvir(iSymB)) =
     &                         Work(ipLJpq+(iA+nPocc(iSymB)-1)+
     &                              (iI-1)*nPorb(iSymB))
                           End Do
                        End Do
                        iType = 2
                        lChoMOia = nPocc(iSymP)*nPvir(iSymB)
                        iAdria = 1 + kOff2ia + kOff1ia(iSymB)
                        Call dDaFile(LuLtra(iType),1,Work(ip_TmpL),
     &                               lChoMOia, iAdria)

*                        Write(6,*) 'ChoVec nr', JVC
*                        Do i = 0, lChoMOia-1
*                           Write(6,*) Work(ip_TmpL+i)
*                        End Do
*                       Put L-vectors ai on disk.
                        Do iA = 1, nPvir(iSymP)
                           Do iI = 1,nPocc(iSymB)
                              Work(ip_TmpL+(iI-1)+(iA-1)*nPocc(iSymB)) =
     &                         Work(ipLJpq+(iI-1)+
     &                              (iA+nPvir(iSymP)-1)*nPorb(iSymB))
                           End Do
                        End Do
                        iType = 3
                        lChoMOai = nPvir(iSymP)*nPocc(iSymB)
                        iAdrai = 1 + kOff2ai + kOff1ai(iSymB)
                        Call dDaFile(LuLtra(iType),1,Work(ip_TmpL),
     &                               lChoMOai, iAdrai)

*                        Write(6,*) 'ChoVec nr', JVC
*                        Do i = 0, lChoMOai-1
*                           Write(6,*) Work(ip_TmpL+i)
*                        End Do

*
*                       Put L-vectors ab on disk.
                        Do iA = 1, nPvir(iSymP)
                           Do iB = 1,nPvir(iSymB)
                              Work(ip_TmpL+(iB-1)+(iA-1)*nPvir(iSymB)) =
     &                         Work(ipLJpq+(iB+nPocc(iSymB)-1)+
     &                              (iA+nPvir(iSymP)-1)*nPorb(iSymB))
                           End Do
                        End Do
                        iType = 4
                        lChoMOab = nPvir(iSymP)*nPvir(iSymB)
                        iAdrab = 1 + kOff2ab + kOff1ab(iSymB)
                        Call dDaFile(LuLtra(iType),1,Work(ip_TmpL),
     &                                  lChoMOab, iAdrab)
*                        Write(6,*) 'ChoVec nr', JVC
*                        Do i = 0, lChoMOab-1
*                           Write(6,*) Work(ip_TmpL+i)
*                        End Do

                        kOff1ab(iSymB) = kOff1ab(iSymB)+
     &                                      nPQ_prodab(jSym)
                        kOff1ij(iSymB) = kOff1ij(iSymB)+
     &                                      nPQ_prodij(jSym)
                        kOff1ia(iSymB) = kOff1ia(iSymB)+
     &                                      nPQ_prodia(jSym)
                        kOff1ai(iSymB) = kOff1ai(iSymB)+
     &                                      nPQ_prodia(jSym)
*
                     End Do
                  EndIf
*
                  CALL CWTIME(TCM4,TWM4)
                  tmotr2(1) = tmotr2(1) + (TCM4 - TCM3)
                  tmotr2(2) = tmotr2(2) + (TWM4 - TWM3)

                  kOff2 = NAp*NAq + kOff2
                  kOff2ij = lChoMOij + kOff2ij
                  kOff2ia = lChoMOia + kOff2ia
                  kOff2ai = lChoMOai + kOff2ai
                  kOff2ab = lChoMOab + kOff2ab
               End Do

C --------------------------------------------------------------------
C --------------------------------------------------------------------

            END DO  ! end batch loop


C --- free memory
            Call GetMem('ChoT','Free','Real',ipChoT,mTvec*nVec)
            Call GetMem('rsL','Free','Real',ipLrs,LREAD)

999         CONTINUE


         END DO   ! loop over red sets

         Call ChoMP2_OpenF(2,1,jSym)
         Do i = 1, nType
            Call DaClos(LuLTra(i))
         End Do

1000     CONTINUE

      END DO   !loop over JSYM

      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1
*
*---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky timing from '//SECNAM
      Write(6,CFmt)'----------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'MO transf. Cholesky vectors     CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ/WRITE VECTORS               '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'1st half-transf.                 '
     &                           //'         ',tmotr1(1),tmotr1(2)
         Write(6,'(2x,A26,2f10.2)')'2nd half-transf.                 '
     &                           //'         ',tmotr2(1),tmotr2(2)
         Write(6,*)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif


      rc  = 0

      Call Cho_X_final(rc)

      Return
      END

**************************************************************
**************************************************************
