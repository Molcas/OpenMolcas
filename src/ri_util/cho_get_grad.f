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
* Copyright (C) 2007, Francesco Aquilante                              *
*               2011, Thomas Bondo Pedersen                            *
************************************************************************
      SUBROUTINE CHO_GET_GRAD(irc,nDen,
     &                        ipDLT,ipDLT2,ipMSQ,
     &                        Txy,nTxy,ipTxy,DoExchange,lSA,
     &                        nChOrb_,ipAorb,nAorb,DoCAS,
     &                        Estimate,Update,
     &                        V_k,U_k,Z_p_k,nnP,npos,nZpk)

************************************************************************
*  Author : F. Aquilante (visiting F. Illas group in Barcelona, Spain, *
*                                                    March-April 2007) *
*                                                                      *
*  Purpose:                                                            *
*         Computation of the relevant quantities for RI                *
*         (and Cholesky) gradient code                                 *
*                                                                      *
*         Coulomb term :  V_k = sum_gd  D_gd L_gd_k                    *
*                                                                      *
*         MP2 Coulomb term : U_k = sum_gd D(MP2)_gd L_gd_k             *
*                                                                      *
*         Active term :   Z_p_k = sum_xy T(xy,p) L_xy_k                *
*                                                                      *
*         Inact. Exchange term: the quantity returned on disk is       *
*                                                                      *
*                     L_ij_k = sum_gd L_gd_k C_gi C_dj                 *
*                                                                      *
*                                                                      *
*  Input:                                                              *
*                                                                      *
*         nDen : is equal to 2 iff Spin Unrestricted                   *
*                4 for SA-CASSCF, otherwise nDen=1                     *
*                                                                      *
*         ipDLT : pointer to the LT-packed and symm. blocked           *
*                 one-body Dmat.                                       *
*                 For spin unrestricted, Dmat = Dalpha + Dbeta         *
*                                                                      *
*         ipDLT2: pointer to the LT-packed and symm. blocked           *
*                 one body MP2 Dmat.                                   *
*                                                                      *
*         ipMSQ : Cholesky MOs stored as C(a,i) symm. blocked          *
*                 with blocks of dim. (nBas,nBas). These are           *
*                 obtained from CD of the 1-particle DMAT.             *
*                 (Two pointers iff alpha and beta spinorbitals)       *
*                                                                      *
*         ipTxy : array (8,8) of pointers to the symm. blocks          *
*                 of the Cholesky decomposed MO-basis (symmetrized)    *
*                 2-body density matrix                                *
*                 T(xy,p) : is stored by compound symmetry JSYM        *
*                            the indices {xy} are stored as PACKED     *
*                            (sym x.le.sym y)                          *
*                                                                      *
*         DoExchange : logical to activate exchange grad. components   *
*                                                                      *
*         nChOrb_ : array of nr. of Cholesky orbitals in each irrep    *
*                                                                      *
*         ipAorb : array of pointers to the active orbitals of each    *
*                  irrep. The storage MUST BE of type C(v,b)           *
*                                                                      *
*         nAorb : array with # of Active orbitals in each irrep        *
*                 (The same orbital basis                              *
*                 in which the 2-body Dmat is expressed)               *
*                                                                      *
*         DoCAS : logical to activate CASSCF grad. components          *
*                                                                      *
*         nScreen : See e.g. LK-screening docum. in SCF                *
*                   or CASSCF read-input routines. Default = 10        *
*                                                                      *
*         dmpK : damping for the LK-screening threshold. Def: 1.0d0    *
*                                                                      *
*         Estimate : logical for LK-screening. Default: .false.        *
*                                                                      *
*         Update : logical for LK-screening. Default: .true.           *
*                                                                      *
*         nnP : array of # of Cholesky vectors for the dec 2-body      *
*               density matrix in each compound symmetry               *
*                                                                      *
*                                                                      *
*  Output:                                                             *
*         irc : return code                                            *
*                                                                      *
*         V_k : array Real*8 for the Coulomb interm. Size=NumCho(1)    *
*                                                                      *
*         U_k : array Real*8 for the mp2 Coulomb interm. Size=NumCho(1)*
*                                                                      *
*         Z_p_k : array Real*8 for the active grad. components.        *
*                  Must be zeroed by the calling routine. Stored       *
*                  according to jSym and blocked after symm. blocks    *
*                  of the active orbitals (square storage).            *
*                                                                      *
*  Modifications:                                                      *
*    August 24, 2011, Thomas Bondo Pedersen:                           *
*       Allow zero vectors on a node.                                  *
*                                                                      *
************************************************************************
      use ChoArr, only: nBasSh, nDimRS
#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: Is_Real_Par
#endif
      Implicit Real*8 (a-h,o-z)

      Logical   timings,DoRead,DoExchange,DoCAS,lSA
      Logical   DoScreen,Estimate,Update,BatchWarn
      Integer   nDen,nChOrb_(8,5),nAorb(8),nnP(8),nIt(5)
      Integer   ipMSQ(nDen),ipAorb(8,*),ipTxy(8,8,2)
      Integer   kOff(8,5), LuRVec(8,3), ipLpq(8), ipLxy(8), iSkip(8)
      Integer   ipDrs(5), ipY, ipYQ, ipML, ipSKsh(5)
      Integer   ipDrs2,ipDLT(5),ipDLT2
      Integer   ipIndx, ipIndik,npos(8,3)
      Integer   iSTSQ(8), iSTLT(8), iSSQ(8,8), nnA(8,8), nInd
      Real*8    tread(2),tcoul(2),tmotr(2),tscrn(2),tcasg(2),tmotr2(2)
      Real*8    Txy(nTxy),V_k(*),Z_p_k(nZpk,*), U_k(*)
      Character*6  Fname
      Character*50 CFmt
      Character*12 SECNAM
      Parameter (SECNAM = 'CHO_GET_GRAD')
      COMMON    /CHOTIME /timings

      parameter (DoRead = .false. )
      parameter (zero = 0.0D0, one = 1.0D0, xone = -1.0D0)
#include "itmax.fh"
#include "Molcas.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "exterm.fh"
#include "chomp2g_alaska.fh"
*#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
#include "print.fh"
      Integer iBDsh(MxShll*8)
      Common /BDshell/ iBDsh

      parameter ( N2 = InfVec_N2 )
      Logical add
      Character*6 mode
      Integer  Cho_F2SP
      External Cho_F2SP
*                                                                      *
************************************************************************
*                                                                      *
      MulD2h(i,j) = iEOR(i-1,j-1) + 1

      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j

      InfVec(i,j,k) = iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)

      IndRed(i,k) = iWork(ip_IndRed-1+nnBstrT(1)*(k-1)+i)

      NNBSTRSH(I,J,K)=IWORK(ip_NNBSTRSH-1+NSYM*NNSHL*(K-1)+NSYM*(J-1)+I)

      ipLab(i) = iWork(ip_Lab+i-1)

      kOffSh(i,j) = iWork(ip_kOffSh+nShell*(j-1)+i-1)

      iShp_rs(i) = iWork(ip_iShp_rs+i-1)

      SvShp(i) = Work(ip_SvShp+i-1)
*
** next is a trick to save memory. Memory in "location 2" is used
** to store this offset array defined later on
*
      iOffShp(i,j) = iWork(ip_iiBstRSh+nSym*nnShl-1+nSym*(j-1)+i)

*
** Jonas 2010
** Thomas Bondo Pedersen, 2013:
**   Usage of these statement functions has been commented out below,
**   hence I comment them out here, too.
*
ctbp  ijList(iSym,jSym,i,j,jDen,jDen2) = ipijList +
ctbp &                        iChOrbR(iSym,jSym,jDen) +
ctbp &                        i + (j-1)*(nChOrb_(iSym,jDen2)+1)
ctbp  ijListTri(iSym,i,j,jDen) = ipijListTri + iChOrbT(iSym,jDen) +
ctbp &                      i + (j-1)*(nChOrb_(iSym,jDen)+1)
*                                                                      *
************************************************************************
*                                                                      *
*     General Initialization                                           *
*                                                                      *
************************************************************************

      iRout = 9
      iPrint = nPrint(iRout)


      CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

      do i=1,2            ! 1 --> CPU   2 --> Wall
         tread(i) = zero  !time read vectors
         tcoul(i) = zero  !time for computing V_k
         tcasg(i) = zero  !time for computing Z_p_k
         tmotr(i) = zero  !time for the MO transf of vectors
         tmotr2(i)= zero  !time for the 2nd MO transf of vectors
         tscrn(i) = zero  !time for screening overhead
      end do

      IREDC = -1  ! unknown reduced set in core

      BatchWarn = .True.
      nInd = 0

      Call set_nnA(nSym,nAorb,nnA)

*
**    Various offsets
*
      MaxB=nBas(1)
      ISTLT(1)=0
      ISTSQ(1)=0
      DO ISYM=2,NSYM
        MaxB=Max(MaxB,nBas(iSym))
        NBB=NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
        NBQ=NBAS(ISYM-1)**2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NBB ! Inactive D matrix
        ISTSQ(ISYM)=ISTSQ(ISYM-1)+NBQ ! Diagonal integrals in full
      END DO

      nnBSQ=0
      DO LSYM=1,NSYM
         DO KSYM=LSYM,NSYM
            ISSQ(KSYM,LSYM) = nnBSQ
            ISSQ(LSYM,KSYM) = nnBSQ ! symmetrization
            nnBSQ = nnBSQ + nBas(kSym)*nBas(lSym)
         END DO
      END DO
*
**     Size for occupied CMO-matrix
*
      Do jDen = 1, nKdens
         lCMOi(jDen) = 0
         iOff_CMOi(1,jDen)=0
         Do i = 1, nSym
            lCMOi(jDen) = lCMOi(jDen) + nBas(i)*nChOrb_(i,jDen)
         End Do
         Do i = 2,nSym
            iOff_CMOi(i,jDen)=iOff_CMOi(i-1,jDen) +
     &                        nBas(i-1)*nChOrb_(i-1,jDen)
         End Do
      End Do

*
**
*
      nI2t=0
      nItmx=0
      Do i=1,5
        nIt(i) = 0
      End Do
      Do jDen=nDen,1,-1
         kOff(1,jDen)=0
         nIt(jDen)=nChOrb_(1,jDen)
         Do i=2,nSym
            kOff(i,jDen)=nIt(jDen)
            nIt(jDen)=nIt(jDen)+nChOrb_(i,jDen)
         End Do
         nI2t=nI2t+nIt(jDen)
         nItmx=Max(nItmx,nIt(jDen))
      End Do
*
**   Initialize pointers to avoid compiler warnings
*
      ipDIAG=ip_Dummy
#if defined (_MOLCAS_MPP_)
      ipjDIAG=ip_Dummy
#endif
      ipDIAH=ip_Dummy
      ipAbs=ip_Dummy
      ipY=ip_Dummy
      ipYQ=ip_Dummy
      ipML=ip_Dummy
      Do i=1,5
        ipSKsh(i)=ip_Dummy
      End Do
      ipIndx=ip_iDummy
      ipIndik=ip_iDummy
      ip_Lab=ip_iDummy
      ip_kOffSh=ip_iDummy
      ip_iShp_rs=ip_iDummy
      ip_SvShp=ip_iDummy
*
      thrv=0.0d0
      xtau=0.0d0
*
**    Construct iBDsh for later use
*
      Do iSyma=1,nSym
         LKsh=0
         Do iaSh=1,nShell
            iSSa=nShell*(iSyma-1)+iaSh
            iBDsh(iSSa) = LKsh
            LKsh = LKsh + nBasSh(iSyma,iaSh)
         End Do
      End Do

!     iShp_rs
      Call GetMem('ip_iShp_rs','Allo','Inte',ip_iShp_rs,nnShl_tot)

************************************************************************
*                                                                      *
*     Initialize a few things for ij-screening //Jonas B               *
*                                                                      *
************************************************************************
      If(DoExchange) Then

         lijList=0
         Do jDen=nKvec,1,-1
            iMOleft=jDen
            iMOright=jDen
            If (DoCAS.and.lSA) iMOright=jDen+2
            iChOrbR(1,1,jDen) = 0
            iChOrbT(1,jDen) = 0
            Do k = 2, nSym
               iChOrbR(k,1,jDen) = iChOrbR(k-1,1,jDen)+
     &                     (nChOrb_(k-1,iMOright)+1)*nChOrb_(1,iMOleft)
            End Do
            Do i=2,nSym
               iChOrbR(1,i,jDen) = iChOrbR(nSym,i-1,jDen)+
     &              (nChOrb_(nSym,iMOright)+1)*nChOrb_(i-1,iMOleft)
               Do k = 2,nSym
                  iChOrbR(k,i,jDen) = iChOrbR(k-1,i,jDen)
     &                 + (nChOrb_(k-1,iMOright)+1)*nChOrb_(i,iMOleft)
               End Do
               iChOrbT(i,jDen) = iChOrbT(i-1,jDen)
     &                 + nChOrb_(i-1,iMOleft)*(nChOrb_(i-1,iMOright)+1)
            End Do
            lijList = max(lijList, iChOrbR(nSym,nSym,jDen)+
     &                nChOrb_(nSym,iMOleft)*(nChOrb_(nSym,iMOright)+1))
         End Do

*
** Define the screening thresholds
*

         Call Get_dScalar('Cholesky Threshold',ThrCom)

         tau = (ThrCom/DBLE(Max(1,nItmx)))*dmpK

         MaxRedT=MaxRed
         Call GAIGOP_SCAL(MaxRedT,'+')

         If (Estimate) tau=tau/DBLE(MaxRedT)
         xtau = sqrt(tau)

         NumVT=NumChT
         Call GAIGOP_SCAL(NumVT,'+')
!        Vector MO transformation screening thresholds
         thrv = ( sqrt(ThrCom/DBLE(Max(1,nItmx)*NumVT)) )*dmpK

#if defined (_MOLCAS_MPP_)
         If (Is_Real_Par() .and. Update) Then
            NNBSTMX=0
            Do i=1,nSym
               NNBSTMX = Max(NNBSTMX,NNBSTR(i,1))
            End Do
            CALL GETMEM('diagJ','Allo','Real',ipjDIAG,NNBSTMX)
            Call FZero(Work(ipjDIAG),NNBSTMX)
         EndIf
#endif

*
** Read the diagonal integrals (stored as 1st red set)
*
         CALL GETMEM('diagI','Allo','Real',ipDIAG,NNBSTRT(1))
         If (Update) CALL CHO_IODIAG(Work(ipDIAG),2) ! 2 means "read"

*
** Allocate memory
*
!        sqrt(D(a,b)) stored in full (squared) dim
         CALL GETMEM('diahI','Allo','Real',ipDIAH,NNBSQ)
         CALL FZERO(Work(ipDIAH),NNBSQ)

         Call GetMem('absc','Allo','Real',ipAbs,MaxB) ! abs(C(l)[k])

         Call GetMem('yc','Allo','Real',ipY,MaxB*nItmx) ! Y(l)[k] vector

         Call GetMem('yq','Allo','Real',ipYQ,nItmx**2) ! Yi[k] vectors

*used to be nShell*something
!        ML[k] lists of largest elements in significant shells
         Call GetMem('MLk1','Allo','Real',ipML,nShell)

!        list of S:= sum_l abs(C(l)[k])
         Call GetMem('SKsh','Allo','Real',ipSKsh(1),nShell*nI2t)
         Do i=2,5
           ipSKsh(i)=ipSKsh(i-1)+nShell*nIt(i-1)    ! for each shell
         End Do

*
** Indx and Indik must be stored for each density, symmetry, etc.
** in case of a batched procedure
*
         Do jDen=1,nKvec
           Do kSym=1,nSym
             nInd = nInd + nChOrb_(kSym,jDen)
           End Do
         End Do

!        Index array
         Call GetMem('Indx','Allo','Inte',ipIndx,(nShell+1)*nInd)

         Call GetMem('Indik','Allo','Inte',ipIndik,
     &               ((nItmx+1)*nItmx+1)*nInd)  !Yi[k] Index array

         Call GetMem('ip_Lab','Allo','Inte',ip_Lab,nShell) ! ipLab

!        kOffSh
         Call GetMem('ip_kOffSh','Allo','Inte',ip_kOffSh,nShell*nSym)

!        shell-pair Frobenius norm of the vectors
         Call GetMem('ip_SvShp','Allo','Real',ip_SvShp,2*nnShl)

*
** Jonas - June 2010:
** allocate memory for rearranged CMO-matrix
*
         Do i=1,nDen
           If (lCMOi(i).gt.0) Then
              Call GetMem('CMO_inv','Allo','Real',ip_CMOi(i), lCMOi(i))
           Else
              ip_CMOi(i)=ip_Dummy
           End If
         End Do
         Call GetMem('ijList','Allo','Inte',ipijList,lijList)
         Call GetMem('ijListTri','Allo','Inte',ipijListTri,lijList)
         Call IZero(iWork(ipijList),lijList)
         Call IZero(iWork(ipijListTri),lijList)

         nQoT = 0
*
** Compute Shell Offsets ( MOs and transformed vectors)
*
         Do iSyma=1,nSym
            LKsh=0
            Do iaSh=1,nShell    ! kOffSh(iSh,iSym)

               iWork(ip_kOffSh+nShell*(iSyma-1)+iaSh-1) = LKsh

               LKsh = LKsh + nBasSh(iSyma,iaSh)
            End Do
         End Do

*
** Determine S:= sum_l C(l)[k]^2  in each shell of C(a,k)
*
         Do jDen=1,nDen
            Do kSym=1,nSym

               Do jK=1,nChOrb_(kSym,jDen)

                  ipMO = ipMSQ(jDen) + ISTSQ(kSym)
     &                 + nBas(kSym)*(jK-1)
*
                  ipSk = ipSKsh(jDen) + nShell*(kOff(kSym,jDen) + jK-1)

                  Do iaSh=1,nShell

                     ipMsh = ipMO + kOffSh(iaSh,kSym)

                     SKsh=zero
                     Do ik=0,nBasSh(kSym,iaSh)-1
                        SKsh = SKsh + Work(ipMsh+ik)**2
                     End Do

                     Work(ipSk+iaSh-1) = SKsh

                  End Do
               End Do
            End Do
         End Do
*
** Reorder CMO-matrix, Needed to construct B-matrix for exchange
** Jonas - June 2010
*
         Do jDen = 1, nKdens
            Do kSym = 1, nSym
*
**If the orbitals come from eigenvalue decomposition, change sign
*
               If (lSA.and.(jDen.ge.3)) Then
                 npos2=npos(ksym,jDen-2)
                 Do jK = 1, nPos2
                    Do jGam = 1, nBas(kSym)
                       ipMO = ipMSQ(jDen) + iSTSQ(kSym)
     &                      + nBas(kSym)*(jK-1) + jGam-1
                       ipMO2 = ip_CMOi(jDen) + iOff_CMOi(kSym,jDen)
     &                      + (jGam-1)*nChOrb_(kSym,jDen) + jK-1
                       Work(ipMO2) = Work(ipMO)
                    End Do
                 End Do
                 Do jK = npos2+1,nChOrb_(kSym,jDen)
                    Do jGam = 1, nBas(kSym)
                       ipMO = ipMSQ(jDen) + iSTSQ(kSym)
     &                      + nBas(kSym)*(jK-1) + jGam-1
                       ipMO2 = ip_CMOi(jDen) + iOff_CMOi(kSym,jDen)
     &                      + (jGam-1)*nChOrb_(kSym,jDen) + jK-1
                       Work(ipMO2) = -Work(ipMO)
                    End Do
                 End Do
               Else
*
                 Do jK = 1, nChOrb_(kSym,jDen)
                    Do jGam = 1, nBas(kSym)
                       ipMO = ipMSQ(jDen) + iSTSQ(kSym)
     &                      + nBas(kSym)*(jK-1) + jGam-1
                       ipMO2 = ip_CMOi(jDen) + iOff_CMOi(kSym,jDen)
     &                      + (jGam-1)*nChOrb_(kSym,jDen) + jK-1
                       Work(ipMO2) = Work(ipMO)
                    End Do
                 End Do
               EndIf
            End Do
         End Do
      End If
*
** Mapping shell pairs from the full to the reduced set
*
      Do iaSh=1,nShell
         Do ibSh=1,iaSh
            iShp = iaSh*(iaSh-1)/2 + ibSh
            iWork(ip_iShp_rs+iShp-1) = Cho_F2SP(iShp)
         End Do
      End Do
************************************************************************
*                                                                      *
*     BIG LOOP OVER VECTORS SYMMETRY                                   *
*                                                                      *
************************************************************************
      DO jSym=1,nSym

         NumCV=NumCho(jSym)
         Call GAIGOP_SCAL(NumCV,'max')
         If (NumCV .lt. 1) GOTO 1000
*
** offsets for active term
*
         iOffZp=0
         Do j=1,jSym-1
            iOffZp = iOffZp + nnP(j)*NumCho(j)
         End Do
*
** Open some files to store exchange auxiliary vectors
*
         If (DoExchange) Then
            iSeed=7
            Do i=1,nSym
               k=muld2h(jSym,i)
               LuRVec(i,1) = IsFreeUnit(iSeed)
               Write(Fname,'(A4,I1,I1)') 'CHTA',i,k
               Call DANAME_MF_WA(LuRVec(i,1),Fname)
               iSeed=iSeed+1
               If (nKvec.ge.2) Then
                  LuRVec(i,2) = IsFreeUnit(iSeed)
                  Write(Fname,'(A4,I1,I1)') 'CHTB',i,k
                  Call DANAME_MF_WA(LuRVec(i,2),Fname)
                  iSeed=iSeed+1
               EndIf
            Enddo
         EndIf
*
** Compute Shell pair Offsets   iOffShp(iSyma,iShp)
*
         LFULL=0

         Do iaSh=1,nShell
          Do ibSh=1,iaSh

           iShp = iaSh*(iaSh-1)/2 + ibSh
           If (iShp_rs(iShp).gt.0) Then
            If (nnBstRSh(Jsym,iShp_rs(iShp),1).gt.0) Then

             Do iSymb=1,nSym
              iSyma=MulD2h(iSymb,Jsym)
              If (iSyma.ge.iSymb) Then

               iWork(ip_iiBstRSh + nSym*nnShl - 1
     &         + nSym*(iShp_rs(iShp)-1) + iSyma) = LFULL

               LFULL = LFULL + nBasSh(iSyma,iaSh)*nBasSh(iSymb,ibSh)
     &        + Min(1,(iaSh-ibSh))*nBasSh(iSyma,ibSh)*nBasSh(iSymb,iaSh)

              EndIf
             End Do
            EndIf
           EndIf
          End Do
         End Do

************************************************************************
*                                                                      *
*     Memory management section                                        *
*                                                                      *
************************************************************************
*
** compute memory needed to store at least 1 vector of JSYM
** and do all the subsequent calculations
*
         mTvec = 0
         MxB=0
         nnOmx=0
         do l=1,nSym
            k=Muld2h(l,JSYM)
            Do jDen=1,nDen
                nnOmx=Max(nnOmx,nChOrb_(l,jDen)*nChOrb_(k,jDen))
                If (nChOrb_(k,jDen).gt.0) Then
                   MxB=Max(MxB,nBas(l)+nChOrb_(l,jDen))
                EndIf
            End Do
            mTvec = mTvec + nAorb(k)*nBas(l)
            If (k.le.l) mTvec = mTvec + nnA(k,l)
         end do

         LFMAX = Max(mTvec,2*LFULL) ! re-use memory for the active vec
         mTvec = nnOmx + Max(MxB,1) ! mem for half transformed + Lik
*
**
*
         iLoc = 3 ! use scratch location in reduced index arrays

         If (NumCho(jSym).lt.1) Then
            JRED1 = 1
            JRED2 = 1
         Else
            JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
            JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last
*                                               !vec
         End If
#if defined (_MOLCAS_MPP_)
         myJRED1=JRED1 ! first red set present on this node
         ntv0=0
#endif

c --- entire red sets range for parallel run
         Call GAIGOP_SCAL(JRED1,'min')
         Call GAIGOP_SCAL(JRED2,'max')

*
** MGD does it need to be so?
*
         DoScreen=.True.
         kscreen=1

         Do JRED=JRED1,JRED2

            If (NumCho(jSym).lt.1) Then
               iVrs=0
               nVrs=0
            Else
               CALL Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)
            End If

            If (nVrs.eq.0) GOTO 999  ! no vectors in that (jred,jsym)

            if (nVrs.lt.0) then
               Write(6,*)SECNAM//
     &          ': Cho_X_nVecRS returned nVrs<0. STOP!'
               call Abend
            endif

            Call Cho_X_SetRed(irc,iLoc,JRED)
c            !set index arrays at iLoc
            if(irc.ne.0)then
              Write(6,*) SECNAM,': cho_X_setred non-zero return code.',
     &                   ' rc= ',irc
              call Abend
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            If(JSYM.eq.1)Then
               Do jden=1,nJdens
                  Call GetMem('rsD','Allo','Real',ipDrs(jden),nRS)
                  Call Fzero(Work(ipDrs(jden)),nRS)
               End Do
               If(iMp2prpt.eq.2) Then
                  Call GetMem('rsD2','Allo','Real',ipDrs2,nRS)
               Else
                  ipDrs2 = ip_Dummy
               End If
            End If

            Call GetMem('MaxM','Max','Real',KDUM,LWORK)

            nVec = Min(LWORK/(nRS+mTvec+LFMAX),nVrs)

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) ' LWORK= ',LWORK
               WRITE(6,*) ' min. mem. need= ',nRS+mTvec+LFMAX
               WRITE(6,*) ' jsym= ',jsym
               WRITE(6,*) ' nRS = ',nRS
               WRITE(6,*) ' mTvec = ',mTvec
               WRITE(6,*) ' LFMAX = ',LFMAX
               irc = 33
               CALL Abend()
               nBatch = -9999  ! dummy assignment
            End If

*                                                                      *
************************************************************************
*                                                                      *
            LREAD = nRS*nVec

            Call GetMem('rsL','Allo','Real',ipLrs,LREAD)
            Call GetMem('ChoT','Allo','Real',ipChoT,mTvec*nVec)
            CALL GETMEM('FullV','Allo','Real',ipLF,LFMAX*nVec)

            If(JSYM.eq.1)Then
C --- Transform the densities to reduced set storage
               mode = 'toreds'
               add  = .false.
               Call play_casscf_sto(irc,iLoc,nJdens,JSYM,ISTLT,ISSQ,
     &                              ipDLT,ipDrs,mode,add)
               If(iMp2prpt .eq. 2) Then
                  Call play_casscf_sto(irc,iLoc,nJdens,JSYM,ISTLT,ISSQ,
     &                              [ipDLT2],[ipDrs2],mode,add)
               End If
            EndIf
*
**  BATCH over the vectors
*

            nBatch = (nVrs-1)/nVec + 1

            If (BatchWarn .and. nBatch.gt.1) Then
               If (iPrint.ge.6) Then
                  Write(6,'(20A3)')('---',I=1,20)
                  Write(6,*)' Batch procedure used.'//
     &                      ' Increase memory if possible!'
                  Write(6,'(20A3)')('---',I=1,20)
                  Write(6,*)
                  Call XFlush(6)
               End If
               BatchWarn = .False.
            EndIf

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

               If (NUMV.le.0 .or.NUMV.ne.JNUM ) then
                  irc=77
                  RETURN
               End If

               CALL CWTIME(TCR2,TWR2)
               tread(1) = tread(1) + (TCR2 - TCR1)
               tread(2) = tread(2) + (TWR2 - TWR1)

************************************************************************
************************************************************************
**                                                                    **
**                                                                    **
**    Coulomb term                                                    **
**           V{#J} = sum_ab  L(ab,{#J}) * D(ab)                       **
**                                                                    **
**                                                                    **
************************************************************************
************************************************************************
               If(JSYM.eq.1)Then

                 CALL CWTIME(TCC1,TWC1)

*
**  Inactive Coulomb term
*
                 Do jden=1,nJdens
                    CALL DGEMV_('T',nRS,JNUM,
     &                         One,Work(ipLrs),nRS,
     &                         Work(ipDrs(jden)),1,
     &                         zero,V_k(jVec+(jDen-1)*NumCho(1)),1)
                 End Do
*
**  MP2 Coulomb term
*
                 If(iMp2prpt .eq. 2) Then
                    CALL DGEMV_('T',nRS,JNUM,
     &                         One,Work(ipLrs),nRS,
     &                         Work(ipDrs2),1,
     &                         zero,U_k(jVec),1)
                 End If
*
                 CALL CWTIME(TCC2,TWC2)
                 tcoul(1) = tcoul(1) + (TCC2 - TCC1)
                 tcoul(2) = tcoul(2) + (TWC2 - TWC1)
               EndIf
************************************************************************
************************************************************************
**                                                                    **
**              E X C H A N G E    T E R M                            **
**                                                                    **
**                                                                    **
************************************************************************
************************************************************************
*
               If (DoExchange) Then

                  CALL CWTIME(TCS1,TWS1)
************************************************************************
*                                                                      *
*   1) Screening                                                       *
*                                                                      *
*      Select only important  ij pairs                                 *
*      For this, one computes the quantity                             *
*         Yik = sum_mu_nu (mu nu | mu nu)^1/2 X_mu_i X_nu_k            *
*      with (mu nu | mu nu) = sum_J (L_mu_nu,J)^2                      *
*                                                                      *
*                                                                      *
*      a) Estimate the diagonals :   D(mu,nu) = sum_J (L_mu_nu,J)^2    *
*                                                                      *
************************************************************************
                  If (Estimate) Then

                     Call Fzero(Work(ipDiag+iiBstR(jSym,1)),
     &                          NNBSTR(jSym,1))

                     Do krs=1,nRS

                        mrs = iiBstR(JSYM,iLoc) + krs
                        jrs = IndRed(mrs,iLoc) ! address in 1st red set

                        Do jvc=1,JNUM

                           ipL = ipLrs + nRS*(jvc-1)

                           Work(ipDiag+jrs-1) = Work(ipDiag+jrs-1)
     &                                     + Work(ipL+krs-1)**2

                        End Do

                     End Do

                  EndIf

                  CALL CWTIME(TCS2,TWS2)
                  tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                  tscrn(2) = tscrn(2) + (TWS2 - TWS1)

                  CALL CWTIME(TCX1,TWX1)

*
** Reorder vectors to Full-dimensions
**
** Vectors are returned in the storage LaJ,b with the restriction:
**    Sym(a).ge.Sym(b)
** and blocked in shell pairs
*

                  CALL FZero(Work(ipLF),LFULL*JNUM)
                  CALL FZero(Work(ip_SvShp),2*nnShl)

                  CALL CHO_getShFull(Work(ipLrs),lread,JNUM,JSYM,
     &                               IREDC,ipLF,Work(ip_SvShp),
     &                               iWork(ip_iShp_rs))

                  CALL CWTIME(TCX2,TWX2)
                  tmotr(1) = tmotr(1) + (TCX2 - TCX1)
                  tmotr(2) = tmotr(2) + (TWX2 - TWX1)

************************************************************************
*                                                                      *
*   1) Screening                                                       *
*                                                                      *
*      b) DH(mu,nu)=sqrt(D(mu,nu))                                     *
*         Only the symmetry blocks with compound symmetry JSYM         *
*         are computed                                                 *
*                                                                      *
************************************************************************
                  IF (DoScreen) THEN

                     CALL CWTIME(TCS1,TWS1)

                     mode = 'tosqrt'
                     ired1 = 1 ! location of the 1st red set
                     add  = .false.
                     nMat = 1
                     Call play_casscf_sto(irc,ired1,nMat,JSYM,ISTLT,
     &                                  ISSQ,[ipDIAH],[ipDIAG],mode,add)

                     CALL CWTIME(TCS2,TWS2)
                     tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                     tscrn(2) = tscrn(2) + (TWS2 - TWS1)

                  ENDIF

************************************************************************
*                                                                      *
*   1) Screening                                                       *
*                                                                      *
*      c) 1st MO transformation of DH(mu,nu)                           *
*            Y(mu)[k] = sum_nu  DH(mu,nu) * |C(nu)[k]|                 *
*                                                                      *
************************************************************************

                  nInd = 0
                  Do jDen=1,nKvec
*
** Choose which MO sets on each side
*
                    iMOleft=jDen
                    iMOright=jDen
                    If (DoCAS.and.lSA) iMOright=jDen+2
*

                    Do kSym=1,nSym

                       lSym=MulD2h(JSYM,kSym)
                       Nik= nChOrb_(kSym,iMOleft)*nChOrb_(lSym,iMOright)
                       nIJR(kSym,lSym,jDen) = Nik
                       nIJ1(kSym,lSym,jDen) = Nik
                       If ((JSYM.eq.1).and.iMOleft.eq.iMOright)
     &                                Nik = nChOrb_(kSym,iMOleft)
     &                                    *(nChOrb_(kSym,iMOleft)+1)/2
                       nIJ1(kSym,lSym,jDen) = Nik
                       If (Nik.eq.0) Go To 98

                       ip_R = ipChoT
     &                      + (nChOrb_(lSym,iMOright)+nBas(lSym))*JNUM

                       Call Fzero(Work(ip_R),Nik*JNUM)

                       Do jK=1,nChOrb_(kSym,iMOleft)

                        CALL FZero(Work(ipChoT),
     &                         (nChOrb_(lSym,iMOright)+nBas(lSym))*JNUM)

                        ipMO = ipMSQ(iMOleft) + ISTSQ(kSym)
     &                       + nBas(kSym)*(jK-1)

                        ipYk = ipY + MaxB*(kOff(kSym,iMOleft)+jK-1)

                        ipYQk = ipYQ + nIt(iMOright)
     &                         *(kOff(kSym,iMOleft)+jK-1)

                        ipMLk = ipML

                        ipIndSh = ipIndx+nInd*(nShell+1)

                        ipIndikk = ipIndik+nInd*((nItmx+1)*nItmx+1)

                        ipSk=ipSKsh(iMOleft)+
     &                       nShell*(kOff(kSym,iMOleft)+jK-1)


                        IF (DoScreen .and. iBatch.eq.1) THEN
                           CALL CWTIME(TCS1,TWS1)
C------------------------------------------------------------------
C --- Setup the screening
C------------------------------------------------------------------
                           ipDIH = ipDIAH + ISSQ(lSym,kSym)

                           Do ik=0,nBas(kSym)-1
                              Work(ipAbs+ik) = abs(Work(ipMO+ik))
                           End Do
                           If (lSym.ge.kSym.and.nBas(lSym).ge.1) Then

                              nBs = Max(1,nBas(lSym))

                              CALL DGEMV_('N',nBas(lSym),nBas(kSym),
     &                                   ONE,Work(ipDIH),nBs,
     &                                       Work(ipAbs),1,
     &                                  ZERO,Work(ipYk),1)

                           Else If (nBas(kSym).ge.1) Then

                              nBs = Max(1,nBas(kSym))

                              CALL DGEMV_('T',nBas(kSym),nBas(lSym),
     &                                   ONE,Work(ipDIH),nBs,
     &                                       Work(ipAbs),1,
     &                                  ZERO,Work(ipYk),1)

                           EndIf

************************************************************************
*                                                                      *
*   1) Screening                                                       *
*                                                                      *
*      d) 2nd MO transformation of DH(mu,nu)                           *
*            Y(i)[k] = sum_mu  |C(mu)[i]| * Y(mu)[k]                   *
*                                                                      *
************************************************************************

*                           iStart=Min(1,lSym-kSym)
*     &                           +jK*(1-Min(1,lSym-kSym)) ! 1 or jK
                           If ((kSym.ne.lSym).or.
     &                        (iMOleft.ne.iMOright)) Then
                               iStart=1
                           Else
                               iStart=jK
                           EndIf

                           nQo=0
                           Do i=iStart,nChOrb_(lSym,iMOright)

                              ipMO_ = ipMSQ(iMOright) + ISTSQ(lSym)
     &                             + nBas(lSym)*(i-1)

                              Do ik=0,nBas(lSym)-1
                                 Work(ipAbs+ik) = abs(Work(ipMO_+ik))
                              End Do
*
                              Work(ipYQk+i-1)=ddot_(nBas(lSym),
     &                             Work(ipAbs),1,Work(ipYk),1)

                              If (Work(ipYQk+i-1).ge.xtau) Then
                                 nQo=nQo+1
                                 If((iBatch.ne.1) .or.
     &                               (JRED.ne.1)) Go To 1111
*                                 ind1 = iWork(ijList(lSym,kSym,
*     &                                  0,jK,jDen,iMOright))+ 1
*                                 iWork(ijList(lSym,kSym,
*     &                                 ind1,jK,jDen,iMOright))=i
*                                 iWork(ijList(lSym,kSym,
*     &                                 0,jK,jDen,iMOright))=ind1
                                 nQoT = nQoT + 1
                                 If((lSym .eq. kSym) .and.
     &                              (i .ne. jK)      .and.
     &                              (iMOright.eq.iMOleft)) Then
*                                    ind2 =iWork(ijList(kSym,kSym,
*     &                                    0,i,jDen,iMOright))+1
*                                    iWork(ijList(kSym,kSym,
*     &                                    ind2,i,jDen,iMOright))=jK
*                                    iWork(ijList(lSym,kSym,
*     &                                    0,i,jDen,iMOright))= ind2
                                    nQoT = nQoT + 1
                                 End If
                                 If((lSym.eq.kSym).and.
     &                              (iMOleft.eq.iMOright)) Then
*                                    ind3=
*     &                                 iWork(ijListTri(kSym,0,jK,jDen))
*     &                                 +1
*                                    iWork(ijListTri(kSym,ind3,jK,jDen))=
*     &                                    i
*                                    iWork(ijListTri(kSym,0,jK,jDen))=
*     &                                     ind3
                                 End If
 1111                            iWork(ipIndikk+nQo)=i
                              Endif

                           End Do
                           iWork(ipIndikk)=nQo
************************************************************************
*                                                                      *
*   1) Screening                                                       *
*                                                                      *
*      e) List the shells present in Y(l)[k] by the largest element    *
*         and sort the list                                            *
*                                                                      *
************************************************************************

                           Do ish=1,nShell
                              YshMax=zero
                              Do ibs=1,nBasSh(lSym,ish)
                                 YshMax = Max(YshMax,
     &                             Work(ipYk+koffSh(ish,lSym)+ibs-1))
                              End Do
                              Work(ipMLk+ish-1) = YshMax
                           End Do

                           Do ish=1,nShell
                              iWork(ipIndSh+ish) = ish
                           End Do

************************************************************************
*                                                                      *
*   1) Screening                                                       *
*                                                                      *
*      f) Screening                                                    *
*                                                                      *
*   Here we use a non-exact bound for the exchange matrix to achieve   *
*   linear scaling. The positive definiteness of the exchange matrix   *
*   combined with the structure of the density matrix makes this       *
*   bound acceptable and likely to be almost exact for what concerns   *
*   the exchange energy                                                *
*                                                                      *
*   The exact bounds (quadratic scaling of the MO transformation)      *
*   would be                                                           *
*      If (Work(ipMLk+jml-1)*Work(ipMLk).ge. tau) then                 *
*                                                                      *
*                                                                      *
************************************************************************

                           numSh=0  ! # of significant shells
                           jml=1
                           Do while (jml.le.nShell)

                              YMax=Work(ipMLk+jml-1)
                              jmlmax=jml

                              Do iml=jml+1,nShell  ! get the max
                                 If (Work(ipMLk+iml-1).gt.YMax) then
                                    YMax = Work(ipMLk+iml-1)
                                    jmlmax = iml
                                 Endif
                              End Do

                              If(jmlmax.ne.jml) then  ! swap positions
                                xTmp = Work(ipMLk+jml-1)
                                iTmp = iWork(ipIndSh+jml)
                                Work(ipMLk+jml-1) = YMax
                                iWork(ipIndSh+jml)=iWork(ipIndSh+jmlmax)
                                Work(ipMLk+jmlmax-1) = xTmp
                                iWork(ipIndSh+jmlmax) = iTmp
                              Endif

                              If ( Work(ipMLk+jml-1) .ge. xtau ) then
                                numSh = numSh + 1
                              else
                                jml=nShell  ! exit the loop
                              endif

                              jml=jml+1

                           End Do

                           iWork(ipIndSh) = numSh

                           CALL CWTIME(TCS2,TWS2)
                           tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                           tscrn(2) = tscrn(2) + (TWS2 - TWS1)
C------------------------------------------------------------------
                        ENDIF    ! Screening setup


                        CALL CWTIME(TCT1,TWT1)

************************************************************************
*                                                                      *
*               E X C H A N G E    T E R M                             *
*                                                                      *
*   2) MO transformation                                               *
*      a) 1st half transformation                                      *
*                                                                      *
*      Transform vectors for shells in the list ML[k]                  *
*                                                                      *
*      Screening based on the Frobenius norm: sqrt(sum_ij  A(i,j)^2)   *
*         || La,J[k] ||  .le.  || Lab,J || * || Cb[k] ||               *
*                                                                      *
************************************************************************

                        IF (lSym.ge.kSym) Then


                            Do iSh=1,iWork(ipIndSh)

                               iaSh = iWork(ipIndSh+iSh)

                               iOffSha = kOffSh(iaSh,lSym)

                               iWork(ip_Lab+iaSh-1)= ipChoT
     &                                         + nChOrb_(lSym,iMOright)
     &                                         * JNUM
     &                                         + iOffSha*JNUM

                               ibcount=0

                               Do ibSh=1,nShell

                                  iOffShb = kOffSh(ibSh,kSym)

                                  iShp = iTri(iaSh,ibSh)

                                If ( iShp_rs(iShp) .gt. 0 ) Then

                                 If ( nnBstRSh(JSym,iShp_rs(iShp),iLoc)*
     &                              nBasSh(lSym,iaSh)*
     &                              nBasSh(kSym,ibSh) .gt. 0
     &                              .and. sqrt(abs(Work(ipSk+ibSh-1)*
     &                                         SvShp(iShp_rs(iShp)) ))
     &                              .ge. thrv )Then

                                    ibcount = ibcount + 1

                                    jOff = iOffShp(lSym,iShp_rs(iShp)) -
     &                                     nBasSh(lSym,ibSh)*
     &                                     nBasSh(kSym,iaSh)*
     &                             Min(0,(iaSh-ibSh))/Max(1,(ibSh-iaSh))

*
**  LaJ,[k] = sum_b  L(aJ,b) * C(b)[k]
** ------------------------------------
*
                                 CALL DGEMV_('N',nBasSh(lSym,iaSh)*JNUM,
     &                                        nBasSh(kSym,ibSh),
     &                                    ONE,Work(ipLF+jOff*JNUM),
     &                                        nBasSh(lSym,iaSh)*JNUM,
     &                                     Work(ipMO+ioffShb),1,
     &                                    ONE,Work(ipLab(iaSh)),1)


                                 EndIf

                                EndIf

                               End Do
*
** The following re-assignement is used later on to check if the
** iaSh vector LaJ[k] can be neglected because identically zero
*

                               iWork(ip_Lab+iaSh-1) = ipLab(iaSh)*
     &                                                Min(1,ibcount)
     &                                    + ipAbs*(1-Min(1,ibcount))

                            End Do


                        Else   ! lSym < kSym


                            Do iSh=1,iWork(ipIndSh)

                               iaSh = iWork(ipIndSh+iSh)

                               iOffSha = kOffSh(iaSh,lSym)

                               iWork(ip_Lab+iaSh-1)= ipChoT
     &                                        + (nChOrb_(lSym,iMOright)
     &                                        + iOffSha)*JNUM
                               ibcount=0

                               Do ibSh=1,nShell

                                 iOffShb = kOffSh(ibSh,kSym)

                                 iShp = iTri(iaSh,ibSh)

                                 If ( iShp_rs(iShp) .gt. 0 ) Then

                                  If (nnBstRSh(JSym,iShp_rs(iShp),iLoc)*
     &                             nBasSh(lSym,iaSh)*
     &                             nBasSh(kSym,ibSh) .gt. 0
     &                             .and. sqrt(abs(Work(ipSk+ibSh-1)*
     &                                        SvShp(iShp_rs(iShp)) ))
     &                             .ge. thrv ) Then

                                   ibcount = ibcount + 1

                                   jOff = iOffShp(kSym,iShp_rs(iShp)) -
     &                                    nBasSh(kSym,iaSh)*
     &                                    nBasSh(lSym,ibSh)*
     &                             Min(0,(ibSh-iaSh))/Max(1,(iaSh-ibSh))

*
**  LJa,[k] = sum_b  L(b,Ja) * C(b)[k]
** ------------------------------------
*
                                   CALL DGEMV_('T',nBasSh(kSym,ibSh),
     &                                       JNUM*nBasSh(lSym,iaSh),
     &                                    ONE,Work(ipLF+jOff*JNUM),
     &                                        nBasSh(kSym,ibSh),
     &                                     Work(ipMO+ioffShb),1,
     &                                    ONE,Work(ipLab(iaSh)),1)

                                  EndIf

                                 Endif

                               End Do

*
** The following re-assignement is used later on to check if the
** iaSh vector LaJ[k] can be neglected because identically zero
*
                               iWork(ip_Lab+iaSh-1) = ipLab(iaSh)*
     &                                               Min(1,ibcount)
     &                                   + ipAbs*(1-Min(1,ibcount))

                            End Do

                        EndIf

************************************************************************
*                                                                      *
*   2) MO transformation                                               *
*      b) 2nd half transformation                                      *
*                                                                      *
************************************************************************

                        nQo=iWork(ipIndikk)

                        Do ir=1,nQo

                          it = iWork(ipIndikk+ir)

                          ipMO = ipMSQ(iMOright) + ISTSQ(lSym)
     &                         + nBas(lSym)*(it-1)

                          ip_Lik = ipChoT + JNUM*(it-1)

                          If (lSym.ge.kSym) Then

                           Do iSh=1,iWork(ipIndSh)

                              iaSh = iWork(ipIndSh+iSh)

                              iaSkip= Min(1,Max(0,
     &                                abs(ipLab(iaSh)-ipAbs)))!=1 or 0

                              Do i=1,iaSkip

                                 ip_Cai = ipMO + kOffsh(iaSh,lSym)
*
**  LJi[k] = sum_a  LaJ[k] * Cai
** ------------------------------
*
                                 CALL DGEMV_('T',nBasSh(lSym,iaSh),JNUM,
     &                                      One,Work(ipLab(iaSh)),
     &                                      nBasSh(lSym,iaSh),
     &                                      Work(ip_Cai),1,
     &                                      one,Work(ip_Lik),1)


                              End Do

                           End Do


                          Else   ! lSym < kSym


                           Do iSh=1,iWork(ipIndSh)

                              iaSh = iWork(ipIndSh+iSh)

                              iaSkip= Min(1,Max(0,
     &                                abs(ipLab(iaSh)-ipAbs)))!=1 or 0

                              Do i=1,iaSkip

                                 ip_Cai = ipMO + kOffsh(iaSh,lSym)
*
**   LJi[k] = sum_a  LJa[k] * Cai
** --------------------------------
*
                                 CALL DGEMV_('N',JNUM,nBasSh(lSym,iaSh),
     &                                      One,Work(ipLab(iaSh)),
     &                                      JNUM,Work(ip_Cai),1,
     &                                      one,Work(ip_Lik),1)


                              End Do

                           End Do

                          EndIf
*
**  Copy LJi[k] in the standard ordered matrix Lik,J
*
                          If ((jSym.eq.1).and.
     &                       (iMOright.eq.iMOleft)) Then
                             itk = it*(it-1)/2 + jK
                          Else
                             itk = nChOrb_(lSym,iMOright)*(jK-1) + it
                          EndIf
                          ip_Rik = ip_R - 1 + itk
                          call dcopy_(JNUM,Work(ip_Lik),1,
     &                                    Work(ip_Rik),Nik)

                        End Do

                        nInd = nInd+1

                        CALL CWTIME(TCT2,TWT2)
                        tmotr2(1) = tmotr(1) + (TCT2 - TCT1)
                        tmotr2(2) = tmotr(2) + (TWT2 - TWT1)


                       End Do  ! loop over k MOs

                       CALL CWTIME(TCT1,TWT1)

************************************************************************
*                                                                      *
*   3) Put to disk                                                     *
*                                                                      *
************************************************************************
                       iAdr = Nik*(JVEC-1)
                       call DDAFILE(LuRVec(lSym,jDen),1,Work(ip_R),
     &                                                Nik*JNUM,iAdr)

                       CALL CWTIME(TCT2,TWT2)
                       tmotr(1) = tmotr(1) + (TCT2 - TCT1)
                       tmotr(2) = tmotr(2) + (TWT2 - TWT1)
 98                    Continue

                    End Do   ! loop over MOs symmetry

                  End Do   ! loop over densities


C ************  END EXCHANGE CONTRIBUTION  ****************

C --- Diagonals updating. It only makes sense if Nscreen > 0

                  If (Update .and. Nscreen .gt. 0) Then

                     CALL CWTIME(TCS1,TWS1)
C ---------------------------------------------------------------------
C --- update the diagonals :   D(a,b) = D(a,b) - sum_J (Lab,J)^2
C
C --- subtraction is done in the 1st reduced set
#if defined (_MOLCAS_MPP_)
                     If (Is_Real_Par()) then

                      Do krs=1,nRS

                        mrs = iiBstR(JSYM,iLoc) + krs
                        jrs = IndRed(mrs,iLoc) - iiBstR(JSYM,1)

                        Do jvc=1,JNUM

                           ipL = ipLrs + nRS*(jvc-1)
                           Work(ipjDiag+jrs-1) = Work(ipjDiag+jrs-1)
     &                                         + Work(ipL+krs-1)**2
                        End Do

                      End Do

                     Else

                      Do krs=1,nRS

                        mrs = iiBstR(JSYM,iLoc) + krs
                        jrs = IndRed(mrs,iLoc) ! address in 1st red set

                        Do jvc=1,JNUM

                           ipL = ipLrs + nRS*(jvc-1)
                           Work(ipDiag+jrs-1) = Work(ipDiag+jrs-1)
     &                                        - Work(ipL+krs-1)**2
                        End Do

                      End Do

                     EndIf

#else
                     Do krs=1,nRS

                        mrs = iiBstR(JSYM,iLoc) + krs
                        jrs = IndRed(mrs,iLoc) ! address in 1st red set

                        Do jvc=1,JNUM

                           ipL = ipLrs + nRS*(jvc-1)
                           Work(ipDiag+jrs-1) = Work(ipDiag+jrs-1)
     &                                        - Work(ipL+krs-1)**2
                        End Do

                     End Do
#endif

                     CALL CWTIME(TCS2,TWS2)
                     tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                     tscrn(2) = tscrn(2) + (TWS2 - TWS1)

                  EndIf

               EndIf ! DoExchange

************************************************************************
************************************************************************
**                                                                    **
**                                                                    **
**    Active term                                                     **
**                                                                    **
**                                                                    **
************************************************************************
************************************************************************
               If (DoCAS) Then

                  CALL CWTIME(TCC1,TWC1)
*
** Set up the skipping flags and the pointers ipLpq
** The memory used before for the full-dimension AO-vectors
**     is now re-used to store half and full transformed
**     vectors in the active space
*
                  lChoa=0
                  Do i=1,nSym

                     k = Muld2h(i,JSYM)
                     iSkip(k) = Min(1,
     &                    nAorb(k)*nBas(i))

                     ipLpq(k) = ipLF + lChoa       ! Lvb,J
                     ipLxy(k) = ipLF + lChoa       ! Lvw,J
     &                        + nAorb(k)*nBas(i)*JNUM

                     lChoa = lChoa + nAorb(k)*nBas(i)*JNUM
                     If (k.le.i) lChoa = lChoa + nnA(k,i)*JNUM

                  End Do

                  iSwap = 0  ! Lvb,J are returned
                  kMOs = 1  !
                  nMOs = 1  ! Active MOs (1st set)
*
                  Do itran=1,nAdens
                     iMO1=1
                     iMO2=1
                     If (itran.eq.2) iMO1=2

************************************************************************
*                                                                      *
*     MO transformation of Cholesky vectors                            *
*                                                                      *
*         1) Lvb,J = sum_a  C(v,a) * Lab,J                             *
*                                                                      *
************************************************************************

                     CALL CHO_X_getVtra(irc,Work(ipLrs),LREAD,jVEC,JNUM,
     &                             JSYM,iSwap,IREDC,nMOs,kMOs,
     &                             ipAorb(1,iMO1),nAorb,ipLpq,iSkip,
     &                             DoRead)

                     if (irc.ne.0) then
                        RETURN
                     endif

************************************************************************
*                                                                      *
*     MO transformation of Cholesky vectors                            *
*                                                                      *
*         2) Lvw,J = sum_b  Lvb,J * C(w,b)                             *
*                                                                      *
************************************************************************
                     If ((JSYM.eq.1).and.(iMO1.eq.iMO2)) Then

                        Do iSymb=1,nSym

                           NAv = nAorb(iSymb)

                           If(NAv.gt.0)Then

                            Do JVC=1,JNUM

                             ipLvb=ipLpq(iSymb)+NAv*NBAS(iSymb)*(JVC-1)
                             ipLvw=ipLxy(iSymb)+nnA(iSymb,iSymb)*(JVC-1)

                             CALL DGEMM_Tri('N','T',NAv,NAv,NBAS(iSymb),
     &                                  One,Work(ipLvb),NAv,
     &                                      Work(ipAorb(iSymb,iMO2)),
     &                                 NAv,Zero,Work(ipLvw),NAv)

                            End Do

                           EndIf

                        End Do

                     Else

                        Do iSymb=1,nSym

                           iSymv = MulD2h(JSYM,iSymb)
                           NAv = nAorb(iSymv)
                           NAw = nAorb(iSymb) ! iSymb=iSymw

                           If(NAv*NAw.ne.0 .and. iSymv.le.iSymb)Then

                            Do JVC=1,JNUM

                             ipLvb=ipLpq(iSymv)+NAv*NBAS(iSymb)*(JVC-1)
                             ipLvw=ipLxy(iSymv)+NAv*NAw*(JVC-1)

                             CALL DGEMM_('N','T',NAv,NAw,NBAS(iSymb),
     &                                  One,Work(ipLvb),NAv,
     &                                      Work(ipAorb(iSymb,iMO2)),
     &                                 NAw,Zero,Work(ipLvw),NAv)

                            End Do

                           EndIf

                        End Do

                     EndIf
************************************************************************
*                                                                      *
*     Evaluation of the Z_p_k                                          *
*                                                                      *
*         Z(p){#J} = sum_xy  T(xy,p) * L(xy,{#J})                      *
*                                                                      *
*     T(xy,p) : is stored by compound symmetry JSYM                    *
*               the indices {xy} are stored as PACKED (sym x.le.sym y) *
*                                                                      *
************************************************************************
                     Do iTxy=itran,nAdens
                       iAvec=itran+iTxy-1
                       Do iSymy=1,nSym

                         iSymx=MulD2h(iSymy,JSYM)

                         If (iSymx.le.iSymy.and.nnA(iSymx,iSymy).ne.0)
     &                      Then

                            ipZp = iOffZp + nnP(jSym)*(JVEC-1) + 1

                            If (iMO1.eq.iMO2) Then
                              CALL DGEMM_('T','N',nnP(jSym),JNUM,
     &                                         nnA(iSymx,iSymy),
     &                               ONE,Txy(ipTxy(iSymx,iSymy,iTxy)),
     &                                   nnP(jSym),
     &                                   Work(ipLxy(iSymx)),
     &                                   nnA(iSymx,iSymy),
     &                               ONE,Z_p_k(ipZp,iAvec),nnP(jSym))
                            Else
*MGD may rearrange the loops
                              Do i=0,nnP(jSym)-1
                                ioff=ipTxy(iSymx,iSymy,iTxy)+
     &                                nnA(iSymx,iSymy)*i
                                Do j=0,JNUM-1
                                  jOff=ipLxy(iSymx)+
     &                                 nAorb(iSymx)*nAorb(iSymy)*j
                                  temp=0.0d0
*MGD don't work with symmetry
                                  Do k=0,nAOrb(iSymx)-1
                                    Do l=0,k
                                       temp=temp+0.5d0*
     &                                     Txy(ioff+k*(k+1)/2+l)*
     &                                    (Work(jOff+k*nAOrb(iSymx)+l)+
     &                                     Work(jOff+l*nAOrb(iSymx)+k))
                                    End Do
                                  End Do
                                  Z_p_k(ipZp+j*nnP(jSym)+i,iAvec)=
     &                                Z_p_k(ipZp+j*nnP(jSym)+i,iAvec)+
     &                                temp
                                End Do
                              End Do
                            EndIf

                         Endif

                       End Do
                     End Do

                  End Do

                  CALL CWTIME(TCC2,TWC2)
                  tcasg(1) = tcasg(1) + (TCC2 - TCC1)
                  tcasg(2) = tcasg(2) + (TWC2 - TWC1)

               EndIf  ! DoCAS

************************************************************************
************************************************************************
**                                                                    **
**    Epilogue                                                        **
**                                                                    **
************************************************************************
************************************************************************

            END DO  ! end batch loop

C --- free memory
            CALL GETMEM('FullV','Free','Real',ipLF,LFMAX*nVec)
            Call GetMem('ChoT','Free','Real',ipChoT,mTvec*nVec)
            Call GetMem('rsL','Free','Real',ipLrs,LREAD)

            If(JSYM.eq.1)Then
              do jden=nJdens,1,-1
               Call GetMem('rsD','Free','Real',ipDrs(jden),nRS)
              end do
              If(iMp2prpt .eq. 2) Then
                 Call GetMem('rsD2','Free','Real',ipDrs2,nRS)
              End If
            EndIf

999         Continue

C --- Screening control section
            DoScreen = kscreen.eq.Nscreen

            if (.not.DoScreen) then
                kscreen = kscreen + 1
            else
                kscreen = 1
            endif

            If (DoExchange) Then
#if defined (_MOLCAS_MPP_)
               If (Is_Real_Par() .and. Update .and. DoScreen) Then
                  Call GaDsum(Work(ipjDiag),nnBSTR(JSYM,1))
                  Call Daxpy_(nnBSTR(JSYM,1),xone,Work(ipjDiag),1,
     &                       Work(ipDiag+iiBstR(JSYM,1)),1)
                  Call Fzero(Work(ipjDiag),nnBSTR(JSYM,1))
               EndIf
C--- Need to activate the screening to setup the contributing shell
C--- indices the first time the loop is entered .OR. whenever other nodes
C--- have performed screening in the meanwhile
               If (Is_Real_Par().and..not.DoScreen.and.nVrs.eq.0) Then
                  ntv0=ntv0+1
                  DoScreen = (JRED.lt.myJRED1 .or. ntv0.ge.Nscreen)
                  if (DoScreen) ntv0=0
               EndIf
#endif
            EndIf

         END DO   ! loop over red sets

         If (DoExchange) Then
           Do jDen=1,nKvec
              Do i=1,nSym
                 Call DACLOS(LuRVec(i,jDen))
              End Do
           End Do
         End If

1000  CONTINUE


      END DO  ! loop over JSYM
*****************************************************************
*              Allocate a field to be used by Compute_A_jk later
*              since allocations cannot be made at that stage
******************************************************************
      If(DoExchange) THen
        nIJMax = 0
        Do jDen = 1, nKvec
           Do iSym1 = 1, nSym
              Do iSym2 = 1, nSym
                 nIJMax = max(nIJMax,nIJR(iSym1,iSym2,jDen))
              End Do
           End Do
        End Do
        ljkVec = 2*nIJMax
        Call GetMem('JKVEC','Allo','Real',ip_VJ,ljkVec)
      End If

      Call GetMem('ip_iShp_rs','Free','Inte',ip_iShp_rs,nnShl_tot)
      If (DoExchange) Then
         Call GetMem('ip_SvShp','Free','Real',ip_SvShp,2*nnShl)
         Call GetMem('ip_kOffSh','Free','Inte',ip_kOffSh,nShell*nSym)
         Call GetMem('ip_Lab','Free','Inte',ip_Lab,nShell)
         Call GetMem('Indik','Free','Inte',ipIndik,(nItmx+1)*nItmx)
         Call GetMem('Indx','Free','Inte',ipIndx,nShell)
         Call GetMem('SKsh','Free','Real',ipSKsh(1),nShell*nI2t)
         Call GetMem('MLk1','Free','Real',ipML,nShell)
         Call GetMem('yq','Free','Real',ipYQ,nItmx**2)
         Call GetMem('yc','Free','Real',ipY,MaxB*nItmx)
         Call GetMem('absc','Free','Real',ipAbs,MaxB)
         CALL GETMEM('diahI','Free','Real',ipDIAH,NNBSQ)
#if defined (_MOLCAS_MPP_)
         If (Is_Real_Par().and.Update)CALL GETMEM('diagJ','Free','Real',
     &                                            ipjDIAG,NNBSTMX)
#endif
         CALL GETMEM('diagI','Free','Real',ipDIAG,NNBSTRT(1))
      EndIf


      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1
#ifdef _CD_TIMING_
      ChoGet_CPU = TOTCPU
      ChoGet_Wall = TOTWALL
#endif
*                                                                      *
*---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky Gradients timing from '//SECNAM
      Write(6,CFmt)'----------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'                                CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ VECTORS                     '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COULOMB CONTRIB.                 '
     &                           //'         ',tcoul(1),tcoul(2)
         Write(6,'(2x,A26,2f10.2)')'SCREENING OVERHEAD               '
     &                           //'         ',tscrn(1),tscrn(2)
         Write(6,'(2x,A26,2f10.2)')'INACT MO-TRANSFORM VECTORS       '
     &                           //'         ',tmotr(1),tmotr(2)
         Write(6,'(2x,A26,2f10.2)')'INACT MO-TRANSFORM VECTORS 2     '
     &                           //'         ',tmotr2(1),tmotr2(2)
         Write(6,'(2x,A26,2f10.2)')'ACTIVE CONTRIB.                  '
     &                           //'         ',tcasg(1),tcasg(2)
         Write(6,*)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif

      irc  = 0


      Return
      END
