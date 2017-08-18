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
* Copyright (C) Valera Veryazov                                        *
*               Luca De Vico                                           *
************************************************************************
*  CHARGE
*
*> @brief
*>   Compute and print Mulliken populations
*> @modified_by V. Veryazov
*> @modified_by L. De Vico
*>
*> @details
*> For a given set of natural orbitals compute the
*> Mulliken population analysis. The final print out
*> reports the populations per center and basis
*> function type as well as the gross atomic charges.
*>
*> - \p iCase = ``2`` used for spin independent case (RHF or RASSCF)
*> - \p iCase = ``3`` same as ``2`` but suitable for rasscf *spin* populations.
*> - \p iCase = ``0`` or ``1`` to calculate alpha and beta populations
*>                    (\p iCase = ``0`` only computes alpha contributions, and
*>                    \p iCase = ``1`` calculates beta, and makes a final print out)
*>
*> The Mulliken population is also used to compute the bond
*> order as \f$ \mathit{BO}_{AB} = \sum_{a\in A} \sum_{b\in B} \mathit{DS}_{ab} \mathit{DS}_{ba} \f$.
*>
*> @param[in] NSYM    Number of irreducible representations
*> @param[in] NBAS    Number of basis functions per irred. rep.
*> @param[in] NAME    Center and function type label per basis function
*> @param[in] CMO     Orbital coefficients
*> @param[in] OCCN    Orbital occupations
*> @param[in] SMAT    Overlap matrix
*> @param[in] iCase   Type of run
*> @param[in] FullMlk Boolean for the type of print
*> @param[in] lSave   Boolean for saving on the Runfile
************************************************************************
      SUBROUTINE CHARGE(NSYM,NBAS,NAME,CMO,OCCN,SMAT,iCase,FullMlk,
     &           lSave)
c a temporary clone for CHARGE util
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "Molcas.fh"
      CHARACTER*(LENIN4) NAME(*)
      DIMENSION NBAS(NSYM),CMO(*),OCCN(*),SMAT(*)
      Logical FullMlk,lSave,Reduce_Prt
      External Reduce_Prt
*
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0
*
      If (iCase.ne.0) Then
         If (iPL.ge.2) Then
            Write (6,*)
            Call CollapseOutput(1,'   Molecular charges:')
            Write (6,'(3X,A)')    '   ------------------'
            Write (6,*)
         End If
      End If
      MXTYP=0
      Do iSym = 1, nSym
         MxTyp=MxTyp+nBas(iSym)
      End Do
      Call Get_iScalar('Unique atoms',nNUC)
      Call GetMem('QQ','ALLO','REAL',ipQQ,MXTYP*nNuc)
      Call FZero(Work(ipQQ),MXTYP*nNuc)
      Call CHARGE_(NSYM,NBAS,NAME,CMO,OCCN,SMAT,iCase,FullMlk,lSave,
     &             MXTYP,Work(ipQQ),nNuc)
      Call GetMem('QQ','FREE','REAL',ipQQ,MXTYP*nNuc)
*
      IF (iCase.ne.0) Then
        If (iPL.ge.2) Then
           Call CollapseOutput(0,'   Molecular charges:')
           Write(6,*)
        End If
      End If

      Return
      End

      SUBROUTINE CHARGE_(NSYM,NBAS,NAME,CMO,OCCN,SMAT,iCase,FullMlk,
     & lSave,MXTYP,QQ,nNuc)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "itmax.fh"
#include "Molcas.fh"
#include "angtp.fh"
#include "real.fh"
#include "WrkSpc.fh"
*
      CHARACTER*(LENIN4) NAME(*)
      DIMENSION NBAS(NSYM),CMO(*),OCCN(*),SMAT(*)
*
c      PARAMETER(MXTYP=maxbfn)
      CHARACTER*(LENIN) CNAME(MXATOM)
      CHARACTER*4 TNAME(MXTYP),TMP
      Character*4 TSwap(MXTYP)
      Character*4 TLbl(MXATOM)
      Character*2 AufBau(19)
      Integer ICNT(MXBAS),ITYP(MXBAS), nStab(MxAtom)
      Integer tNUC, NPBonds, AtomA, AtomB, nBas2
      Real*8 QQ(MXTYP,nNuc),QSUM(MXATOM)
      Real*8 Q2(MXATOM), QSUM_TOT(MXATOM)
      Logical FullMlk,lSave
c      Character*(LENIN) LblCnt(MxAtom)
      Character*(LENIN4) LblCnt4(MxAtom)
      Character*(LENIN4) Atom_A, Atom_B
      Save ipqswap
      Save ipDSswap
      Logical DMN_SpinAV
      Logical Do_SpinAV
      COMMON  / SPAVE_L  / Do_SpinAV
      COMMON  / SPAVE_I  / ip_DSc
      External Get_ProgName
      Character*100 ProgName, Get_ProgName
      Logical DoBond,Reduce_Prt
      External Reduce_Prt
      Data AufBau/'1s',
     &            '2s',          '2p',
     &            '3s',          '3p',
     &            '4s',     '3d','4p',
     &            '5s',     '4d','5p',
     &            '6s','4f','5d','6p',
     &            '7s','5f','6d','7p'/
*                                                                      *
************************************************************************
*                                                                      *
*---- Statement function
*
      Fac(i) = DBLE(nStab(i))/DBLE(nSym)
*                                                                      *
************************************************************************
*                                                                      *
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0

      xsg=-1.0d0
      DMN_SpinAV=.false.
      If (iCase.eq.0 .or. iCase.eq.1) Then
         If (Do_SpinAV) DMN_SpinAV=.true.
         If (iCase.eq.1) xsg=1.0d0
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
      Do i = 1, mxTyp
         TName(i)='    '
      End Do
*
*----------------------------------------------------------------------*
*     Get the name of the calling module.                              *
*     If CPFMCPF no bond analysis is done.                             *
*----------------------------------------------------------------------*
*
      ProgName=Get_ProgName()
      Call Upcase(ProgName)
      Call LeftAd(ProgName)
      iEnd = 1
 93   If (ProgName(iEnd:iEnd).ne.' ') Then
         iEnd=iEnd+1
         Go To 93
      End If

      DoBond = .True.
      If (ProgName(1:iEnd) .eq. 'CPF') DoBond = .False.
*
*----------------------------------------------------------------------*
*     Set the Mulliken Bond Order threshold for printout               *
*----------------------------------------------------------------------*
*
      BOThrs = 0.5D0
*
*----------------------------------------------------------------------*
*     GET THE TOTAL NUMBER OF BASIS FUNCTIONS AND CHECK LIMITS         *
*----------------------------------------------------------------------*
*
      NBAST=0
      Do I=1,NSYM
        NBAST=NBAST+NBAS(I)
      End Do
      IF(NBAST.GT.MXBAS) GOTO 991
*
*----------------------------------------------------------------------*
*     Find the list of unique center labels                            *
*----------------------------------------------------------------------*
*
      Call Get_cArray('Unique Atom Names',CNAME,LENIN*nNuc)
      Call Get_iArray('nStab',nStab,nNuc)

*
*----------------------------------------------------------------------*
*     Find the center label for each basis function                    *
*----------------------------------------------------------------------*
*
      Do I=1,NBAST
         ICNT(I)=-1
         Do J=1,NNUC
            If (NAME(I)(1:LENIN).EQ.CNAME(J)) ICNT(I)=J
         End Do
      End Do

*
*----------------------------------------------------------------------*
*     Find the type label for each basis function                      *
*----------------------------------------------------------------------*
*
      NXTYP=0
      Call ICopy(nBAST,0,0,ITYP,1)
      Do I=1,NBAST
        If (ICNT(I).lt.0) Go To 99  ! skip pseudo center
        Do J=1,NXTYP
          If (J.gt.MxTyp) Then
             Write (6,*) 'Charge: J.gt.MxTyp'
             Write (6,*) 'J=',J
             Write (6,*) 'MxTyp=',MxTyp
             Write (6,*) 'Increase MxType and recompile!'
             Call Abend()
          End If
          If (NAME(I)(LENIN1:LENIN4).EQ.TNAME(J)) Then
            ITYP(I)=J
            Go To 99
          End If
        End Do
        NXTYP=NXTYP+1
        TNAME(NXTYP)=NAME(I)(LENIN1:LENIN4)

*
        ITYP(I)=NXTYP
 99     Continue
      End Do

       lqSwap=NNUC+NNUC*NXTYP

       if(iCase.eq.0) then
c instead of printing charges we dump everything into a memory
c same with DS matrix

         Call GetMem('CHRG_SWP','ALLO','REAL',ipqSwap,lqSwap)

         If (DoBond) Then
           Call Allocate_Work(ipDSswap, (NBAST*NBAST))
         End If

       endif
*
*----------------------------------------------------------------------*
*     Do some trivial sorting of the type labels                       *
*----------------------------------------------------------------------*
*
*     Sort with respect to radial index
*
      ix=0
      jx=0
      Do i = 1, NxTyp-1
         ix = iChar(TNAME(i)(1:1))-iChar('1')+1
*        Put polarization and diffuse functions last
         if (tName(i)(1:1).eq.'*') ix = 99
         Do j = i+1, NxTyp
            jx = iChar(TNAME(j)(1:1))-iChar('1')+1
            if (tName(j)(1:1).eq.'*') jx = 99
            If (ix.gt.jx) Then
               iSwap=ix
               ix=jx
               jx=iSwap
               TMP=TNAME(i)
               TNAME(i)=TNAME(j)
               TNAME(j)=TMP
            End If
         End Do
      End Do
*
*     Sort with respect to angular index
*
      iAng = 0
      jAng = 0
      ix = 1
 666  iix = iChar(tName(ix)(1:1))
      jx = ix
      Do i = Min(ix+1,NxTyp), NxTyp
         If (iChar(tName(i)(1:1)).eq.iix) jx = i
      End Do
*
      Do i = ix, jx-1
         Do k = 0, iTabMx
            If (AngTp(k).eq.tName(i)(2:2)) iAng=k
         End Do
         Do j = i+1, jx
            Do l = 0, iTabMx
               If (AngTp(l).eq.tName(j)(2:2)) jAng=l
            End Do
            If (iAng.gt.jAng) Then
               iSwap=iAng
               iAng=jAng
               jAng=iSwap
               TMP=TNAME(i)
               TNAME(i)=TNAME(j)
               TNAME(j)=TMP
            End If
         End Do
      End Do
C     Write (*,*) ' Sorted n subrange'
C     Do i = ix, jx
C        Write (*,*) TName(i)
C     End Do
*
*     Now sort with respect to the magnetic index
*
      If (iAng.gt.9) Go To 888 ! skip sorting!
      iEnd = jx
      iStart = ix
 777  Do k = 0, iTabMx
         If (AngTp(k).eq.tName(iStart)(2:2)) iAng=k
      End Do
      jEnd = iStart
      Do i = Min(iStart+1,iEnd),iEnd
         If (tName(i)(2:2).eq.AngTp(iAng)) jEnd=i
      End Do
      if (jEnd.gt.mxtyp) then
        call abend
      end if
*
      i0 = iChar('1') - 1
*
      iM = 0
      jM = 0
      If (iAng.eq.1) Then
         Do i = iStart, jEnd-1
            If (tName(i)(3:3).eq.'x') iM = 1
            If (tName(i)(3:3).eq.'z') iM = 0
            If (tName(i)(3:3).eq.'y') iM =-1
            Do j = i+1, jEnd
               If (tName(j)(3:3).eq.'x') jM = 1
               If (tName(j)(3:3).eq.'z') jM = 0
               If (tName(j)(3:3).eq.'y') jM =-1
               If (jM.gt.iM) Then
                  iSwap=iM
                  iM=jM
                  jM=iSwap
                  TMP=TNAME(i)
                  TNAME(i)=TNAME(j)
                  TNAME(j)=TMP
               End If
            End Do
         End Do
      Else If (iAng.ge.2) Then
         Do i = iStart, jEnd-1
            iM = iChar(tName(i)(3:3)) - i0
            If (tName(i)(4:4).eq.'-') iM = -iM
            Do j = i+1, jEnd
               jM = iChar(tName(j)(3:3)) - i0
               If (tName(j)(4:4).eq.'-') jM = -jM
               If (jM.gt.iM) Then
                  iSwap=iM
                  iM=jM
                  jM=iSwap
                  TMP=TNAME(i)
                  TNAME(i)=TNAME(j)
                  TNAME(j)=TMP
               End If
            End Do
         End Do
      End If
*
      If (jEnd.ne.iEnd) Then
          iStart = jEnd + 1
          Go To 777
      End If
*
 888  Continue
      If (jx.ne.NxTyp) Then
         ix = jx + 1
         Go To 666
      End If
*
*     Sort according to AufBau
*
      iStart = 1
      Do iAB = 1, 19
         Do i = 1, NxTyp
            If (TName(i)(1:2).eq.AufBau(iAB)) Then
               TSwap(iStart) = TName(i)
               TName(i)='    '
               iStart = iStart + 1
            End If
         End Do
      End Do
      Do i = 1, NxTyp
         If (TName(i).ne.'    ') Then
            TSwap(iStart) = TName(i)
            TName(i)='    '
            iStart = iStart + 1
         End If
      End Do
      Do i = 1, NxTyp
         TName(i) = TSwap(i)
      End Do

*
      Do I=1,NBAST
         If (ICNT(I).lt.0) Go To 98  ! skip pseudo center
         Do J=1,NXTYP
            If (NAME(I)(LENIN1:LENIN4).eq.TNAME(J)) Then
               ITYP(I)=J
               Go To 98
             End If
        End Do
 98     Continue
      End Do
*
*----------------------------------------------------------------------*
*     Get the total number of atoms tNUC, regardless of symmetry       *
*----------------------------------------------------------------------*
*
      Call Get_iScalar('LP_nCenter', tNUC)
*                                                                      *
*----------------------------------------------------------------------*
*     Bond analysis initialization                                     *
*----------------------------------------------------------------------*
*
      If (DoBond) Then
*                                                                      *
*----------------------------------------------------------------------*
*     In case of symmetry we need the desymmetrization matrix,         *
*     for the bond order calculation only.                             *
*----------------------------------------------------------------------*
*
      If (nSym.gt.1) then
*
         Call Allocate_Work(ipP,NBAST**2)
         Call Allocate_Work(ipPInv,NBAST**2)
         Call Get_dArray('SM',Work(ipP),NBAST**2)
#ifdef _DEBUG_
         Call RecPrt('SM',' ',Work(ipP),NBAST,NBAST)
#endif
         Call MINV(Work(ipP),Work(ipPInv),ISING,DET,NBAST)
#ifdef _DEBUG_
         Call RecPrt('SMInv',' ',Work(ipPInv),NBAST,NBAST)
#endif
         Call DGeTMi(Work(ipPInv),NBAST,NBAST)
      End If
*
*     Pick up index array of which center a basis function belongs to.
*     If no symmetry, it is the same as ICNT(I).
*
      Call Allocate_iWork(ip_center,NBAST)
      Call Get_iArray('Center Index',iWork(ip_center),NBAST)
*                                                                      *
************************************************************************
*
*----------------------------------------------------------------------*
*     Initialize symmetric density D_tmp and overlap S_tmp matrices,   *
*     block D_blo and S_blo matrices (if symmetry)                     *
*     plus asymmetric D, S and DS matrices                             *
*----------------------------------------------------------------------*
*
      Call Allocate_Work(ipD_tmp, (NBAST*NBAST))
      Call Allocate_Work(ipS_tmp, (NBAST*NBAST))
      Call Allocate_Work(ipD    , (NBAST*NBAST))
      Call Allocate_Work(ipS    , (NBAST*NBAST))
      Call Allocate_Work(ipDS   , (NBAST*NBAST))
      Do I=1,(NBAST*NBAST)
          Work(ipD_tmp +I-1)=Zero
          Work(ipS_tmp +I-1)=Zero
          Work(ipD     +I-1)=Zero
          Work(ipS     +I-1)=Zero
          Work(ipDS    +I-1)=Zero
      End Do

      If (nSym.gt.1) then
          nBas2=0
          Do I=1, nsym
              nBas2=nBas2+nBas(i)*nBas(i)
          End Do
          Call Allocate_Work(ipD_blo, nBas2)
          Call Allocate_Work(ipS_blo, nBas2)
          Do I =1, nBas2
              Work(ipD_blo + I -1) = Zero
              Work(ipS_blo + I -1) = Zero
          End Do
      End If

*
*----------------------------------------------------------------------*
*     Find the center label for each atom, regardless of symmetry      *
*----------------------------------------------------------------------*
*

*     Just atom label. It's a double of the next one,
*     but someone could find it usefull in future

      Call Get_LblCnt_All(TLbl)

*     Atom labels plus symmetry generator

      Call Get_cArray('LP_L',LblCnt4,(LENIN4)*tNUC)
c      Do i=1,tNUC
c       LblCnt(i)(1:LENIN)=LblCnt4(i)(1:LENIN)
c      EndDo
*
*----------------------------------------------------------------------*
*     Initialize bond order vector                                     *
*----------------------------------------------------------------------*
*
      NPBonds = tNUC*(tNUC-1)/2
      Call Allocate_Work(ipBonds, NPBonds)
      Call FZero(Work(ipBonds),nPBonds)
*
*----------------------------------------------------------------------*
*     End of Bond analysis initialization                              *
*----------------------------------------------------------------------*
*
      End If
*                                                                      *
*
*----------------------------------------------------------------------*
*     Compute Mulliken atomic charges for each center and basis        *
*     function type                                                    *
*----------------------------------------------------------------------*
*

c      NDIM=NNUC*MXTYP
      NDIM=NXTYP*NNUC
      Call FZero(QQ,nDim)
      IB=0
      IS=0
      IMO=0
      Do ISYM=1,NSYM
        NB=NBAS(ISYM)
        IF ( NB.NE.0 ) THEN
          IMN=0
          Do MY=1,NB
            Do NY=1,MY
              IMN=IMN+1
              DMN=Zero
              ISMO=IMO
              Do IO=1,NB
                DMN=DMN+OCCN(IO+IB)*CMO(ISMO+MY)*CMO(ISMO+NY)
                ISMO=ISMO+NB
              End Do

              If (DMN_SpinAV) Then
                 DMN = DMN + xsg*Work(ip_DSc+(NY+IB-1)*NBAST+MY+IB-1)
              EndIf

              If (DoBond) then
*  Save the Density matrix element (my.ny) and (ny,my) in work(ipD_tmp)
*  Save the Overlap matrix element (my.ny) and (ny,my) in work(ipS_tmp)
              Work(ipD_tmp + (NY+IB-1) * NBAST + MY+IB -1)=DMN
              Work(ipD_tmp + (MY+IB-1) * NBAST + NY+IB -1)=DMN
              Work(ipS_tmp + (NY+IB-1) * NBAST + MY+IB -1)=SMAT(IMN+IS)
              Work(ipS_tmp + (MY+IB-1) * NBAST + NY+IB -1)=SMAT(IMN+IS)
              End If

              If ( MY.NE.NY ) DMN=Two*DMN
              MYNUC=ICNT(MY+IB)
              NYNUC=ICNT(NY+IB)
              MYTYP=ITYP(MY+IB)
              NYTYP=ITYP(NY+IB)
              If ( MY.EQ.NY ) Then
                TERM=SMAT(IMN+IS)*DMN
                If (MYNUC.gt.0) QQ(MYTYP,MYNUC)=QQ(MYTYP,MYNUC)+TERM
              Else
                TERM=Half*SMAT(IMN+IS)*DMN
                If (MYNUC.gt.0) QQ(MYTYP,MYNUC)=QQ(MYTYP,MYNUC)+TERM
                If (NYNUC.gt.0) QQ(NYTYP,NYNUC)=QQ(NYTYP,NYNUC)+TERM
              End If
           End Do
          End Do
          IB=IB+NB
          IS=IS+(NB+NB**2)/2
          IMO=IMO+NB**2
        End If
      End Do

*
*----------------------------------------------------------------------*
*     Density and overlap matrix handling for bond order               *
*----------------------------------------------------------------------*
*
      If (DoBond) Then

#ifdef _DEBUG_
      Call RecPrt('Density Matrix = ', ' ',Work(ipD_tmp), NBAST, NBAST)
      Call RecPrt('Overlap Matrix = ', ' ',Work(ipS_tmp), NBAST, NBAST)
      E=Zero
      Do I=1, NBAST
          Do J=1, NBAST
              E=E+ Work(ipD_tmp + (J-1) * NBAST + I - 1) *
     &             Work(ipS_tmp + (J-1) * NBAST + I - 1)
          End Do
      End Do
      Write(6,*)
      Write(6,*) 'Number of electrons as sum of D and S elements = ', E
#endif

*
*     In case of symmetry, we desymmetrize D and S through D_blo and S_blo
*
      If (nSym.gt.1) then
        iBlo = 0
        iSum = 0
        Do i = 1, NSYM
            If (nbas(i).ne.0) then
                Do j = 0, nbas(i) - 1
                    Do k = 0, nbas(i) - 1
                        Work(ipD_blo + iBlo) = Work(ipD_tmp +
     &                      (j+iSum)*NBAST + iSum + k)
                        Work(ipS_blo + iBlo) = Work(ipS_tmp +
     &                      (j+iSum)*NBAST + iSum + k)
                        iBlo = iBlo +1
                    End Do
                End Do
                iSum = iSum + nbas(i)
            End IF
        End Do

#ifdef _DEBUG_
        Write(6,*) 'D_blo = '
        Do i=1,nBas2
            Write(6,*) (Work(ipD_blo +I -1))
        End Do
        Write(6,*) 'S_blo = '
        Do i=1,nBas2
            Write(6,*) (Work(ipS_blo +I -1))
        End Do
#endif

*SVC: no need to have maximum number of basis sets here?
*       nScr=MXBAS*NBAST
        nScr=NBAST**2
        iSyLbl=1
        Call Allocate_Work(ipScr,nScr)
        Call Desymmetrize(Work(ipD_blo),nBas2,Work(ipScr),nScr,
     &                    Work(ipD),nBas,NBAST,Work(ipP),nSym,
     &                    iSyLbl)
        Call Free_Work(ipScr)

        Call Allocate_Work(ipScr,nScr)
        Call Desymmetrize(Work(ipS_blo),nBas2,Work(ipScr),nScr,
     &                    Work(ipS),nBas,NBAST,Work(ipPInv),nSym,
     &                    iSyLbl)
        Call Free_Work(ipScr)
*
*     Otherwise we simply copy D and S tmp into D and S
*
      Else
         call dcopy_(nBasT**2,Work(ipD_tmp),1,Work(ipD),1)
         call dcopy_(nBasT**2,Work(ipS_tmp),1,Work(ipS),1)
C        Do I=1,NBAST*NBAST
C           Work(ipD+I-1)=Work(ipD_tmp+I-1)
C           Work(ipS+I-1)=Work(ipS_tmp+I-1)
C        End Do
      End If

#ifdef _DEBUG_
      Write(6,*)'After Desymmetrization'
C     Call RecPrt('Density Matrix = ', ' ', Work(ipD), NBAST, NBAST)
C     Call RecPrt('Overlap Matrix = ', ' ', Work(ipS), NBAST, NBAST)
      Write (6,*) 'Dens=',DDot_(nBast**2,Work(ipD),1,Work(ipD),1),
     &                    DDot_(nBast**2,Work(ipD),1,One,0)
      Write (6,*) 'Ovrl=',DDot_(nBast**2,Work(ipS),1,Work(ipS),1),
     &                    DDot_(nBast**2,Work(ipS),1,One,0)
      Write (6,*) 'DO  =',DDot_(nBast**2,Work(ipS),1,Work(ipD),1)
      E=Zero
      Do I=1, NBAST
          Do J=1, NBAST
              E=E+ Work(ipD + (J-1) * NBAST + I - 1) *
     &             Work(ipS + (J-1) * NBAST + I - 1)
          End Do
      End Do
      Write(6,*)
      Write(6,*) 'Number of electrons as sum of D by S elements = ', E
#endif

*
*  Finally, we compute the DS matrix as product of D and S
*
      Call DGEMM_('N','N',
     &            NBAST,NBAST,NBAST,
     &            1.0d0,Work(ipD),NBAST,
     &            Work(ipS),NBAST,
     &            0.0d0,Work(ipDS),NBAST)
*
#ifdef _DEBUG_
      Call RecPrt('DS Matrix = ',' ', Work(ipDS),    NBAST, NBAST)
      E=Zero
      Do I=1,NBAST
          E=E+Work(ipDS + (I-1) * NBAST + I -1)
      End Do
      Write(6,*)
      Write(6,*) 'Number of electrons as sum of the DS diagonal = ', E
#endif
*
* in case of first call for UHF we dump everything only
*
      If (iCase.eq.0) Then
         Do I=1,NBAST
           Work(ipDSswap + I - 1) = Work(ipDS + I - 1)
         End Do
      End If
*
* in case of second call for UHF we add what dumped before
* and release swap memory
*
      If (iCase.eq.1) Then
        Do I=1,NBAST
          Work(ipDS + I - 1) = Work(ipDS + I - 1) +
     &                         Work(ipDSswap + I - 1)
        End Do
        Call Free_Work(ipDSswap)
*
#ifdef _DEBUG_
         Call RecPrt('DS Matrix = ',' ', Work(ipDS),    NBAST, NBAST)
         E=Zero
         Do I=1,NBAST
            E=E+Work(ipDS + (I-1) * NBAST + I -1)
         End Do
         Write(6,*)
         Write(6,*)'Number of electrons as sum of the DS diagonal = ', E
#endif
*
         End If
*
      End If
*
*----------------------------------------------------------------------*
*     Compute gross atomic charges                                     *
*----------------------------------------------------------------------*
*
      Do I=1,NNUC
        QSUMI=Zero
        Do J=1,NXTYP
          QSUMI=QSUMI+QQ(J,I)
        End Do
        QSUM(I)=QSUMI
      End Do
c if iCase=0, or 1 we need to put/get QSUM
        Do i=1,NNUC
         if(iCase.eq.0) Work(ipqSwap+i-1)=QSUM(I)
         if(iCase.eq.1) QSUM_TOT(I)=QSUM(I)+Work(ipqSwap+i-1)
         if(iCase.ge.2) QSUM_TOT(I)=QSUM(I)
        Enddo

*
*----------------------------------------------------------------------*
*     Pick up the nuclear charge                                       *
*----------------------------------------------------------------------*
*
      If (iCase.ne.0) then
         Call Allocate_Work(ip_Charge,nNuc)
         Call Get_dArray('Effective nuclear charge',
     &      Work(ip_Charge),nNuc)
         Do iNuc = 0, nNuc-1
           Work(ip_Charge+iNuc) = Work(ip_Charge+iNuc)
     &                          * DBLE(nSym / nStab(iNuc+1))
         End Do
         Call DaXpY_(nNuc,-One,QSUM_TOT,1,Work(ip_Charge),1)
         If(lSave) Call Put_dArray('Mulliken Charge',
     &      Work(ip_Charge),nNuc)
      End If
*
*----------------------------------------------------------------------*
*     Compute the 'Mulliken' Bond Order                                *
*----------------------------------------------------------------------*
*
      If ((DoBond) .AND. (tNUC.gt.1) .AND. (iCase.ge.1)) Then
*
#ifdef _DEBUG_
        Write (6,*) 'nPBonds,tNuc=',nPBonds,tNuc
        Do MY=1,NBAST
            AtomA=iWork(ip_center+MY-1)
            Write (6,*) 'AtomA,My=',AtomA,My
        End Do
#endif
        Do MY=1,NBAST
            AtomA=iWork(ip_center+MY-1)
            If (ICNT(MY).le.0) Go To 95    ! skip pseudo center
            Do NY=1,MY
              AtomB=iWork(ip_center+NY-1)
              If (ICNT(NY).le.0)      Go To 94  ! skip pseudo center
              If (AtomA.eq.AtomB)  Go To 94  ! same atom
*
              iPair = (Max(AtomA,AtomB)-1)
     &              * (Max(AtomA,AtomB)-2)/2
     &              +  Min(AtomA,AtomB)
              jPair=ipBonds-1+iPair
*
              Work(jPair)=Work(jPair) +
     &                    Work(ipDS+ (NY-1)*NBAST + MY-1) *
     &                    Work(ipDS+ (MY-1)*NBAST + NY-1)
*
#ifdef _DEBUG_
              Write(6,*)'Bond Number=',iPair
              Write(6,*)'Atom numbers = ',AtomA,AtomB
              Write(6,*)'Bond number = ', iPair,
     &                  'bond order = ',Work(jPair)
              Write(6,*) 'Work(ipDS+ (NY-1) * NBAST + MY-1) =',
     &                    Work(ipDS+ (NY-1) * NBAST + MY-1)
              Write(6,*) 'Work(ipDS+ (MY-1) * NBAST + NY-1) =',
     &                    Work(ipDS+ (MY-1) * NBAST + NY-1)
#endif
 94           Continue
           End Do
 95        Continue
        End Do

*     distant atoms could have negative bond order, set to zero

        Do I=1, NPBonds
            If (Work(ipBonds+I-1).lt.Zero) Work(ipBonds+I-1)=Zero
        End Do

#ifdef _DEBUG_
        Write(6,*)'Bond order vector'
        Call TriPrt('Bonds','(10F10.5)',Work(ipBonds),tNUC-1)
#endif

      End If
*
*----------------------------------------------------------------------*
*     Printout section                                                 *
*----------------------------------------------------------------------*
*
*
      If (iCase.eq.0) Then
c        first call for UHF, so just dump numbers to swap
         IEND=0
         ik=0
         Do IST=1,NNUC,6
            IEND=MIN(IEND+6,NNUC)
            Do IT=1,NXTYP
               Do j=IST,IEND
                  Work(ipqSwap+NNUC+ik)=QQ(IT,J)
                  ik=ik+1
               End Do
            End Do
         End Do
      End If
*
      If (iCase.eq.1.and.iPL.ge.2) Then
c second call, make a real print out
         If (FullMlk) Then
            Write(6,'(6X,A)')
     &      'Mulliken charges per centre and basis function type'
            Write(6,'(6X,A)')
     &      '---------------------------------------------------'
         Else
            Write(6,'(6X,A)')
     &      'Mulliken charges per centre'
            Write(6,'(6X,A)')
     &      '---------------------------'
         End If
         IEND=0
         ik=0
         ikk=0
         Do IST=1,nNuc,6
            Fact = DBLE(nStab(ist))/DBLE(nSym)
            IEND=MIN(IEND+6,nNuc)
            Write(6,*)
            Write(6,'(11X,6(14X,A,4X))')
     &         (CNAME(I),I=IST,IEND)
            Write(6,'(11X,6(A12,A12))')
     &         (' alpha','  beta',I=IST,IEND)
            Do IT=1,NXTYP
               Do J=IST,IEND
                  Q2(J)=Work(ipqSwap+NNUC+ik)
                  ik=ik+1
               End Do
               If (FullMlk) then
                  Write(6,'(6X,A4,1X,12F12.4)')TNAME(IT),
     &                 (Fac(j)*Q2(J),Fac(j)*QQ(IT,J), J=IST,IEND)
               End If
            End Do
*
            Do J=IST,IEND
               Q2(J)=Work(ipqSwap+ikk)
               ikk=ikk+1
            End Do
*
            Write(6,'(6X,A,12F12.4)')'Total',
     &         (Fac(i)*Q2(I),Fac(i)*QSUM(I),I=IST,IEND)
            Write(6,'(6X,A,6(6X,F12.4,6X))')'Total',
     &         (Fac(i)*(Q2(I)+QSUM(I)),I=IST,IEND)
            Write(6,*)
            Write(6,'(6X,A,6(5X,F12.4,7X))')'Charge',
     &         (Fac(i)*Work(ip_Charge+I-1),I=IST,IEND)
         End Do
         Write(6,*)
         Write(6,'(6X,A,F12.6)') 'Total electronic charge=',
     &                 DDot_(nNuc,One,0,QSum_TOT,1)
         Write(6,*)
         Write(6,'(6X,A,F12.6)') 'Total            charge=',
     &                 DDot_(nNuc,One,0,Work(ip_Charge),1)
*
      End If
      If (iCase.eq.1) Then
         Call Free_Work(ip_Charge)
         Call GetMem('CHRG_SWP','FREE','REAL',ipqSwap,lqSwap)
      End If
*
      If ((iCase.eq.2.and.iPL.ge.2) .or. (iCase.eq.3.and.iPL.ge.2)) Then
c icase=2 for usual mulliken, =2 for spin population.
*
         If (FullMlk) Then
            If (iCase.eq.2) then
             Write(6,'(6X,A)')
     &   'Mulliken charges per centre and basis function type'
            else
             Write(6,'(6X,A)')
     &   'Mulliken spin population per centre and basis function type'
            end if
            Write(6,'(6X,A)')
     &         '---------------------------------------------------'
         Else
            If (iCase.eq.2) then
             Write(6,'(6X,A)')'Mulliken charges per centre'
            else
             Write(6,'(6X,A)')
     &         'Mulliken spin population per centre'
            end if
            Write(6,'(6X,A)')
     &         '---------------------------'
         End If
*
         IEND=0
         Do IST=1,nNuc,12
            Fact = DBLE(nStab(ist))/DBLE(nSym)
            IEND=MIN(IEND+12,nNuc)
            Write(6,*)
            Write(6,'(11X,12(2X,A))') (CNAME(I),I=IST,IEND)
            If (FullMlk) then
               Do IT=1,NXTYP
                  Write(6,'(6X,A4,1X,12F8.4)')TNAME(IT),
     &              (Fac(j)*QQ(IT,J),J=IST,IEND)
               End Do
            endIf
            Write(6,'(6X,A,12F8.4)')'Total',(Fac(i)*QSUM(I),I=IST,IEND)
            If (iCase.ne.3) Then
              Write(6,*)
              Write(6,'(6X,A,12F8.4)')'N-E  ',
     &             (Fac(i)*Work(ip_Charge+I-1),I=IST,IEND)
            End If
         End Do
         if(iCase.eq.3) then
           Write(6,*)
           Write(6,'(6X,A,F12.6)') 'Total electronic spin=',
     &                 DDot_(nNuc,One,0,QSum,1)
         else
           Write(6,*)
           Write(6,'(6X,A,F12.6)') 'Total electronic charge=',
     &                 DDot_(nNuc,One,0,QSum,1)
           Write(6,*)
           TCh=DDot_(nNuc,One,0,Work(ip_Charge),1)
           Write(6,'(6X,A,F12.6)') 'Total            charge=',
     &                    DDot_(nNuc,One,0,Work(ip_Charge),1)
         Call xml_dDump('FormalCharge','Total charge','a.u',0,TCh,1,1)
         End If
      End If
      If (iCase.ge.2) Then
         Call Free_Work(ip_Charge)
      EndIf

*  Mulliken bond order print

      If (iPL.le.2) Go To 9999
      If((iCase.ge.1) .AND. (iCase.le.2) .AND. (tNUC.gt.1)
     &    .AND. (DoBond)) then
        Write(6,'(6X,A)')
     &   'Mulliken Bond Order analysis'
        Write(6,'(6X,A)')
     &   '----------------------------'
        Write(6,'(6X,A,F5.3,A)')
     &   'Only bonds with order larger than ',BOThrs,' are printed'
        Write(6,*)
        If (nSym.gt.1) then
            Write(6,'(6X,A)')
     &      'Atom A -   Generator  Atom B -   Generator  Bond Order'
            Do I=1, tNUC-1
             Do J=I+1, tNUC
                iPair = (J-1)*(J-2)/2 + I
                BO = Work(ipBonds -1 + iPair)
              If (BO .ge. BOThrs) then
               Atom_A =  LblCnt4(I)
*    &          TLbl(I)//' -  '//LblCnt(I)(Index(LblCnt(I),':')+1:)
               Atom_B =   LblCnt4(J)
*    &          TLbl(J)//' -  '//LblCnt(J)(Index(LblCnt(J),':')+1:)
               Write(6,'(8X,A,9X,A,9X,F7.3)')
     &          Atom_A, Atom_B, BO
              End If
             End Do
            End Do
        Else
            Write(6,'(6X,A)')'Atom A  Atom B  Bond Order'
            Do I=1, tNUC-1
              Do J=I+1, tNUC
                 iPair = (J-1)*(J-2)/2 + I
                 BO = Work(ipBonds -1 + iPair)
                If (BO .ge. BOThrs) then
                    Write(6,'(8X,A,4X,A,4X,F7.3)')
     &               LblCnt4(I),LblCnt4(J),BO
*    &               TLbl(I),TLbl(J),BO
                End If
              End Do
            End Do
        End If
        Write(6,*)
      End If

 9999 Continue
      If (DoBond) Then
          If (nSym.gt.1) then
              Call Free_Work(ipP)
              Call Free_Work(ipPInv)
              Call Free_Work(ipD_blo)
              Call Free_Work(ipS_blo)
          End If
          Call Free_iWork(ip_center)
          Call Free_Work(ipD_tmp)
          Call Free_Work(ipS_tmp)
          Call Free_Work(ipD)
          Call Free_Work(ipS)
          Call Free_Work(ipDS)
          Call Free_Work(ipBonds)
      End If
*
      Return
*
*----------------------------------------------------------------------*
*     Error Exits                                                      *
*----------------------------------------------------------------------*
*
991   Write(6,'(/6X,A)')
     &'The number of basis functions exceeds the present limit'
      Call Abend
*992   Write(6,'(/6X,A)')
*     &'The number of basis functions exceeds the present limit'
*      Call Abend
*993   Write(6,'(/6X,A)')
*     &'Warning: Total charge is not equal to number of electrons'
*      Call Abend
      End
