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
      SubRoutine RdAB()
      use Arrays, only: CMO
      use stdalloc, only: mma_allocate
      use MCLR_Data, only: ChDisp, lDisp
      use input_mclr, only: nSym,nIsh,nBas,nOrb,Perturbation,
     &                      McKinley,iMethod,nIsh,ntITri,ntISqr,ntBSqr,
     &                      nDisp,PT2,ESCF,nDel,ntPert,ntIsh
      Implicit None

      Character(LEN=8) Label
      Integer iRC, iOpt, Length, iSym, iS, iDum
*
*
      Perturbation='NONE'
*
      If (MCKINLEY) Then
       LABEL='TDISP'
       iRc=-1
       iOpt=0
       Call RdMck(iRC,iOpt,Label,idum,ntpert,idum)
       If (iRC.ne.0) Then
          Write (6,*) 'RdAB: Error reading MCKINT'
          Write (6,'(A,A)') 'Label=',Label
          Call Abend()
       End If
       LABEL='PERT'
       iRc=-1
       iOpt=0
       Call cRdMck(iRC,iOpt,Label,idum,Perturbation,idum)
       If (iRC.ne.0) Then
          Write (6,*) 'RdAB: Error reading MCKINT'
          Write (6,'(A,A)') 'Label=',Label
          Call Abend()
       End If
      End If
*

      If (iMethod.eq.1) Then
         Call Get_dScalar('Last energy',ESCF)
         Call Get_iArray('nIsh',nIsh,nSym)
         Call Get_iArray('nDel',nDel,nSym)
*----------------------------------------------------------------------*
*     Precompute the total sum of variables and size of matrices       *
*----------------------------------------------------------------------*
         ntIsh=0
         ntItri=0
         ntIsqr=0
         ntBsqr=0
         Length=0
         Do iSym=1,nSym
            ntIsh=ntIsh+nIsh(iSym)
            ntItri=ntItri+nIsh(iSym)*(nIsh(iSym)+1)/2
            ntIsqr=ntIsqr+nIsh(iSym)*nIsh(iSym)
            ntbSQR=ntbsqr+nbas(isym)**2
            norb(isym)=nbas(isym)-ndel(isym)
            Length=Length+nBas(iSym)*nOrb(iSym)
         End Do
         Call mma_allocate(CMO,Length,Label='CMO')
         Call Get_dArray_chk('Last orbitals',CMO,Length)
      End If

*
      If (MCKINLEY) Then
         Label='ldisp   '
         iRc=-1
         iOpt=0
         Call RdMck(iRc,iOpt,label,idum,ldisp,idum)
         If (iRC.ne.0) Then
            Write (6,*) 'RdAB: Error reading MCKINT'
            Write (6,'(A,A)') 'Label=',Label
            Call Abend
         End If
         nDisp=0
         Do iS=1,nSym
            nDisp=nDisp+lDisp(iS)
         End Do
         If (ndisp.ne.0) Then
            Label='Chdisp  '
            iRc=-1
            iOpt=0
            Call cRdmck(iRc,iOpt,label,idum,ChDisp(1),idum)
            If (iRC.ne.0) Then
               Write (6,*) 'RdAB: Error reading MCKINT'
               Write (6,'(A,A)') 'Label=',Label
               Call Abend
            End If
         End if
      End If
*
      If (PT2) Then
         Call icopy(nsym,[0],0,ldisp,1)
         ldisp(1)=1
      End If
*
      Return
      End
