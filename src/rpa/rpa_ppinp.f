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
* Copyright (C) 2013, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine RPA_PPInp()
C
C     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
C
C     Input postprocessing.
C
      Implicit None
#include "rpa_config.fh"
#include "rpa_data.fh"
#include "WrkSpc.fh"

      Character*9 SecNam
      Parameter (SecNam='RPA_PPInp')

      ! set RPAModel
      If (dRPA) Then
         If (SOSEX) Then
            RPAModel='SOSX@'//Reference(1:3)
         Else
            RPAModel='dRPA@'//Reference(1:3)
         End If
      Else
         ! this should never happen
         Call RPA_Warn(3,SecNam//': internal error [RPAModel]')
         RPAModel='None@Non'
      End If

      ! freeze orbitals
      Call RPA_Freezer()

      ! print config after input processing
      Call RPA_PrInp()

      End
************************************************************************
      Subroutine RPA_PrInp()
C
C     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
C
C     Print RPA configuration after input processing.
C
      Implicit None
#include "rpa_config.fh"
#include "rpa_data.fh"
#include "WrkSpc.fh"

      Character*9 SecNam
      Parameter (SecNam='RPA_PrInp')
      Integer lPaper
      Parameter (lPaper=132)

      Integer  RPA_iUHF, RPA_LENIN8
      External RPA_iUHF, RPA_LENIN8

      Character*3 lIrrep(8)
      Character*7 spin(2)
      Character*8 Fmt1, Fmt2
      Character*13 orbitals
      Character*120 Line,BlLine,StLine

      Integer iUHF
      Integer lLine
      Integer nLine
      Integer l_orbitals
      Integer i, j, k
      Integer left
      Integer iSym
      Integer LENIN8
      Integer nB
      Integer ip_Name, l_Name
      Integer iCount

      Real*8 Dummy(1)

      Integer p, q
      Real*8 epsi, epsa
      epsi(p,q)=Work(ip_OccEn(q)-1+p)
      epsa(p,q)=Work(ip_VirEn(q)-1+p)

      ! set restricted(1)/unrestricted(2)
      iUHF=RPA_iUHF()

      ! set labels
      If (iUHF.eq.1) Then
         orbitals='orbitals'
         l_orbitals=8
         spin(1)=' '
         spin(2)=' '
      Else If (iUHF.eq.2) Then
         orbitals='spin-orbitals'
         l_orbitals=13
         spin(1)='(alpha)'
         spin(2)='(beta)'
      Else
         Write(6,'(A,I6)') 'iUHF=',iUHF
         Call RPA_Warn(3,SecNam//': iUHF error')
         orbitals=' '
         l_orbitals=1
         spin(1)=' '
         spin(2)=' '
      End If

      ! set irrep names
      Call Get_cArray('Irreps',lIrrep,24)
      Do iSym=1,nSym
         Call RightAd(lIrrep(iSym))
      End Do

      ! init blank line and "star" line
      lLine=len(Line)
      Do i=1,lLine
        BlLine(i:i)=' '
        StLine(i:i)='*'
      End Do
      left=(lPaper-lLine)/2
      Write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'

      ! print title from input
      If (nTitle.gt.0) Then
         Write(6,*)
         nLine=nTitle+5
         Do i=1,nLine
           Line=BlLine
           If (i.eq.1 .or. i.eq.nLine) Line=StLine
           If (i.eq.3) Line='Title:'
           If (i.ge.4 .and. i.le.nLine-2) Line=Title(i-3)
           Call Center(Line)
           Write(6,Fmt1) '*'//Line//'*'
         End Do
         Write(6,*)
      End If

      ! print coordinates of the molecule
      If (iPrint.ge.2) Then
         Call PrCoor()
      End If

      ! print orbital info
      If (iPrint.ge.2) Then
         Write(6,*)
         Write(6,Fmt2//'A,2(1X,A))') Reference,'reference',orbitals
         j=len(Reference)+11+l_orbitals
         Write(6,Fmt2//'80A1)') ('-',i=1,j)
         If (Reference(2:3).eq.'KS') Then
            Write(6,Fmt2//'A,A)') 'DFT functional: ',DFTFunctional
         End If
         Write(6,*)
         Write(6,Fmt2//'A,T47,8I4)') 'Symmetry species',
     *                               (iSym,iSym=1,nSym)
         Write(6,Fmt2//'A,T47,8(1X,A))') '                ',
     *                               (lIrrep(iSym),iSym=1,nSym)
         Write(6,Fmt2//'A,T47,8I4)') 'Number of basis functions',
     *                               (nBas(iSym),iSym=1,nSym)
         Write(6,Fmt2//'A,T47,8I4)') 'Number of orbitals',
     *                               (nOrb(iSym),iSym=1,nSym)
         Do k=1,iUHF
            Write(6,Fmt2//'A,2(1X,A),T47,8I4)')
     *      'Frozen occupied',orbitals,spin(k),
     *      (nFro(iSym,k),iSym=1,nSym)
         End Do
         Do k=1,iUHF
            Write(6,Fmt2//'A,2(1X,A),T47,8I4)')
     *      'Active occupied',orbitals,spin(k),
     *      (nOcc(iSym,k),iSym=1,nSym)
         End Do
         Do k=1,iUHF
            Write(6,Fmt2//'A,2(1X,A),T47,8I4)')
     *      'Active virtual',orbitals,spin(k),
     *      (nVir(iSym,k),iSym=1,nSym)
         End Do
         Do k=1,iUHF
            Write(6,Fmt2//'A,2(1X,A),T47,8I4)')
     *      'Deleted virtual',orbitals,spin(k),
     *      (nDel(iSym,k),iSym=1,nSym)
         End Do
      End If

      ! print orbital energies
      If (iPrint.ge.2) Then
         iCount=0
         Do k=1,iUHF
            Do iSym=1,nSym
               iCount=iCount+nFro(iSym,k)
            End Do
         End Do
         If (iCount.gt.0) Then
            Write(6,*)
            Write(6,*)
            Write(6,Fmt2//'A,1X,A,T47)')
     *      'Energies of the frozen occupied',orbitals
            Do k=1,iUHF
               i=0
               Do iSym=1,nSym
                  If (nFro(iSym,k).gt.0) then
                     Write(6,*)
                     Write(6,Fmt2//'A,I2,2(1X,A),(T40,5F14.6))')
     *               'symmetry species',iSym,lIrrep(iSym),spin(k),
     *               (epsi(i+j,k),j=1,nFro(iSym,k))
                     i=i+nFro(iSym,k)+nOcc(iSym,k)
                  End If
               End Do
            End Do
         End If
         Write(6,*)
         Write(6,*)
         Write(6,Fmt2//'A,1X,A,T47)')
     *   'Energies of the active occupied',orbitals
         Do k=1,iUHF
            i=0
            Do iSym=1,nSym
               If (nOcc(iSym,k).gt.0) then
                  Write(6,*)
                  Write(6,Fmt2//'A,I2,2(1X,A),(T40,5F14.6))')
     *            'symmetry species',iSym,lIrrep(iSym),spin(k),
     *            (epsi(i+nFro(iSym,k)+j,k),j=1,nOcc(iSym,k))
                  i=i+nFro(iSym,k)+nOcc(iSym,k)
               End If
            End Do
         End Do
         Write(6,*)
         Write(6,*)
         Write(6,Fmt2//'A,1X,A,T47)')
     *   'Energies of the active virtual',orbitals
         Do k=1,iUHF
            i=0
            Do iSym=1,nSym
               If (nVir(iSym,k).gt.0 ) then
                  Write(6,*)
                  Write(6,Fmt2//'A,I2,2(1X,A),(T40,5F14.6))')
     *            'symmetry species',iSym,lIrrep(iSym),spin(k),
     *            (epsa(i+j,k),j=1,nVir(iSym,k))
                  i=i+nVir(iSym,k)+nDel(iSym,k)
               End If
            End Do
         End Do
         iCount=0
         Do k=1,iUHF
            Do iSym=1,nSym
               iCount=iCount+nDel(iSym,k)
            End Do
         End Do
         If (iCount.gt.0) Then
            Write(6,*)
            Write(6,*)
            Write(6,Fmt2//'A,1X,A,T47)')
     *      'Energies of the deleted virtual',orbitals
            Do k=1,iUHF
               i=0
               Do iSym=1,nSym
                  If (nDel(iSym,k).gt.0 ) then
                     Write(6,*)
                     Write(6,Fmt2//'A,I2,2(1X,A),(T40,5F14.6))')
     *               'symmetry species',iSym,lIrrep(iSym),spin(k),
     *               (epsa(i+nVir(iSym,k)+j,k),j=1,nDel(iSym,k))
                     i=i+nVir(iSym,k)+nDel(iSym,k)
                  End If
               End Do
            End Do
         End If
      End If

      ! print orbitals
      If (iPrint.ge.2) Then
         LENIN8=RPA_LENIN8()
         nB=nBas(1)
         Do iSym=2,nSym
            nB=nB+nBas(iSym)
         End Do
         l_Name=LENIN8*nB
         Call GetMem('Name','Allo','Char',ip_Name,l_Name)
         Call Get_cArray('Unique Basis Names',cWork(ip_Name),LENIN8*nB)
         Do k=1,iUHF
            Call PriMO(Reference//' '//orbitals//' '//spin(k),
     *                 .false.,.true.,-9.9d9,9.9d9,
     *                 nSym,nBas,nOrb,cWork(ip_Name),
     *                 Work(ip_EMO(k)),Dummy,Work(ip_CMO(k)),-1)
         End Do
         Call GetMem('Name','Free','Char',ip_Name,l_Name)
      End If

      ! flush output buffer
      Call xFlush(6)

      End
