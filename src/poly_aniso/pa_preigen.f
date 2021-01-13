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
      Subroutine pa_preigen( exch, lmax, ibas, Dipol, AnisoLines1,
     &                       AnisoLines3, AnisoLines9, KE, DM_exchange,
     &                       JITO_exchange,
     &                       WLIN, WDIP, WKEX, WDMO, WITO, W, Z, iPrint)
c this function prints the energies and eigenvectors of the interaction Hamiltonians:
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)          :: exch
      Integer, intent(in)          :: lmax
      Integer, intent(in)          :: iPrint
      Integer, intent(in)          :: ibas(exch,lmax)
!     eigenstates of all Lines interactions: 1+3+9
      Real(kind=8), intent(in)    :: WLIN(exch)
      Real(kind=8), intent(in)    :: WDIP(exch)
      Real(kind=8), intent(in)    :: WKEX(exch)
      Real(kind=8), intent(in)    :: WDMO(exch)
      Real(kind=8), intent(in)    :: WITO(exch)
      Real(kind=8), intent(in)    :: W(exch)
      Complex(kind=8), intent(in) :: Z(exch,exch)
      Logical, intent(in)          :: AnisoLines1
      Logical, intent(in)          :: AnisoLines3
      Logical, intent(in)          :: AnisoLines9
      Logical, intent(in)          :: Dipol
      Logical, intent(in)          :: DM_exchange
      Logical, intent(in)          :: KE
      Logical, intent(in)          :: JITO_exchange
c local variables
      Integer           :: i,m,ipar,nf,j,jEnd,iss
      Character(len=80) :: fmtline
      Logical           :: dbg

      Write(6,*)
      Write(6,'(100a)') (('%'),i=1,100)
      Write(6,'(10x,a)') 'EigenValues of the Magnetic Interaction'
      Write(6,'(100a)') (('%'),i=1,100)

      dbg=.false.

      If(dbg) Write(6,*) 'pa_preigen:   AnisoLines1 =',AnisoLines1
      If(dbg) Write(6,*) 'pa_preigen:   AnisoLines3 =',AnisoLines3
      If(dbg) Write(6,*) 'pa_preigen:   AnisoLines9 =',AnisoLines9
      If(dbg) Write(6,*) 'pa_preigen:   DM_exchange =',DM_exchange
      If(dbg) Write(6,*) 'pa_preigen:         Dipol =',Dipol
      If(dbg) Write(6,*) 'pa_preigen:            KE =',KE
      If(dbg) Write(6,*) 'pa_preigen: JITO_exchange =',JITO_exchange
      If(dbg) Write(6,*) WDMO
      If(dbg) Write(6,*) WITO

!----------------------------------------------------------------------!
      If ((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &        Dipol.and.KE.and..not.JITO_exchange) Then
        Write(6,'(6a)') 'Coupled|',
     &                  '     Lines Model      |',
     &                  ' Dipole-Dipole Inter. |',
     &                  '   Kinetic Exchange   |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'

        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WLIN(i),'|',WDIP(i),'|',
     &                  WKEX(i),'|',W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'



!----------------------------------------------------------------------!
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &         Dipol.and..not.KE.and..not.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  '     Lines Model      |',
     &                  ' Dipole-Dipole Inter. |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WLIN(i),'|',WDIP(i),'|',
     &                     W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'



!----------------------------------------------------------------------!
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &        .not.Dipol.and.KE.and..not.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  '     Lines Model      |',
     &                  '   Kinetic Exchange   |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WLIN(i),'|',WKEX(i),'|',
     &                     W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'


!----------------------------------------------------------------------!
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &             Dipol.and.KE.and..not.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  ' Dipole-Dipole Inter. |',
     &                  '   Kinetic Exchange   |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WDIP(i),'|',WKEX(i),'|',
     &                     W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'


!----------------------------------------------------------------------!
      Else If( (AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and..not.Dipol.and..not.KE.and..not.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  '     Lines Model      |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WLIN(i),'|',W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'


!----------------------------------------------------------------------!
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &             .and. Dipol.and..not.KE.and..not.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  ' Dipole-Dipole Inter. |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WDIP(i),'|',W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'


!----------------------------------------------------------------------!
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and..not.Dipol.and.KE.and..not.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  '   Kinetic Exchange   |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WKEX(i),'|',W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'







! JITO_exchange active below:
!----------------------------------------------------------------------!
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and..not.Dipol.and.KE.and.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  '   Kinetic Exchange   |',
     &                  '     ITO Exchange     |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WKEX(i),'|',WITO(i),'|', W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'



!----------------------------------------------------------------------!
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and..not.Dipol.and..not.KE.and.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  '     ITO Exchange     |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WITO(i),'|', W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'



!----------------------------------------------------------------------!
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and.Dipol.and..not.KE.and.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  ' Dipole-Dipole Inter. |',
     &                  '     ITO Exchange     |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WDIP(i),'|',WITO(i),'|',
     &                     W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'



!----------------------------------------------------------------------!
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and.Dipol.and.KE.and.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  ' Dipole-Dipole Inter. |',
     &                  '   Kinetic Exchange   |',
     &                  '     ITO Exchange     |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WDIP(i),'|',WKEX(i),'|',WITO(i),'|',
     &                     W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'




!----------------------------------------------------------------------!
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and..not.Dipol.and..not.KE.and.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  '     Lines Model      |',
     &                  '     ITO Exchange     |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WLIN(i),'|',WITO(i),'|',
     &                     W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'



!----------------------------------------------------------------------!
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and.Dipol.and..not.KE.and.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  '     Lines Model      |',
     &                  ' Dipole-Dipole Inter. |',
     &                  '     ITO Exchange     |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WLIN(i),'|',WDIP(i),'|',WITO(i),'|',
     &                     W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'



!----------------------------------------------------------------------!
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and..not.Dipol.and.KE.and.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  '     Lines Model      |',
     &                  '   Kinetic Exchange   |',
     &                  '     ITO Exchange     |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  WLIN(i),'|',WKEX(i),'|',WITO(i),'|',
     &                     W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'



!----------------------------------------------------------------------!
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and.Dipol.and.KE.and.JITO_exchange ) Then
        Write(6,'(6a)') 'Coupled|',
     &                  '     Lines Model      |',
     &                  ' Dipole-Dipole Inter. |',
     &                  '   Kinetic Exchange   |',
     &                  '     ITO Exchange     |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,6(f20.12,2x,A))') i,'|',
     &                  WLIN(i),'|',WDIP(i),'|',WKEX(i),'|',WITO(i),'|',
     &                     W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'


!----------------------------------------------------------------------!
      Else

        Write(6,'(6a)') 'Coupled|',
     &                  '  Magnetic Exchange   |',
     &                  '         Total Magnetic Interaction          |'
        Write(6,'(6a)') 'State  |',
     &                  '  absolute val., cm-1 |',
     &                  '  absolute val., cm-1 |  relative val., cm-1 |'
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
        Do i=1,exch
          Write(6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',
     &                  W(i),'|',W(i),'|',W(i)-W(1),'|'
        End Do
        Write(6,'(6a)') '-------|',
     &                  '----------------------|',
     &                  '---------------------------------------------|'
      End If











ccc  eigenvectors
      If(iPrint>=3) Then

      Write(6,*)
      Write(6,'(100a)') (('%'),i=1,100)
      Write(6,'(10x,a)') 'EigenVectors of the Total Magnetic '//
     &                   'Interaction'
      Write(6,'(100a)') (('%'),i=1,100)


      If(exch>=4) Then
         Write(6,'(A,16A)') '-',('---',m=1,lmax),'|',
     &                    ('---------------------------|',i=1,4)
         Write(fmtline,'(A,i2,A)') '(1x,',lmax,'A,39x,A,38x,A)'
         Write(6,fmtline) ('   ',m=1,lmax),
     &                 'eigenvectors of the exchange matrix','|'
      Else If(exch==3) Then
         Write(6,'(A,16A)') '-',('---',m=1,lmax),'|',
     &                    ('---------------------------|',i=1,3)
         Write(fmtline,'(A,i2,A)') '(1x,',lmax,'A,25x,A,24x,A)'
         Write(6,fmtline) ('   ',m=1,lmax),
     &                 'eigenvectors of the exchange matrix','|'
      Else If(exch==2) Then
         Write(6,'(A,16A)') '-',('---',m=1,lmax),'|',
     &                    ('---------------------------|',i=1,2)
         Write(fmtline,'(A,i2,A)') '(1x,',lmax,'A,11x,A,10x,A)'
         Write(6,fmtline) ('   ',m=1,lmax),
     &                 'eigenvectors of the exchange matrix','|'
      Else If(exch==1) Then
         Write(6,'(A,16A)') '-',('---',m=1,lmax),'|',
     &                     '---------------------------|'
         Write(fmtline,'(A,i2,A)') '(1x,',lmax,'A,5x,A,3x,A)'
         Write(6,fmtline) ('   ',m=1,lmax),
     &                 'exchange eigenvector','|'
      End If


      Do J=1,exch,4
        jEnd=MIN(exch,J+3)
        Write(6,'(A,16A)') '-',('---',m=1,lmax),'|',
     &       ('---------------------------|',i=j,jEnd)


        If (lmax.eq.1) Then

          Write(6,'(A,5A)') 'Exch|',
     &                     ('      exchange state       |',i=j,jEnd)

          Write(6,'(A,6A)') ' bas|',
     &                     ('                           |',i=j,jEnd)
          Write(6,'(A,6(a,i5,a))') '    |',
     &                     ('         ',i,'             |',i=j,jEnd)
          Write(fmtline,'(A,i2,A)')
     &                     '(A,',lmax,'i3,a,6(a,f19.12,2x,a))'
          Write(6,fmtline) ' ',(i,i=1,lmax),'|',
     &                     ('  E =',w(i)-w(1),' |',i=j,jEnd)

        Else If(lmax.eq.2) Then

          Write(6,'(A,6A)') 'Exch.  |',
     &                     ('      exchange state       |',i=j,jEnd)

          Write(6,'(A,6A)') 'basis  |',
     &                     ('                           |',i=j,jEnd)
          Write(6,'(A,6(a,i5,a))') 'on site|',
     &                     ('         ',i,'             |',i=j,jEnd)
          Write(fmtline,'(A,i2,A)')
     &                    '(A,',lmax,'i3,a,6(a,f19.12,2x,a))'
          Write(6,fmtline) ' ',(i,i=1,lmax),'|',
     &                     ('  E =',w(i)-w(1),' |',i=j,jEnd)

        Else If(lmax.eq.3) Then

          Write(fmtline,'(A,i2,A,i2,A)') '(A,A,4A)'
          Write(6,fmtline) ' Exchange |',
     &                     ('      exchange state       |',i=j,jEnd)

          Write(6,'(A,6A)') ' basis on |',
     &                     ('                           |',i=j,jEnd)
          Write(6,'(A,6(a,i5,a))') '    sites |',
     &                     ('         ',i,'             |',i=j,jEnd)
          Write(fmtline,'(A,i2,A)')
     &                     '(A,',lmax,'i3,a,6(a,f19.12,2x,a))'
          Write(6,fmtline) ' ',(i,i=1,lmax),'|',
     &                     ('  E =',w(i)-w(1),' |',i=j,jEnd)

        Else

          ipar=mod(lmax,2)
            nf=(3*lmax-8+ipar)/2
          Write(fmtline,'(A,i2,A,i2,A)') '(',nf,'x,A,',
     &                               nf,'x,A,4A)'
          Write(6,fmtline) ' Exchange','|',
     &                     ('      exchange state       |',i=j,jEnd)
          Write(6,fmtline) ' basis on','|',
     &                     ('                           |',i=j,jEnd)
          Write(fmtline,'(A,i2,A,i2,A)') '(',nf,'x,A,',
     &                      nf,'x,A,4(a,i5,a))'
          Write(6,fmtline) '    sites','|',
     &                     ('         ',i,'             |',i=j,jEnd)
          Write(fmtline,'(A,i2,A)')
     &                     '(A,',lmax,'i3,a,6(a,f19.12,2x,a))'
          Write(6,fmtline) ' ',(i,i=1,lmax),'|',
     &                     ('  E =',w(i)-w(1),' |',i=j,jEnd)

        End If


          Write(6,'(A,16A)') '-',('---',m=1,lmax),'|',
     &                     ('----- Real ------- Imag ---|',i=j,jEnd)

          Write(fmtline,'(A,i2,A)')
     &                         '(A,',lmax,'i3,A,5(2F13.9,1x,a))'
          Do iss=1,exch
             Write(6,fmtline) '<',
     &            (ibas(iss,m)+1,m=1,lmax),'|',
     &            ( Z(iss,i),'|',i=j,jEnd)
          End Do  !iss
          Write(6,'(A,16A)') '-',('---',m=1,lmax),'|',
     &                     ('---------------------------|',i=j,jEnd)
          Write(6,*)
      End Do ! j

      Else
         Write(6,'(10x,a)') 'Printing of EigenVectors of the '//
     &                      'Total Magnetic Interaction is Suppressed.'
         Write(6,'(10x,a)') 'Set PRLV >=3 to activate printing of '//
     &                      'this information.'
      End If


      Return
      End
