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
      Subroutine VDer_PCM(nAt,nTs,nS,AtmC,AtmChg,EF_n,EF_e,
     &     Tessera,iSphe,DerTes,DerPunt,DerRad,DerCentr,VDer)
      Implicit Real*8 (a-h,o-z)
      Dimension AtmC(3,nAt),AtmChg(nAt)
      Dimension EF_n(3,*),EF_e(3,*),VDer(nTs,*)
      Dimension Tessera(4,*),iSphe(*)
      Dimension DerTes(nTs,nAt,3),DerPunt(nTs,nAt,3,3)
      Dimension DerRad(nS,nAt,3),DerCentr(nS,nAt,3,3)
cpcm_solvent very temporary! read the potential derivatives from file
      open(1,file='DerPot.dat',status='old',form='formatted')
      do 1134 iAt = 1, nAt
        do 1134 iCoord = 1, 3
          Index = 3 * (iAt-1) + iCoord
          do 1134 iTs = 1, nTs
            read(1,*)VDer(iTs,Index)
 1134 continue
      close(1)
cpcm_solvent end
c     Loop on atoms and coordinates
      Do 100 iAt = 1, nAt
        Do 100 iCoord = 1, 3
          Index = 3 * (iAt-1) + iCoord
*
*---- (Total) derivative of the electronic + nuclear potential
*
          Do 200 iTs = 1, nTs
            L = iSphe(iTs)
c           Derivative of the representative point
            dX = DerPunt(iTs,iAt,iCoord,1)+DerCentr(L,iAt,iCoord,1)
            dY = DerPunt(iTs,iAt,iCoord,2)+DerCentr(L,iAt,iCoord,2)
            dZ = DerPunt(iTs,iAt,iCoord,3)+DerCentr(L,iAt,iCoord,3)
c           Distance tessera - nucleus
            DTessNuc = Sqrt( (Tessera(1,iTs) - AtmC(1,iAt))**2 +
     &                       (Tessera(2,iTs) - AtmC(2,iAt))**2 +
     &                       (Tessera(3,iTs) - AtmC(3,iAt))**2 )
            dCoord = Tessera(iCoord,iTs) - AtmC(iCoord,iAt)
c           Deriv. of the nuclear potential (with fixed repres. point)
            DVNuc = - AtmChg(iAt) * dCoord / DTessNuc**3
c           Electric field (electrons and nuclei) times the derivative of
c           the repres. point
            Fld_e = EF_e(1,iTs)*dX + EF_e(2,iTs)*dY + EF_e(3,iTs)*dZ
            Fld_n = EF_n(1,iTs)*dX + EF_n(2,iTs)*dY + EF_n(3,iTs)*dZ
c           Total deriv. of the potential
            VDer(iTs,Index) = VDer(iTs,Index) - Fld_e + DVNuc + Fld_n
cpcm_solvent
C     In MCLR the electron electric field is not computed, probably because of
C     DoRys set to False in mclr. Setting DoRys to True causes the program
C     to stop because the abdata.ascii file is not found.
          if(iat.eq.1.and.icoord.eq.1.and.its.eq.1)
     &      write(6,'(''In the loop VDer,Fld_e,DVNuc,Fld_n'',4f12.6)')
     &      VDer(1,Index),Fld_e,DVNuc,Fld_n
          if(iat.eq.nat.and.icoord.eq.3.and.its.eq.1)
     &      write(6,'(''In the loop VDer,Fld_e,DVNuc,Fld_n'',4f12.6)')
     &      VDer(1,Index),Fld_e,DVNuc,Fld_n
cpcm_solvent end

  200     Continue
  100 Continue
cpcm_solvent
      write(6,'(''At the end of DerPot,VDer(1,ind),VDer(nTs,ind)'')')
      Do 101 iAt = 1, nAt
        Do 101 iCoord = 1, 3
          Index = 3 * (iAt-1) + iCoord
          write(6,'(2f12.6)') VDer(1,Index), VDer(nTs,Index)
  101 continue
cpcm_solvent end
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real_array(DerTes)
        Call Unused_real_array(DerRad)
      End If
      End
