C=======================================================================
      SUBROUTINE QUINIT
C-----------------------------------------------------------------------
C!  - User initialisation
C!======================================================================
#include "qcde.h"
      common /conta/nevent, ndect,nhad
#include "qmacro.h"
      nevent=0
      ndect = 0
      nhad = 0
  999 return
      end
****
C=======================================================================
      SUBROUTINE QUEVNT (QT,KT,QV,KV)
C-----------------------------------------------------------------------
C!  - User event analysis routine
C!======================================================================
#include "qcde.h"
      common /conta/nevent, ndect,nhad
#include "qhac.h"
#include "parcut.h" 
      DIMENSION QT(KCQVEC,1), KT(KCQVEC,1), QV(KCQVRT,1), KV(KCQVEC,1)
      EXTERNAL XALEPHOK,XHADRSELEC
      LOGICAL  XALEPHOK,XHADRSELEC
      EXTERNAL IPWFLAG
#include "qmacro.h"
* First before any cut for any event
      if (xmcev) then
         write(25,212)krun,kevt,qelep
 212     format("MC_TRUE_BEFORE_CUT RUN = ",
     .        i5," EVENT",i6," ECM = ",F6.3," GEV")
         
         do i=KFMCT,KLMCT
            if((KSTABC(i).eq. 1) .or.
     >           (KSTABC(i).eq. 2) .or.
     >           (KSTABC(i).eq.-3)     ) then
               IPWF = ipwflag(ktpcod(i))
               if (ipwf.lt. 0  .AND. ipwf .ge. -10) then
                  write (6,*) " Unknown code ", ktpcod(i),"   ",cqtpn(i)
               else if (ipwf .gt. 5) then
                  write (6,*) " ??? ", ipwf,ktpcod(i),"   ",cqtpn(i)
                  
               endif
               write (25,213)qx(i),qy(i),qz(i),
     .              qm(i),qch(i),ipwf,ktpcod(i),cqtpn(i)
 213           format("px= ",f7.3," py= ",f7.3," pz= ",f7.3,
     .              " m= ",f7.3,
     .              " charge ",f4.1, " pwflag ",i4,
     .              " pid ",i5," pname ",A)
               
            endif
         enddo
         write(25,201)
      endif
      if (xmcev .OR. XALEPHOK(DUMMY)) then
         nevent = nevent + 1 
      else
         return
      endif
      if (qelep .gt. 100.) return 
      if(knevt/1000*1000 .eq. knevt)
     .     write(6,*)'run/event --> ',knevt,krun,kevt,nevent,ndect,nhad

C -- compute thrust axis with the preselected objects
      call QJOPTR('EF',' ')
      call QJTHRU(THRUST,'THRUST',KRECO)
      if (THRUST.eq.0.) then
        write(6,*) 'No thrust',KRUN,KEVT
        RETURN                ! ------------>
      endif
      ITHRU = KPDIR('THRUST',KRECO)
      if (ITHRU.le.0.or.QP(ITHRU).eq.0.) THEN
        CALL QWMESE('Problems in calculating thrust axis')
        RETURN                ! ------------>
      endif
      ndect=ndect+1
      if (XHADRSELEC(DUMMY)) then
         nhad = nhad + 1
      else
         return
      endif


      if (xmcev) then
         write(26,214)krun,kevt,qelep
 214     format("MC_TRUE_AFTER_CUT RUN = ",
     .        i5," EVENT",i6," ECM = ",F6.3," GEV")

         do i=KFMCT,KLMCT
            if((KSTABC(i).eq. 1) .or.
     >           (KSTABC(i).eq. 2) .or.
     >           (KSTABC(i).eq.-3)     ) then
               write (26,213)qx(i),qy(i),qz(i),
     .              qm(i),qch(i),ipwflag(ktpcod(i)),ktpcod(i),cqtpn(i)
            endif
         enddo
         write(26,201)
      endif

      if (XMCEV) then
         write(27,111)krun,kevt,qelep
 111  format("MC_RECO RUN = ",i5," EVENT",i6," ECM = ",F6.3," GEV")
      else
         write(27,112)krun,kevt,qelep
 112  format("ALEPH_DATA RUN = ",i5," EVENT",i6," ECM = ",F6.3," GEV")
      endif

      do ief = kfeft, kleft
         if (xlock(ief)) then
            if (xefo(ief)) then
               smass = qm(ief)
               if (smass.lt.0) smass=0.0
               write (27,113)qx(ief),qy(ief),qz(ief),
     .              smass,qch(ief),kefoty(ief)
 113           format("px= ",f7.3," py= ",f7.3," pz= ",f7.3,
     .              " m= ",f7.3,
     .              " charge ",f4.1, " pwflag ",i4)
            endif
         endif
      enddo
      write(27,201)
 201  format("END_EVENT")
      End
*********************************************************************
      SUBROUTINE QUTERM
C-----------------------------------------------------------------------
C!  - User termination
C!======================================================================
#include "qcde.h"
      common /conta/nevent, ndect,nhad
#include "qhac.h"
#include "qmacro.h"
      write(25,202)
      write(26,202)
      write(27,202)
 202  format("END_FILE")

      write(6,*)'-------------------------------------'
      write(6,*)'---- statisticazzi ------------------'
      write(6,*)'-------------------------------------'
      write(6,*)'analizzati: ',knevt,' eventi'
      write(6,*)'di cui: ',nevent, ' clas16 '
      write(6,*)'        ',ndect,' aleph ok'
      write(6,*)'        ',nhad,' adronici'
C-----------------------------------------------------------------------
  999 RETURN
      END
*=======================================================================
      SUBROUTINE FLLEPT(ITK,LFLAG)
* FLag LEPTons                                 Paolo Azzurri 8.12.97
*   -- lepton ID for isolated charged tracks --
* Input: ITK charged track
* Output: LFLAG = +1 good muon id from qmuido 
*                 -1 other muon in hcal crack
*                 +2 good electron in ecal
*                 -2 electron in ecal crack
*                  0 No lepton identification
*----------------------------------------------------------------------
#include "qcde.h"
      DATA EECMAX  / 1.0/ ! Max ECAL energy associated to mu in hcal crack 
      DATA FECMAX  / 0.1/ ! Max en(ECAL)/P  associated to mu in hcal crack 
      DATA EHCMAX  / 5.0/ ! Max HCAL energy associated to mu in hcal crack 
      DATA FHCMAX  / 0.1/ ! Max en(HCAL)/P  associated to mu in hcal crack 
      DATA RIMIN  /-2.5/   ! Min RI (dE/dx) for e
      DATA RIMAX  /999./   ! Max RI (dE/dx) for e
      DATA RLMIN  /-2.5/   ! Min R3 (RL) for e
      DATA RLMAX  /100./   ! Max R3 (RL) for e
      DATA RTMIN  /-8.0/   ! Min R2 (RT) for e 
      DATA RTMAX  /100./   ! Max R2 (RT) for e
      DATA NHPMAX  /10 /   ! Max hcal fired planes ass to e in ecal crack
      EXTERNAL MGEO_KHMAIG
#include "qmacro.h"
      LFLAG = 0
      IF (ITK.LE.0) GOTO 999
**  muon identification criteria
*    Get qmuido flags
**      CALL QMUIDO(ITK,IRUN,IBE,IBT,IM1,IM2,NEXP,NFIR,N10,
**     &     N03,XMULT,RAPP,ANG,ISHAD,SUDNT,IDF,IMCF,IER)
      IDF = KMUIIF(ITK)
      IF (IDF.GE.1) THEN
        LFLAG = 1
        GOTO 999
*    Recover hcal cracks 
      ELSE
        IFLG=MGEO_KHMAIG(ITK)
        IF (IFLG.EQ.1) THEN
* find Ecal and Hcal associated deposits
          CALL ENEREAD(ITK,IEC,IHC)
          EEC = 0
          IF (IEC.GT.0) EEC = QE(IEC)
          EHC = 0 
          IF (IHC.GT.0) EHC = QE(IHC)
          IF (EEC.GT.0 .AND. EEC.LT.EECMAX .AND. EEC.LT.FECMAX*QE(ITK)
     &         .AND. EHC.LT.EHCMAX .AND. EHC.LT.FHCMAX*QE(ITK))THEN
            LFLAG = -1
            GOTO 999
          ENDIF
        ENDIF 
      ENDIF
** electron identification estimators and criteria 
      RTE=RTCOR(ITK)
      RLE=RLCOR(ITK)
      RIE=RICOR(ITK)
      IF(XEID(ITK)) THEN
        IF (RIE.GT.RIMIN .AND.
     &       RLE.GT.RLMIN .AND. RLE.LT.RLMAX .AND.
     &       RTE.GT.RTMIN .AND. RTE.LT.RTMAX) THEN
          LFLAG = 2
          GOTO 999
        ENDIF
      ENDIF
*    Get helix distance from ECAL cracks
      CALL ECCRACK (ITK,DCRACK,IEC)
      IF (DCRACK.LE.0) THEN
*   Only Dedx if trough ecal crack
        IF (RIE.GT.RIMIN .AND. RIE.LT.RIMAX) THEN
* Look also at the HCAL Hadronic Digital Pattern
          HDE = -400
          KPPOB = NLINK('PPOB',0)
          KPCRL = NLINK('PCRL',0)
          IF (KPPOB.NE.0 .AND. KPCRL.NE.0) THEN
            NPCRL = LROWS(KPCRL)
            DO I = 1, NPCRL
              IF (ITABL(KPCRL,I,3).EQ.KTN(ITK) ) THEN
                IPOINT = ITABL(KPCRL,I,5)
                IF (IPOINT.GT.0) THEN
*  Values from the associated ppob obj with max energy 
                  ENO = RTABL(KPPOB,IPOINT,2)
                  IF (ENO.GT.HDE) THEN
                    HDE = ENO   ! Digital Energy
                    HIP = ITABL(KPPOB,IPOINT,5) ! Interaction Plane
                    HFP = ITABL(KPPOB,IPOINT,6) ! First Plane
                    HLP = ITABL(KPPOB,IPOINT,7) ! Last Plane
                    HPD = RTABL(KPPOB,IPOINT,9) ! Pattern Density
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF
          IF (HDE.GT.0 .AND. HFP.GT.0 .AND. HLP-HFP.LE.NHPMAX ) THEN
            LFLAG = -2
            GOTO 999
          ENDIF
        ENDIF
      ENDIF
*
 999  CONTINUE
      RETURN
      END
*============================================================================
      subroutine eneread(itk,iec,ihc)
C
C     Author : A.Venturi   26/05/94
C     Modified: P.Azzurri  17/05/01
C
C     Evaluate the ECAL and HCAL energies of a track
C
C     INPUT parameter    itk = track number
C
C     OUTPUT parameters  iec = ECAL associated object
C                        ihc = HCAL associated object
C
C------------------------------------------------------------------
#include "qcde.h"
C
#include "qmacro.h"
C
C...Sum the ECAL objects associated to the track
      iec=-1.
      if(knecal(itk).eq.0) goto 11
      iec = kvnew(dummy)
      call qvzero(iec)
      do 10 j=1,knecal(itk)
        jec=kecal(itk,j)
        call qvaddn(iec,jec)
 10   continue
 11   continue
C...Sum the HCAL energy associated to the track
      ihc=-1.
      if(knhcal(itk).eq.0) goto 16
      ihc = kvnew(dummy)
      call qvzero(ihc)
      do 15 j=1,knhcal(itk)
        jhc=khcal(itk,j)
        call qvaddn(ihc,jhc)
 15   continue
 16   continue
c
      RETURN
      END
C
C------------------------------------------------------------------------

      Subroutine ECCRACK (itr,dist,iflag)
*
* Input : itr : charged track
* Output: dist: to closest ecal crack at itr helix entrance into ecal
*        iflag: (eqedge flag)      0     =  OK
C                     bit 0 set = point in the phi crack
C                     bit 1 set =   "    "  "  Overlap BL-EC
C                     bit 2 set =   "    "  "  EC edge towards beam hole
C                     bit 3 set =   "  outside EMCAL theta,phi range
*
      Save rmass

#include "qdecl.h"
#include "qcde.h"
#include "qhac.h"

      Integer itr, iflag

c     Radius of Ecal barrel (roughly between stacks 1 and 2):
      Real radcal
      Parameter (radcal=197.)

c     abs(z) at Ecal endcaps (between stacks 1 and 2):
      Real zcal
      Parameter (zcal=263.)

      Real dist, phiedg, pxy, pxyz, par(7), parsw(7)

      Real rmass(1) / 0.000511 /

      Integer eqedge, ibits, efndsc
      External eqedge, ibits, efndsc, vdot

#include "qmacro.h"

      dist = 10000.
      iflag = -1

c     Compute track position at entrance to Ecal.
c     Fill vector of track starting parameters:
      pxy = qpt(itr)
      pxyz = qp(itr)
      par(1) = qfrfd0(itr)*qy(itr)/pxy
      par(2) = -qfrfd0(itr)*qx(itr)/pxy
      par(3) = qfrfz0(itr)
      par(4) = qx(itr)/pxyz
      par(5) = qy(itr)/pxyz
      par(6) = qz(itr)/pxyz
      par(7) = pxyz

c  Swim to the Ecal. 
      Call Auhcyl(radcal,zcal,15.,qch(itr),par,parsw,irc)
      if (irc .eq. 0) go to 1000
c  Check whether the track went into a crack:
      iflag = eqedge(parsw,dist,phiedg)
C
 1000 continue
      return
      end
*====================================================================
      INTEGER FUNCTION MGEO_KHMAIG(ITK)
C----------------------------------------------------------------------
C!  - Compute the HMAD geometry flag even if HMAD is not present
C!
C!     Author:    A. Venturi 
C!  
C!            Extracted from HMROAD
C!
C!            Input: ITK (I) : ALPHA track number
C!  
C!            Output: 
C!                        = 0  track in active zone
C!                        = 1  dead zone between two barrel modules
C!                        = 2  dead zone between two endcap modules
C!                        = 3  barrel module notch
C!                        = 4  endcap module spacer
C!                        = 5  not active zone of tube
C!                        = 6  first plane in thin endcap
C!                        = -1 error condition
C!  
C!            To use this function the following statements should appear
C!            at the very beginning of QUEVNT:
C! 
C!            CALL BDROP(IW,'TREX')
C!            CALL FTRACK
C!
C!            and this statement should be added in QUNEWR:
C!
C!            IF(XMINI) CALL MUNEWR
C!
C!======================================================================
#include "qcde.h"
#include "trexjj.h"
C
      INTEGER ALGTDB
      PARAMETER(NPLH=23)
      COMMON/HMGCON/NBEGTU(NPLH),SPWID1(NPLH),SPWID2(NPLH)
     +              ,NOSSP1(NPLH),NOSSP2(NPLH)
      DIMENSION BORD(NPLH,3),BORDS(NPLH,4),ZMX(4),ZMXS(5)
      DATA ICALL/0/
#include "qmacro.h"
C
      MGEO_KHMAIG=0
      ITRAC=KTN(ITK)
C
      JTREX=NLINK('TREX',ITRAC)
      NRPT=LROWS(JTREX)
C
C               load geometrical data
      IF(ICALL.EQ.0) THEN
         JDBAS=ALGTDB(KUCONS,'HBGEHBDEHEDEHEGEHETC',
     +         KRUN)
         IF(JDBAS.EQ.0) THEN
            IGFL=-1
            GO TO  99
         ENDIF
         NHBGE=NAMIND('HBGE')
         JHBGE=IW(NHBGE)
         NHBDE=NAMIND('HBDE')
         JHBDE=IW(NHBDE)
         NHEDE=NAMIND('HEDE')
         JHEDE=IW(NHEDE)
         NHEGE=NAMIND('HEGE')
         JHEGE=IW(NHEGE)
         NHETC=NAMIND('HETC')
         JHETC=IW(NHETC)
C                 tube lengths for ordinary barrel modules
         ZMX(1)=RTABL(JHETC,1,12)/2.
         ZMX(2)=RTABL(JHETC,1,11)/2.
         ZMX(3)=RTABL(JHETC,1,10)/2.
         ZMX(4)=ZMX(1)
C                  tube lengths for special barrel modules
         ZMXS(1)=RTABL(JHETC,1,13)/2.
         ZMXS(2)=RTABL(JHETC,1,12)/2.
         ZMXS(3)=RTABL(JHETC,1,11)/2.
         ZMXS(4)=RTABL(JHETC,1,10)/2.
         ZMXS(5)=RTABL(JHETC,1,12)/2.
C                  barrel radius and half lenght
         RBAR=RTABL(JHBGE,1,4)
         HBLG=RTABL(JHBGE,1,6)
C                  barrel dead zones
C                  width between two barrel modules
         DWID=RTABL(JHBDE,1,10)
         DANG=ATAN2(2.*DWID,RBAR)
C                  number and position of spacers
         DO 3 IG=1,NPLH
         NBEGTU(IG)=ITABL(JHBDE,1,13+IG)
         SPWID1(IG)=RTABL(JHBDE,1,82+IG)
         SPWID2(IG)=RTABL(JHBDE,1,105+IG)
         NOSSP1(IG)=ITABL(JHBDE,1,36+IG)
         NOSSP2(IG)=ITABL(JHBDE,1,59+IG)
 3       CONTINUE
C                  end cap dead zones
C                  width between endcaps modules
         DWEC=RTABL(JHEDE,1,4)
         REEC=RTABL(JHEGE,1,11)
         DAEC=ATAN2(2.*DWEC,REEC)
C                  radius of the pole
         REPO=RTABL(JHEGE,1,10)
C                  width and distance between spacers
         SPWD=RTABL(JHEDE,1,5)
         SPDT=RTABL(JHEDE,1,6)
         XLIM = DWEC + 2.*(SPWD + SPDT)
         XMAEC=DWEC + 4.*(SPWD + SPDT)
         W8TU=RTABL(JHETC,1,5)
C           limits of the regions with different tube length
C                 in ordinary barrel modules
         DO 4 I=1,NPLH
         BORD(I,1) = DWID + 4.*W8TU
         IF(NOSSP1(I).GE.0.AND.NOSSP2(I).GE.0) THEN
             BORD(I,2) = BORD(I,1) + 4.*W8TU + SPWID1(I)
             ICK = 0
             IF((NBEGTU(I)-NOSSP2(I)).GT.1) ICK = 1
             BORD(I,3) = BORD(I,2) + FLOAT(NBEGTU(I) - 9)*W8TU +
     +                   FLOAT(ICK)*SPWID2(I)
         ELSE
             IF(NOSSP1(I).GE.0) THEN
                 IF(NOSSP1(I).LE.8) THEN
                   BORD(I,2) = BORD(I,1) + 4.*W8TU + SPWID1(I)
                   BORD(I,3) = BORD(I,2) +FLOAT(NBEGTU(I)-9)*W8TU
                 ELSEIF((NBEGTU(I)-NOSSP1(I)).GT.1) THEN
                   BORD(I,2) = BORD(I,1) + 4.*W8TU
                   BORD(I,3)=BORD(I,2)+FLOAT(NBEGTU(I)-9)*W8TU+SPWID1(I)
                 ELSE
                   BORD(I,2) = BORD(I,1) + 4.*W8TU
                   BORD(I,3)=BORD(I,2)+FLOAT(NBEGTU(I)-9)*W8TU
                 ENDIF
             ELSE
                 BORD(I,2) = BORD(I,1) + 4.*W8TU
                 BORD(I,3) = BORD(I,2) + FLOAT(NBEGTU(I) - 9)*W8TU
             END IF
         END IF
 4       CONTINUE
C               limits of regions with different tube length
C                   in barrel special modules
         DO 6 I=1,NPLH
         BORDS(I,1) = BORD(I,1) - W8TU
         BORDS(I,2) = BORD(I,1)
         BORDS(I,3) = BORD(I,2)
         BORDS(I,4) = BORD(I,3)
 6       CONTINUE
C
         ICALL=1
      ENDIF
C
      IGFL=0
C               loop over extrapol. points
      ICB = 0
      DO 60 I=1,NRPT
C               region.module  and  plane
C                of the point
      IRGM1=ITABL(JTREX,I,JTRERG)
      IREGM=IRGM1/100
C                   hcal points
      IF(IRGM1.LE.2500) GO TO 60
      ICB = ICB + 1
      IF (ICB.GT.NPLH) GOTO 60
C
      IPLA=ITABL(JTREX,I,JTRELA)
      IF(IPLA.EQ.24) GO TO 60
C                expected coordin. and multiple
C                scattering r.m.s.
      XASS=RTABL(JTREX,I,JTREX3)
      YASS=RTABL(JTREX,I,JTREX3+1)
      ZASS=RTABL(JTREX,I,JTREX3+2)
      SIG=RTABL(JTREX,I,JTRERM)
C
      IARM=IABS(IREGM)
      MOD=IABS(IRGM1)-IARM*100
C
C                    end caps
C                    finds coordinate in module r.f.
      IF(IARM.EQ.33 .OR. IARM.EQ.34) THEN
         IF(IPLA.EQ.1) THEN
               GO TO 60
         END IF
         PHIP=ATAN2(YASS,XASS)
         IF(PHIP.LT.0.)PHIP=TWOPI+PHIP
         ISES=PHIP/PIBY3+1.
         XSRL=((-1)**(ISES+1))*YASS
         IF(ISES.EQ.2.OR.ISES.EQ.5)THEN
            PHIR=PIBY3
            IF(ZASS.GT.0.)THEN
               XSRL=-XASS*SIN(PHIR)+YASS*COS(PHIR)
            ELSE
               XSRL=XASS*SIN(PHIR)+YASS*COS(PHIR)
            ENDIF
            IF(ISES.EQ.5)XSRL=-XSRL
         ENDIF
         IF((FLOAT(ISES)*PIBY3-PHIP).LT.DAEC.OR.
     +   (PHIP-FLOAT((ISES-1))*PIBY3).LT.DAEC) THEN
               IGFL=2
               GO TO 40
         ENDIF
         IF(XSRL.GT.XMAEC) THEN
               IGFL=2
               GOTO 40
         ENDIF
         IF(XSRL.LT.DWEC) THEN
               IGFL=2
               GOTO 40
         ENDIF
         DO 10 ISP=1,4
         IF(XSRL.LE.(FLOAT(ISP)*(SPDT+SPWD)+DWEC).AND.
     +      XSRL.GE.(FLOAT(ISP)*(SPDT+SPWD)-SPWD+DWEC)) THEN
               IGFL=4
               GO TO 40
         ENDIF
 10      CONTINUE
         IF((XASS**2+YASS**2).GT.(REPO**2)) THEN
           IF(IPLA.LE.8) THEN
             IGFL=6
             GOTO 40
           ENDIF
         ENDIF
C
C                    barrel
C                    finds coordinate in module r.f.
      ELSEIF(IARM.EQ.31 .OR. IARM.EQ.32) THEN
         PHIM=FLOAT(MOD-1)*PIBY6
         PHIR=PHIM-PIBY2
         XSRL=(XASS-RBAR*COS(PHIM))*COS(PHIR)
         XSRL=XSRL+(YASS-RBAR*SIN(PHIM))*SIN(PHIR)
         IF(XSRL.GE.0.) THEN
            NMTX = MOD*2 - 2
            IF(NMTX.EQ.0) NMTX = 24
         ELSE
            NMTX = MOD*2 - 1
         END IF
         XSRL=ABS(XSRL)
C
         IF(XSRL.LE.DWID) THEN
                          IGFL=1
                          GO TO 40
                          END IF
         XMAX=FLOAT(NBEGTU(IPLA))*W8TU+SPWID1(IPLA)+SPWID2(IPLA)+DWID
         IF(XSRL.GE.XMAX)THEN
               IGFL=1
               GO TO 40
         END IF
         IF(IARM.EQ.32)THEN
              IGFL=3
              GO TO 40
         ENDIF
C
         IF(NMTX.EQ.6.OR.NMTX.EQ.7) THEN
C              check if the extr. point is in tube active zone
C                     of barrel special module
           ZTU = ZMXS(5)
           DO 37 IA = 4,1,-1
           IF(XSRL.LT.BORDS(IPLA,IA)) ZTU = ZMXS(IA)
 37        CONTINUE
           IF(ABS(ZASS).GT.ZTU) THEN
              IGFL=3
              GO TO 40
           END IF
         ELSE
C              check if the extr. point is in tube active zone
C                  of barrel ordinary module
           ZTU = ZMX(4)
           DO 38 IA = 3,1,-1
           IF(XSRL.LT.BORD(IPLA,IA)) ZTU = ZMX(IA)
 38        CONTINUE
           IF(ABS(ZASS).GT.ZTU) THEN
              IGFL=3
              GO TO 40
           END IF
         ENDIF
      ELSE
         GO TO 60
      ENDIF
C
 40   CONTINUE
 60   CONTINUE
C
 99   CONTINUE
C
      MGEO_KHMAIG=IGFL
C
      RETURN
      END
C
*====================================================================
CCCC The following is provided by Paolo Azzurri (paolo.azzurri@cern.ch)
C
      SUBROUTINE TRBINS(ITR,IPB,ICB)
C---------------------------------------------------------------------
c input ITR : track number
c output IPB : momentum bin number
c output ICB : cos(theta) bin number
c not output IDY : data year bin number (1990=1, 2000=11)
C---------------------------------------------------------------------
#include "qcde.h"
#include "qmacro.h"
      parameter (npb=4,ncb=19,ndy=11)
      dimension xpb(npb), xcb(ncb), ifr(ndy)! vector of first data run
      data xpb / 0., 20., 40., 70./ 
      data xcb /0.0,0.2,0.3,0.4,0.5,0.6,0.65,0.70,0.73,0.75,
     &     0.77,0.79,0.81,0.83,0.85,0.87,0.90,0.93,0.95/ 
      data ifr / 6898, 10564, 14883, 20165, 25911, 35734, 
     &     41275, 43364, 45196, 49217, 53072/ 
      PTR =  QP(ITR)
      do ipb = npb, 1, -1
        if (ptr.gt.xpb(ipb)) goto 10
      enddo      
 10   CTR = ABS(QCT(ITR))
      do icb = ncb, 1, -1
        if (ctr.gt.xcb(icb)) goto 30 ! for data year goto 20
      enddo      
 20   continue
      if (.not.xmcev) then
        do idy = ndy, 1, -1
          if (krun.ge.ifr(idy)) goto 30
        enddo 
      endif
      idy = 0
 30   continue
      RETURN
      END
C
C---------------------------------------------------------------------
      FUNCTION RICOR(ITR)
C---------------------------------------------------------------------
c This function returns the corrected RI estimator for data and mc 
c  from 1994 Z0 and 189GeV bhabha events
c input ITR : track number
C---------------------------------------------------------------------
#include "qcde.h"
      PARAMETER (NORI = 1000)
      DIMENSION AMC(4,19),RMC(4,19)
      DIMENSION ADT(4,19),RDT(4,19)
      DIMENSION AADT(11)
C
c mean and sigma RI distributions in momentum,cos(th) bins for DATA.
      DATA ADT/-0.102, 0.033,0.078,-0.006,
     &      0.016,0.112,0.138,0.074,
     &     -0.032,0.084,0.113,0.067,
     &     -0.032,0.137,0.105,0.121,
     &     -0.207,0.151,0.103,0.040,
     &     -0.219,0.143,0.097,0.008,
     &     -0.214,0.034,0.051,-0.027,
     &     -0.193,0.025,0.020,-0.064,
     &     -0.157,-0.007,0.036,-0.088,
     &     -0.272,-0.078,-0.039,-0.085,
     &     -0.326,-0.064,-0.018,-0.128,
     &     -0.270,-0.103,-0.126,-0.135,
     &     -0.456,-0.113,-0.090,-0.179,
     &     -0.523,-0.186,-0.162,-0.196,
     &     -0.594,-0.267,-0.242,-0.225,
     &     -0.656,-0.333,-0.361,-0.317,
     &     -0.495,-0.348,-0.400,-0.334,
     &     -0.470,-0.387,-0.474,-0.411,
     &     -0.584,-0.542,-0.602,-0.629/
      DATA RDT/ 1.270,0.934,0.910,0.868,
     &      1.158,0.886,0.895,0.902,
     &      1.080,0.919,0.889,0.939,
     &      1.160,0.907,0.889,0.878,
     &      1.208,0.903,0.880,0.932,
     &      1.070,0.864,0.888,0.918,
     &      1.249,0.910,0.865,0.888,
     &      1.125,0.921,0.859,0.854,
     &      0.910,0.889,0.867,0.876,
     &      0.983,0.848,0.877,0.880,
     &      1.042,0.876,0.855,0.888,
     &      1.142,0.863,0.864,0.852,
     &      1.073,0.907,0.867,0.857,
     &      0.975,0.921,0.874,0.854,
     &      0.923,0.918,0.879,0.888,
     &      0.959,0.954,0.874,0.865,
     &      1.003,0.927,0.893,0.913,
     &      1.008,0.973,0.899,0.883,
     &      1.008,0.958,0.890,0.851/
c mean and sigma RI distributions in momentum,cos(th) bins for MC.
      DATA AMC/-0.049,0.342,0.426,0.463,
     &       0.014,0.410,0.505,0.540,
     &       0.023,0.467,0.569,0.579,
     &       0.100,0.561,0.661,0.695,
     &       0.228,0.675,0.758,0.819,
     &       0.240,0.761,0.848,0.868,
     &       0.302,0.804,0.902,0.929,
     &       0.270,0.871,0.964,1.007,
     &       0.368,0.966,1.025,1.054,
     &       0.448,1.000,1.052,1.071,
     &       0.435,1.005,1.068,1.112,
     &       0.448,0.992,1.036,1.100,
     &       0.450,0.976,1.035,1.051,
     &       0.336,0.908,0.982,1.010,
     &       0.322,0.903,0.938,0.966,
     &       0.319,0.807,0.858,0.872,
     &       0.330,0.793,0.818,0.843,
     &       0.245,0.515,0.551,0.576,
     &       0.189,0.285,0.314,0.333/
      DATA RMC/1.276,1.128,1.061,1.063,
     &        1.385,1.151,1.071,1.066,
     &        1.394,1.132,1.086,1.083,
     &        1.336,1.170,1.103,1.105,
     &        1.373,1.192,1.114,1.100,
     &        1.463,1.212,1.133,1.124,
     &        1.333,1.217,1.144,1.134,
     &        1.338,1.227,1.161,1.154,
     &        1.216,1.205,1.155,1.164,
     &        1.219,1.198,1.174,1.185,
     &        1.218,1.217,1.186,1.178,
     &        1.248,1.220,1.168,1.165,
     &        1.202,1.174,1.148,1.146,
     &        1.234,1.150,1.134,1.126,
     &        1.178,1.150,1.122,1.106,
     &        1.141,1.106,1.084,1.087,
     &        1.145,1.063,1.059,1.058,
     &        1.048,1.004,1.005,1.006,
     &        1.075,0.975,0.965,0.976/
* additional overall year by year shifts in Z0 data 
      DATA AADT /0.00, 0.00, 0.05, 0.12, 0.00,
     &     0.01, 0.00, 0.11, 0.27, 0.00, 0.20 /
#include "qmacro.h"
      CALL TRBINS(ITR,IP,IC)
      RMASS = QPMASS('e+')
      RICOR = NORI
      CALL QDEDX(ITR,1,RMASS,1.,RI,NS,TL,RIEXP,SIGM,JER)
      IF (JER.LE.0) THEN
       RIRAW = (RI - RIEXP)/SIGM
       IF (XMCEV) RICOR = (RIRAW - AMC(IP,IC))/RMC(IP,IC)
       IF (.NOT.XMCEV) RICOR = (RIRAW - ADT(IP,IC))/RDT(IP,IC)
      ELSE
C no de/dx, check pad de/dx
        XNS = 0.
        CALL QPADX(ITR,1,RMASS,1.,RI,XNS,TL,RIEXP,SIGM,JER)
        IF (JER.LE.0) RICOR = (RI - RIEXP)/SIGM
      ENDIF
C..
      RETURN
      END
C
C---------------------------------------------------------------------
      FUNCTION RLCOR(ITR)
C---------------------------------------------------------------------
c This function returns the corrected RL estimator for data and mc 
c  from 1994 Z0 and 189GeV bhabha events
** TO BE USED WITH THE FIXPRL CORRECTION TO THE MONTE CARLO RL VALUES  
**              (activated by the UPRL alpha card)  
c input ITR : track number
C---------------------------------------------------------------------
#include "qcde.h"
      PARAMETER (NORL = 1000)
      DIMENSION AMC(4,19),RMC(4,19)
      DIMENSION ADT(4,19),RDT(4,19)
c
c mean and sigma RL distributions in momentum,cos(th) bins for DT and MC
      DATA (ADT(1,J),j=1,19)/-0.125,-0.090,-0.084,-0.154,-0.185,
     &     -0.163,-0.135,0.023,0.571,0.235,0.223,0.053,-0.159,
     &     -0.221,-0.165,-0.196,-0.177,-0.300,-0.461/
      DATA (RDT(1,J),j=1,19)/1.182,1.224,1.130,1.279,1.268,1.170,
     &     1.295,1.256,0.894,1.279,1.075,1.385,1.303,1.253,1.099,
     &     1.108,1.184,1.180,1.090/
      DATA (AMC(1,J),j=1,19)/-0.155,-0.180,-0.226,-0.191,-0.266,
     &     -0.175,-0.183,-0.178,0.023,-0.120,-0.219,-0.390,-0.530,
     &     -0.396,-0.324,-0.302,-0.337,-0.270,-0.601/
      DATA (RMC(1,J),j=1,19)/1.081,1.182,1.150,1.119,1.207,1.201,
     &     1.223,1.173,1.007,1.231,1.195,1.233,1.189,1.158,1.122,
     &     1.107,1.109,1.087,1.064/
      DATA (ADT(2,J),j=1,19)/-0.007,0.040,-0.003,-0.068,-0.025,
     &     -0.025,-0.006,0.048,0.265,0.091,0.180,0.364,-0.064,
     &     -0.193,-0.116,-0.153,-0.067,-0.142,-0.243/
      DATA (RDT(2,J),j=1,19)/1.072,1.047,1.097,1.104,1.105,1.082,
     &     1.142,1.118,0.913,1.194,0.997,1.329,1.192,1.109,1.044,
     &     1.043,1.097,1.062,1.010/
      DATA (AMC(2,J),j=1,19)/0.065,0.024,0.020,0.008,-0.018,-0.034,
     &     0.006,0.064,-0.354,-0.060,-0.042,0.044,-0.209,-0.165,
     &     -0.154,-0.125,-0.164,-0.097,-0.084/
      DATA (RMC(2,J),j=1,19)/1.005,0.996,1.023,1.039,1.072,1.138,
     &     1.125,1.104,0.960,1.244,1.163,1.122,1.018,1.011,0.966,
     &     0.935,0.924,0.913,0.937/
      DATA (ADT(3,J),j=1,19)/-0.021,-0.037,-0.034,-0.030,-0.046,
     &     -0.055,-0.026,0.068,0.177,0.039,0.099,0.417,-0.053,-0.143,
     &     -0.150,-0.122,0.080,0.158,0.126/
      DATA (RDT(3,J),j=1,19)/1.017,1.043,1.040,1.038,1.038,1.018,
     &     1.046,0.969,0.835,1.092,0.903,1.259,1.100,1.053,1.013,
     &     0.995,1.104,1.096,1.021/
      DATA (AMC(3,J),j=1,19)/0.078,0.050,0.047,0.030,0.012,0.018,
     &     0.037,0.066,-0.412,-0.075,-0.072,0.142,-0.210,-0.176,
     &     -0.135,-0.121,-0.143,-0.002,0.125/
      DATA (RMC(3,J),j=1,19)/1.017,1.024,1.025,1.036,1.046,1.038,
     &     1.050,1.082,0.997,1.291,1.220,1.067,0.969,0.969,0.956,
     &     0.959,0.959,0.964,0.979/
      DATA (ADT(4,J),j=1,19)/-0.091,-0.122,-0.139,-0.192,-0.125,
     &     -0.121,-0.111,-0.056,-0.030,-0.237,-0.072,0.235,-0.172,
     &     -0.257,-0.208,-0.184,-0.019,-0.041,-0.113/
      DATA (RDT(4,J),j=1,19)/0.862,0.824,0.892,0.871,0.814,0.830,
     &     0.838,0.773,0.728,0.862,0.773,0.970,0.867,0.873,0.813,
     &     0.831,0.919,0.886,0.900/
      DATA (AMC(4,J),j=1,19)/0.046,0.062,0.021,0.015,-0.007,-0.009,
     &     0.007,0.103,-0.707,-0.349,-0.085,0.175,-0.298,-0.250,
     &     -0.213,-0.164,-0.185,-0.091,0.027/
      DATA (RMC(4,J),j=1,19)/1.026,1.042,1.032,1.039,1.063,1.059,
     &     1.046,1.163,1.108,1.150,1.503,1.167,0.952,0.941,0.945,
     &     0.946,0.946,0.956,0.954/
#include "qmacro.h"
C...check if UPRL cards is activated for MC events
      if(xmcev.and.iw(namind('UPRL')).le.0) then
        nwarn=nwarn+1
        if(nwarn.le.10) 
     +       call qwmese('RLCOR Warning - MC events need the UPRL'//
     +       ' card to correct properly RL !')
      endif
C
      CALL TRBINS(ITR,IP,IC)
C
      RLCOR = NORL
      IF (.NOT.XEID(ITR)) GOTO 999
      RLRAW = QEIDRI(ITR,3)
      IF (XMCEV) RLCOR = (RLRAW - AMC(IP,IC))/RMC(IP,IC)
      IF (.NOT.XMCEV) RLCOR = (RLRAW - ADT(IP,IC))/RDT(IP,IC)
 999  RETURN
      END
C
C---------------------------------------------------------------------
      FUNCTION RTCOR(ITR)
C---------------------------------------------------------------------
c This function returns the corrected RT estimator for data and mc 
c  from 1994 Z0 and 189GeV bhabha events
c input ITR : track number
C---------------------------------------------------------------------
#include "qcde.h"
      PARAMETER (NORT = 1000)
      DIMENSION AMC(4,19),RMC(4,19)
      DIMENSION ADT(4,19),RDT(4,19)
c
c mean and sigma RT distributions in momentum,cos(th) bins for DT and MC
      DATA (ADT(1,J),j=1,19)/0.877,0.800,0.798,0.700,0.631,0.440,
     &     0.331,0.016,-0.269,-1.048,-0.942,-0.243,0.246,0.524,
     &     0.563,0.461,0.540,0.650,0.611/
      DATA (RDT(1,J),j=1,19)/1.609,1.706,1.583,1.529,1.644,1.506,
     &     1.507,1.359,1.219,1.298,1.486,1.595,1.606,1.678,1.648,
     &     1.700,1.874,1.804,2.114/
      DATA (AMC(1,J),j=1,19)/0.795,0.753,0.653,0.644,0.541,0.413,
     &     0.265,0.100,0.113,-0.853,-0.545,0.269,0.565,0.571,0.670,
     &     0.733,0.880,0.711,0.514/
      DATA (RMC(1,J),j=1,19)/1.716,1.748,1.666,1.684,1.605,1.591,
     &     1.490,1.479,1.331,1.349,1.537,1.468,1.638,1.640,1.577,
     &     1.687,1.806,1.856,2.192/
      DATA (ADT(2,J),j=1,19)/2.189,2.259,2.223,2.075,2.001,1.974,
     &     1.801,1.642,1.218,-0.349,0.927,1.126,1.856,2.030,2.021,
     &     2.061,2.282,2.415,2.417/
      DATA (RDT(2,J),j=1,19)/1.774,1.765,1.746,1.760,1.811,1.755,
     &     1.763,1.860,2.136,1.976,1.949,1.866,1.854,1.752,1.858,
     &     1.721,1.726,1.712,1.840/
      DATA (AMC(2,J),j=1,19)/2.173,2.155,2.104,1.968,1.842,1.709,
     &     1.606,1.560,1.897,0.294,1.022,1.192,1.718,1.796,1.845,
     &     1.880,2.080,1.943,1.888/
      DATA (RMC(2,J),j=1,19)/1.743,1.765,1.760,1.779,1.793,1.806,
     &     1.842,1.835,1.808,2.025,1.878,1.737,1.814,1.817,1.756,
     &     1.763,1.721,1.764,1.878/
      DATA (ADT(3,J),j=1,19)/0.078,0.102,0.161,0.099,0.082,0.017,
     &     -0.047,-0.126,-0.211,-1.686,-0.190,-0.358,0.182,0.207,
     &     0.220,0.053,-0.156,-0.555,-0.829/
      DATA (RDT(3,J),j=1,19)/1.131,1.085,1.086,1.105,1.115,1.174,
     &     1.184,1.259,1.575,1.447,1.443,1.231,1.211,1.180,1.165,
     &     1.124,1.245,1.276,1.319/
      DATA (AMC(3,J),j=1,19)/0.020,0.000,-0.026,-0.074,-0.152,-0.195,
     &     -0.233,-0.050,1.330,-0.848,0.112,-0.340,-0.053,-0.135,
     &     -0.148,-0.240,-0.163,-0.468,-0.863/
      DATA (RMC(3,J),j=1,19)/1.110,1.096,1.112,1.111,1.127,1.127,
     &     1.155,1.345,1.645,1.569,1.592,1.095,1.114,1.147,1.079,
     &     1.074,1.119,1.183,1.237/
      DATA (ADT(4,J),j=1,19)/0.247,0.407,0.395,0.334,0.263,0.244,
     &     0.217,0.168,-0.075,-1.261,0.470,0.218,0.564,0.580,0.517,
     &     0.522,0.559,0.603,0.751/
      DATA (RDT(4,J),j=1,19)/1.130,1.080,1.144,1.120,1.115,1.295,
     &     1.212,1.346,1.729,1.729,1.507,1.218,1.173,1.261,1.118,
     &     1.095,1.016,1.002,1.150/
      DATA (AMC(4,J),j=1,19)/0.325,0.310,0.253,0.225,0.189,0.136,
     &     0.118,0.315,1.450,-0.210,0.695,0.056,0.207,0.156,0.169,
     &     0.239,0.343,0.302,0.363/
      DATA (RMC(4,J),j=1,19)/1.077,1.103,1.058,1.105,1.141,1.174,
     &     1.189,1.333,1.672,1.704,1.522,1.182,1.140,1.096,0.960,
     &     0.950,0.921,0.889,0.990/
#include "qmacro.h"
      CALL TRBINS(ITR,IP,IC)

      RTCOR = NORT
      IF (.NOT.XEID(ITR)) GOTO 999
      RTRAW = QEIDRI(ITR,2)
      IF (XMCEV) RTCOR = (RTRAW - AMC(IP,IC))/RMC(IP,IC)
      IF (.NOT.XMCEV) RTCOR = (RTRAW - ADT(IP,IC))/RDT(IP,IC)
 999  RETURN
      END
************************************************************************
      subroutine genereco(itk,mc_cod,nmoth)
************************************************************************
* routine written by Gian Musolino 26/9/96 (from GRDM by David Jaffe)
* input: itk alpha reconstructed MC track number
* output: mc_cod vector of MC codes for track ITK and all the
*         mothers till 7th generation
************************************************************************
#include "qcde.h"
#include "qhac.h"
      dimension mc_cod(10)
#include "qmacro.h"
************************************************************************
* Find the best match for each reconstructed track. 'Best' means that
* the reconstructed track shares the most hits with the monte carlo
* track, if two or more tracks share the same number of hits then take
* the first track.
************************************************************************
      hits = 0
      ibest = -1
      do i=1,7
        mc_cod(i) = -1
      enddo
      do 20 j = 1,knmtch(itk)
        if(ksmtch(itk,j) .gt. hits) then
          hits = ksmtch(itk,j)
          ibest = j
        endif
 20   continue
************************************************************************
* found best match, now get its code
************************************************************************
      if(ibest .gt. 0) then
        j = kmtch(itk,ibest)
        if(xmc(j)) then
          mc_cod(1) = ktpcod(j)
          mc_cod(10) = j
          mtk = j
        else
          goto 999
        endif
************************************************************************
* no match found for this track
************************************************************************
      else
        goto 999
      endif
************************************************************************
* mc_cod(1) is the code of the MC particle matching ITK
* mtk is the number of the MC particle matching ITK
************************************************************************
      moth = mtk
************************************************************************
* Get the codes for all the mothers of ITK if available.
************************************************************************
      igen = 1
 50   continue
      igen = igen + 1
      if(igen .gt. 9) goto 999
      idau = moth
      nmoth = knmoth(idau)
      if(idau .gt .0) then
        if(knmoth(idau) .eq. 1) then
          moth = kmoth(idau,1)
          momc = ktpcod(moth)
          if(igen .le. 9) mc_cod(igen) = momc
          if(momc .gt. 350) goto 999
          goto 50
        endif
      endif
 999  continue
      return
      end
************************************************************************
C
*====================================================================
      subroutine mctrue(mchz,pt1,pt2,einv,evis,xmasz,csz0,
     &                  csmis,csel,cs1,cs2)
************************************************************************
* routine written by Paul Spanish on 08/05/02 
* input: no input
* output: mchz --> monte carlo true on the HZ event
* output: pt1/2   --> muon pt
* output: einv  --> invisible energy of the missing ele
* output: evis  --> energy of the visible ele
* output: xmasz  --> dimuon invariant mass
* output: csz0  --> cos theta Z0
* output: csmis  --> cos theta invis el
* output: csel  --> cos theta visib el
************************************************************************
#include "qcde.h"
#include "qhac.h"
#include "qmacro.h"
      nmu = 0
      ntau = 0
      nz0  = 0
      itau = 0
      ipos = 0
      idau = 0
      imoth = 0
      nmoth = 0
      mchz = 0
      pt2 = 0 
      pt1 = 0 
      einv = 0
      evis = 0 
      xmasz = 0
      csz0 = 0
      csmis = 0
      csel = 0
      cs1 = 0
      cs2 = 0
      iz  = kvnew(dummy)
      if (xmcev) then
         do imc = kfmct, klmct 
            idau = imc
            kod = ktpcod(idau)
            if (kod.eq.33.or.kod.eq.34) then
              nmoth = knmoth(idau)
              imoth = kmoth(idau,1)  
*              if (nmoth.eq.0.or.ktpcod(imoth).eq.44) then
                 ntau = ntau+1
*              endif
            endif
            if (kod.eq.5.or.kod.eq.6) then
              nmoth = knmoth(idau)
              imoth = kmoth(idau,1)
              if (nmoth.eq.0.or.ktpcod(imoth).eq.44) then
                 nmu = nmu+1  
                 if(qch(imc) .gt. 0) then 
                    imup = imc                  
                 else 
                    imum = imc 
                 endif
              endif 
            endif 
         enddo
      endif  
*
* Signal Definition
*
      if ((nmu.eq.2).and.(ntau.eq.4)) then 
         mchz = 1 
      endif
      if (abs(mchz) .eq. 1) then
            pt1 = sqrt(qx(imum)**2+qy(imum)**2)
            pt2 = sqrt(qx(imup)**2+qy(imup)**2)
            xmasz = qm2(imum,imup)
            csz0 = xcosz0
            cs1 = xcs1
            cs2 = xcs2
      endif
 999  continue
      return
      end

      subroutine bremadd(ich,ibrem)
c-----------------------------------------------------------
c  Add bremstrahlung photons to leptons
c     Input:    ich,   the lepton index
c     Output:   ibrem, a new track including all brems
c Patrick Janot
c-----------------------------------------------------------
#include "parcut.h"
#include "qcde.h"
#include "qmacro.h"
****
      ibrem = kvnew(dummy)
      call qvcopy(ibrem,ich)
      theti=acos(min(1.,max(-1.,qct(ibrem))))*180./qqpi
      phii=qph(ibrem)*180./qqpi
      do kga=kfrtot(3),klstot(4)
        thetnt=acos(min(1.,max(-1.,qct(kga))))*180./qqpi
        phint=qph(kga)*180./qqpi
        delthi=abs(theti-thetnt)
        delphi=qch(ibrem)*(phii-phint)
        if (delphi.gt.180.) delphi=delphi-360.
        if (delphi.lt.-180.) delphi=delphi+360.
        if ((delthi.lt.2.).and.(delphi.gt.-3.).and.
     .       (delphi.lt.6.)) then
          call qvaddn(ibrem,kga)
        endif
      enddo
      return
      end
************************************************************************
      SUBROUTINE QUIHIS
CKEY HIST INIT /USER
C----------------------------------------------------------------------
C! initialize histogram package
C called from QMINIT
C                                                   H.Albrecht 20.09.88
C                                       modified    J.Boucrot  03.05.93
C----------------------------------------------------------------------
C
#include "qcde.h"
      PARAMETER (LQPAW = 1000000)
      COMMON /PAWC/ PAWC(LQPAW)
C
C       HBOOK initialisation
C
      IF (KDEBUG.NE.0)  THEN
         WRITE ( KUPRNT , 1001 ) LQPAW
         IF (KUPTER.NE.0) WRITE ( KUPTER , 1001 ) LQPAW
      ENDIF
      CALL HLIMIT (LQPAW)
C
C       read HISTOGRAM AND HTITLE card and open histogram file
C
      CALL QMIHIS
C
 1001 FORMAT ('0_QUIHIS_ Init HBOOK4 with ',I8,' words working space')
      END

      LOGICAL function XALEPHOK(DUMMY)
#include "qdecl.h"
#include "qcde.h"
#include "qmacro.h"
      XALEPHOK = .true. 
      if(.not.xvitc)then
         XALEPHOK =  .false.
      endif
      if(.not.xvtpc)then
         XALEPHOK =  .false.
      endif
      if(.not.xvtpcd)then
         XALEPHOK =  .false.
      endif
      if(.not.xvecal)then
         XALEPHOK =  .false.
      endif
      if(.not.xvhcal)then
         XALEPHOK =  .false.
      endif
      if(.not.xvdeok(dummy))then
         XALEPHOK =  .false.
      endif
      if(KEVEES.NE.1)then
         XALEPHOK =  .false.
      endif
      return 
      end
      
      LOGICAL function XHADRSELEC(DUMMY)
#include "qdecl.h"
#include "qcde.h"
      integer nch,nne
      integer i
      integer nt,type
      integer ntr
      logical xok
      real ech
      real ct,st,th,ph,pt,ee,d0,z0
      real cts
#include "qmacro.h"
      xhadrselec = .false.
      call QJOPTR('CH',' ')
      call QLZER(KRECO)
      nch = 0
      ech = 0.
      nne = 0
      do i=KFEFT,KLEFT
         if(.not.XEFO(i)) then
            write(6,*)
     >      'QCD_HADSEL: ENFLW-object without EFOL info => ',i
         else
            ct   = QCT(i)
            st   = sqrt(max(1.e-10,1.-ct*ct))
            th   = acos(max(min(ct,1.),-1.)) * 180./QQPI
            ph   = QPH(i) * 180./QQPI
            pt   = QPT(i)
            ee   = QE(i)
            d0   = QDB(i)
            z0   = QZB(i)
            nt   = KFRTNT(i)
            type = KEFOTY(i)
            if(     type.lt.0) then
               write(6,*)
     >              'QCD_EFLOW: unknown EFLOW object type = ',type
            else if(type.le.3) then ! charged tracks       
               xok = .true.
               if(abs(d0).gt.     2.   ) xok = .false.
               if(abs(z0).gt.    10.   ) xok = .false.
               if(abs(ct).gt.     0.94 ) xok = .false.
               if(    pt .lt.     0.1  ) xok = .false.
               if(    nt .lt.     4    ) xok = .false.
               if(xok) then
                  nch = nch + 1
                  ech = ech + ee
                  call QLOCK(i)
               endif
            else if(type.eq.4) then   ! neutral ECAL objects  
               xok = .true.
               if(abs(ct).gt.     0.98 ) xok = .false.
               if(    ee .lt.     0.80 ) xok = .false.
               if(xok) then 
                  nne = nne + 1
                  call QLOCK(i)
               endif
            else if(type.eq.5) then   ! neutral HCAL objects
               xok = .true.
               if(abs(ct).gt.     0.98 ) xok = .false.
               if(    ee .lt.     0.80 ) xok = .false.
               if(xok) then
                  call QLOCK(i)
                  nne = nne + 1
               endif
            else if(type.gt.8) then
               write(6,*)
     >         'QCD_HADSEL: unknown EFLOW object type = ',type
            endif
         endif
      enddo

      call QLREV(KRECO)
      call QJSPHE(sph,'sph_axis',KRECO)
      if(nch.gt.0) then
         cts = abs(QCT(KPDIR('sph_axis',KRECO)))
      else
         cts = 2.
      endif
      ntr = nch + nne
      call QLREV(KRECO)
      if(nch.lt.5) RETURN
      if(ech.lt. 15. ) RETURN
      if(cts.gt.  0.82 ) RETURN
      if(ntr.lt. 0) RETURN
      xhadrselec = .true.

      return
      end
      FUNCTION IPWFLAG(KCOD)
      IPWFLAG = -3
      IF (KCOD .EQ. 8 .OR. KCOD .EQ. 9 .OR. KCOD .EQ. 11 .OR.
     .    KCOD .EQ. 12 .OR. KCOD .EQ. 14 .OR. KCOD .EQ. 15) THEN
         IPWFLAG = 0
      ELSE IF( KCOD .EQ. 5 .OR. KCOD .EQ.6) THEN
         IPWFLAG = 1
      ELSE IF( KCOD .EQ. 2 .OR. KCOD .EQ.3) THEN
         IPWFLAG = 2
      ELSE IF (KCOD .EQ. 16 .OR. KCOD .EQ. 18 .OR. KCOD .EQ. 26 .OR. 
     .         KCOD .EQ. 23 .OR. KCOD .EQ. 31 .OR. KCOD .EQ. 19 .OR.
     .         KCOD .EQ. 21 .OR. KCOD .EQ. 27 .OR. KCOD .EQ. 30 .OR.
     .         KCOD .EQ. 29 .OR. KCOD .EQ. 22) THEN
         IPWFLAG = 3
      ELSE IF( KCOD .EQ. 1 .OR. KCOD .EQ. 7) THEN
         IPWFLAG = 4
      ELSE IF( KCOD .EQ. 10 .or. KCOD .EQ. 13 .OR. KCOD .EQ. 25) THEN
         IPWFLAG = 5
      ELSE IF (KCOD .EQ. 4 .OR. (KCOD .GE. 58 .AND. KCOD .LE. 63) ) THEN
         IPWFLAG = -11
      ELSE 
        IPWFLAG = -2 
      ENDIF
      RETURN
      END
