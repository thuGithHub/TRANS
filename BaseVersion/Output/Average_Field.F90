!******for turbomacinery simulation, P0 should be modified!

SUBROUTINE Average_Fields
    USE Global
    IMPLICIT NONE
    include "../tecio.F90"
    INTEGER:: iBlock
    INTEGER:: NI,NJ,NK
    INTEGER:: NI1,NJ1,NK1
    INTEGER:: I,J,K,L
    INTEGER:: KZ
    REAL:: TIME
    CHARACTER(LEN=6):: cc
    character(LEN=1):: NULLCHR
    Integer*4::   Debug,II, III,VIsDouble
    POINTER (NullPtr,Nulp)
    INTEGER*4:: Nulp(*)
    REAL*4,ALLOCATABLE:: VARs1(:,:,:)
    REAL*4,ALLOCATABLE:: VARs2(:,:,:)
    REAL*4,ALLOCATABLE:: VARs3(:,:,:)
    REAL*4,ALLOCATABLE:: VARs4(:,:,:)
    REAL*4,ALLOCATABLE:: VARs5(:,:,:)
    REAL*4,ALLOCATABLE:: VARs6(:,:,:)
    REAL*4,ALLOCATABLE:: VAR_value(:)
    INTEGER*4,ALLOCATABLE:: varlocation(:)
    INTEGER*4::   VAR_num
    REAL:: Den,RmiuT,vK,vOmega,vab_01,vab_02,vab_03,vab_04,vab_05,vab_06,vab_07,vab_08,vab_09,vab_10,vab_11,vab_12,&
        &   vab_13,vab_14,vab_15,vab_16,vab_17,vab_18,vab_19,vab_20,vab_21,vab_22,vab_23,vab_24,vab_25,vab_26,QQ,AA
    REAL:: DuDx,DuDy,DuDz,DvDx,DvDy,DvDz,DwDx,DwDy,DwDz,S0
    REAL:: upup,vpvp,wpwp,upvp,upwp,vpwp
    !REAL:: uurms,vvrms,wwrms,uvrms,uwrms,vwrms
    
    DO iBlock = 1, Max_Block
        CALL GLOBAL_SetPointers(iBlock)
    
        NI=ThisBlock%NI
        NJ=ThisBlock%NJ
        NK=ThisBlock%NK
        NI1=ThisBlock%NI1
        NJ1=ThisBlock%NJ1
        NK1=ThisBlock%NK1

        IF(Ki == KI_aver)THEN
            DO J=1,NJ1
            DO K=1,NK1
            DO I=1,NI1
!	        for instantaneous values
!!!!!!  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
!     >'X Y Z U V W R M Mu k P0 P  T  uu vv ww uv uw vw'
                Den   =  V(1,I,J,K)
                RmiuT = Rmiu(I,J,K)
                vK    =  V(6,I,J,K) / Den
                vOmega=  V(7,I,J,K) / Den
                vOmega= max(vOmega, Rof)
!!!!!!
                vab_04= V(2,I,J,K) / Den                        !! u
                vab_05= V(3,I,J,K) / Den                        !! v
                vab_06= V(4,I,J,K) / Den                        !! w
                vab_07= Den

                QQ    = sqrt(vab_04*vab_04+vab_05*vab_05+vab_06*vab_06)
                AA    = SQRT(1.4_8*PP(I,J,K)/ Den )             !! velocity of sound
                vab_08= QQ/AA                                   !! local Mach number
                vab_09= Rmiu(I,J,K)                             !! eddy viscosity
                vab_10=  V(6,I,J,K) / Den                       !! kinetic energy
                vab_11= PP(I,J,K)*(1.0_8+0.2_8*vab_08*vab_08)**3.5_8/PPF    !! p0 total pressure
                vab_12= PP(I,J,K)                                           !! p  static pressure
                vab_13=  T(I,J,K)                                           !! temperature
!			stress
                DuDx= dQdxyz( 1,I,J,K)
                DuDy= dQdxyz( 2,I,J,K)
                DuDz= dQdxyz( 3,I,J,K)
                DvDx= dQdxyz( 4,I,J,K)
                DvDy= dQdxyz( 5,I,J,K)
                DvDz= dQdxyz( 6,I,J,K)
                DwDx= dQdxyz( 7,I,J,K)
                DwDy= dQdxyz( 8,I,J,K)
                DwDz= dQdxyz( 9,I,J,K)
!   
                S0   = -2.0_8/3.0_8*(dUdx+dVdy+dWdz)

                upup = RmiuT/Ref*(S0+2.0_8*dUdx)-2.0_8/3.0_8*Den*vK
                vpvp = RmiuT/Ref*(S0+2.0_8*dVdy)-2.0_8/3.0_8*Den*vK
                wpwp = RmiuT/Ref*(S0+2.0_8*dWdz)-2.0_8/3.0_8*Den*vK
                upvp = RmiuT/Ref*(      dUdy+dVdx)
                upwp = RmiuT/Ref*(      dUdz+dWdx)
                vpwp = RmiuT/Ref*(      dVdz+dWdy)
    
                upup = -upup 
                vpvp = -vpvp 
                wpwp = -wpwp 
                upvp = -upvp 
                upwp = -upwp 
                vpwp = -vpwp 
!!!!!!
                vab_14  = upup
                vab_15  = vpvp
                vab_16  = wpwp
                vab_17  = upvp
                vab_18  = upwp
                vab_19  = vpwp
!			for average
                Ub(I,J,K)=vab_04
                Vb(I,J,K)=vab_05
                Wb(I,J,K)=vab_06

                 Rb(I,J,K)=vab_07
                aMb(I,J,K)=vab_08
                rMb(I,J,K)=vab_09
                vKb(I,J,K)=vab_10
    
                Ptb(I,J,K)=vab_11
                Pb(I,J,K)=vab_12
                Tb(I,J,K)=vab_13

                uub(I,J,K)=vab_14
                vvb(I,J,K)=vab_15
                wwb(I,J,K)=vab_16
                uvb(I,J,K)=vab_17
                uwb(I,J,K)=vab_18
                vwb(I,J,K)=vab_19

           !uuaver(I,J,K)=vab_04*vab_04
           !vvaver(I,J,K)=vab_05*vab_05
           !wwaver(I,J,K)=vab_06*vab_06
           !uvaver(I,J,K)=vab_04*vab_05
           !vwaver(I,J,K)=vab_05*vab_06
           !uwaver(I,J,K)=vab_04*vab_06
           !ppaver(I,J,K)=vab_12*vab_12


            ENDDO
            ENDDO
            ENDDO
        ELSE
            Time=float(Ki-KI_aver)
            DO J=1,NJ1
            DO K=1,NK1
            DO I=1,NI1
!	        for instantaneous values
            Den   =  V(1,I,J,K)
            vK    =  V(6,I,J,K) / Den
            RmiuT = Rmiu(I,J,K)
            vOmega=  V(7,I,J,K) / Den
            vOmega= max(vOmega, Rof)
!!!!!!
            vab_04= V(2,I,J,K) / Den                    !! u
            vab_05= V(3,I,J,K) / Den                    !! v
            vab_06= V(4,I,J,K) / Den                    !! w
            vab_07= Den

            QQ    = sqrt(vab_04*vab_04+vab_05*vab_05+vab_06*vab_06)
            AA    = SQRT(1.4_8*PP(I,J,K)/ Den )                     !! velocity of sound
            vab_08= QQ/AA!! local Mach number
            vab_09= Rmiu(I,J,K)         !! eddy viscosity

            vab_11= PP(I,J,K)*(1.0_8+0.2_8*vab_08*vab_08)**3.5_8/PPF    !! p0 total pressure
            vab_12= PP(I,J,K)                           !! p  static pressure
            vab_13=  T(I,J,K)                           !! temperature
!!!!!!
            DuDx= dQdxyz( 1,I,J,K)
            DuDy= dQdxyz( 2,I,J,K)
            DuDz= dQdxyz( 3,I,J,K)
            DvDx= dQdxyz( 4,I,J,K)
            DvDy= dQdxyz( 5,I,J,K)
            DvDz= dQdxyz( 6,I,J,K)
            DwDx= dQdxyz( 7,I,J,K)
            DwDy= dQdxyz( 8,I,J,K)
            DwDz= dQdxyz( 9,I,J,K)

            S0   = -2.0_8/3.0_8*(dUdx+dVdy+dWdz)
            upup = RmiuT/Ref*(S0+2.0_8*dUdx)-2.0_8/3.0_8*Den*vK
            vpvp = RmiuT/Ref*(S0+2.0_8*dVdy)-2.0_8/3.0_8*Den*vK
            wpwp = RmiuT/Ref*(S0+2.0_8*dWdz)-2.0_8/3.0_8*Den*vK
            upvp = RmiuT/Ref*(      dUdy+dVdx)
            upwp = RmiuT/Ref*(      dUdz+dWdx)
            vpwp = RmiuT/Ref*(      dVdz+dWdy)

            upup = -upup 
            vpvp = -vpvp 
            wpwp = -wpwp 
            upvp = -upvp 
            upwp = -upwp 
            vpwp = -vpwp 

            vab_10 = V(6,I,J,K) / Den           !! kinetic energy

            vab_14  = upup
            vab_15  = vpvp
            vab_16  = wpwp
            vab_17  = upvp
            vab_18  = upwp
            vab_19  = vpwp

            vab_20  = vab_04*vab_04
            vab_21  = vab_05*vab_05
            vab_22  = vab_06*vab_06
            vab_23  = vab_04*vab_05
            vab_24  = vab_05*vab_06
            vab_25  = vab_04*vab_06
            vab_26  = vab_12*vab_12
!!!!!!
             Ub(I,J,K) = ( Ub(I,J,K)*Time + vab_04)/(Time+1.)
             Vb(I,J,K) = ( Vb(I,J,K)*Time + vab_05)/(Time+1.)
             Wb(I,J,K) = ( Wb(I,J,K)*Time + vab_06)/(Time+1.)
             Rb(I,J,K) = ( Rb(I,J,K)*Time + vab_07)/(Time+1.)
            aMb(I,J,K) = (aMb(I,J,K)*Time + vab_08)/(Time+1.)
            rMb(I,J,K) = (rMb(I,J,K)*Time + vab_09)/(Time+1.)
            vKb(I,J,K) = (vKb(I,J,K)*Time + vab_10)/(Time+1.)
            Ptb(I,J,K) = (Ptb(I,J,K)*Time + vab_11)/(Time+1.)
             Pb(I,J,K) = ( Pb(I,J,K)*Time + vab_12)/(Time+1.)
             Tb(I,J,K) = ( Tb(I,J,K)*Time + vab_13)/(Time+1.)
            uub(I,J,K) = (uub(I,J,K)*Time + vab_14)/(Time+1.)
            vvb(I,J,K) = (vvb(I,J,K)*Time + vab_15)/(Time+1.)
            wwb(I,J,K) = (wwb(I,J,K)*Time + vab_16)/(Time+1.)
            uvb(I,J,K) = (uvb(I,J,K)*Time + vab_17)/(Time+1.)
            uwb(I,J,K) = (uwb(I,J,K)*Time + vab_18)/(Time+1.)
            vwb(I,J,K) = (vwb(I,J,K)*Time + vab_19)/(Time+1.)

         !uuaver(I,J,K) = (uuaver(I,J,K)*Time + vab_20)/(Time+1.)
         !vvaver(I,J,K) = (vvaver(I,J,K)*Time + vab_21)/(Time+1.)
         !wwaver(I,J,K) = (wwaver(I,J,K)*Time + vab_22)/(Time+1.)
         !uvaver(I,J,K) = (uvaver(I,J,K)*Time + vab_23)/(Time+1.)
         !vwaver(I,J,K) = (vwaver(I,J,K)*Time + vab_24)/(Time+1.)
         !uwaver(I,J,K) = (uwaver(I,J,K)*Time + vab_25)/(Time+1.)
         !ppaver(I,J,K) = (ppaver(I,J,K)*Time + vab_26)/(Time+1.)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

       IF(Ki >= KI_statis)THEN
        
         Time=float(Ki-KI_statis)
       
         DO K=1,NK1
         DO J=1,NJ1
         DO I=1,NI1
            Den   =  V(1,I,J,K)
            vab_01= V(2,I,J,K) / Den        !! u
            vab_02= V(3,I,J,K) / Den        !! v
            vab_03= V(4,I,J,K) / Den        !! w
            
            vab_05= PP(I,J,K)               !! P
            vab_06= T(I,J,K)                !! T
            vab_07= Den                     !! den
            
            uurms(I,J,K)= (uurms(I,J,K)*Time+ (vab_01-Ub(I,J,K))*(vab_01-Ub(I,J,K)))/(Time+1.)
            vvrms(I,J,K)= (vvrms(I,J,K)*Time+ (vab_02-Vb(I,J,K))*(vab_02-Vb(I,J,K)))/(Time+1.)
            wwrms(I,J,K)= (wwrms(I,J,K)*Time+ (vab_03-Wb(I,J,K))*(vab_03-Wb(I,J,K)))/(Time+1.)
            uvrms(I,J,K)= (uvrms(I,J,K)*Time+ (vab_01-Ub(I,J,K))*(vab_02-Vb(I,J,K)))/(Time+1.)
            uwrms(I,J,K)= (uwrms(I,J,K)*Time+ (vab_01-Ub(I,J,K))*(vab_03-Wb(I,J,K)))/(Time+1.)
            vwrms(I,J,K)= (vwrms(I,J,K)*Time+ (vab_02-Vb(I,J,K))*(vab_03-Wb(I,J,K)))/(Time+1.)
             
            pprms(I,J,K)= (pprms(I,J,K)*Time+ (vab_05-pb(I,J,K))*(vab_05-pb(I,J,K)))/(Time+1.)

         ENDDO
         ENDDO
         ENDDO
      ENDIF
!!!!!!
     
      KZ=MOD(KI,KQ)
      IF(Kz  == 0 )THEN
        OPEN(7,FILE=ThisBlock%FLfld,form="unformatted")
        DO K=1,NK1
        DO J=1,NJ1
        DO I=1,NI1
              WRITE(7)Xc(I,J,K), Yc(I,J,K), Zc(I,J,K),&
     &                Ub(I,J,K), Vb(I,J,K), Wb(I,J,K),&
     &                Rb(I,J,K),aMb(I,J,K),rMb(I,J,K),vKb(I,J,K),&
     &                Ptb(I,J,K), Pb(I,J,K), Tb(I,J,K),&
     &                uub(I,J,K),vvb(I,J,K),wwb(I,J,K),&
     &                uvb(I,J,K),uwb(I,J,K),vwb(I,J,K),&
     &                uurms(I,J,K), vvrms(I,J,K), wwrms(I,J,K),&
     &                uvrms(I,J,K), vwrms(I,J,K), uwrms(I,J,K), pprms(I,J,K)
         ENDDO
         ENDDO
         ENDDO
         CLOSE(7)
        VAR_num =26
        ALLOCATE(varlocation(VAR_num))
        
        varlocation(1:3)=1
        varlocation(4:VAR_num)=0       
        
         NullPtr=0
         NULLCHR   = CHAR(0)
         Debug     = 0
         VIsDouble = 0
!    >'X Y Z U V W R M Mu k P0 P T uu vv ww uv uw vw'
        I = TecIni100('SIMPLE DATASET'//NULLCHR,&
         &'X Y Z U V W R M Mu k P0 P T uu vv ww uv uw vw  uusol vvsol wwsol uvsol vwsol uwsol ppsol'&
         &//NULLCHR,&
         & trim(ThisBlock%FLaver)//NULLCHR,&
         &'.'  //NULLCHR,&
         &             Debug,&
         &             VIsDouble)
        
        I = TecZne100('Simple Zone'//NULLCHR,&
         &             0, &
         &             NI, &
         &             NJ, &
         &             NK, &
         &             0, &
         &             0, &
         &             0, &
         &             1, &
         &   0,0,varlocation,nulp,0)
!
        ALLOCATE(VARs1(NI,NJ,NK))
        ALLOCATE(VARs2(NI,NJ,NK))
        ALLOCATE(VARs3(NI,NJ,NK))
        ALLOCATE(VARs4(NI,NJ,NK))
        ALLOCATE(VARs5(NI,NJ,NK))
        ALLOCATE(VARs6(NI,NJ,NK))
        III=NI*NJ*NK
        ALLOCATE(var_value(III))
      !XYZc
        VARs1(1:NI,1:NJ,1:NK)=XX(1:NI,1:NJ,1:NK)
        VARs2(1:NI,1:NJ,1:NK)=YY(1:NI,1:NJ,1:NK)
        VARs3(1:NI,1:NJ,1:NK)=ZZ(1:NI,1:NJ,1:NK)
!1
        var_value(1:III)=RESHAPE(VARs1(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!2
        var_value(1:III)=RESHAPE(VARs2(1:NI,1:NJ,1:NK),(/III/))
         II      = TecDat100(III,var_value,VIsDouble)
!3
        var_value(1:III)=RESHAPE(VARs3(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
      
        III=NI1*NJ1*NK1
      !Ub Vb Wb  Denb  
        VARs1(1:NI1,1:NJ1,1:NK1)=Ub(1:NI1,1:NJ1,1:NK1)
        VARs2(1:NI1,1:NJ1,1:NK1)=Vb(1:NI1,1:NJ1,1:NK1)
        VARs3(1:NI1,1:NJ1,1:NK1)=Wb(1:NI1,1:NJ1,1:NK1)
        VARs4(1:NI1,1:NJ1,1:NK1)=Rb(1:NI1,1:NJ1,1:NK1)
!4
        var_value(1:III)=RESHAPE(VARs1(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!5
        var_value(1:III)=RESHAPE(VARs2(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!6
        var_value(1:III)=RESHAPE(VARs3(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!7
        var_value(1:III)=RESHAPE(VARs4(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
      !M Mu K P0 PP T
         !var(8)  = aM
         !var(9)  = RmiuT
         !var(10) = vK
         !var(11) = P0
         !var(12) = PP(I,J,K)
         !var(13) =  T(I,J,K)
        VARs1(1:NI1,1:NJ1,1:NK1)=aMb(1:NI1,1:NJ1,1:NK1)
        VARs2(1:NI1,1:NJ1,1:NK1)=rMb(1:NI1,1:NJ1,1:NK1)
        VARs3(1:NI1,1:NJ1,1:NK1)=vKb(1:NI1,1:NJ1,1:NK1)
        VARs4(1:NI1,1:NJ1,1:NK1)=Ptb(1:NI1,1:NJ1,1:NK1)
        VARs5(1:NI1,1:NJ1,1:NK1)=pb(1:NI1,1:NJ1,1:NK1)
        VARs6(1:NI1,1:NJ1,1:NK1)=Tb(1:NI1,1:NJ1,1:NK1)
!8
        var_value(1:III)=RESHAPE(VARs1(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!9
        var_value(1:III)=RESHAPE(VARs2(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!10
        var_value(1:III)=RESHAPE(VARs3(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!11
        var_value(1:III)=RESHAPE(VARs4(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!12
        var_value(1:III)=RESHAPE(VARs5(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!13     
        var_value(1:III)=RESHAPE(VARs6(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
        !uu vv ww
        VARs1(1:NI1,1:NJ1,1:NK1)=uub(1:NI1,1:NJ1,1:NK1)
        VARs2(1:NI1,1:NJ1,1:NK1)=vvb(1:NI1,1:NJ1,1:NK1)
        VARs3(1:NI1,1:NJ1,1:NK1)=wwb(1:NI1,1:NJ1,1:NK1)
        VARs4(1:NI1,1:NJ1,1:NK1)=uvb(1:NI1,1:NJ1,1:NK1)
        VARs5(1:NI1,1:NJ1,1:NK1)=uwb(1:NI1,1:NJ1,1:NK1)
        VARs6(1:NI1,1:NJ1,1:NK1)=vwb(1:NI1,1:NJ1,1:NK1)
!14
        var_value(1:III)=RESHAPE(VARs1(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!15
        var_value(1:III)=RESHAPE(VARs2(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!16
        var_value(1:III)=RESHAPE(VARs3(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!17
        var_value(1:III)=RESHAPE(VARs4(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!18
        var_value(1:III)=RESHAPE(VARs5(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!19     
        var_value(1:III)=RESHAPE(VARs6(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
        !rms
        VARs1(1:NI1,1:NJ1,1:NK1)=uurms(1:NI1,1:NJ1,1:NK1)
        VARs2(1:NI1,1:NJ1,1:NK1)=vvrms(1:NI1,1:NJ1,1:NK1)
        VARs3(1:NI1,1:NJ1,1:NK1)=wwrms(1:NI1,1:NJ1,1:NK1)
        VARs4(1:NI1,1:NJ1,1:NK1)=uvrms(1:NI1,1:NJ1,1:NK1)
        VARs5(1:NI1,1:NJ1,1:NK1)=vwrms(1:NI1,1:NJ1,1:NK1)
        VARs6(1:NI1,1:NJ1,1:NK1)=uwrms(1:NI1,1:NJ1,1:NK1)
!20
        var_value(1:III)=RESHAPE(VARs1(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!21
        var_value(1:III)=RESHAPE(VARs2(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!22
        var_value(1:III)=RESHAPE(VARs3(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!23
        var_value(1:III)=RESHAPE(VARs4(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!24
        var_value(1:III)=RESHAPE(VARs5(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!25     
        var_value(1:III)=RESHAPE(VARs6(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)

        VARs1(1:NI1,1:NJ1,1:NK1)=pprms(1:NI1,1:NJ1,1:NK1)
!26
        var_value(1:III)=RESHAPE(VARs1(1:NI1,1:NJ1,1:NK1),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
            
        I = TecEnd100()

        DEALLOCATE(VARs1)
        DEALLOCATE(VARs2)
        DEALLOCATE(VARs3)
        DEALLOCATE(VARs4)
        DEALLOCATE(VARs5)
        DEALLOCATE(VARs6)
        DEALLOCATE(var_value)
        DEALLOCATE(Varlocation)
    ENDIF
    enddo
END SUBROUTINE
