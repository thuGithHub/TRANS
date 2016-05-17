SUBROUTINE INPUT_Inititializes
    USE Global
    IMPLICIT NONE
    INTEGER:: iBlock
    INTEGER:: NI1,NJ1,NK1
    INTEGER:: I,J,K,L,II,LL
    INTEGER:: MLTEMP
    REAL:: Coff_dst
    REAL:: Rat
    REAL:: C1,C2,Smax,tmax,aKin,v3d,aOin
    REAL:: VV
    !REAL:: uurms,vvrms,wwrms,uvrms,uwrms,vwrms
    REAL::aaa,bbb,ccc
    real::ptemp    
    
    DO iBlock = 1, Max_Block
        CALL GLOBAL_SetPointers(iBlock)
        NI1=ThisBlock%NI1
        NJ1=ThisBlock%NJ1
        NK1=ThisBlock%NK1

    IF(KR == 1) THEN !read inits
        IF(Num_dat == 0)THEN
            OPEN(UNIT=15, FILE=ThisBlock%FLDAT, MODE='READ', FORM="UNFORMATTED") !, SHARE='DENYWR')
        ELSE !Datb
            OPEN(UNIT=15, FILE=ThisBlock%FLDATb, MODE='READ', FORM="UNFORMATTED") !, SHARE='DENYWR')
        ENDIF
      
        IF(KR == 0 .and. Kind_Hybrid/=0) THEN
            IF(Num_dat == 0)THEN
                OPEN(UNIT=15, FILE=ThisBlock%FLDAT_RANS, MODE='READ', FORM="UNFORMATTED") !, SHARE='DENYWR')
            ELSE !Datb
                OPEN(UNIT=15, FILE=ThisBlock%FLDATb_RANS, MODE='READ', FORM="UNFORMATTED") !, SHARE='DENYWR')
            ENDIF
        ENDIF
      
        read(15) KI_ini   
        read(15) ((((V(LL,I,J,K),K=1,NK1),J=1,NJ1),I=1,NI1),LL=1,ML)
        read(15) (((Rmiu(i,j,k),K=1,NK1),j=1,NJ1),I=1,NI1)
        if(Kind_dual==1)then
            read(15) ((((Wt0(LL,i,j,k),k=1,Nk1),j=1,NJ1),i=1,NI1),LL=1,ML)
        elseif(Kind_dual==2)then
            read(15) ((((Wt0(LL,i,j,k),k=1,Nk1),j=1,NJ1),i=1,NI1),LL=1,ML)
            read(15) ((((Wt1(LL,i,j,k),k=1,Nk1),j=1,NJ1),i=1,NI1),LL=1,ML)
        endif     
        close(15)
 
        KI = KI_ini
        IF (KR == 0 .and. Kind_Hybrid/=0) THEN
            KI = 0 !not affected by RANS method
            KI_ini =0
        ENDIF
    
        IF (Kind_Hybrid /= 0) THEN !time order upgrade
            if (Kind_dual == 1 .and. Kind_dual_RANS == 0 ) then
                DO L=1,ML
                    Wt0(L,I,J,K)=V(L,I,J,K)
                ENDDO
            endif
            if (Kind_dual == 2 .and. Kind_dual_RANS == 0 ) then
                DO L=1,ML
                  Wt0(L,I,J,K)=V(L,I,J,K)
                  Wt1(L,I,J,K)=V(L,I,J,K)
                ENDDO
            endif
            if (Kind_dual == 2 .and. Kind_dual_RANS == 1 ) then
                DO L=1,ML
                    Wt0(L,I,J,K)=Wt0(L,I,J,K)
                    Wt1(L,I,J,K)=Wt0(L,I,J,K)
                ENDDO
            endif
        ENDIF

        IF( Krr == 1 ) THEN
            OPEN(UNIT=21, FILE=ThisBlock%FLfld, MODE='READ', FORM="UNFORMATTED" )!, SHARE='DENYWR')
            DO K=1,NK1
            DO J=1,NJ1
            DO I=1,NI1
                READ(21)aaa,bbb,ccc, Ub(I,J,K), Vb(I,J,K), Wb(I,J,K), &
                    &    Rb(I,J,K),aMb(I,J,K),rMb(I,J,K),vKb(I,J,K), &
                    &  Ptb(I,J,K), Pb(I,J,K), Tb(I,J,K),             &
                    &   uub(I,J,K),vvb(I,J,K),wwb(I,J,K),            &
                    &   uvb(I,J,K),uwb(I,J,K),vwb(I,J,K),            &
                    &   uurms(I,J,K), vvrms(I,J,K), wwrms(I,J,K),    &
                       uvrms(I,J,K), vwrms(I,J,K), uwrms(I,J,K), pprms(I,J,K)                   
            ENDDO
            ENDDO 
            ENDDO
            CLOSE(21)
        ENDIF  !Krr
        !!!!!!
        KI=KI+1

    ELSE !init from beginning

        Coff_dst = 4.6 / Q_dst
        write(*,*) Coff_dst
        DO K=1,NK1
        DO J=1,NJ1
        DO I=1,NI1
            V(1,I,J,K)=1.0
            V(2,I,J,K)=1.0!1.0!vibn(i,j,k)/normalS(1)!gridV(1,i,j,k)       !Vxf*0.5
            V(3,I,J,K)=omega(1)*Zc(i,j,k)!vjbn(i,j,k)/normalS(2)!gridv(2,i,j,k)       !Vyf
            V(4,I,J,K)=-omega(1)*Yc(i,j,k)!vkbn(i,j,k)/normalS(3)!gridv(3,i,j,k)       !Vzf
            ptemp=1.0/(1.4*XM*XM)
            PP(i,j,k)=Ptemp
               
            T(I,J,K)=RXM2*PP(I,J,K)   !1.0
            Rmiu(I,J,K)=Rmiuf
            F1(I,J,K)=1.0
            if(If_NS_uniform == 0)then
               if(Dst(I,J,K) < Q_dst )then
                  Rat=1.-exp(-Coff_dst*Dst(I,J,K))
                  V(2,I,J,K)=V(2,i,j,k)*Rat !Vxf*Rat                      
               endif
            endif
            V(5,I,J,K)=PP(i,j,k)/0.4+0.5*(V(2,i,j,k)**2.0+V(3,i,j,k)**2.0+V(4,i,j,k)**2.0)-  &   !ENF
            &       0.5*(omega(1)*rad(i,j,k))**2.0
            if(Kind_model == 1 .or. Kind_model == 2 .or. Kind_model == 5 .or. Kind_model == 6 .or.Kind_model == 7  )then
                if(If_Turb_uniform == 1)then
                    V(6,I,J,K)=Rkf
                    V(7,I,J,K)=Rof
                endif
                if(If_Turb_uniform == 0)then
                    C1=45.8
                    C2=1.68
                    Smax=C2/(2.*C1)
                    tmax=-C1*Smax*Smax+C2*Smax

                    aKin=amax1(Rkf,-C1*Dst(I,J,K)**2+C2*Dst(I,J,K))
                    v3d =100.*aKin/tmax
                    aOin=amax1(-12444.*Dst(I,J,K)+0.54,aKin/v3d)
                    V(6,I,J,K)=aKin
                    V(7,I,J,K)=aOin
                endif
            endif
            !			dual1->2
            MLtemp = ML

            if (Kind_dual == 1) then
               DO L=1,MLtemp
                  Wt0(L,I,J,K)=V(L,I,J,K)
               ENDDO

            elseif(Kind_dual == 2 ) then
               DO L=1,MLtemp
                  Wt0(L,I,J,K)=V(L,I,J,K)
                  Wt1(L,I,J,K)=V(L,I,J,K)
               ENDDO
            endif
         ENDDO
         ENDDO
         ENDDO
    ENDIF !KR=1/=1
    
        XM2 =Xm*Xm
        RXM2=XM2*1.4    !by ydd
        DO K=1,NK1
        DO J=1,NJ1
        DO I=1,NI1
!            gam(I,J,K)=1.4
            VV=V(2,I,J,K)*V(2,I,J,K)+V(3,I,J,K)*V(3,I,J,K)+V(4,I,J,K)*V(4,I,J,K)
            !PP(I,J,K)=(gam(I,J,K)-1.0)*V(5,I,J,K)-0.5*(gam(I,J,K)-1.0)*(VV/V(1,I,J,K)+ &
            PP(I,J,K)=(1.4-1.0)*V(5,I,J,K)-0.5*(1.4-1.0)*(VV/V(1,I,J,K)+ &  !by ydd
            &   -(omega(1)*rad(i,j,k))**2.0*V(1,I,J,K))
            T(I,J,K)=RXM2*PP(I,J,K)/V(1,I,J,K)
        ENDDO
        ENDDO
        ENDDO
    enddo !do block
END SUBROUTINE
