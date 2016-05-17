SUBROUTINE OUTPUT_InstFields
    USE Global
    IMPLICIT NONE
    INTEGER:: iBlock,i,j,k
    character(len=100)::filename
  
!    call CornerPoints
    if(if_VTK==0)then 
        DO iBlock = 1, Max_Block
            CALL GLOBAL_SetPointers(iBlock)
            filename=trim(ThisBlock%FLplt)//".plt"
            CALL OUTPUT_InstField_Tec(filename)
        ENDDO
    elseif(if_VTK==1)then
        DO iBlock = 1, Max_Block
            CALL GLOBAL_SetPointers(iBlock)
            filename=trim(ThisBlock%FLplt)//".vtk"
            CALL OUTPUT_InstField_VTK(filename)
        ENDDO
!    elseif(if_VTK==2)then
!        DO iBlock = 1, Max_Block
!            CALL GLOBAL_SetPointers(iBlock)
!            call output_tec_point
!        ENDDO
    endif
        
        

END SUBROUTINE OUTPUT_InstFields


SUBROUTINE OUTPUT_InstField_Tec(filename)
    USE Global
    IMPLICIT NONE
    include "../tecio.F90"
 
    integer::iBlock   
    INTEGER*4:: I,J,K,L,M,N,iqq,jqq,kqq
    INTEGER*4:: NI, NJ, NK
    INTEGER*4:: NI1, NJ1, NK1
    INTEGER*4:: NI2, NJ2, NK2
    integer::KIout
      CHARACTER(LEN=100):: OutputFile,filename
      CHARACTER(LEN=6):: cc
      character(LEN=1):: NULLCHR
      Integer*4::   Debug,II, III,VIsDouble
      !real*4:: var(100)
      Integer*4::Is

      REAL*4::Den,UU,VV,WW,vK,Uabs,Vabs,Wabs
      REAL*4::asonic,comm,Ps,Ts,P0rel,P0,T0rel,T0,Mach,Mabs
      REAL*4::Ox,Oy,Oz,Omg,Yita4
      REAL*4::DuDx,DuDy,DuDz
      REAL*4::DvDx,DvDy,DvDz
      REAL*4::DwDx,DwDy,DwDz
      REAL*4::W12,W13,W23
      REAL*4::S11,S22,S33,S12,S13,S23
      REAL*4::Vort,Skl
      REAL*4::Qv
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Integer*4:: Iav,Jav,Kav
        real::count
      REAL*4::Den_av,UU_av,VV_av,WW_av,vK_av,Uabs_av,Vabs_av,Wabs_av
      real*4::ma_av,Mabs_av,Ps_av,P0rel_av,P0_av,Ts_av,T0rel_av,T0_av
      REAL*4::Ox_av,Oy_av,Oz_av,Omg_av,Qv_av,Rmiu_av,shock_av,yita4_av
      REAL*4::angular,Vol_Wei
      real::Sij(3,3),Wij(3,3)
      real::Tuu,Tuv,Tuw,Tvw,Tvv,Tww,Tuu_av,Tuv_av,Tuw_av,Tvw_av,Tvv_av,Tww_av,div,TmiuT,frr,frr_av
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

    angular=omega(1)
    NullPtr=0
    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1
         
    outputfile=trim(filename)
         !IF(Kind_Model==1.and.Kind_hybrid==0 )THEN
             !IF(KI >= KI_c)THEN
                !DO Is=0,999
                  !IF( KI == INT(KI_c+Is*Ks) )THEN
                        !WRITE(cc, '(I3.3)') Is
                        !WRITE(cc, '(I6.6)') KI

                        !WRITE(ThisBlock%FLPLT, '(A,A,A,A)') trim(FilePathPrefix), 'resu/RASST_', cc,".PLT"
                        !WRITE(ThisBlock%FLPLT, '(A,A,A,A)') trim(FilePathPrefix), 'resu/RASST_', cc,".PLT"
                  !ENDIF
                !ENDDO
             !ENDIF
         !ENDIF
!        OutputFile = trim(ThisBlock%FLPLT)
!         IF(Kind_hybrid==0 )THEN
!             IF(KI >= KI_c)THEN
                !DO Is=0,999
                  !IF( KI == INT(KI_c+Is*Ks) )THEN
                        !WRITE(cc, '(I3.3)') Is
                        !WRITE(cc, '(I6.6)') KI
!                        WRITE(cc, '(I6.6)') (KI-KI_c)/Ksss

                        !WRITE(ThisBlock%FLPLT, '(A,A,A,A)') trim(FilePathPrefix), 'resu/RASST_', cc,".PLT"
!                        WRITE(OutputFile, '(A,A,A,A)'), trim(ThisBlock%FLPLT), '_C', cc, ".PLT"
                  !ENDIF
                !ENDDO
!             ENDIF
!         ENDIF
        if(KI>KI_c)then
            KIout=KI-KI_c
!            if(mod(KIout,Ksss)==0)then
                write(cc,'(I6.6)')KIout/Ksss
!            endif
            write(OutputFile, '(A,A,A,A)'), trim(ThisBlock%FLPLT), '_C', cc, ".PLT"
        endif
!    WRITE(*,*) 'Output:', OutputFile
        ALLOCATE(VARs1(NI,NJ,NK))
        ALLOCATE(VARs2(NI,NJ,NK))
        ALLOCATE(VARs3(NI,NJ,NK))
        ALLOCATE(VARs4(NI,NJ,NK))
        ALLOCATE(VARs5(NI,NJ,NK))
        ALLOCATE(VARs6(NI,NJ,NK))
        VAR_num =34
        ALLOCATE(varlocation(VAR_num))
        III=NI*NJ*NK
        ALLOCATE(var_value( III))
        
        NULLCHR   = CHAR(0)
        Debug     = 0
        VIsDouble = 0   !!!
        varlocation(1:3)=1
        varlocation(4:VAR_num)=0
        if(kind_field ==1)  varlocation(4:VAR_num)=1
!        I = TecIni100('SIMPLE DATASET'//NULLCHR, &
!     &  'X Y Z  U V W R M Mu k P0 P T Ox Oy Oz Om Q gm Rmiudis Shock gam' &
!     &  'X Y Z U V W R M Mu k P0 P T Ox Oy Oz Om Q gm Rmiudis Shock gam R1 R2 R3 R4 R5 & 
!     &  R6 Tve Mabs Pt Tt' !        

        I = TecIni100('SIMPLE DATASET'//NULLCHR, &
        &   'X Y Z density U V W Uabs Vabs Wabs Ma Mabs Ps P0rel P0 Ts T0rel T0 Ox Oy Oz Om &
        &   Q Yita4 Rmiu k shock fr uu vv ww uv uw vw' &
        &                  //NULLCHR, &
        &             trim(OutputFile)//NULLCHR, &
        &             '.'  //NULLCHR, &
        &             Debug, &
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
!       &             'POINT'//NULLCHR, &
!       &             NULLCHR) 
!1
        var_value(1:III)=RESHAPE(XX(1:NI,1:NJ,1:NK),(/III/))    !XX
        II      = TecDat100(III,var_value,VIsDouble)
!2    
        var_value(1:III)=RESHAPE(YY(1:NI,1:NJ,1:NK),(/III/))    !YY
        II      = TecDat100(III,var_value,VIsDouble)
!3
        var_value(1:III)=RESHAPE(ZZ(1:NI,1:NJ,1:NK),(/III/))    !ZZ
        II      = TecDat100(III,var_value,VIsDouble)
      
        if(kind_field ==0) then        ! center dataset  
            III=NI1*NJ1*NK1
        do K=1,NK1
        do J=1,NJ1
        do I=1,NI1
            Den=V(1,i,j,k)
            UU=V(2,i,j,k)/Den
            VV=V(3,i,j,k)/Den
            WW=V(4,i,j,k)/Den
            Uabs=UU
            Vabs=VV-angular*Zc(i,j,k)
            Wabs=WW+angular*Yc(i,j,k)
            
            Vars1(i,j,k)=UU
            Vars2(i,j,k)=VV
            Vars3(i,j,k)=WW
            Vars4(i,j,k)=Uabs
            Vars5(i,j,k)=Vabs
            Vars6(i,j,k)=Wabs
        enddo
        enddo
        enddo
!4
        var_value(1:III)=RESHAPE(V(1,1:NI1,1:NJ1,1:NK1),(/III/))    !density
        var_value=var_value*Roref
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!5
        var_value(1:III)=RESHAPE(VARs1(1:NI1,1:NJ1,1:NK1),(/III/))  !relative U
        var_value=var_value*Vref
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!6
        var_value(1:III)=RESHAPE(VARs2(1:NI1,1:NJ1,1:NK1),(/III/))  !relative V
        var_value=var_value*Vref
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!7
        var_value(1:III)=RESHAPE(VARs3(1:NI1,1:NJ1,1:NK1),(/III/))  !relative w
        var_value=var_value*Vref
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!8
        var_value(1:III)=RESHAPE(VARs4(1:NI1,1:NJ1,1:NK1),(/III/))  !absolute V
        var_value=var_value*Vref
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!9      
        var_value(1:III)=RESHAPE(VARs5(1:NI1,1:NJ1,1:NK1),(/III/))  !absolute V
        var_value=var_value*Vref
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!10        
        var_value(1:III)=RESHAPE(VARs6(1:NI1,1:NJ1,1:NK1),(/III/))  !absolute V
        var_value=var_value*Vref
        II      = TecDat100(III,var_value(1:III),VIsDouble)
      
        do K=1,NK1
        do J=1,NJ1
        do I=1,NI1
            Den    = V(1,I,J,K)
            UU=Vars1(i,j,k)
            VV=Vars2(i,j,k)
            WW=Vars3(i,j,k)
            Uabs=Vars4(i,j,k)
            Vabs=Vars5(i,j,k)
            Wabs=Vars6(i,j,k)
            asonic=sqrt(1.4*PP(i,j,k)/Den)
            Vars1(i,j,k)=sqrt(UU*UU+VV*VV+WW*WW)/asonic
            Vars2(i,j,k)=sqrt(Uabs**2.0+Vabs**2.0+Wabs**2.0)
        enddo
        enddo
        enddo

!11
        var_value(1:III)=RESHAPE(VARs1(1:NI1,1:NJ1,1:NK1),(/III/))  !relative Ma
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!12
        var_value(1:III)=RESHAPE(VARs2(1:NI1,1:NJ1,1:NK1),(/III/))  !absoulute Ma
        II      = TecDat100(III,var_value(1:III),VIsDouble)
            
        do k=1,NK1
        do j=1,NJ1
        do i=1,NI1
            Mach=Vars1(i,j,k)
            Mabs=Vars2(i,j,k)
            comm=(1.0+0.2*Ma*Ma)
            Vars3(i,j,k)=PP(i,j,k)*(comm**3.5)
            Vars5(i,j,k)=T(i,j,k)*comm
            comm=1.0+Mabs*Mabs
            Vars4(i,j,k)=PP(i,j,k)*(comm*3.5)
            Vars6(i,j,k)=T(i,j,k)*comm
        enddo
        enddo
        enddo

!13
        var_value(1:III)=RESHAPE(PP(1:NI1,1:NJ1,1:NK1),(/III/))  ! Ps
        var_value=var_value*Pref
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!14
        var_value(1:III)=RESHAPE(VARs3(1:NI1,1:NJ1,1:NK1),(/III/))  !relative total P
        var_value=var_value*Pref
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!15
        var_value(1:III)=RESHAPE(VARs4(1:NI1,1:NJ1,1:NK1),(/III/))  !absolute total P
        var_value=var_value*Pref
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!16
        var_value(1:III)=RESHAPE(T(1:NI1,1:NJ1,1:NK1),(/III/))  !Ts
        var_value=var_value*Tinf
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!17
        var_value(1:III)=RESHAPE(VARs5(1:NI1,1:NJ1,1:NK1),(/III/))  !relative total T
        var_value=var_value*Tinf
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!18
        var_value(1:III)=RESHAPE(VARs6(1:NI1,1:NJ1,1:NK1),(/III/))  !absolute total T
        var_value=var_value*Tinf
        II      = TecDat100(III,var_value(1:III),VIsDouble)
      
        do K=1,NK1
        do J=1,NJ1
        do I=1,NI1
            Dudx = DqDxyz(1,I,J,K)
            Dudy = DqDxyz(2,I,J,K)
            Dudz = DqDxyz(3,I,J,K)
            Dvdx = DqDxyz(4,I,J,K)
            Dvdy = DqDxyz(5,I,J,K)
            Dvdz = DqDxyz(6,I,J,K)
            Dwdx = DqDxyz(7,I,J,K)
            Dwdy = DqDxyz(8,I,J,K)
            Dwdz = DqDxyz(9,I,J,K)
    
            Sij(1,1)=DuDx
            Sij(2,2)=DvDy
            Sij(3,3)=DwDz
            Sij(1,2)=0.5*(DuDy+DvDx)
            Sij(1,3)=0.5*(DwDx+DuDz)
            Sij(2,3)=0.5*(DvDz+DwDy)
            Sij(2,1)=Sij(1,2)
            Sij(3,1)=Sij(1,3)
            Sij(3,2)=Sij(2,3)
            Wij(1,1)=0.0
            Wij(2,2)=0.0
            Wij(3,3)=0.0
            Wij(1,2)=0.5*(DuDy-DvDx)
            Wij(1,3)=0.5*(DuDz-DwDx)
            Wij(2,3)=0.5*(DvDz-DwDy)
            Wij(2,1)=-Wij(1,2)
            Wij(3,1)=-Wij(1,3)
            Wij(3,2)=-Wij(2,3)
            Ox=Wij(3,2)
            Oy=Wij(1,3)
            Oz=Wij(2,1)
            Omg    =2.0*sqrt(Ox*Ox+Oy*Oy+Oz*Oz)
            ! Qv=0.5*(|WijWij|-|SijSij|)
                Qv=0.0
                do jqq=1,3
                do iqq=1,3
                    Qv=Qv+Wij(iqq,jqq)*Wij(iqq,jqq)-Sij(iqq,jqq)*Sij(iqq,jqq)
                enddo
                enddo
                Qv=0.5*Qv
                ! Yita4=SijWjkWki
                Yita4=0.0
                do kqq=1,3
                do jqq=1,3
                do iqq=1,3
                    Yita4=Yita4+Sij(iqq,jqq)*Wij(jqq,kqq)*Wij(kqq,iqq)
                enddo
                enddo
                enddo
            
            VARs1(I,J,K) = Ox
            VARs2(I,J,K) = Oy
            VARs3(I,J,K) = Oz
            VARs4(I,J,K) = Omg
            VARs5(I,J,K) = Qv
            VARs6(I,J,K) =Yita4
        enddo
        enddo
        enddo
        
!19
        var_value(1:III)=RESHAPE(VARs1(1:NI1,1:NJ1,1:NK1),(/III/))  !Ox
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!20
        var_value(1:III)=RESHAPE(VARs2(1:NI1,1:NJ1,1:NK1),(/III/))  !Oy
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!21
        var_value(1:III)=RESHAPE(VARs3(1:NI1,1:NJ1,1:NK1),(/III/))  !Oz
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!22
        var_value(1:III)=RESHAPE(VARs4(1:NI1,1:NJ1,1:NK1),(/III/))  !Omg
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!23
        var_value(1:III)=RESHAPE(VARs5(1:NI1,1:NJ1,1:NK1),(/III/))  !Q
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!24      
        var_value(1:III)=RESHAPE(VARs6(1:NI1,1:NJ1,1:NK1),(/III/))  !Q
        II      = TecDat100(III,var_value(1:III),VIsDouble)
        do K=1,NK1
        do J=1,NJ1
        do I=1,NI1
            Den    = V(1,I,J,K)
            VARs1(I,J,K) = Rmiu(i,j,k)
            Vars2(i,j,k)=V(6,i,j,k)/Den
            VARs3(I,J,K) = shock(I,J,K)
            VARs4(I,J,K)=f_r1(I,J,K)
        enddo
        enddo
        enddo
!25
        var_value(1:III)=RESHAPE(VARs1(1:NI1,1:NJ1,1:NK1),(/III/))  !Rmiu
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!26
        var_value(1:III)=RESHAPE(VARs2(1:NI1,1:NJ1,1:NK1),(/III/))  !k
        II      = TecDat100(III,var_value(1:III),VIsDouble)
!27     
        var_value(1:III)=RESHAPE(VARs3(1:NI1,1:NJ1,1:NK1),(/III/))  !shock
        II      = TecDat100(III,var_value(1:III),VIsDouble)

        var_value(1:III)=RESHAPE(VARs4(1:NI1,1:NJ1,1:NK1),(/III/))  !shock
        II      = TecDat100(III,var_value(1:III),VIsDouble)

        do K=1,NK1
        do J=1,NJ1
        do I=1,NI1
            Dudx = DqDxyz(1,I,J,K)
            Dudy = DqDxyz(2,I,J,K)
            Dudz = DqDxyz(3,I,J,K)
            Dvdx = DqDxyz(4,I,J,K)
            Dvdy = DqDxyz(5,I,J,K)
            Dvdz = DqDxyz(6,I,J,K)
            Dwdx = DqDxyz(7,I,J,K)
            Dwdy = DqDxyz(8,I,J,K)
            Dwdz = DqDxyz(9,I,J,K)
            TmiuT=Rmiu(i,j,k)
            vK=V(6,i,j,k)
            div=dudx+dvdy+dwdz
    
            Tuu=(2.0*dudx-2.0*div/3.0)*TmiuT-2.0*vK/3.0
            Tvv=(2.0*dvdy-2.0*div/3.0)*TmiuT-2.0*vK/3.0
            Tww=(2.0*dwdz-2.0*div/3.0)*TmiuT-2.0*vk/3.0
            Tuv=(Dudy+dVdx)*TmiuT
            Tuw=(Dudz+dwdx)*TmiuT
            Tvw=(Dvdz+dwdy)*TmiuT
            VARs1(i,j,k)=Tuu
            VARs2(i,j,k)=Tvv
            VARs3(i,j,k)=Tww
            VARs4(i,j,k)=Tuv
            VARs5(i,j,k)=Tuw
            VARs6(i,j,k)=Tvw
        enddo
        enddo
        enddo
        var_value(1:III)=RESHAPE(VARs1(1:NI1,1:NJ1,1:NK1),(/III/))  !shock
        II      = TecDat100(III,var_value(1:III),VIsDouble)
            
        var_value(1:III)=RESHAPE(VARs2(1:NI1,1:NJ1,1:NK1),(/III/))  !shock
        II      = TecDat100(III,var_value(1:III),VIsDouble)

        var_value(1:III)=RESHAPE(VARs3(1:NI1,1:NJ1,1:NK1),(/III/))  !shock
        II      = TecDat100(III,var_value(1:III),VIsDouble)

        var_value(1:III)=RESHAPE(VARs4(1:NI1,1:NJ1,1:NK1),(/III/))  !shock
        II      = TecDat100(III,var_value(1:III),VIsDouble)

        var_value(1:III)=RESHAPE(VARs5(1:NI1,1:NJ1,1:NK1),(/III/))  !shock
        II      = TecDat100(III,var_value(1:III),VIsDouble)

        var_value(1:III)=RESHAPE(VARs6(1:NI1,1:NJ1,1:NK1),(/III/))  !shock
        II      = TecDat100(III,var_value(1:III),VIsDouble)

    else              ! vertex dataset 
        III=NI*NJ*NK    

        do k=1,NK
        do j=1,NJ
        do i=1,NI
            den_av=0.0
            count=0.0
            do kav=k-1,k
            do jav=j-1,j
            do iav=i-1,i
                den_av=den_av+V(1,iav,jav,kav)*Vol(Iav,Jav,Kav)
                count=count+Vol(Iav,Jav,Kav)
            enddo
            enddo
            enddo
            Vars1(i,j,k)=den_av/count*Roref
        enddo
        enddo
        enddo
!4
        var_value(1:III)=RESHAPE(VARs1(1:NI,1:NJ,1:NK),(/III/))  !density
        II      = TecDat100(III,var_value,VIsDouble)

        do K=1,NK
        do J=1,NJ
        do I=1,NI
            UU_av = 0.0
            VV_av = 0.0
            WW_av = 0.0
            Uabs_av=0.0
            Wabs_av=0.0
            Vabs_av=0.0
            count = 0.0
            do kav = k-1,k
            do jav = j-1,j
            do iav = i-1,i
                Den=V(1,Iav,Jav,Kav)
        if(den>0.0)then
                vol_wei=Vol(iav,jav,kav)

                UU=V(2,Iav,Jav,Kav)/Den
                VV=V(3,Iav,Jav,Kav)/Den
                WW=V(4,Iav,Jav,Kav)/Den
                UU_av=UU_av+UU*Vol_wei
                VV_av=VV_av+VV*Vol_wei
                WW_av=WW_av+WW*Vol_wei
                Uabs_av=Uabs_av+UU*Vol_wei
                Vabs_av=Vabs_av+(VV-angular*Zc(iav,jav,kav))*Vol_wei
                Wabs_av=Wabs_av+(WW+angular*Yc(iav,jav,kav))*Vol_wei
                count=count+Vol_wei
        endif
            enddo
            enddo
            enddo
            Vars1(i,j,k)=UU_av/count
            Vars2(i,j,k)=VV_av/count
            Vars3(i,j,k)=WW_av/count
            Vars4(i,j,k)=Uabs_av/count
            Vars5(i,j,k)=Vabs_av/count
            Vars6(i,j,k)=Wabs_av/count
        enddo
        enddo
        enddo
!5
        var_value(1:III)=RESHAPE(VARs1(1:NI,1:NJ,1:NK),(/III/))  !relative U
        var_value=var_value*Vref
        II      = TecDat100(III,var_value,VIsDouble)
!6
        var_value(1:III)=RESHAPE(VARs2(1:NI,1:NJ,1:NK),(/III/))  !relative V
        var_value=var_value*Vref
        II      = TecDat100(III,var_value,VIsDouble)
!7
        var_value(1:III)=RESHAPE(VARs3(1:NI,1:NJ,1:NK),(/III/))  !relative w
        var_value=var_value*Vref
        II      = TecDat100(III,var_value,VIsDouble)
!8
        var_value(1:III)=RESHAPE(VARs4(1:NI,1:NJ,1:NK),(/III/))  !absolute V
        var_value=var_value*Vref
        II      = TecDat100(III,var_value,VIsDouble)
!9      
        var_value(1:III)=RESHAPE(VARs5(1:NI,1:NJ,1:NK),(/III/))  !absolute V
        var_value=var_value*Vref
        II      = TecDat100(III,var_value,VIsDouble)
!10        
        var_value(1:III)=RESHAPE(VARs6(1:NI,1:NJ,1:NK),(/III/))  !absolute V
        var_value=var_value*Vref
        II      = TecDat100(III,var_value,VIsDouble)

        do k=1,NK
        do j=1,NJ
        do i=1,NI
            Ma_av=0.0
            Mabs_av=0.0
            P0rel_av=0.0
            P0_av=0.0
            T0rel_av=0.0
            T0_av=0.0
            count=0.0
            do kav=k-1,k
            do jav=j-1,j
            do iav=i-1,i
                den=V(1,iav,jav,kav)
        if(den>0.00001)then
                Vol_wei=Vol(iav,jav,kav)

                UU=V(2,iav,jav,kav)/den
                VV=V(3,iav,jav,kav)/den
                WW=V(4,iav,jav,kav)/den
                Uabs=UU
                Vabs=VV-angular*Zc(iav,jav,kav)
                Wabs=WW+angular*Yc(iav,jav,kav)
                asonic=sqrt(1.4*PP(iav,jav,kav)/den)

                Mach=sqrt(UU*UU+VV*VV+WW*WW)/asonic
                comm=1.0+0.2*Mach*Mach
                Ma_av=Ma_av+Mach*Vol_wei
                P0rel_av=P0rel_av+PP(iav,jav,kav)*comm**3.5*Vol_wei
                T0rel_av=T0rel_av+T(iav,jav,kav)*comm*Vol_wei

                Mabs=sqrt(Uabs*Uabs+Vabs*Vabs+Wabs*Wabs)/asonic
                Mabs_av=Mabs_av+Mabs*Vol_wei
                comm=1.0+0.2*Mabs*Mabs
                P0_av=P0_av+PP(iav,jav,kav)*comm**3.5*Vol_wei
                T0_av=T0_av+T(iav,jav,kav)*comm*Vol_wei
    
                count=count+Vol_wei
        endif
            enddo
            enddo
            enddo
                Vars1(i,j,k)=Ma_av/count
                Vars2(i,j,k)=Mabs_av/count
                Vars3(i,j,k)=P0rel_av/count*Pref
                Vars4(i,j,k)=P0_av/count*Pref
                Vars5(i,j,k)=T0rel_av/count*Tinf
                Vars6(i,j,k)=T0_av/count*Tinf
        enddo
        enddo
        enddo
!11
        var_value(1:III)=RESHAPE(VARs1(1:NI,1:NJ,1:NK),(/III/))  !relative Ma
        II      = TecDat100(III,var_value,VIsDouble)
!12
        var_value(1:III)=RESHAPE(VARs2(1:NI,1:NJ,1:NK),(/III/))  !absolute Ma
        II      = TecDat100(III,var_value,VIsDouble)

        do k=1,NK
        do j=1,NJ
        do i=1,NI
            Ps_av=0.0
            Ts_av=0.0
            count=0.0
            do kav=k-1,k
            do jav=j-1,j
            do iav=i-1,i
                vol_wei=vol(iav,jav,kav)
                Ps_av=Ps_av+PP(iav,jav,kav)*Vol_wei
                Ts_av=Ts_av+T(iav,jav,kav)*Vol_wei
                count=count+Vol_wei
            enddo
            enddo
            enddo
            Vars1(i,j,k)=Ps_av/count*Pref
            Vars2(i,j,k)=Ts_av/count*Tinf
        enddo
        enddo
        enddo
!13
        var_value(1:III)=RESHAPE(Vars1(1:NI,1:NJ,1:NK),(/III/))  ! Ps
        II      = TecDat100(III,var_value,VIsDouble)
!14
        var_value(1:III)=RESHAPE(VARs3(1:NI,1:NJ,1:NK),(/III/))  !relative total P
        II      = TecDat100(III,var_value,VIsDouble)
!15
        var_value(1:III)=RESHAPE(VARs4(1:NI,1:NJ,1:NK),(/III/))  !absolute total P
        II      = TecDat100(III,var_value,VIsDouble)
!16
        var_value(1:III)=RESHAPE(Vars2(1:NI,1:NJ,1:NK),(/III/))  !Ts
        II      = TecDat100(III,var_value,VIsDouble)
!17
        var_value(1:III)=RESHAPE(VARs5(1:NI,1:NJ,1:NK),(/III/))  !relative total T
        II      = TecDat100(III,var_value,VIsDouble)
!18
        var_value(1:III)=RESHAPE(VARs6(1:NI,1:NJ,1:NK),(/III/))  !absolute total T
        II      = TecDat100(III,var_value,VIsDouble)

        !Qs
        do K=1,NK
        do J=1,NJ
        do I=1,NI
            Ox_av  =0.0
            Oy_av  =0.0
            Oz_av  =0.0
            Omg_av =0.0
            Qv_av  =0.0
            Yita4_av=0.0
            Tuu_av=0.0
            Tuv_av=0.0
            Tuw_av=0.0
            Tvw_av=0.0
            Tvv_av=0.0
            Tww_av=0.0

            count  =0.0
            do kav = k-1,k
            do jav = j-1,j
            do iav = i-1,i
                vol_wei=Vol(iav,jav,kav)
                Den =    V(1,Iav,Jav,Kav)
                Dudx = DqDxyz(1,Iav,Jav,Kav)
                Dudy = DqDxyz(2,Iav,Jav,Kav)
                Dudz = DqDxyz(3,Iav,Jav,Kav)
                Dvdx = DqDxyz(4,Iav,Jav,Kav)
                Dvdy = DqDxyz(5,Iav,Jav,Kav)
                Dvdz = DqDxyz(6,Iav,Jav,Kav)
                Dwdx = DqDxyz(7,Iav,Jav,Kav)
                Dwdy = DqDxyz(8,Iav,Jav,Kav)
                Dwdz = DqDxyz(9,Iav,Jav,Kav)
                Sij(1,1)=DuDx
                Sij(2,2)=DvDy
                Sij(3,3)=DwDz
                Sij(1,2)=0.5*(DuDy+DvDx)
                Sij(1,3)=0.5*(DwDx+DuDz)
                Sij(2,3)=0.5*(DvDz+DwDy)
                Sij(2,1)=Sij(1,2)
                Sij(3,1)=Sij(1,3)
                Sij(3,2)=Sij(2,3)
                Wij(1,1)=0.0
                Wij(2,2)=0.0
                Wij(3,3)=0.0
                Wij(1,2)=0.5*(DuDy-DvDx)
                Wij(1,3)=0.5*(DuDz-DwDx)
                Wij(2,3)=0.5*(DvDz-DwDy)
                Wij(2,1)=-Wij(1,2)
                Wij(3,1)=-Wij(1,3)
                Wij(3,2)=-Wij(2,3)
                Ox=Wij(3,2)
                Oy=Wij(1,3)
                Oz=Wij(2,1)
                Vort=2.0*sqrt(Ox*Ox+Oy*Oy+Oz*Oz)
                ! Qv=0.5*(|WijWij|-|SijSij|)
                Qv=0.0
                do jqq=1,3
                do iqq=1,3
                    Qv=Qv+Wij(iqq,jqq)*Wij(iqq,jqq)-Sij(iqq,jqq)*Sij(iqq,jqq)
                enddo
                enddo
                Qv=0.5*Qv
                ! Yita4=SijWjkWki
                Yita4=0.0
                do kqq=1,3
                do jqq=1,3
                do iqq=1,3
                    Yita4=Yita4+Sij(iqq,jqq)*Wij(jqq,kqq)*Wij(kqq,iqq)
                enddo
                enddo
                enddo
                ! lamda2 criterion
                ! lamada2=?
                Ox_av=Ox_av+Ox*Vol_wei
                Oy_av=Oy_av+Oy*Vol_wei
                Oz_av=Oz_av+Oz*Vol_wei
                Omg_av=Omg_av+Vort*Vol_wei
                Qv_av=Qv_av+Qv*Vol_wei
                Yita4_av=Yita4_av+Yita4*Vol_wei
                count  = count   +Vol_wei
            enddo
            enddo
            enddo
            VARs1(I,J,K) = Ox_av /count
            VARs2(I,J,K) = Oy_av /count
            VARs3(I,J,K) = Oz_av /count
            VARs4(I,J,K) = Omg_av/count
            VARs5(I,J,K) = Qv_av /count
            VARs6(i,j,k) = Yita4_av/count
        enddo
        enddo
        enddo
!19
        var_value(1:III)=RESHAPE(VARs1(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!20
        var_value(1:III)=RESHAPE(VARs2(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!21
        var_value(1:III)=RESHAPE(VARs3(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!22
        var_value(1:III)=RESHAPE(VARs4(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!23
        var_value(1:III)=RESHAPE(VARs5(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!24
        var_value(1:III)=RESHAPE(VARs6(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
        do K=1,NK
        do J=1,NJ
        do I=1,NI
            Rmiu_av=0.0
            vK_av=0.0
            shock_av  =0.0
            frr_av=0.0
            count     =0
            do kav = k-1,k
            do jav = j-1,j
            do iav = i-1,i
                den=V(1,iav,jav,kav)
        if(den>0.0000)then
                Vol_wei=Vol(iav,jav,kav)
                Rmiu_av=Rmiu_av+Rmiu(iav,jav,kav)*Vol_wei
                vK_av=vK_av+V(6,iav,jav,kav)/den*Vol_wei
                shock_av  = shock_av+shock(Iav,Jav,Kav)*Vol_wei
                frr_av=frr_av+f_r1(Iav,Jav,Kav)*Vol_wei
                count  = count   +   Vol_wei
        endif
            enddo
            enddo
            enddo
            VARs1(I,J,K) = Rmiu_av/count
            VARs2(I,J,K) = vK_av/count
            VARs3(I,J,K) = shock_av /count
            VARs4(i,j,k)=frr_av/count
        enddo
        enddo
        enddo
!25
        var_value(1:III)=RESHAPE(VARs1(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!26
        var_value(1:III)=RESHAPE(VARs2(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
!27    
        var_value(1:III)=RESHAPE(VARs3(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
        
        var_value(1:III)=RESHAPE(VARs4(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
        
        do K=1,NK
        do J=1,NJ
        do I=1,NI
            Tuu_av=0.0
            Tuv_av=0.0
            Tuw_av=0.0
            Tvw_av=0.0
            Tvv_av=0.0
            Tww_av=0.0
            count  =0.0
            
            do kav=k-1,k    
            do jav=j-1,j
            do iav=i-1,i
                den=V(1,Iav,Jav,Kav)
                Dudx = DqDxyz(1,Iav,Jav,Kav)
                Dudy = DqDxyz(2,Iav,Jav,Kav)
                Dudz = DqDxyz(3,Iav,Jav,Kav)
                Dvdx = DqDxyz(4,Iav,Jav,Kav)
                Dvdy = DqDxyz(5,Iav,Jav,Kav)
                Dvdz = DqDxyz(6,Iav,Jav,Kav)
                Dwdx = DqDxyz(7,Iav,Jav,Kav)
                Dwdy = DqDxyz(8,Iav,Jav,Kav)
                Dwdz = DqDxyz(9,Iav,Jav,Kav)
                
                Vol_wei=Vol(iav,jav,kav)
                div=dudx+dvdy+dwdz
                TmiuT=Rmiu(iav,jav,kav)
                vK=V(6,iav,jav,kav)

            !    Tuu=TmiuT*(2.0*DuDx-2.0*div/3.0)*Vref*miuref-2.0*vK/3.0*roref*Vref*Vref
                Tuu=TmiuT*(2.0*DuDx-2.0*div/3.0)/Ref-2.0*vK/3.0
                Tuu_av=Tuu_av+Tuu*Vol_wei/den
            !    Tuu_av=Tuu_av+TmiuT*(2.0*Dudx-2.0*div/3.0)-2.0*vK
            !    Tvv=TmiuT*(2.0*Dvdy-2.0*div/3.0)*Vref*miuref-2.0*vK/3.0*roref*Vref*Vref
                Tvv=TmiuT*(2.0*Dvdy-2.0*div/3.0)/Ref-2.0*vK/3.0
                Tvv_av=Tvv_av+Tvv*Vol_wei/den
            !    Tvv_av=Tvv_av+TmiuT*(2.0*DvDy-2.0*div/3.0)-2.0*vK
            !    Tww=TmiuT*(2.0*DwDz-2.0*div/3.0)*Vref*miuref-2.0*vK/3.0*roref*Vref*Vref
                Tww=TmiuT*(2.0*DwDz-2.0*div/3.0)/Ref-2.0*vK/3.0
                Tww_av=Tww_av+Tww*Vol_wei/den
            !    Tww_av=Tww_av+TmiuT*(2.0*DwDz-2.0*div/3.0)-2.0*vK
                Tuv_av=Tuv_av+(DuDy+DvDx)*TmiuT*Vol_wei/den
                Tuw_av=Tuw_av+(DuDz+DwDx)*TmiuT*Vol_wei/den
                Tvw_av=Tvw_av+(DvDz+DwDy)*TmiuT*Vol_wei/den
                
                count=count+Vol_wei
            enddo
            enddo
            enddo
            VARs1(i,j,k)=-Tuu_av/count
            VARs2(i,j,k)=-Tvv_av/count
            VARs3(i,j,k)=-Tww_av/count
            VARs4(i,j,k)=-Tuv_av/count/Ref
            VARs5(i,j,k)=-Tuw_av/count/Ref
   !         VARs6(i,j,k)=Tvw_av/count/Ref
            VARs6(i,j,k)=-0.5*(Tuu_av+Tvv_av+Tww_av)/count
        enddo
        enddo
        enddo
! for reynolds stress output, by ydd
        var_value(1:III)=RESHAPE(VARs1(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
        
        var_value(1:III)=RESHAPE(VARs2(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
        
        var_value(1:III)=RESHAPE(VARs3(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
        
        var_value(1:III)=RESHAPE(VARs4(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
        
        var_value(1:III)=RESHAPE(VARs5(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
        
        var_value(1:III)=RESHAPE(VARs6(1:NI,1:NJ,1:NK),(/III/))
        II      = TecDat100(III,var_value,VIsDouble)
        

    endif
    I = TecEnd100()
    DEALLOCATE(VARs1)
    DEALLOCATE(VARs2)
    DEALLOCATE(VARs3)
    DEALLOCATE(VARs4)
    DEALLOCATE(VARs5)
    DEALLOCATE(VARs6)
    DEALLOCATE(var_value)
    DEALLOCATE(Varlocation)

END SUBROUTINE OUTPUT_InstField_Tec
