subroutine ParameterSet
    use global
    implicit none

    PI      = 4*ATAN(1.)
    Alfa    = AoA/180.*PI
    Beta_Yaw= AoYaw/180.*PI

    VxF = COS(Beta_Yaw) * COS(Alfa)
    VyF = COS(Beta_Yaw) * SIN(Alfa)
    VzF = SIN(Beta_Yaw)

    gam0=1.4

    XM2=Xm*Xm;         RXM2=XM2*gam0
    PPF=1./gam0/XM2;    ENF=PPF/(gam0-1.0)+0.5
    !GamaTT=1.4
    Csthlnd = 110.0/ Tinf
    Csth1=Tinf/273.15

    PrT=0.9
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! CENTERED ONLY
    !PARAS AV
    E2I=1./E2;        E2J=1./E2;         E2K=1./E2
    E4I=1./E4;        E4J=1./E4;         E4K=1./E4

    !TURB PARAS AV
    E2Imd=1./E2md;    E2Jmd=1./E2md;     E2Kmd=1./E2md
    E4Imd=1./E4md;    E4Jmd=1./E4md;     E4Kmd=1./E4md
    
    IF (IF_Turb_org) THEN
        a1=0.31;  aKapa=0.41;  Beta_star=0.09
        IF(Kind_Model == 1)THEN   !!!!!SST
            sigk1=0.85;   sigo1=0.500;    Beta1=0.0750;     Gama1=5./9.
            sigk2=1.00;   sigo2=0.856;   Beta2=0.0828;     Gama2=0.440
        ENDIF

        IF(Kind_Model == 2)THEN    !!!WD+2006     
            sigk1=0.6;    sigo1=0.5;      Beta1=0.0708;     Gama1=0.52
            sigk2=0.6;    sigo2=0.5;      Beta2=0.0708;     Gama2=0.52
            PrT  =8./9.
        ENDIF
        !     Far field
      !!!!!! Rmiu
        IF(Kind_model == 1 .or.Kind_model ==2) THEN
            Rmiuf=1.8e-5/1.2    !Rkf/(Rof+Tiny)
        ENDIF
        IF(Kind_Model == 1 .or. Kind_model == 2  )THEN
            Rkf=1.5*FSTI**2.0   !9.E-9
            Rof=1.2*Rkf/Rmiuf   !R1.E-6
!        Rkf=0.0087 !1.05d-7 !9.E-9
!        Rof=0.00087 !2.10d-8 !1.E-6
!        Rgf =0.01_8
!        Rgf=sqrt(rkf)  !by xu
!        Rkf =1.05d-7  !4.86E-6 !1.05E-7
!        Rof =2.10d-8  !3.2076/3.3E6
!        Rgf =0.01_8
        ENDIF
        IF(IF_setRkof == 1) THEN
            Rkf = Rkfset
            Rof = Rofset
            Rgf = 0.01
        ENDIF


    ENDIF
endsubroutine

subroutine NonDimensionlization
    use global
    implicit none
    real::asonic,comm

        T0in=T0in/Tinf
    asonic=sqrt(1.4*Tinf*287)
    Vref=Xm*asonic
    Lref=1.0
    miuref=(1.0+Csthlnd)/(T0in*Csth1+Csthlnd)*(Csth1*T0in)**1.5*1.7e-5
!    miuref=1.8E-5
    Roref=1.2
    Ref=Roref*Vref*Lref/miuref
!    Roref=Ref*miuref/(Vref*Lref)
    Pref=Roref*Vref**2.0
    P0in=P0in/Pref
    Pouthub=Pouthub/Pref
    omega(1:3)=omega(1:3)*2.0*pi*Lref/(60.0*Vref)
    comm=1.0+0.2*XM*XM
    Ps00=P0in/(Comm**3.5)
    Ts00=T0in/Comm
!    if(myid==0)then
!        write(*,*)"Pref,Tref,Ref,Vref=",Pref,Tinf,Ref,Vref
!        write(*,*)"P0in,T0in,Pouthub,omega=",P0in,T0in,Pouthub,omega(1)
!   endif

endsubroutine
