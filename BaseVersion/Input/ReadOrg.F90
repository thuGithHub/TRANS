!!!!!!!!!!!!Read Org File!!!!!!!!!!!!!!!11
SUBROUTINE INPUT_ReadOrg
    USE Global
    IMPLICIT NONE
    real::refOrg
    
    INTEGER:: TEMP
    
        
    OPEN(UNIT=1,FILE= 'original.org', MODE='READ')  !, SHARE='DENYNONE') !single file
    
    READ(1,*)
    READ(1,*) KR,MA,Kind_RotMod                 !
    READ(1,*)
    READ(1,*)
    READ(1,*) IF_parallel,KO,KQ,Num_dat,kind_field,if_VTK
    READ(1,*)
    READ(1,*)
    READ(1,*) IF_Dst,Alscale,XWc,YWc,ZWc,Area_Ref,C_Ref,B_Ref 
    READ(1,*)
    READ(1,*)
    READ(1,*) XM,AoA,AoYaw,RefOrg,Tinf,If_EquT,Wall_Temp,Delta_BL   !replace Ref with RefOrg, by ydd
    READ(1,*)
    READ(1,*)
    READ(1,*) IF_NS_Uniform,Q_Dst,VisRatio,FSTI
    READ(1,*)
    READ(1,*)
    READ(1,*) Kind_Spatial,Kind_Roe,E2,E4,F_limiter1,F_limiter2,gama_max,gama_min,IF_Precondition 
    READ(1,*)
    READ(1,*)
    READ(1,*) DT_LUSGS_Main,CFL_Sub_start,CFL_Sub_end,CFL_Sub_step
    READ(1,*)
    READ(1,*)
    READ(1,*) Kind_Dual,Kchld,  IF_RK 
    READ(1,*)
    READ(1,*) 
    READ(1,*) IF_Euler,IF_Turb_Org,IF_setRkof, Rkfset, Rofset
    READ(1,*)
    READ(1,*)
    READ(1,*) KItm,Kind_Model,Prod_Lmt,If_wallF
    READ(1,*)
    READ(1,*)
    READ(1,*) IF_Comppress,Kind_CC
    READ(1,*)
    READ(1,*)
    READ(1,*) Kind_Hybrid,Kind_SubHY
    READ(1,*)
    READ(1,*)
    READ(1,*) IF_TRANSITION_org,Tuinf,RmiuLmax
    READ(1,*)
    READ(1,*)
    READ(1,*) If_fix_Transition,Trans_location
    READ(1,*)
    READ(1,*)
    READ(1,*) Krr,Ki_Aver,KI_statis,KI_C,KSss
    
    ClOSE(1)

    IF_Local_Main =0
    IF_Local_Sub  =1
    CFL_Main      =1.
    Kind_Dual_RANS = 2

END SUBROUTINE INPUT_ReadOrg
