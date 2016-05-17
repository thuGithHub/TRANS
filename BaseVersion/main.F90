PROGRAM TRANSF95
    USE Global
    IMPLICIT NONE

    INTEGER:: KI_org, KI_residual, KI_surf, KZ
    INTEGER:: KKtm ,II,JJ,KK,i,j,k
    INTEGER:: Num_subIte
    INTEGER:: iBlock,mainLoop
    integer::Num_BC,NI,NJ,NK
    integer,pointer::IJKBC(:,:)
   
    Tiny=1.E-20
    KI_ini=0

    CALL INPUT_ReadOrg
    CALL PARALLEL_Initialize
    call ParameterSet
    CALL COMMON_Define_BCinput_parameter 
    CALL INPUT_Read_GlobalInfo          
    call PLOT3D_Read_And_Init
    CALL INPUT_Read_BlockDimensions
    call NonDimensionlization 
    
    CALL OUTPUT_Define_OutputFiles      
    CALL PARALLEL_Dimensions            
    CALL INPUT_ReadGrids                 
    CALL INPUT_ReadSurf                 
    CALL GEOM_CalcDSTs                  
    CALL INPUT_Inititializes            
    call GeomTransform
    call GeomTransform1
           DO iBlock = 1, Max_Block
                CALL GLOBAL_SetPointers(iBlock)
                    if(ThisBlock%NBlockGlobal==2)then
                        do k=0,0
                        do j=1,10
                        do i=1,10
!                            write(*,*)"Geo",i,j,k,thtc(i,j,k),rad(i,j,k),Dst(i,j,k)
!                            write(*,*)"SD",SD(1,1:3,i,j,k)
                        enddo
                        enddo
                        enddo
                    endif
            ENDDO
   
 
    IF(KR == 0) KI=1
    do MainLoop=KI,MA
        KI_org=MOD(KI,10)
        IF(KI_org == 0)THEN         !org file reread every 10 steps
            CALL INPUT_ReadOrg
        ENDIF

        if (Kind_dual /= 0) Num_subIte = Kchld
        if (Kind_dual == 0) Num_subIte = 1
        if (IF_RK == 1)   Num_subIte = 3       

        !Reset when IF_TURB/IF_TRANS changed
        IF_turb=0;    IF_transition=0
        IF (IF_euler ==0 .and. IF_turb_org == 1 .and. (KI >= KItm .or. Kind_hybrid > 0)) THEN
            IF_turb = 1
        ENDIF


        DO Itr_Chld=1,Num_subIte
            DO iBlock = 1, Max_Block           
                CALL GLOBAL_SetPointers(iBlock) 
                CALL TIME_SpectralRadius
            ENDDO
             !BC transfer
             ! vars of 0-order:  U V W P T k omega muT
            CALL BC_Boundarys(-1)    ! boundary in cells center
            CALL PARALLEL_BCs(-1)
            DO iBlock = 1, Max_Block
                CALL GLOBAL_SetPointers(iBlock) 
                CALL SPACEINTERP_GetSurfCenterVars  ! primitive values on the surface
            ENDDO
            DO iBlock = 1, Max_Block
                CALL GLOBAL_SetPointers(iBlock)
                if(If_euler == 0)then
                    CALL DIFF_DqDxyz    ! gradient, including gradients values in ghost cells center
                endif
            ENDDO
            CALL BC_Boundarys(0)     
            CALL PARALLEL_BCs(0)    

            do iblock=1,Max_block
                call global_SetPointers(iblock)
                call RCCorrection
            enddo

            DO iBlock = 1, Max_Block
                CALL GLOBAL_SetPointers(iBlock)
                D=0.
                if(If_euler == 0) then
                    CALL DIFFUSION_Centered   ! 
                endif
                if(If_TURB) then
                    CALL TURBTWO_DIFFUSION_Centered   ! 
                endif
            ENDDO       

            do iblock=1,Max_Block   !spatial discretization
                call Global_SetPointers(iblock)
                select case(Kind_Spatial)
                case (1)
                    call CONVECT_MUSCL2
                case (2)
                    call CONVECT_MUSCL
                case (3)
                    call CONVECT_MDCD
                case (4)
                    call CONVECT_MDCDHYB
                case default
                    write(*,*)MyId,"Unsupported Scheme!"
                endselect 
                if(IF_TURB) call TURBTWO_CONVECT_MUSCL
            enddo

            DO iBlock = 1, Max_Block
               CALL GLOBAL_SetPointers(iBlock)
                SELECT CASE (Kind_Roe)
                CASE (1)
                    CALL CONVECT_Roe_Original
                CASE (2)
                    CALL CONVECT_Roe_Rotated
                case(3)
                    call CONVECT_AUSMPW_Plus   !AUSMPW+
                END SELECT
                IF (IF_Turb) THEN
                    CALL TURBTWO_CONVECT_Upwind
                ENDIF
            ENDDO

                
            IF(IF_RK == 0)then   !Time Marching, LUSGS        
                DO iBlock = 1, Max_Block
                    CALL GLOBAL_SetPointers(iBlock)
                    CALL TIME_LUSGS             
                    IF (IF_turb) THEN
                        CALL TURBTWO_TIME_LUSGS_SST
                    ENDIF
                    CALL TIME_Update_V          
                    KI_residual = MOD(KI,KO)          
                    if (KI_residual == 0) then
                        CALL OUTPUT_Residual_NS  ! 
                        if (IF_Turb) then
                            call TURBTWO_OUTPUT_Residual  !
                        endif
                    endif
                ENDDO 
            ELSEIF(IF_RK==1) then   !Rk marching
                DO iBlock = 1, Max_Block
                    CALL GLOBAL_SetPointers(iBlock)
                    CALL TIME_RK_Save
                    CALL TIME_RK 
                    IF (IF_turb) THEN
                        CALL TURBTWO_TIME_LUSGS_SST  !
                    ENDIF
                    CALL TIME_RK_V_Update_New               
               
                    KI_residual = MOD(KI,KO)          
                    if (KI_residual == 0) then
                            CALL OUTPUT_Residual_NS  
                        if (IF_Turb) then
                            call TURBTWO_OUTPUT_Residual  
                        endif
                    endif   
                ENDDO
            ENDIF
        ENDDO 
    

        KI_surf=MOD(KI,KQ)      !OutPut
        IF (KI_surf == 0 ) THEN
            CALL OUTPUT_Surfaces  !OUTPUT_Surfaces.F90
        ELSEIF(KI >= KI_c) THEN
            IF (MOD(KI, Ksss) == 0) THEN
                CALL OUTPUT_Surfaces
            ENDIF
        ENDIF
        CALL output_pressure_samples    
        CALL OUTPUT_ClCdCms             !
        
        if(KI >= Ki_aver)THEN   
            CALL Average_Fields       
            CALL Average_Surfaces      
        endif
        
        KI=KI+1       ! main iteration
    !!!!!!
        if(Kind_dual /= 0 ) then   
            DO iBlock = 1, Max_Block
                CALL GLOBAL_SetPointers(iBlock)
                CALL TIME_Update_U  
            ENDDO
        endif
      
        KZ=MOD(KI,KQ)
        IF(KZ == 0) THEN
            WRITE(*,*)"SAVING..."
            CALL OUTPUT_InstFields  
            !call UVWp_blk001
            CALL OUTPUT_SaveFields  
            WRITE(*,*)"SAVE DONE...!"
        ELSEIF(KI >= KI_c) THEN
            IF (MOD(KI, Ksss) == 0) THEN
                CALL OUTPUT_InstFields
            ENDIF
        ENDIF
    enddo  
        
    CALL PARALLEL_Finish
END PROGRAM TRANSF95
