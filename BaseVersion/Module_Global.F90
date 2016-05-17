!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Global Vars & Parameters
MODULE Global
    USE Block       !in Module_Block.f90
    ! Block,每个计算块的变量，定义在Module_Block.f90文件中。
    IMPLICIT NONE
    SAVE   !!RIGHT?
    INTEGER, PARAMETER:: NULL=10001 !!WHERE to use?
    !!!!!! Storage Dimension，存储的主要变量个数。
  !  INTEGER, PARAMETER:: ML=8
  !  by luoj  Num of chemical reaction equations:7
     !INTEGER, PARAMETER:: ML=15
     INTEGER, PARAMETER:: ML=7
        CHARACTER(LEN=6):: FilePathPrefix !MPI idxxx/

        !!!!!! Wall Surface Dimension，壁面变量
        CHARACTER(LEN=80):: FLsur
    !    character(len=100),allocatable::FL_Res(:)
        !Num_Surface，壁面数目
        INTEGER:: Num_Surface
        
        TYPE SurfStruct
        INTEGER:: NSur, MSur
        REAL, ALLOCATABLE:: X(:,:), Y(:,:), Z(:,:) !Xsur
    END TYPE SurfStruct
    TYPE(SurfStruct), ALLOCATABLE:: AllSurfs(:)
            

   
    
    !!!!!! Block Definitions，网格块的定义
    !Max_Block:最多网格块数； Max_Block_BC：每块的最多边界数；Max_BC_Item：最多边界类型数
    INTEGER::  Max_Block, Max_Block_BC, Max_BC_Item !!MAY NOT BE USED
        
    !!!!!! Block BC Definitions /BC_whole/，块边界定义
    INTEGER::  I_Suf, I_Ibgn, I_Iend, I_Jbgn, I_Jend, I_Kbgn, I_Kend, I_Kind, I_Kindsub, I_Blk_Ex, I_Blkbd_Ex, I_Lcross, I_Lrev_A, I_Lrev_B,I_Period

    !!!!!! Far Physics，远场边界参数
    REAL::    XM,Vxf,Vyf,Vzf,Ref,Xm2,RXM2,PPf,Enf,Rkf,Rof,Rmiuf,Csthlnd,Csth1,Ps00,Ts00 !设定远场参数，马赫数，速度，压力，湍动能，粘性系数等
    REAL::    AoA,Alfa,Aoyaw,Beta_Yaw,PI !设定攻角、侧滑角等
    REAL::    Tiny, Tinf,Wall_Temp !设定来流温度、壁面温度等
    INTEGER:: If_EquT             !是否等温壁
    INTEGER:: IF_SetRkof          !是否自行设定来流湍流度
    REAL::     Rkfset, Rofset,Tuinf,RmiuLmax  !设定的来流湍动能、耗散率、湍流度等
    
    
    
    INTEGER:: IF_RealGas 
    !!!!!! Geometry，几何信息
    REAL:: Xwc,Ywc,Zwc,Area_Ref,C_Ref,B_Ref !网格中心坐标，网格面积，参考弦长，参考展长
    REAL:: Alscale, Q_Dst, Delta_BL         !网格缩放因子，给定边界层厚度


    INTEGER:: IF_WALLF !是否采用壁函数

 !!added by ydd, rotating velocity and coordinate origin
    real::omega(3),Rot_Ori(3)
    integer::nblade,npassage
    real::Pref,Vref,Lref,miuref,Roref
    real::VisRatio,FSTI     !inlet turbulence information, miuT/miu, free stream turbulence intensity
    INTEGER:: Kind_Spatial, Kind_Roe, Kind_Hybrid_Scheme ! 离散格式种类，Roe求解器类型，混合格式种类
    REAL:: F_limiter1,F_limiter2,gama_max,gama_min       ! 混合格式里的参数设置； gama_max(min)，MDCD格式中的控制参数

    REAL, PARAMETER:: cmu=0.09

    INTEGER:: Kind_Dual, Kchld !
    INTEGER:: IF_RK            ! if Runge-Kutta, 1-yes, 0-no-LUSGS
    INTEGER:: IF_Local_Main    !
    REAL:: CFL_Main, DT_LUSGS_Main !
    INTEGER:: IF_Local_Sub         !
    REAL:: CFL_Sub, DT_LUSGS_Sub,CFL_Sub_start,CFL_Sub_end   
    Integer::CFL_Sub_step                                   
    
    REAL:: Prod_Lmt                                 
    INTEGER:: If_compress, Kind_subCC           !if compressible modification and modification type

    !!!!!! Output
    
    INTEGER:: KR,Krr,MA,KI,KI_ini,Itr_Chld,Kind_Dual_RANS !set /KYE/
    INTEGER:: If_Parallel, N_Proc, KO, KQ, Num_Dat,kind_field,if_VTK!set /Sav/
    !!!!!! NS
    INTEGER:: If_Euler, IF_Turb, IF_Precondition  !, Kind_NS_Diff !set /EuNS/
    
    REAL:: PrT !set /PrantlT/
    INTEGER:: If_Favre !set /Favre/
    !INTEGER:: IF_ThinLayer

    !!	  switches for reading
    INTEGER:: If_turb_org, If_transition_org
    !!!!!! Turb
    INTEGER:: KItm, Kind_Model, IF_Dst, Kind_Turbequ_Diff
    INTEGER:: IOrder_MD
    
    !!!!!! Init
    INTEGER:: If_NS_Uniform, If_Turb_Uniform !set /INIFLD/
    
    !!!!!! Aritficial Viscosity
    REAL:: E2, E4
    REAL:: E2I, E2J, E2K, OMGns !set /E22/
    REAL:: E4I, E4J, E4K !set /E44/
    REAL:: E2md, E4md
    REAL:: E2Imd,E2Jmd,E2Kmd,OMGmd
    INTEGER:: IF_ArtiV_Md !set /E22md/
   	REAL:: E4Imd,E4Jmd,E4Kmd !set /E44md/
   	REAL:: CH1,CH2,CH3
    parameter(CH1=3.,CH2=1.,CH3=2.)
   	
   	!!!!!! Turb Equation
   	INTEGER:: Kind_Prod !set /TwoEqu_Prodction/
	INTEGER:: IF_comppress,Kind_CC !set /compressibleCorrection/
	INTEGER:: IF_FixTrans
	REAL:: XFix_TranS !set /fix_transition/
	
	!!!!!! Transition Equation
	INTEGER:: If_Transition,If_fix_Transition
    REAL::    Trans_location
	REAL::     rgf,sigg0,GMXC1,GMXC2,GMXC3,GMXC4,GMXC5,GMXC6,GMXC7,GMXC8,GMXC9

    !!!!!! Hybrid Methods
    INTEGER:: Kind_Hybrid, Kind_SubHY,Ki_Aver, KI_statis,Ki_C, Ksss, Ki_Surfout !surfout?? !set /Hybrd/

	REAL:: Cdes_Omeg, Cdes_Epsi !set /para_DES_hybrid/
	REAL:: Ciddes_Omeg, Ciddes_Epsi !set /para_IDDES_hybrid/
    parameter( Cdes_omeg=0.78, Cdes_epsi=0.61, Ciddes_omeg=0.78, Ciddes_epsi=0.61)
	
    REAL, PARAMETER:: Akapa_Iddes=0.41, Cw=0.15, Ct=1.87, Cl=5.00 !constant ct and cl should be redefined by computing channel flow but it is not carried out in our codes  
	REAL:: Cdd, Css !set /para_zonal_hybrid                                                      

    !!!!!! SAS
    REAL:: Csmag, Cwale !set /para_sas/
    
    
    !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!Real Gas!!!! by luoj
    REAL, PARAMETER:: Rconst=287.0,B1=989.0504 ,B2=-0.009595592 ,B3=0.0001041469 ,B4=-4.433065E-8, B5=5.879263E-12 
    REAL:: TTD,gam0
    !  Chemical reations Vars by luoj
    INTEGER:: If_Chemical, IF_Catalysis,IF_TwoTemp,IF_Heatgas
    REAL, PARAMETER:: R=8314                         ! 摩尔气体常数

    !!!!!!!!!!!!!!!!!!!!!!!!!Turbulence constants
    REAL:: A1,AKapa,Beta_Star  !,Iorder_Md   !/TwoE/
    REAL:: Sigk1,Sigo1,Beta1,Gama1,Cmuc1   !/TwoE1/
    REAL:: Sigk2,Sigo2,Beta2,Gama2,Cmuc2   !/TwoE2/
!added by ydd for turbomachinery
    real:: P0in,T0in,ma_Income,Pouthub
    real:: Inxpos,Outxpos
    integer::IF_Rot,Kind_RotMod
     
     INCLUDE "mpif.h"
    !!!!!!!!!????????????????paras?
    !!!!!! Parallel
      INTEGER, ALLOCATABLE:: M_blockids(:),PBlock(:),M_BlockIDs1(:)  !(Max_block)
      INTEGER:: Num_Block_All
      INTEGER, ALLOCATABLE:: ProcNo_All(:), BlockNo_All(:)    

      INTEGER, PARAMETER:: NVtot=63,MIJK0=600000*NVtot  !200000*NVtot)
	  REAL(kind=4):: A(MIJK0), B(MIJK0)

	  !common /MPItransfer/ A(MIJK0), B(MIJK0)
	  INTEGER:: Status(MPI_STATUS_SIZE)
	  INTEGER:: Myid,Ierr,Mpi_Rc,NumProcs
	  !common /parallel/ myid,ierr,mpi_rc,numprocs,status
      INTEGER*4:: MPI_Temp_Displ(1000), Mpi_Temp_C(1000)
      INTEGER, PARAMETER:: MaxNmpibcpairs = 5000
      INTEGER:: Nmpibcpairs, Mpi_BCpairs(MaxNmpibcpairs, 7)
      INTEGER, PARAMETER:: MaxNmpibcpairs_this = 1000, MaxNumProcs = 1000
      INTEGER:: Nmpibcpairs_this, Nmpibcpairs_all(MaxNumProcs), Mpi_BCpairs_this(MaxNmpibcpairs_this, 7)
      INTEGER:: Mpi_buffer(7*MaxNmpibcpairs),Mpi_buffer_this(7*MaxNmpibcpairs_this)
      INTEGER, PARAMETER:: MaxNmpiTAGs = 5000
      INTEGER:: NmpiTAGs, Mpi_TAGS(MaxNmpiTAGs)
      integer,parameter :: bufferSize=10000000
      real :: bufferForISend(bufferSize)  !! it's better to dynamics allocate this 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!! Block VARs
    type NameAllBlockConn
        integer::NDim(3),Transform(3),revMatrix(3)
        integer :: NCon   !ydd_number of connection
        integer,allocatable :: Prange(:,:),CPrange(:,:)  ! ydd_point range?
        integer :: NBC    !ydd_number of bc
        !! the order of conn and bc is the same as that in ABlock
        integer,allocatable::IJKBC(:,:)
      end type
      type(NameAllBlockConn),allocatable :: AllBlockName(:)

    integer::NAllBlock,NBlock,NumProc,MyProc
    TYPE(BlockStruct), ALLOCATABLE, TARGET:: AllBlocks(:) !(NBlocks?)
    TYPE(BlockStruct), POINTER:: ThisBlock
    
    !!Pointers (Main)
    REAL, POINTER:: XX(:,:,:), YY(:,:,:), ZZ(:,:,:) !/XYZ/(MIJK)
    REAL, POINTER:: Xc(:,:,:), Yc(:,:,:), Zc(:,:,:) !/XXYYZZc/(0MIJK)
    REAL, POINTER:: Dst(:,:,:), ALeng(:,:,:), ASeng(:,:,:) !/Dist/ Dst(MIJK), ALeng(0MIJK), ASeng(MIJK)
    REAL, POINTER:: Vol(:,:,:) !/VLH/(MIJK)
    REAL, POINTER:: SD(:,:,:,:,:) !/CCI/(3,3,0MIJK)
    REAL, POINTER:: Grad(:,:,:,:) !/Grd/(3,0MIJK)
!added by ydd
    real,pointer::thtf(:,:,:,:),thtc(:,:,:)
    real,pointer::rad(:,:,:)
    real,pointer::f_r1(:,:,:),f_r2(:,:,:)
    
    REAL, POINTER:: Alagm(:,:,:,:) !/AGM/(3,0MIJK)
    
    REAL, POINTER:: Dtm(:,:,:) !/DTM/(MIJK)
    REAL, POINTER:: F(:,:,:,:), Q(:,:,:,:) !/FQQ/(ML,0MIJK) MAYBE ML should be at last
    REAL, POINTER:: D(:,:,:,:) !/Fdd/(ML,0MIJK)
    REAL, POINTER:: Src(:,:,:,:), Dsrc(:,:,:,:) !/Source/(2,MIJK)
    
    REAL, POINTER:: V(:,:,:,:) !(ML,-2MIJK)

    REAL, POINTER:: VL(:,:,:,:,:) !(3,15,-2MIJK)
    REAL, POINTER:: VR(:,:,:,:,:) !(3,15,-2MIJK)
    real,pointer::rccl(:,:,:,:),rccr(:,:,:,:)      !by ydd

    !added by ydd
    real,pointer::vibn(:,:,:),vjbn(:,:,:),vkbn(:,:,:),gridV(:,:,:,:,:),radSurf(:,:,:,:)

    REAL, POINTER:: PP(:,:,:) !(-2MIJK)
    REAL, POINTER:: T(:,:,:) !(-1MIJK)
    REAL, POINTER:: Rmiu(:,:,:) !(-1MIJK)
    
    REAL, POINTER:: Rmiudis(:,:,:) !(-1MIJK)
    REAL, POINTER:: DQdxyz(:,:,:,:)  !(18,-1MIJK) 
  
    REAL,POINTER::Ub(:,:,:)           !add by dzw05
    REAL,POINTER::Vb(:,:,:)
    REAL,POINTER::Wb(:,:,:)

    REAL,POINTER::Rb(:,:,:)
    REAL,POINTER::aMb(:,:,:)
    REAL,POINTER::rMb(:,:,:)
    REAL,POINTER::vKb(:,:,:)


    REAL,POINTER::Pb(:,:,:)
    REAL,POINTER::Ptb(:,:,:)
    REAL,POINTER::Tb(:,:,:)

    REAL,POINTER:: uub(:,:,:)
    REAL,POINTER:: vvb(:,:,:)
    REAL,POINTER:: wwb(:,:,:)
    REAL,POINTER:: uvb(:,:,:)
    REAL,POINTER:: uwb(:,:,:)
    REAL,POINTER:: vwb(:,:,:)

    REAL,POINTER:: pprms(:,:,:)
    
    REAL,POINTER:: uurms(:,:,:)
    REAL,POINTER:: vvrms(:,:,:)
    REAL,POINTER:: wwrms(:,:,:)
    REAL,POINTER:: uvrms(:,:,:)
    REAL,POINTER:: vwrms(:,:,:)
    REAL,POINTER:: uwrms(:,:,:)
    !
    REAL, POINTER:: Rds(:,:,:,:) !(3,0MIJK)
    REAL, POINTER:: Wt0(:,:,:,:) !(ML,MIJK)
    REAL, POINTER:: Wt1(:,:,:,:) !(ML,MIJK)
    REAL, POINTER:: Wt2(:,:,:,:) !(ML,MIJK)

    INTEGER, POINTER:: MARKWALLi0(:,:), MARKWALLim(:,:)
    INTEGER, POINTER:: MARKWALLj0(:,:), MARKWALLjm(:,:)
    INTEGER, POINTER:: MARKWALLk0(:,:), MARKWALLkm(:,:)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Turbulence
    REAL, POINTER:: F1(:,:,:), F2(:,:,:)  
    REAL, POINTER:: FunDES(:,:,:)
    REAL, POINTER:: Fscheme(:,:,:)
    REAL, POINTER:: shock(:,:,:)
    

CONTAINS

    SUBROUTINE GLOBAL_SetPointers(Block)
        IMPLICIT NONE
        INTEGER:: Block
        
        !TYPE(BlockStruct), POINTER:: ThisBlock
        
        ThisBlock => AllBlocks(Block)
        XX => ThisBlock%XX
        YY => ThisBlock%YY
        ZZ => ThisBlock%ZZ
        
        Xc => ThisBlock%Xc
        Yc => ThisBlock%Yc
        Zc => ThisBlock%Zc
        
        Dst => ThisBlock%Dst
        Aleng => ThisBlock%Aleng
        ASeng => ThisBlock%ASeng
        Vol => ThisBlock%Vol
        SD => ThisBlock%SD
        Grad => ThisBlock%Grad
!added by ydd
        thtf=>ThisBlock%thtf
        thtc=>ThisBlock%thtc
        rad=>ThisBlock%rad
        f_r1=>ThisBlock%f_r1
        f_r2=>ThisBlock%f_r2

        Alagm => ThisBlock%Alagm
        
        Dtm => ThisBlock%Dtm

        F => ThisBlock%F
        Q => ThisBlock%Q
        D => ThisBlock%D

        Src => ThisBlock%Src
        DSrc => ThisBlock%DSrc
        
        !U => ThisBlock%U
        V => ThisBlock%V
        VL => ThisBlock%VL
        VR => ThisBlock%VR
        rccl=>ThisBLock%rccl  !by ydd
        rccr=>ThisBLock%rccr  !by ydd
!added by ydd
        vibn=>ThisBlock%vibn
        vjbn=>ThisBlock%vjbn
        vkbn=>ThisBlock%vkbn
        gridV=>ThisBlock%gridV
        radSurf=>ThisBlock%radSurf
        PP => ThisBlock%PP
        T => ThisBlock%T
        Rmiu => ThisBlock%Rmiu
        
        Rmiudis => ThisBlock%Rmiudis
        DQdxyz => ThisBlock%DQdxyz

    Ub  =>  ThisBlock%Ub         !add by dzw05
    Vb  =>  ThisBlock%Vb  
    Wb  =>  ThisBlock%Wb  

    Rb  =>  ThisBlock%Rb  
    aMb =>  ThisBlock%aMb 
    rMb =>  ThisBlock%rMb 
    vKb =>  ThisBlock%vKb 


    Pb  =>  ThisBlock%Pb  
    Ptb =>  ThisBlock%Ptb 
    Tb  =>  ThisBlock%Tb  

    uub =>  ThisBlock%uub 
    vvb =>  ThisBlock%vvb 
    wwb =>  ThisBlock%wwb 
    uvb =>  ThisBlock%uvb 
    uwb =>  ThisBlock%uwb 
    vwb =>  ThisBlock%vwb 


            
 !uuaver => ThisBlock%uuaver 
 !vvaver => ThisBlock%vvaver 
 !wwaver => ThisBlock%wwaver 
 !uvaver => ThisBlock%uvaver 
 !uwaver => ThisBlock%uwaver 
 !vwaver => ThisBlock%vwaver 
 !ppaver => ThisBlock%ppaver 
 pprms =>  ThisBlock%pprms
 
 uurms =>  ThisBlock%uurms
 vvrms =>  ThisBlock%vvrms
 wwrms =>  ThisBlock%wwrms
 uvrms =>  ThisBlock%uvrms
 vwrms =>  ThisBlock%vwrms
 uwrms =>  ThisBlock%uwrms 

        Rds => ThisBlock%Rds
        Wt0 => ThisBlock%Wt0
        Wt1 => ThisBlock%Wt1
        Wt2 => ThisBlock%Wt2

        !is this useful?
        MARKWALLi0 => ThisBlock%MARKWALLi0
        MARKWALLim => ThisBlock%MARKWALLim
        MARKWALLj0 => ThisBlock%MARKWALLj0
        MARKWALLjm => ThisBlock%MARKWALLjm
        MARKWALLk0 => ThisBlock%MARKWALLk0
        MARKWALLkm => ThisBlock%MARKWALLkm
        !!! Turbulence
        F1 => ThisBlock%F1
        F2 => ThisBlock%F2
        
        FunDES => ThisBlock%FunDES
        Fscheme => ThisBlock%Fscheme
        shock   => ThisBlock%shock

        
    END SUBROUTINE GLOBAL_SetPointers

END MODULE Global
