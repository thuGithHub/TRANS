MODULE Block
    IMPLICIT NONE

    type :: ConnectivityStruct
      integer::lcross, lreva,lrevb
      integer :: ConType, PerType,BCType,SubBCType
      real  ::    property(10)
      integer :: Transform(3),revMatrix(3),TransformT(3)
      integer :: pointStart(3),pointEnd(3),CPointStart(3),CPointEnd(3),IJKDirect, LeftOrRight
      integer :: IJKAdd(3),  CIJKAdd(3), CEIJKAdd(3),CIJKAdd1(3),CIJKAdd2(3), &
            &   CEIJKAdd1(3),CEIJKAdd2(3),IJKAdd1(3),IJKAdd2(3) 
      integer :: TargetBlock, TargetConn
      integer :: TargetPStart(3),TargetPEnd(3), TPLijk(3),  &
              &   TargetPLD(3),TargetPRU(3), Offset(3),COffset(3), &
              &   CTargetPStart(3),CTargetPEnd(3), TCLijk(3)
                 !! CLijk: side length i,j,k direct, for mpi_recv bufr length
!                !! targetBlock side start and end, necessary for data transform between cells and sides
!                !! because pointStart<pointEnd, (PointStart,PointEnd) also represent cell range,
!                !! so no need for CpointStart,CPointEnd
      !! periodic connectivity
      real :: rotateCenter(3),rotateAngle(3),Translation(3)

      end type ConnectivityStruct

    TYPE BlockStruct

        !!!!!! Dimensions
        INTEGER:: MI, MJ, MK
        !INTEGER:: MIJK !MAY not needed
        integer::NBlockGlobal 
        !!!!!! Dimension Indicators (USED in geom?)
        INTEGER:: NIp2,NJp2,NKp2
        INTEGER:: NIp1,NJp1,NKp1
        INTEGER:: NI  ,NJ  ,NK, NDim(3)  
        INTEGER:: NI1, NJ1, NK1
        INTEGER:: NI2, NJ2, NK2
        INTEGER:: NI3, NJ3, NK3

        !!!!!! Output File Names
        CHARACTER (100) :: FLGRD, FLdst, FLgeo  !/FLgeo/
        CHARACTER (100) :: FLHIS,FLHISk,FLHISo,FLHISg,FLres,FLpre !/FLproc/
        CHARACTER (100) :: FLHISm !/FLprsa/
        CHARACTER (100) :: FLDAT,FLDATb,FLDAT_rans,FLDATb_rans !/FLdat/
        CHARACTER (100) :: FLPLT ,FLPLT_point !/FLins/ instant
        CHARACTER (100) :: FLaver,FLfld             !/FLins/ aver
        
        !!!!!! Surface Paras and File Names
        INTEGER:: Num_Wall,Num_sample
        INTEGER:: ID_samp
        INTEGER, ALLOCATABLE:: Id_Sur(:), IJK_Sur(:) !(Nblk_wall?)
        INTEGER, ALLOCATABLE:: IJK_BE(:,:) !(Nblk_wall?,6)
        INTEGER, ALLOCATABLE:: Isam(:),  Jsam(:), Ksam(:) !sample point 
        CHARACTER (LEN=100), ALLOCATABLE :: FL_Sur(:), FL_Cldm(:),FL_Suraver(:),FL_sam(:)       !(Nblk_wall?) !!FORMAT right? 
        REAL,ALLOCATABLE::CLDM(:,:,:)
        REAL,ALLOCATABLE:: SAMPLE(:,:,:)  !!!den, vk,pre,temp

        !!!!!! Surf paras for MPI?
        INTEGER:: ID_Present_Blk, Num_Block_BC
        INTEGER, ALLOCATABLE:: IJKBC(:,:) !(Max_Block_BC?, Max_BC_item?)
        
    !    !!!!!! Geometry
        REAL, ALLOCATABLE:: XX(:,:,:), YY(:,:,:), ZZ(:,:,:) !/XYZ/(MIJK)
        REAL, ALLOCATABLE:: Xc(:,:,:), Yc(:,:,:), Zc(:,:,:) !/XXYYZZc/(0MIJK)
        REAL, ALLOCATABLE:: Dst(:,:,:), ALeng(:,:,:), ASeng(:,:,:) !/Dist/ Dst(MIJK), ALeng(0MIJK), ASeng(MIJK)
        REAL, ALLOCATABLE:: Vol(:,:,:) !/VLH/(MIJK)
    !added by ydd, angles information
        real,allocatable::thtf(:,:,:,:),thtc(:,:,:)
        real,allocatable::rad(:,:,:)
        real,allocatable::f_r1(:,:,:),f_r2(:,:,:)   !for roatation modification, by ydd

        REAL, ALLOCATABLE:: SD(:,:,:,:,:) !/CCI/(3,3,0MIJK)
        REAL, ALLOCATABLE:: Grad(:,:,:,:) !/Grd/(3,0MIJK)
        REAL, ALLOCATABLE:: Alagm(:,:,:,:) !/AGM/(3,0MIJK)
        
!        REAL, ALLOCATABLE:: S1_I0(:,:), S2_I0(:,:), S3_I0(:,:) !/RNImin/(MJK)
!        REAL, ALLOCATABLE:: S1_Im(:,:), S2_Im(:,:), S3_Im(:,:) !/RNImax/(MJK)
!        REAL, ALLOCATABLE:: S1_J0(:,:), S2_J0(:,:), S3_J0(:,:) !/RNJmin/(MIK)
!        REAL, ALLOCATABLE:: S1_Jm(:,:), S2_Jm(:,:), S3_Jm(:,:) !/RNJmax/(MIK)
!        REAL, ALLOCATABLE:: S1_K0(:,:), S2_K0(:,:), S3_K0(:,:) !/RNKmin/(MIJ)
!        REAL, ALLOCATABLE:: S1_Km(:,:), S2_Km(:,:), S3_Km(:,:) !/RNKmax/(MIJ)
!        REAL, ALLOCATABLE:: Ti0(:,:), Tim(:,:) !/RNNi/(MJK)
!        REAL, ALLOCATABLE:: Tj0(:,:), Tjm(:,:) !/RNNj/(MIK)
!        REAL, ALLOCATABLE:: Tk0(:,:), Tkm(:,:) !/RNNk/(MIJ)
        !	wall surface marks
        INTEGER, ALLOCATABLE:: MARKWALLi0(:,:), MARKWALLim(:,:)
        INTEGER, ALLOCATABLE:: MARKWALLj0(:,:), MARKWALLjm(:,:)
        INTEGER, ALLOCATABLE:: MARKWALLk0(:,:), MARKWALLkm(:,:)
	  
	  !COMMON /WALLBCMARKi_blk001/ MARKWALLi0(MJ,MK), MARKWALLim(MJ,MK)
	  !COMMON /WALLBCMARKj_blk001/ MARKWALLj0(MI,MK), MARKWALLjm(MI,MK)
	  !COMMON /WALLBCMARKk_blk001/ MARKWALLk0(MI,MJ), MARKWALLkm(MI,MJ)
        
        !!!!!! Time Step
        REAL, ALLOCATABLE:: Dtm(:,:,:) !/DTM/(MIJK)

        !!!!!! Main Working (Flux?) VARs
        REAL, ALLOCATABLE:: F(:,:,:,:), Q(:,:,:,:) !/FQQ/(ML,0MIJK) MAYBE ML should be at last
        REAL, ALLOCATABLE:: D(:,:,:,:) !/Fdd/(ML,0MIJK)
        
        REAL, ALLOCATABLE:: Src(:,:,:,:), Dsrc(:,:,:,:) !/Source/(2,MIJK)
        
        !!!!!! Main Physics/Gradient VARs
!        REAL, ALLOCATABLE:: U(:,:,:,:) !(ML,-2MIJK)
        REAL, ALLOCATABLE:: V(:,:,:,:) !(ML,-2MIJK)

    !added by ydd grdi velocity
    real,allocatable::vibn(:,:,:),vjbn(:,:,:),vkbn(:,:,:),gridV(:,:,:,:,:),radSurf(:,:,:,:)
	    REAL, ALLOCATABLE:: VL(:,:,:,:,:)
        REAL, ALLOCATABLE:: VR(:,:,:,:,:) 
        real,allocatable:: rccl(:,:,:,:),rccr(:,:,:,:) !by ydd, for (wr)^2       
 
        REAL, ALLOCATABLE:: PP(:,:,:) !(-2MIJK)
        
        REAL, ALLOCATABLE:: T(:,:,:) !(-1MIJK)
        REAL, ALLOCATABLE:: Rmiu(:,:,:) !(-1MIJK)
        
        
        REAL, ALLOCATABLE:: Rmiudis(:,:,:) !(-1MIJK)

    REAL, ALLOCATABLE:: DQdxyz(:,:,:,:) !(21,-1MIJK)

    REAL,ALLOCATABLE::Ub(:,:,:)           !add by dzw05
    REAL,ALLOCATABLE::Vb(:,:,:)
    REAL,ALLOCATABLE::Wb(:,:,:)

    REAL,ALLOCATABLE::Rb(:,:,:)
    REAL,ALLOCATABLE::aMb(:,:,:)
    REAL,ALLOCATABLE::rMb(:,:,:)
    REAL,ALLOCATABLE::vKb(:,:,:)


    REAL,ALLOCATABLE::Pb(:,:,:)
    REAL,ALLOCATABLE::Ptb(:,:,:)
    REAL,ALLOCATABLE::Tb(:,:,:)

    REAL,ALLOCATABLE:: uub(:,:,:)
    REAL,ALLOCATABLE:: vvb(:,:,:)
    REAL,ALLOCATABLE:: wwb(:,:,:)
    REAL,ALLOCATABLE:: uvb(:,:,:)
    REAL,ALLOCATABLE:: uwb(:,:,:)
    REAL,ALLOCATABLE:: vwb(:,:,:)


    	     
    !REAL,ALLOCATABLE:: uuaver(:,:,:)
    !REAL,ALLOCATABLE:: vvaver(:,:,:)
    !REAL,ALLOCATABLE:: wwaver(:,:,:)
    !REAL,ALLOCATABLE:: uvaver(:,:,:)
    !REAL,ALLOCATABLE:: uwaver(:,:,:)
    !REAL,ALLOCATABLE:: vwaver(:,:,:)
    !REAL,ALLOCATABLE:: ppaver(:,:,:)
    REAL,ALLOCATABLE:: pprms(:,:,:)
    
    REAL,ALLOCATABLE:: uurms(:,:,:)
    REAL,ALLOCATABLE:: vvrms(:,:,:)
    REAL,ALLOCATABLE:: wwrms(:,:,:)
    REAL,ALLOCATABLE:: uvrms(:,:,:)
    REAL,ALLOCATABLE:: vwrms(:,:,:)
    REAL,ALLOCATABLE:: uwrms(:,:,:)

    REAL, ALLOCATABLE:: Rds(:,:,:,:) !(3,0MIJK)

        !!!!!! Residuals
        INTEGER:: Ima,Jma,Kma
        REAL:: DUa,DRSa
        !INTEGER:: Ima1,Jma1,Kma1
        !INTEGER:: Ima2,Jma2,Kma2
        !INTEGER:: Ima3,Jma3,Kma3
        !INTEGER:: Ima4,Jma4,Kma4
        !INTEGER:: Ima5,Jma5,Kma5
        !REAL:: DUa1,DRSa1
        !REAL:: DUa2,DRSa2
        !REAL:: DUa3,DRSa3
        !REAL:: DUa4,DRSa4
        !REAL:: DUa5,DRSa5

        !!!!!! Dual Time
        REAL, ALLOCATABLE:: Wt0(:,:,:,:) !(ML,MIJK)
        REAL, ALLOCATABLE:: Wt1(:,:,:,:) !(ML,MIJK)
        REAL, ALLOCATABLE:: Wt2(:,:,:,:) !(ML,MIJK)
!        REAL, ALLOCATABLE:: Wt2(:,:,:,:) !(ML,MIJK)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! Turbulence

        REAL, ALLOCATABLE:: F1(:,:,:), F2(:,:,:) ! COMMON /BlendF_blk001/ F1( 0:MI, 0:MJ, 0:MK), F2( 0:MI, 0:MJ, 0:MK) 
        !hybrid

        !hybrid roe scheme
        REAL, ALLOCATABLE:: Fscheme(:,:,:),FunDES(:,:,:)
        REAL, ALLOCATABLE:: shock(:,:,:)

        INTEGER:: Imk,Jmk,Kmk  !,DUk,DRSk   !/Residual_k_blk001/
        REAL:: DUk,DRSk
        INTEGER:: Imo,Jmo,Kmo  !,DUo,DRSo   !/Residual_o_blk001/
        REAL:: DUo,DRSo
        INTEGER:: Imm,Jmm,Kmm  !,DUm,DRSm   !/Residual_m_blk001/
        REAL:: DUm,DRSm
        
        INTEGER:: Img,Jmg,Kmg  !,DUg,DRSg   !/Residual_g_blk001/
        REAL:: DUg,DRSg
        !COMMON /Residual_gm_blk001/ Img,Jmg,Kmg,DUg,DRS
        integer :: NConnect
        type(ConnectivityStruct),allocatable :: AConnectivity(:)

    END TYPE BlockStruct    
    
END MODULE Block
