SUBROUTINE INPUT_Read_BlockDimensions
    USE Global
    IMPLICIT NONE
    INTEGER:: iBlock
    INTEGER:: IJK_wall(6)
    CHARACTER(LEN=80):: well
    CHARACTER(LEN=100):: dimensionFile
    TYPE(BlockStruct), POINTER:: Block
    INTEGER:: NBlkBC, No_BC
    INTEGER, POINTER:: Id_Sur(:), IJK_Sur(:)
    INTEGER, POINTER:: IJK_BE(:,:)
    INTEGER, POINTER:: IJKBC(:,:)
    INTEGER:: i_idxa0,i_idxa1,i_idxb0,i_idxb1
    INTEGER:: N,NN
    INTEGER:: I,J,K

    Max_BC_item=15
    DO iBlock = 1, Max_Block
        Block => AllBlocks(iBlock)
        WRITE (dimensionFile, '(A, A, I3.3, A)') FilePathPrefix, "grd/dimension_blk", M_blockids(iBlock), "_m.dat"
        OPEN (UNIT=13, FILE=dimensionFile, MODE='READ')
        READ(13,*)
        READ(13,*)
        READ(13,*) well, Block%NI, Block%NJ, Block%NK

        Block%NI1=Block%NI-1
        Block%NJ1=Block%NJ-1
        Block%NK1=Block%NK-1
        
        Block%NI2=Block%NI-2
        Block%NJ2=Block%NJ-2
        Block%NK2=Block%NK-2

        Block%NI3=Block%NI-3
        Block%NJ3=Block%NJ-3
        Block%NK3=Block%NK-3

        Block%NIp1=Block%NI+1
        Block%NJp1=Block%NJ+1
        Block%NKp1=Block%NK+1

        Block%NIp2=Block%NI+2
        Block%NJp2=Block%NJ+2
        Block%NKp2=Block%NK+2

        !may not NEEDED
        Block%MI=Block%NIp2
        Block%MJ=Block%NJp2
        Block%MK=Block%NKp2

        !!!!!!
        READ(13,*)
        READ(13,*)
        READ(13,*) well, Block%id_present_blk, Block%Num_block_BC

        ALLOCATE(Block%IJKBC(Block%Num_block_BC, Max_BC_item))
        IJKBC => Block%IJKBC

        read(13,*)
        read(13,*)
        DO NBlkbc = 1, Block%Num_Block_BC
            READ(13,*) well, No_BC, IJKBC(Nblkbc, i_suf),                   &
             &                                          i_idxa0    ,        &
             &                                          i_idxa1    ,        &
             &                                          i_idxb0    ,        &
             &                                          i_idxb1    ,        &
             &                             IJKBC(Nblkbc,i_kind    ),        &
             &                             IJKBC(Nblkbc,i_kindsub ),        &
             &                             IJKBC(Nblkbc,i_blk_ex  ),        &
             &                             IJKBC(Nblkbc,i_blkbd_ex),        &
             &                             IJKBC(Nblkbc,i_lcross  ),        &
             &                             IJKBC(Nblkbc,i_lrev_a  ),        &
             &                             IJKBC(Nblkbc,i_lrev_b  ),        &
             &                             IJKBC(NBlkbc,i_period)

            !IJK From to
            !		 I direction
            IF( IJKBC(Nblkbc,i_suf     ) == 1 .or. IJKBC(Nblkbc,i_suf     ) == 2     ) THEN
              if ( i_idxa0 == 0 ) i_idxa0 = 1
              if ( i_idxa1 == 0 ) i_idxa1 = Block%NJ1
              if ( i_idxb0 == 0 ) i_idxb0 = 1
              if ( i_idxb1 == 0 ) i_idxb1 = Block%NK1

              if (IJKBC(Nblkbc,i_suf     ) == 1      ) then
                  IJKBC(Nblkbc,i_Ibgn    ) =  1
                  IJKBC(Nblkbc,i_Iend    ) =  1
              else                                          
                  IJKBC(Nblkbc,i_Ibgn    ) =  Block%NI1
                  IJKBC(Nblkbc,i_Iend    ) =  Block%NI1
              endif
                    IJKBC(Nblkbc,i_Jbgn) =  i_idxa0
                    IJKBC(Nblkbc,i_Jend) =  i_idxa1
                    IJKBC(Nblkbc,i_Kbgn) =  i_idxb0
                    IJKBC(Nblkbc,i_Kend) =  i_idxb1
            ENDIF

            !		 J direction
            IF( IJKBC(Nblkbc,i_suf     ) == 3 .or. IJKBC(Nblkbc,i_suf     ) == 4     ) THEN
              if ( i_idxa0 == 0 ) i_idxa0 = 1
              if ( i_idxa1 == 0 ) i_idxa1 = Block%NI1
              if ( i_idxb0 == 0 ) i_idxb0 = 1
              if ( i_idxb1 == 0 ) i_idxb1 = Block%NK1

              if (IJKBC(Nblkbc,i_suf     ) == 3      ) then 
                  IJKBC(Nblkbc,i_Jbgn    ) =  1
                  IJKBC(Nblkbc,i_Jend    ) =  1
              else                                          
                  IJKBC(Nblkbc,i_Jbgn    ) =  Block%NJ1
                  IJKBC(Nblkbc,i_Jend    ) =  Block%NJ1
              endif
                    IJKBC(Nblkbc,i_Ibgn    ) =  i_idxa0
                    IJKBC(Nblkbc,i_Iend    ) =  i_idxa1
                    IJKBC(Nblkbc,i_Kbgn    ) =  i_idxb0
                    IJKBC(Nblkbc,i_Kend    ) =  i_idxb1
            ENDIF

            !		 K direction
            IF( IJKBC(Nblkbc,i_suf     ) == 5 .or. IJKBC(Nblkbc,i_suf     ) == 6     ) THEN
              if ( i_idxa0 == 0 ) i_idxa0 = 1
              if ( i_idxa1 == 0 ) i_idxa1 = Block%NI1
              if ( i_idxb0 == 0 ) i_idxb0 = 1
              if ( i_idxb1 == 0 ) i_idxb1 = Block%NJ1

              if (IJKBC(Nblkbc,i_suf     ) == 5      ) then 
                  IJKBC(Nblkbc,i_Kbgn    ) =  1
                  IJKBC(Nblkbc,i_Kend    ) =  1
              else                                              
                  IJKBC(Nblkbc,i_Kbgn    ) =  Block%NK1
                  IJKBC(Nblkbc,i_Kend    ) =  Block%NK1
              endif
                IJKBC(Nblkbc,i_Ibgn    ) =  i_idxa0
                IJKBC(Nblkbc,i_Iend    ) =  i_idxa1
                IJKBC(Nblkbc,i_Jbgn    ) =  i_idxb0
                IJKBC(Nblkbc,i_Jend    ) =  i_idxb1
            ENDIF
 
             if(ijkbc(Nblkbc,i_blk_ex) == 0 )ijkbc(Nblkbc,i_blk_ex) =  Block%id_present_blk !is cut BC
       enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
    !																			!
    !	  MARK of Wall Boundary													!
    !	  useful now? (diffusion, surface integration or ?)						!
    !																			!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
      ALLOCATE(Block%MARKWALLi0(Block%NJ1,Block%NK1))
      ALLOCATE(Block%MARKWALLim(Block%NJ1,Block%NK1))
      ALLOCATE(Block%MARKWALLj0(Block%NI1,Block%NK1))
      ALLOCATE(Block%MARKWALLjm(Block%NI1,Block%NK1))
      ALLOCATE(Block%MARKWALLk0(Block%NI1,Block%NJ1))
      ALLOCATE(Block%MARKWALLkm(Block%NI1,Block%NJ1))
      Block%MARKWALLi0 = 0
      Block%MARKWALLim = 0
      Block%MARKWALLj0 = 0
      Block%MARKWALLjm = 0
      Block%MARKWALLk0 = 0
      Block%MARKWALLkm = 0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
    !	  WALL SURFACE
      read(13,*) 
      read(13,*) 
      read(13,*) well, Block%Num_wall
      
      !Allocate
      ALLOCATE(Block%Id_Sur(Block%Num_wall))
      ALLOCATE(Block%IJK_Sur(Block%Num_wall))
      ALLOCATE(Block%IJK_BE(Block%Num_wall,6))
      Id_Sur => Block%Id_Sur
      IJK_Sur => Block%IJK_Sur
      IJK_BE => Block%IJK_BE
      
      ALLOCATE(Block%FL_Sur(Block%Num_wall))
      ALLOCATE(Block%FL_Cldm(Block%Num_wall))
      ALLOCATE(Block%FL_Suraver(Block%Num_wall))
      ALLOCATE(Block%ClDM(Block%Num_wall,8,0:KQ-1))

      do N=1,Block%Num_wall
        read(13,*) well,(IJK_wall(NN),NN=1,6),Block%FL_sur(N),Block%FL_cldm(N),Block%FL_suraver(N)
        Id_sur(N)=IJK_wall(1)
        IJK_sur(N)=IJK_wall(2)

        !I direction
        if (IJK_sur (N) == 1 .or. IJK_sur (N) == 2 ) then
            if (IJK_sur (N) == 1 )then              
                IJK_BE(N,1) = 1
                IJK_BE(N,2) = 1
            else                                        
                IJK_BE(N,1) = Block%NI1
                IJK_BE(N,2) = Block%NI1
            endif
                IJK_BE(N,3) = IJK_wall (3)
                IJK_BE(N,4) = IJK_wall (4)
                IJK_BE(N,5) = IJK_wall (5)
                IJK_BE(N,6) = IJK_wall (6)
            if(IJK_wall (3) == 0 ) IJK_BE(N,3) = 1
            if(IJK_wall (4) == 0 ) IJK_BE(N,4) = Block%NJ1
            if(IJK_wall (5) == 0 ) IJK_BE(N,5) = 1
            if(IJK_wall (6) == 0 ) IJK_BE(N,6) = Block%NK1
         endif

        !J direction
         if (IJK_sur (N) == 3 .or. IJK_sur (N) == 4 ) then
            if (IJK_sur (N) == 3 )then              
                IJK_BE(N,3) = 1
                IJK_BE(N,4) = 1
            else                            
                IJK_BE(N,3) = Block%NJ1
                IJK_BE(N,4) = Block%NJ1
            endif
                IJK_BE(N,1) = IJK_wall (3)
                IJK_BE(N,2) = IJK_wall (4)
                IJK_BE(N,5) = IJK_wall (5)
                IJK_BE(N,6) = IJK_wall (6)

            if(IJK_wall (3) == 0 ) IJK_BE(N,1) = 1
            if(IJK_wall (4) == 0 ) IJK_BE(N,2) = Block%NI1
            if(IJK_wall (5) == 0 ) IJK_BE(N,5) = 1
            if(IJK_wall (6) == 0 ) IJK_BE(N,6) = Block%NK1
         endif

         !K direction
         if (IJK_sur (N) == 5 .or. IJK_sur (N) == 6 ) then
            if (IJK_sur (N) == 5 )then                  
                IJK_BE(N,5) = 1
                IJK_BE(N,6) = 1
            else                                    
                IJK_BE(N,5) = Block%NK1
                IJK_BE(N,6) = Block%NK1
            endif
                IJK_BE(N,1) = IJK_wall (3)
                IJK_BE(N,2) = IJK_wall (4)
                IJK_BE(N,3) = IJK_wall (5)
                IJK_BE(N,4) = IJK_wall (6)
            if(IJK_wall (3) == 0 ) IJK_BE(N,1) = 1
            if(IJK_wall (4) == 0 ) IJK_BE(N,2) = Block%NI1
            if(IJK_wall (5) == 0 ) IJK_BE(N,3) = 1
            if(IJK_wall (6) == 0 ) IJK_BE(N,4) = Block%NJ1
         endif
  
    ! WALL MARKS !!!!!!!!!!!!!!!!!!!!
        if (IJK_sur(N)==1) then
            do J=IJK_BE(N,3),IJK_BE(N,4)
            do K=IJK_BE(N,5),IJK_BE(N,6)
                Block%MARKWALLi0(J,K)=1
            ENDDO
            ENDDO
        endif
        if (IJK_sur(N)==2) then
            do J=IJK_BE(N,3),IJK_BE(N,4)
            do K=IJK_BE(N,5),IJK_BE(N,6)
                Block%MARKWALLim(J,K)=1
            ENDDO
            ENDDO
        endif
        if (IJK_sur(N)==3) then
            do I=IJK_BE(N,1),IJK_BE(N,2)
            do K=IJK_BE(N,5),IJK_BE(N,6)
                Block%MARKWALLj0(I,K)=1
            ENDDO
            ENDDO
        endif
        if (IJK_sur(N)==4) then
            do I=IJK_BE(N,1),IJK_BE(N,2)
            do K=IJK_BE(N,5),IJK_BE(N,6)
                Block%MARKWALLjm(I,K)=1
            ENDDO
            ENDDO
        endif
        if (IJK_sur(N)==5) then
            do I=IJK_BE(N,1),IJK_BE(N,2)
            do J=IJK_BE(N,3),IJK_BE(N,4)
                Block%MARKWALLk0(I,J)=1
            ENDDO
            ENDDO
        endif
        if (IJK_sur(N)==6) then
            do I=IJK_BE(N,1),IJK_BE(N,2)
            do J=IJK_BE(N,3),IJK_BE(N,4)
                Block%MARKWALLkm(I,J)=1
            ENDDO
            ENDDO
        endif
    enddo !Num_Wall
    !	  sample points     !dzw,20121115
    read(13,*) 
    read(13,*) 
    read(13,*) well, Block%Num_sample
          
          ALLOCATE(Block%Isam(Block%Num_sample))
          ALLOCATE(Block%Jsam(Block%Num_sample)) 
          ALLOCATE(Block%Ksam(Block%Num_sample))
          ALLOCATE(Block%FL_sam(Block%Num_sample))   
          ALLOCATE(Block%SAMPLE(6,Block%Num_sample,0:KQ-1)) 
    !!
!     	  if(Block%Num_sample > Max_sample ) then
 !   	     write(*,*) "Please increase the Max_sample as ",Block%Num_sample, " in blk: ", Block%id_present_blk
  ! 	      pause
  ! 	      endif
          do N=1,Block%Num_sample
             read(13,*)well,Block%ID_samp,Block%Isam(N),Block%Jsam(N),Block%Ksam(N),Block%FL_sam(N)
          enddo
        close(13)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        !!!Allocate VARS
        ALLOCATE(Block%XX(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        ALLOCATE(Block%YY(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        ALLOCATE(Block%ZZ(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        
        ALLOCATE(Block%Xc(0:Block%NIp2, 0:Block%NJp2, 0:Block%NKp2))
        ALLOCATE(Block%Yc(0:Block%NIp2, 0:Block%NJp2, 0:Block%NKp2))
        ALLOCATE(Block%Zc(0:Block%NIp2, 0:Block%NJp2, 0:Block%NKp2))
        
        ALLOCATE(Block%Dst(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        ALLOCATE(Block%Aleng(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        ALLOCATE(Block%ASeng(Block%NIp2, Block%NJp2, Block%NKp2))

        ALLOCATE(Block%Vol(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        ALLOCATE(Block%SD(3, 3, -2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        ALLOCATE(Block%Grad(3, -2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        ALLOCATE(Block%Alagm(3, -2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        ALLOCATE(Block%Dtm(Block%NIp2, Block%NJp2, Block%NKp2))
!added by ydd
        allocate(Block%thtf(3,-2:Block%NIp2,-2:Block%NJp2,-2:Block%NKp2))
        allocate(Block%thtc(-2:Block%NIp2,-2:Block%NJp2,-2:Block%NKp2))
        allocate(Block%rad(-2:Block%NIp2,-2:Block%NJp2,-2:Block%NKp2))
        allocate(Block%f_r1(-2:Block%NIp2,-2:Block%NJp2,-2:Block%NKp2))
        allocate(Block%f_r2(-2:Block%NIp2,-2:Block%NJp2,-2:Block%NKp2))

        ALLOCATE(Block%F(ML, 0:Block%NIp2, 0:Block%NJp2, 0:Block%NKp2))
        ALLOCATE(Block%Q(ML, 0:Block%NIp2, 0:Block%NJp2, 0:Block%NKp2))
        ALLOCATE(Block%D(ML, 0:Block%NIp2, 0:Block%NJp2, 0:Block%NKp2))

        ALLOCATE(Block%Src(2, Block%NIp2, Block%NJp2, Block%NKp2))
        ALLOCATE(Block%Dsrc(2, Block%NIp2, Block%NJp2, Block%NKp2))

        ALLOCATE(Block%V(ML, -2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))

        ALLOCATE(Block%VL(3,8, -2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))   !by ydd
        ALLOCATE(Block%VR(3,8, -2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))   !by ydd
        ALLOCATE(Block%rccl(3, -2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))   !by ydd
        ALLOCATE(Block%rccr(3, -2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))   !by ydd

!added by ydd
        allocate(Block%vibn(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        allocate(Block%vjbn(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        allocate(Block%vkbn(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        allocate(Block%gridV(3,3,-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        allocate(Block%radSurf(3,-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))

        ALLOCATE(Block%PP(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        ALLOCATE(Block%T(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        ALLOCATE(Block%Rmiu(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
        
        ALLOCATE(Block%Rmiudis(-1:Block%NIp2, -1:Block%NJp2, -1:Block%NKp2))
       ALLOCATE(Block%Dqdxyz(18, -1:Block%NIp2, -1:Block%NJp2, -1:Block%NKp2))  !by ydd
     
           ALLOCATE(BLOCK%Ub(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))         !add by dzw05
       ALLOCATE(BLOCK%Vb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE(BLOCK%Wb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))

       ALLOCATE(BLOCK%Rb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE(BLOCK%aMb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE(BLOCK%rMb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE(BLOCK%vKb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))


       ALLOCATE(BLOCK%Pb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE(BLOCK%Ptb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE(BLOCK%Tb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))

       ALLOCATE (BLOCK%uub(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE (BLOCK% vvb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE (BLOCK%wwb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE (BLOCK%uvb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE (BLOCK%uwb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE (BLOCK%vwb(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))


      
       !ALLOCATE (BLOCK%uuaver(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       !ALLOCATE (BLOCK%vvaver(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       !ALLOCATE (BLOCK%wwaver(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       !ALLOCATE (BLOCK%uvaver(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       !ALLOCATE (BLOCK%uwaver(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       !ALLOCATE (BLOCK%vwaver(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       !ALLOCATE (BLOCK% ppaver(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE (BLOCK% pprms(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
   
       ALLOCATE (BLOCK% uurms(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE (BLOCK% vvrms(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE (BLOCK% wwrms(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE (BLOCK% uvrms(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE (BLOCK% uwrms(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       ALLOCATE (BLOCK% vwrms(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))


        ALLOCATE(Block%Rds(3, 0:Block%NIp2, 0:Block%NJp2, 0:Block%NKp2))
        ALLOCATE(Block%Wt0(ML, Block%NIp2, Block%NJp2, Block%NKp2))
        ALLOCATE(Block%Wt1(ML, Block%NIp2, Block%NJp2, Block%NKp2))
        ALLOCATE(Block%Wt2(ML, Block%NIp2, Block%NJp2, Block%NKp2))


        Block%XX=0.;  Block%YY=0.;  Block%ZZ=0.
        Block%Xc=0.;  Block%Yc=0.;  Block%Zc=0.
        Block%Dst=0.; Block%Aleng=0.; Block%Aseng=0.
        Block%Vol=0.
        Block%SD=0.;  Block%Grad=0.
        
        Block%Alagm=0.
        
        Block%Dtm=0.
        Block%F=0.;  Block%Q=0.;  Block%D=0.
        Block%Src=0.; Block%Dsrc=0.
        !Block%U=0.; 
        Block%V=0.
        
        Block%VL=0.; Block%VR=0.
        !added by ydd
        Block%vibn=0.0; Block%vjbn=0.0; Block%vkbn=0.0
        Block%gridV=0.0
        Block%radSurf=0.0
        Block%thtf=0.0
        Block%thtc=0.0
        Block%rad=0.0
        Block%f_r1=1.0
        Block%f_r2=1.0
        block%rccl=0.0;block%rccr=0.0

        Block%PP=0.;  Block%T=0.
        Block%Rmiu=0.
            
        Block%Rmiudis=0.
        Block%Dqdxyz=0.
        Block%Ub =0.        !add by dzw05
        Block%Vb =0.  
        Block%Wb =0. 

        Block%Rb  =0.
        Block%aMb =0.
        Block%rMb =0.
        Block%vKb =0.

        Block%Pb =0. 
        Block%Ptb =0.
        Block%Tb  =0.

        Block%uub =0.
        Block%vvb =0.
        Block%wwb =0.
        Block%uvb =0.
        Block%uwb =0.
        Block%vwb =0.
        Block%pprms =0.
        Block%uurms =0.
        Block%vvrms =0.
        Block%wwrms =0.
        Block%uvrms =0.
        Block%uwrms =0.
        Block%vwrms =0.
        
        Block%Rds=0.
        Block%Wt0=0.; Block%Wt1=0.
 
        
        ALLOCATE(Block%F1(0:Block%NIp2, 0:Block%NJp2, 0:Block%NKp2))
        ALLOCATE(Block%F2(0:Block%NIp2, 0:Block%NJp2, 0:Block%NKp2))
        
        Block%F1=0.; Block%F2=0.
        
        !Hybrid
        ALLOCATE(Block%FunDES(0:Block%NIp2, 0:Block%NJp2, 0:Block%NKp2))
        ALLOCATE(Block%Fscheme(0:Block%NIp2, 0:Block%NJp2, 0:Block%NKp2))
        ALLOCATE(Block%shock(-2:Block%NIp2, -2:Block%NJp2, -2:Block%NKp2))
       
        Block%FunDES=0.
        Block%Fscheme=0.
        Block%shock=0.
    
    ENDDO
    
END SUBROUTINE INPUT_Read_BlockDimensions

