SUBROUTINE BC_Boundarys(sweep)
    USE Global
    IMPLICIT NONE
    
    INTEGER:: sweep
      !0: VL in side surf (as BC)
      !1: VL(Dq is in center, do in Boundary_blk001) in ghost cell(VL outside surf) and connect
      !2: VL,VR in side(use external to replace internal, external is original VL on surface) 
      !   and  Connect in Vexternal?

      ! need to know in which dir: I,J,K, only one side neded
      ! all VL, VR:  I dir, Jend&Kend+1 but what about cut?

    
    INTEGER:: Nblkbc
    INTEGER, POINTER:: IJKBC(:,:)
    INTEGER:: iBlock
    INTEGER:: Kind_BC
    INTEGER:: iblk_ex, iblkbd_ex
    INTEGER:: iblknum_ex

    DO iBlock = 1, Max_Block
        CALL GLOBAL_SetPointers(iBlock)
        IJKBC => ThisBlock%IJKBC
        DO Nblkbc = 1,ThisBlock%Num_block_BC
            Kind_BC = ijkBC(Nblkbc, i_kind)
!            iblk_ex = ijkBC(Nblkbc,     i_blk_ex  )
!         if (Kind_BC == 1 .and.  Kindsub_BC == 0 ) call BC_inviscidwall(Ibgn,Iend,Jbgn,Jend,Kbgn,Kend  ,sweep,IJK,minormax)
            if(Kind_BC==1) call BC_Wall(NBlkbc,sweep)       
            if(Kind_BC==2) call BC_far(NBlkbc,sweep)
            if(Kind_BC==3) call BC_symm(NBlkbc,sweep)
            if(Kind_BC==5) call BC_Inlet(NBlkbc,sweep)
            if(Kind_BC==7) call BC_Outlet(Nblkbc,sweep)
            if(Kind_BC==9) then
                iblk_ex=ijkBC(Nblkbc,i_blk_ex)
                do iblknum_ex = 1, Max_Block
!                    if(AllBlocks(iblknum_ex)%ID_Present_Blk==iblk_ex) THEN
                    if(M_blockids(iblknum_ex)==iblk_ex) THEN
!                    if(ProcNo_All(iblk_ex)==MyID) THEN
!                        iblknum_ex=BlockNo_All(iblk_ex)
                        call BC_P2P(Nblkbc,iblknum_ex,sweep)
                        EXIT !only once is enough
                    endif
                enddo
            endif
        ENDDO
    ENDDO
END SUBROUTINE BC_Boundarys
