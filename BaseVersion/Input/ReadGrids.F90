SUBROUTINE INPUT_ReadGrids
    USE Global
    IMPLICIT NONE
    INTEGER iBlock
    TYPE(BlockStruct), POINTER:: Block
    INTEGER NI,NJ,NK !SHOULD be Block%MIJK-2
    INTEGER I,J,K

    DO iBlock = 1, Max_Block
        Block => AllBlocks(iBlock)
        CALL GLOBAL_SetPointers(iBlock)
        NI=Block%NI
        NJ=Block%NJ
        NK=Block%NK
            WRITE(Block%FLgrd, '(A, A, I3.3, A)') FilePathPrefix, "grd/grid_blk", Block%ID_Present_Blk, ".dat" 
            OPEN(UNIT=13, FILE=Block%FLgrd, MODE='READ')
            read(13,*); read(13,*); read(13,*); read(13,*);
            read(13,*); read(13,*); read(13,*); read(13,*);
            read(13,*); read(13,*); read(13,*)
            DO K=1,NK
            DO J=1,NJ
            DO I=1,NI
                READ(13,*) XX(I,J,K),YY(I,J,K),ZZ(I,J,K)
            ENDDO
            ENDDO 
            ENDDO
            CLOSE(13)
            do k=1,nk
            do j=1,nj
            do i=1,ni
                xx(i,j,k) = xx(i,j,k) / alscale
                YY(i,j,k) = YY(i,j,k) / alscale
                ZZ(i,j,k) = ZZ(i,j,k) / alscale
            ENDDO
            ENDDO
            ENDDO
    ENDDO
    
END SUBROUTINE INPUT_ReadGrids

    
SUBROUTINE INPUT_ReadSurf
    USE Global
    IMPLICIT NONE
    CHARACTER(LEN=80):: well
    INTEGER:: I,J,K,N,index_pre
    INTEGER:: Ns, Ms
    REAL:: X,Y,Z

    FLsur=trim(FilePathPrefix) // "grd/Surf.dat"
    OPEN(UNIT=2, FILE=FLsur, MODE='READ')
    read(2,*) well, Num_Surface
    
    ALLOCATE(AllSurfs(Num_Surface))
    
    do N=1,Num_surface
        read(2,*)well, index_pre, Ns, Ms
        AllSurfs(N)%NSur = Ns
        AllSurfs(N)%MSur = Ms
        ALLOCATE(AllSurfs(N)%X(Ns, Ms))
        ALLOCATE(AllSurfs(N)%Y(Ns, Ms))
        ALLOCATE(AllSurfs(N)%Z(Ns, Ms))
    enddo

    read(2,*);  read(2,*);  read(2,*)
    read(2,*);  read(2,*);  read(2,*)

    do N=1,Num_surface
        read(2,*);  read(2,*);  read(2,*)
        read(2,*);  read(2,*)
        do J=1,AllSurfs(N)%MSur  !Index2(N)
        do I=1,AllSurfs(N)%NSur  !index1(N)
            read(2,*) X, Y, Z !Xsur(N,I,J),Ysur(N,I,J),Zsur(N,I,J)
            AllSurfs(N)%X(I,J) = X/alscale
            AllSurfs(N)%Y(I,J) = Y/alscale
            AllSurfs(N)%Z(I,J) = Z/alscale
        enddo
        enddo
    enddo
    close(2)
END SUBROUTINE INPUT_ReadSurf
