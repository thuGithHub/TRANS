SUBROUTINE OUTPUT_Define_OutputFiles
    USE Global
    IMPLICIT NONE
    CHARACTER(LEN=3) strBlkno
    INTEGER:: iBlock
    CHARACTER(LEN=14):: strPrefix, strPrefixRANS
    TYPE(BlockStruct), POINTER:: Block

    if(Kind_Hybrid==0)then
        if(Kind_Model==1)then   !SST, RANS
            strPrefix = 'resu/RASST_blk'
        elseif(Kind_Model==2)then   !WD+2006, RANS
            strPrefix = 'resu/W2006_blk'
        else
            write(*,*)Myid,"Wrong Output file names for RANS!"
        endif
    else
        if(Kind_Model==1)then   !SST, Hybrid
            strPrefix = 'resu/HYBRIDSST_blk'
            strPrefixRANS = 'resu/RASST_blk'
        else
            write(*,*)Myid,"Wrong Output file names for Hybrid!"
        endif
    endif

    DO iBlock = 1, Max_Block
        Block => AllBlocks(iBlock)
        WRITE(strBlkno, '(I3.3)') M_blockids(iBlock) !iBlock
        WRITE(Block%FLHIS  , '(A,A,A,A)') trim(FilePathPrefix), strPrefix, strBlkno, ".his"
        WRITE(Block%FLHISk , '(A,A,A,A)') trim(FilePathPrefix), strPrefix, strBlkno, ".kk"
        WRITE(Block%FLHISo , '(A,A,A,A)') trim(FilePathPrefix), strPrefix, strBlkno, ".oo"
        WRITE(Block%FLHISm , '(A,A,A,A)') trim(FilePathPrefix), strPrefix, strBlkno, ".mm"
        WRITE(Block%FLHISg , '(A,A,A,A)') trim(FilePathPrefix), strPrefix, strBlkno, ".gg"

        WRITE(Block%FLDat  , '(A,A,A,A)') trim(FilePathPrefix), strPrefix, strBlkno, ".dat"
        WRITE(Block%FLDatb , '(A,A,A,A)') trim(FilePathPrefix), strPrefix, strBlkno, ".dtb"

!       WRITE(Block%FLplt  , '(A,A,A,A)') trim(FilePathPrefix), strPrefix, strBlkno, ".plt"
        WRITE(Block%FLplt  , '(A,A,A)') trim(FilePathPrefix), strPrefix, strBlkno

        WRITE(Block%FLaver  , '(A,A,A,A)') trim(FilePathPrefix), strPrefix, strBlkno, ".aver"
        WRITE(Block%FLfld  , '(A,A,A,A)') trim(FilePathPrefix), strPrefix, strBlkno, ".FLfld"

        !WRITE(Block%FLDAT_rans  , *) strPrefix, iBlock, ".dat"
        WRITE(Block%FLDat_RANS  , '(A,A,A,A)') trim(FilePathPrefix), strPrefixRANS, strBlkno, ".dat" !for RANS, Will this be invalid?

    ENDDO   

end subroutine 
