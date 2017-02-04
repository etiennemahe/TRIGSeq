function setup_TRIGSeq_Util
    try
        cd([pwd,'/TRIGSeq_Util/bwa-0.7.12']);
        unix('make');
        here = pwd;
        unix([here,'/bwa index ',here,'/IGVGenes.fasta']);
        unix([here,'/bwa index ',here,'/IGJGenes.fasta']);
        unix([here,'/bwa index ',here,'/TRVGenes.fasta']);
        unix([here,'/bwa index ',here,'/TRJGenes.fasta']);
    catch e
        rethrow(e);
    end
end