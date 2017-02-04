%% TRIGSeq (v 0.1.10)
%
% A BWA-based tool for IG & TR clonotyping.
%
% TRIGSeqAlign: Performs parallel V and J gene alignment, sorts alignment
% data, identifies IG/TR rearrangements (if present) as well as CDR3
% regions (if present). TRIGSeqAnalyze also disambiguates BWA alignments by
% selecting the longest x highest MAPQ score alignments, and defaults to
% the lowest alphabetically ordered gene (or allele) based on unique
% headers. TRIGSeqAnalyze outputs a SAM-formatted file containing two 
% entries per unique input FASTQ-formatted header, one corresponding to the
% V gene and other other corresponding to the J gene alignment result.
%
% Included in the appertaining optional alignment section one to five
% data fields:
%
% cl = clonotype (defined by V-gene:J-gene:CDR3, as applicable)
% cs = cdr3 sequence startindex relative to input sequence
% ce = cdr3 sequence end index relative to input sequence 
% ts = TRIG Score
% mq = combined MAPQ score
%
% The TRIGSeq TRIGScore is a quality metric:
% TRIGScore = min(qual) * ((vloal/lenV)*nminV + (jloal/lenJ)*nminJ)
% where qual := PHRED quality
% vloal := length of the v-gene alignment 
% lenV := length of corresponding v-gene reference sequence
% nminV := (pre-calculated) the smallest number of base substitutions
% to the given reference v-gene sequence needed to produce the next most
% similar reference v-gene sequence
% jloal := length of the j-gene alignment
% lenJ := length of corresponding j-gene reference sequence
% nminJ := (pre-calculated) the smallest number of base substitutions
% to the given reference j-gene sequence needed to produce the next most
% similar reference j-gene sequence
%
% WARNING: A properly structured (and unaltered) TRIGSeq_Util folder is
% required. This folder contains the curated human V and J gene
% lists and reference sequences, the curated list of valid V-J
% gene rearrangement pairings and the curated reference data for
% calculation of the TRIGScore. Advanced users, or users wishing to apply
% this software to other species, will need to properly format their 
% respective reference files for these purposes. Contact
% etienne.mahe@medportal.ca for details.
%
% Usage:
% TRIGSeqAnalyze(<path/to/TRIGSeq_Util/folder>,'InputFile',
% <input>,'MoleculeType',<'TR' or 'IG'>,'Verbose',<optional
% 'no'>,'OutputPath',<optional path for output
% files>,'ReadLengthMax',<integer>,'ReadLengthMin',<integer>,
% 'OffTargetExclusion',<'yes','no'>)
%
% Output: 3 sam-formatted files: <filename>.vgenes.sam;
% <filename>.jgenes.sam, <filename>.<igseq/trseq>.sam
%
% The .vgenes.sam and .jgenes.sam are direct BWA outputs, whereas the 
% .<igseq/trseq>.sam is the output of the TRIGSeqAnalyze algorithm,
% containing the sam-formatted V/J gene on-target reads (in duplicate) as
% well as rearrangement data (as applicable)
%
% Written by Etienne Mahe, University of Calgary, (c) 2014-2017
% etienne.mahe@medportal.ca
% BWA is used in conjunction with this software, per open source licensing;
% usage outside of the research setting may require express permission from
% various other authors
%%
function TRIGSeqAlign(PATH_TO_TRIGSEQ_UTIL,varargin)
    
    % gather & validate input
    tic;
    p = inputParser;
    p.FunctionName = 'TRIGSeq_v0.1.10';
    defaultFile = [PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/IMGTSampleFASTQ.fastq'];
    defaultMolecule = 'IG';
    outputPath = '.';
    verbNo = {0,'no',false,'0','false','n','NO','No','N','FALSE','False'};
    validMoleculeType = {'TR','IG'};
    addParameter(p,'InputFile',defaultFile,@(x) isFileExtCorrect(x));
    addParameter(p,'MoleculeType',defaultMolecule,@(x) ...
        any(validatestring(x,validMoleculeType)));
    addParameter(p,'Verbose', 1);
    addParameter(p,'OutputPath',outputPath);
    addParameter(p,'ReadLengthMax',1000);
    addParameter(p,'ReadLengthMin',200);
    addParameter(p,'ExcludeOffTargets','No');
    parse(p,varargin{:});
    verb = ~sum(cellfun(@(y) isequal(p.Results.Verbose,y),verbNo));
    offt = ~sum(cellfun(@(y) isequal(p.Results.ExcludeOffTargets,y),verbNo));
    % prepare output files & initiate reference data
    [~,nm,~]=fileparts(p.Results.InputFile);
    if isequal(exist([p.Results.OutputPath,'/.',nm,'data'],'dir'),7)
        rmdir([p.Results.OutputPath,'/.',nm,'data'],'s');
    end
    mkdir([p.Results.OutputPath,'/.',nm,'data']);
    tempFile = [p.Results.OutputPath,'/.',nm,'data','/tempFile.txt'];
    tempHeader = [p.Results.OutputPath,'/.',nm,'data','/tempDictionary.txt'];
    if isequal(p.Results.MoleculeType,'IG')
        rearfile = [p.Results.OutputPath,'/',nm,'.igseq.sam'];
    else
        rearfile = [p.Results.OutputPath,'/',nm,'.trseq.sam'];
    end
    if strcmp(p.Results.MoleculeType,'TR')
        RefTable1 = readtable([PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/RefTable1TR.dat'],'Delimiter',',','ReadVariableNames',true,'ReadRowNames',true);
    else
        RefTable1 = readtable([PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/RefTable1IG.dat'],'Delimiter',',','ReadVariableNames',true,'ReadRowNames',true);
    end
    RefTable2 = readtable([PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/RefTable2.dat'],'Delimiter',',','ReadVariableNames',true,'ReadRowNames',true);
    
    % confirm input
    if verb
        disp(['Welcome to ',p.FunctionName]);
        disp([datestr(now),': Input file: ',p.Results.InputFile]);
        disp([datestr(now),': Molecule type: ',p.Results.MoleculeType]);
        disp([datestr(now),': Output file path: ',p.Results.OutputPath]);
        disp([datestr(now),': Temporary file path: ',p.Results.OutputPath,'/.',nm,'data']);
        disp([datestr(now),': Maximum allowed read length: ',num2str(p.Results.ReadLengthMax)]);
        disp([datestr(now),': Minimum allowed read length: ',num2str(p.Results.ReadLengthMin)]);
        disp([datestr(now),': Off target exclusion: ',p.Results.ExcludeOffTargets]);
    end
    
    % warnings
    if verb
        if p.Results.ReadLengthMax > 1000
            warning([datestr(now),': Allowance of read length > 1kb may result in incorrect alignment-based V or J gene assignment']);
        end
        if p.Results.ReadLengthMin < 200
            warning([datestr(now),': Allowance of ultra-short read length may result in low hit rate']);
        end
    end
    
    % perform alignment
    TRIGSeq_align(PATH_TO_TRIGSEQ_UTIL,p.Results.InputFile,p.Results.MoleculeType,[p.Results.OutputPath,'/',nm,'.vgenes.sam'],[p.Results.OutputPath,'/',nm,'.jgenes.sam'],verb);
    if verb
        disp([datestr(now),': Filtering alignments...']);
    end
    
    % initiate MapReduce: SLOW starts here
    dsvj = datastore({[p.Results.OutputPath,'/',nm,'.vgenes.sam'],[p.Results.OutputPath,'/',nm,'.jgenes.sam']},'Type','tabulartext','Delimiter','\t','ReadVariableNames',false,'TextscanFormats',{'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s'},'FileExtensions',{'.sam'});
    mapreduce(dsvj,@vjMapper,@vjReducer,'Display','off','OutputFolder',[p.Results.OutputPath,'/.',nm,'data']);
    
    % prepend dictionary & add PG field data
    f_rearFile = fopen(rearfile,'a+');
    f_tempHeader = fopen(tempHeader,'r');
    f_tempFile = fopen(tempFile,'r');
    while ~feof(f_tempHeader)
          tl = fgetl(f_tempHeader);
          if ~strncmp(tl,'@PG',3)
            fwrite(f_rearFile, sprintf('%s\n',tl));
          end
    end
    fwrite(f_rearFile,sprintf('%s\n',['@PG',char(9),'ID:',p.FunctionName]));
    fclose(f_tempHeader);
    while ~feof(f_tempFile)
          tl = fgetl(f_tempFile);
          fwrite(f_rearFile, sprintf('%s\n',tl));
    end
    fclose(f_tempFile);
    fclose(f_rearFile); 
    if verb
        disp([datestr(now),': Done. Total elapsed time ',num2str(toc/60),' mins.']);
    end
    % end TRIGSeqAnalyze
    
    % Begin nested functions
    % map function
    function vjMapper(data,~,intermKVStore)
        
        % place sam headers in temp header file
        ind = strncmp(data.Var1,'@',1);
        if ~isequal(sum(ind),0)
            header = data(ind,1:3);
            fid = fopen(tempHeader,'a+');
            for ii=1:sum(ind)
                t = table2array(header(ii,:));
                fprintf(fid,'%s',[t{1},char(9),t{2},char(9),t{3},char(10)]);
            end
            fclose(fid);
            data = data(~ind,:);
        end
        
        % filter out off-target reads if set
        if offt
            ind = logical(cell2mat(cellfun(@(x) isequal(4,str2num(x)),data.Var2,'UniformOutput',false)));
            data = data(~ind,:);
        end
        
        % consider only reads < ReadLengthMax & > ReadLengthMin
        l = cell2mat(cellfun(@(x) length(x),data.Var10,'UniformOutput',false));
        ind = (l<p.Results.ReadLengthMax) & (p.Results.ReadLengthMin<l);
        data = data(ind,:);
        
        % work on FASTQ headers
        intermVals = rowfun(@(x1,x2,x3,x4,x5,x6,x7) [x1,x2,x3,x4,x5,x6,x7],data,'InputVariables',{'Var2','Var3','Var4','Var5','Var6','Var10','Var11'},'OutputFormat','cell');
        addmulti(intermKVStore,data.Var1,intermVals);
    end

    % Reduce function
    function vjReducer(intermKey,intermValIter,~)
        vFlagList = {}; 
        jFlagList = {}; 
        vGeneList = {}; 
        jGeneList = {};
        vPosnList = {}; 
        jPosnList = {};
        vMapQList = {}; 
        jMapQList = {};
        vCigarList = {}; 
        jCigarList = {};
        seq = ''; 
        qual = '';
        while hasnext(intermValIter)
            value = getnext(intermValIter);
            flag = value{1};
            gene = value{2};
            posn = str2num(value{3});
            mapQ = str2num(value{4});
            cigar = value{5};
            seq = value{6};
            qual = value{7};
            
            % separate V from J
            if strfind(gene,'V')
                vFlagList = [vFlagList;flag];
                vGeneList = [vGeneList;gene];
                vPosnList = [vPosnList;posn];
                vMapQList = [vMapQList;mapQ];
                vCigarList = [vCigarList;cigar];
            else
                jFlagList = [jFlagList;flag];
                jGeneList = [jGeneList;gene];
                jPosnList = [jPosnList;posn];
                jMapQList = [jMapQList;mapQ];
                jCigarList = [jCigarList;cigar];
            end
        end
        
        % sort default to lowest numbered gene and allele
        [vGeneList,indv] = sort(vGeneList); 
        vFlagList = vFlagList(indv); 
        vPosnList = vPosnList(indv); 
        vMapQList = vMapQList(indv); 
        vCigarList = vCigarList(indv); 
        vLOAL = cellfun(@lengthOfAlignment,vCigarList,'UniformOutput',false);
        
        [jGeneList,indj] = sort(jGeneList); 
        jFlagList = jFlagList(indj); 
        jPosnList = jPosnList(indj); 
        jMapQList = jMapQList(indj); 
        jCigarList = jCigarList(indj); 
        jLOAL = cellfun(@lengthOfAlignment,jCigarList,'UniformOutput',false);
        
        % select best alignments
        [~,inxV] = max(cell2mat(vMapQList).* cell2mat(vLOAL));
        [~,inxJ] = max(cell2mat(jMapQList).* cell2mat(jLOAL));
        
        % pull out best alignment data
        vFlag = cell2mat(vFlagList(inxV));
        jFlag = cell2mat(jFlagList(inxJ));
        vGene = cell2mat(vGeneList(inxV));
        jGene = cell2mat(jGeneList(inxJ)); 
        vPosn = cell2mat(vPosnList(inxV));
        jPosn = cell2mat(jPosnList(inxJ)); 
        vMapQ = cell2mat(vMapQList(inxV)); 
        jMapQ = cell2mat(jMapQList(inxJ));
        vCigar = cell2mat(vCigarList(inxV)); 
        jCigar = cell2mat(jCigarList(inxJ)); 
        vloal = cell2mat(vLOAL(inxV)); 
        jloal = cell2mat(jLOAL(inxJ));
        
        % check for valid rearrangement
        try
            t = RefTable1(vGene,jGene);
            RearrangementStatus = logical(table2array(t));
        catch
            RearrangementStatus = 0;
        end
        
        % if rearranged, identify CDR3 (if any) and compute MAPQ & TRIG scores
        cdr3seq = ''; 
        cdr3start = []; 
        cdr3end = []; 
        vjMapQ = []; 
        TRIGscore = [];
        if RearrangementStatus
            cdr3 =  findCDR3(seq);
            cdr3seq = cdr3.Sequence;
            cdr3start = cdr3.Start;
            cdr3end = cdr3.End;
            vjMapQ = vMapQ + jMapQ;
            lenV = double(table2array(RefTable2(vGene,'Lengths')));
            nminV = double(table2array(RefTable2(vGene,'MinNValues')));
            lenJ = double(table2array(RefTable2(jGene,'Lengths')));
            nminJ = double(table2array(RefTable2(jGene,'MinNValues')));
            TRIGscore = double(min(qual))*((vloal/lenV)*nminV + (jloal/lenJ)*nminJ);
        end
        
        % proper clonotyping requires AA of CDR3:
        cdr3seq = nt2aa(cdr3seq);
        % adjust output to conform appropriately to sam format
        if length(vGene) > 3
            vgeneonly = vGene(1:end-3);
        else
            vgeneonly = vGene;
        end
        if length(jGene) > 3
            jgeneonly = jGene(1:end-3);
        else
            jgeneonly = jGene;
        end
        if isempty(vFlag)
            vFlag = 4;
        end
        if isempty(jFlag)
            jFlag = 4;
        end
        if isempty(vPosn)
            vPosn = 0;
        end
        if isempty(jPosn)
            jPosn = 0;
        end
        if isempty(vMapQ)
            vMapQ = 0;
        end
        if isempty(jMapQ)
            jMapQ = 0;
        end
        if isempty(vGene)
            vGene = '*';
        end
        if isempty(jGene)
            jGene = '*';
        end
        if isempty(vCigar)
            vCigar = '*';
        end
        if isempty(jCigar)
            jCigar = '*';
        end
        
        % write sam & include custom sam fields:
        % cl = clonotype (v & j genes only, not alleles)
        % cs = cdr3 sequence startindex relative to input sequence
        % ce = cdr3 sequence end index relative to input sequence 
        % ts = TRIG Score
        % mq = combined MAPQ score
        
        write2sam(tempFile,intermKey,num2str(vFlag),vGene,num2str(vPosn),num2str(vMapQ),vCigar,seq,qual,[vgeneonly,':',jgeneonly,':',cdr3seq],cdr3start,cdr3end,TRIGscore,vjMapQ);
        write2sam(tempFile,intermKey,num2str(jFlag),jGene,num2str(jPosn),num2str(jMapQ),jCigar,seq,qual,[vgeneonly,':',jgeneonly,':',cdr3seq],cdr3start,cdr3end,TRIGscore,vjMapQ);
    end
    
    % bwa alignment function
    function TRIGSeq_align(PATH_TO_TRIGSEQ_UTIL,INPUT_FASTQ,TRorIG_INPUT,OUTPUT_V_SAM,OUTPUT_J_SAM,v)
        tic;
        if ~v
            warning('off','all');
        end
        if ispc
            error([datestr(now),': UNIX-based Platform Required.']);
        end
        if nargin ~= 6
            error([datestr(now),': Insufficient Input Arguments']);
        end
        if ~strcmpi(TRorIG_INPUT,'TR') & ~strcmpi(TRorIG_INPUT,'IG')
            error([datestr(now),': Define TR/IG Molecule Type']);
        end
        if strcmp(TRorIG_INPUT,'TR')
            cmdV = [PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/bwa-0.7.12/bwa mem ',PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/bwa-0.7.12/TRVGenes.fasta ',INPUT_FASTQ,' > ',OUTPUT_V_SAM];
            cmdJ = [PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/bwa-0.7.12/bwa mem ',PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/bwa-0.7.12/TRJGenes.fasta ',INPUT_FASTQ,' > ',OUTPUT_J_SAM];
        else
            cmdV = [PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/bwa-0.7.12/bwa mem ',PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/bwa-0.7.12/IGVGenes.fasta ',INPUT_FASTQ,' > ',OUTPUT_V_SAM];
            cmdJ = [PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/bwa-0.7.12/bwa mem ',PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/bwa-0.7.12/IGJGenes.fasta ',INPUT_FASTQ,' > ',OUTPUT_J_SAM];
        end
        if v
            disp([datestr(now),': Running BWA...']);
        end
        [statusV,coutV] = unix(cmdV);
        [statusJ,coutJ] = unix(cmdJ);
        if statusV ~= 0 | statusJ ~= 0
            error([datestr(now),': BWA error.',' ',coutV,' ',coutJ]);
        end
        if v
            disp([datestr(now),': Done. Total Time ',num2str(toc/60),' mins.']);
        end    
    end

    % calculate length of alignment
    function o=lengthOfAlignment(INPUT_CIGAR)
        delim = {'M','I','D','N','S','H','P','=','X'};
        vars = {'M','I','D','=','X'};
        cigar = INPUT_CIGAR;
        tally = 0;
        while ~isempty(cigar)
            [token,remain]=strtok(cigar,delim);
            if ~isempty(remain)
                let = remain(1);
                if ismember(let,vars)
                    tally = tally + str2num(token);
                end
                if length(remain)>1
                    cigar = remain(2:end);
                else
                    cigar = '';
                end
            else
                cigar = '';
            end
        end
        o = tally;
    end

    % CDR3 finder
    function outstr=findCDR3(sequence)
        % C-X{5,21}-F/W
        delimTR = 'TG(C|T)(((A|C|T|G){3}){5,21})(TGG|TTC|TTT)';
        %sequence
        [strt,e] = regexpi(sequence,delimTR);
        if isequal(length(strt),0)
            sequence = seqreverse(sequence);
            [strt,e] = regexpi(sequence,delimTR);
        end
        outstr(1).Sequence = sequence(strt:e);
        outstr(1).Start = strt;
        outstr(1).End = e;
    end

    % check file extension
    function tf=isFileExtCorrect(input)
        [~,w,x]=fileparts(input);
        if isequal(x,'.fastq') | isequal(x,'.fq')
            tf = 1;
        elseif isequal(x,'.gz')
            [~,~,x]=fileparts(w);
            if isequal(x,'.fastq') | isequal(x,'.fq')
                tf = 1;
            else
                tf = 0;
            end
        else 
            tf = 0;
            
        end
    end

    % write sam
    function write2sam(file,header,flag,gene,posn,mapq,cigar,seq,qual,cl,cs,ce,ts,mq)
        fid = fopen(file,'a+');
        fprintf(fid,'%s',[header,char(9),flag,char(9),gene,char(9),posn,char(9),mapq,char(9),cigar,char(9),'*',char(9),'0',char(9),'0',char(9),seq,char(9),qual,char(9),'cl:Z:',cl,char(9),'cs:Z:',num2str(cs),char(9),'ce:Z:',num2str(ce),char(9),'ts:Z:',num2str(ts),char(9),'mq:Z:',num2str(mq),char(10)]);
        fclose(fid);
    end
end
