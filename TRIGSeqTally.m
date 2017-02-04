%% TRIGSeq (v 0.1.10)
%
% TRIGSeqTally analyzes .igseq.sam or .trseq.sam (TRIGSeqAnalyzed produced)
% files to tally unique IG or TR rearrangements by clonotype, as well as 
% computing count-based IG or TR gene coverage. From these data,
% TRIGSeqTally also computes an estimated VAF (eVAF) for each clonotype,
% based on likelihood ratio modeling (see doi: 10.1186/1471-2105-15-299). 
% TRIGSeqAnalyze will filter results by input TRIGScore or MAPQ scores. 
% This tool will also compute the fraction of clonotype-specific reads 
% showing evidence of somatic hypermutation for input .igseq.sam files
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
% <input>,'OutputPath',<optional path for output
% files>,'Verbose',<optional 'no'>,'MAPQScoreCutOff',<optional MAPQ score
% cut-off value>,'TRIGScoreCutOff',<optional TRIGscore cut-off value>,
% 'SHRatioCutOff',<optional somatic hypermutation cut-off value>)
%
% Output: <input>.<igseg/trseq>.coverage.txt & 
% <input>.<igseg/trseq>.tally.txt
%
% Written by Etienne Mahe, University of Calgary, (c) 2014-2017
% etienne.mahe@medportal.ca
%%
function TRIGSeqTally(PATH_TO_TRIGSEQ_UTIL,varargin)
    tic;
    p = inputParser;
    p.FunctionName = 'TRIGSeq_v0.1.10';
    defaultInputFile = [PATH_TO_TRIGSEQ_UTIL,'/IMGTSampleFASTQ.igseq.sam'];
    noOpts = {0,'no',false,'0','false','n','NO','No','N','FALSE','False'};
    addParameter(p,'OutputPath','.');
    addParameter(p,'InputFile',defaultInputFile);
    addParameter(p,'Verbose',1);
    addParameter(p,'MAPQScoreCutOff',0);
    addParameter(p,'TRIGScoreCutOff',0);
    addParameter(p,'SHRatioCutOff',0.02); % based on guidelines
    parse(p,varargin{:});
    verb = ~sum(cellfun(@(y) isequal(p.Results.Verbose,y),noOpts));
    if verb
        disp(['Welcome to ',p.FunctionName]);
        disp([datestr(now),': Loading file ',p.Results.InputFile]);
    end
    
    % validate input
    [~,nm,ext]=fileparts(p.Results.InputFile);
    if strcmp(ext,'.sam')
        [~,nm,ext]=fileparts(nm);
        if strcmp(ext,'.igseq')
            if verb
                disp([datestr(now),': Analysis molecule type: IG']);
            end
            vgeneref = fastaread(strcat(PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/bwa-0.7.12/IGVGenes.fasta'));
            jgeneref = fastaread(strcat(PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/bwa-0.7.12/IGJGenes.fasta'));
            outfile = [p.Results.OutputPath,'/',nm,'.igseq.tally.txt'];
            covfile = [p.Results.OutputPath,'/',nm,'.igseq.coverage.txt'];
            mol = 1;
        elseif strcmp(ext,'.trseq')
            if verb
                disp([datestr(now),': Analysis molecule type: TR']);
            end
            vgeneref = fastaread(strcat(PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/bwa-0.7.12/TRVGenes.fasta'));
            jgeneref = fastaread(strcat(PATH_TO_TRIGSEQ_UTIL,'/TRIGSeq_Util/bwa-0.7.12/TRJGenes.fasta'));
            outfile = [p.Results.OutputPath,'/',nm,'.trseq.tally.txt'];
            covfile = [p.Results.OutputPath,'/',nm,'.trseq.coverage.txt'];
            mol = 0;
        else
            error([datestr(now),': Reference to file with invalid sub-extension type; igseq or trseq file required']);
        end
    else
        error([datestr(now),': Reference to file with invalid extension type; sam file required']);
    end
    if isdir([p.Results.OutputPath,'/.',nm,'data'])
        rmdir([p.Results.OutputPath,'/.',nm,'data'],'s');
    end
    mkdir([p.Results.OutputPath,'/.',nm,'data'])
    
    % index to sam file & extract headers to hidden header file
    if verb
        disp([datestr(now),': Extracting sam data...']);
    end
    sam = BioIndexedFile('SAM',p.Results.InputFile,[p.Results.OutputPath,'/.',nm,'data'],'Verbose',verb);
    fid = fopen([p.Results.OutputPath,'/.',nm,'data/headers.txt'],'w+');
        cellfun(@(x) fprintf(fid,'%s',x,char(10)),getKeys(sam),'UniformOutput',false);
    fclose(fid);
    
    % map reduce for coverage, clonotype tally, eVAF & SH
    if verb
        disp([datestr(now),': Calculating...']);
    end
    vcats = unique(cellfun(@(x) x(1:end-3),{vgeneref(:).Header},'UniformOutput',false));
    jcats = unique(cellfun(@(x) x(1:end-3),{jgeneref(:).Header},'UniformOutput',false));   
    ds = datastore([p.Results.OutputPath,'/.',nm,'data/headers.txt'],'Type','tabulartext','Delimiter','\n','ReadVariableNames',false);
    
    % bring relevant data into memory
    dataTally = readall(mapreduce(ds,@hMapper,@hReducer,'Display','off','OutputFolder',[p.Results.OutputPath,'/.',nm,'data']));
    
    if height(dataTally)<= 2
        indcov = cell2mat(cellfun(@(x) isequal(x,'subset_cov'),dataTally.Key,'UniformOutput',false));
        coverageTable = dataTally.Value{indcov};
        writetable(coverageTable,covfile,'WriteRowNames',true,'Delimiter','\t');
            
        if verb
            disp([datestr(now),': Done: No identifiable rearrangements. Total elapsed time = ',num2str(toc/60),' mins.']);
        else
            error([datestr(now),': Done: No identifiable rearrangements.']);
        end
    else
        % pull out data indices
        indrl = cell2mat(cellfun(@(x) isequal(x,'read_length'),dataTally.Key,'UniformOutput',false));
        indcov = cell2mat(cellfun(@(x) isequal(x,'subset_cov'),dataTally.Key,'UniformOutput',false));
        
        % coverage
        coverageTable = dataTally.Value{indcov};
        
        % Reads-on-target
        ROT = sum(sum(table2array(coverageTable))) - table2array(coverageTable(1,1));
        
        % Number rearranged
        Vonly = sum(table2array(coverageTable(2:end,1)));
        Jonly = sum(table2array(coverageTable(1,2:end)));
        numRearranged = ROT - Vonly - Jonly;
        
        % eVAF
        A = 1-numel(vcats)-numel(jcats);
        dataTally = dataTally(~(indrl | indcov),:);
        numericdata = cell2mat(dataTally.Value);
        clonotypelist = dataTally.Key;
        svar = cell2mat(cellfun(@(x) table2array(coverageTable(x{1},x{2})),cellfun(@(x) strsplit(x,':'),dataTally.Key,'UniformOutput',false),'UniformOutput',false));
        B = numel(vcats) + numel(jcats) - ROT - svar;
        eVAF = (-B-sqrt(B.^2-4*A*svar)/(2*A)) .* (numericdata(:,1)./svar);
        
        % output fraction of total rearranged per clonotype
        frac = numericdata(:,1)/numRearranged;
        
        if mol
            T = array2table([eVAF,frac,numericdata(:,2:end)],'VariableNames',{'eVAF','FractionOfTotalRearranged','TRIGScore','MAPQ','SHFraction'},'RowNames',clonotypelist);
            T.Properties.DimensionNames = {'Clonotype','Data'};
            writetable(T,outfile,'WriteRowNames',true,'Delimiter','\t');
            writetable(coverageTable,covfile,'WriteRowNames',true,'Delimiter','\t');
            if verb
                disp([datestr(now),': Done. Total elapsed time = ',num2str(toc/60),' mins.']);
                disp([datestr(now),': Refer to results in file ',outfile]);
            end
        else
            T = array2table([eVAF,frac,numericdata(:,2:end-1)],'VariableNames',{'eVAF','FractionOfTotalRearranged','TRIGScore','MAPQ'},'RowNames',clonotypelist);
            T.Properties.DimensionNames = {'Clonotype','Data'};
            writetable(T,outfile,'WriteRowNames',true,'Delimiter','\t');
            writetable(coverageTable,covfile,'WriteRowNames',true,'Delimiter','\t');
            if verb
                disp([datestr(now),': Done. Total elapsed time = ',num2str(toc/60),' mins.']);
                disp([datestr(now),': Refer to results in file ',outfile]);
            end
        end
    end
    
    % begin nested functions
    % map
    function hMapper(data,~,intermKVStore)
        
        % collate sam data
        str = cell2mat(cellfun(@(x) trigseqSamDataExtractor(sam,x),data.Var1,'UniformOutput',false));
        sub_read_length_sum = sum(cell2mat({str(:).ReadLength}));
        
        % coverage
        [~,locv]=ismember({str(:).VGene},vcats);
        [~,locj]=ismember({str(:).JGene},jcats);
        locv = locv + 1; locj = locj + 1;
        subset_cov = accumarray([locv',locj'],1,[numel(vcats)+1 numel(jcats)+1])/2;
        
        % consider only rearranged entries & filter by MAPQ & TRIGScore
        % cut-offs
        indr = cell2mat(transpose(cellfun(@(x) isempty(x),strfind({str(:).clonotype},'*:'),'UniformOutput',false)));
        indt = transpose((cell2mat({str(:).trigscore}) > p.Results.TRIGScoreCutOff));
        indm = transpose((cell2mat({str(:).mapq}) > p.Results.MAPQScoreCutOff));
        ind = indr & indt & indm;
        [clonotypes,~,idx] = unique({str(ind).clonotype});
        counts = countcats(categorical(idx))/2;
        trigscore = accumarray(idx,cell2mat({str(ind).trigscore}),[numel(clonotypes) 1],@mean);
        mapq = accumarray(idx,cell2mat({str(ind).mapq}),[numel(clonotypes) 1],@mean);
        
        % somatic hypermutation assessment
        [vsht,vshl] = arrayfun(@(x) variation_from_reference(cell2mat(x)),{str(ind).VCigar},'UniformOutput',false);
        [jsht,jshl] = arrayfun(@(x) variation_from_reference(cell2mat(x)),{str(ind).JCigar},'UniformOutput',false);
        sh = (accumarray(idx,transpose((cell2mat(vsht)+cell2mat(jsht))./(cell2mat(vshl)+cell2mat(jshl))))>p.Results.SHRatioCutOff);
        
        intermKeys = [{'read_length'},{'subset_cov'},clonotypes];
        intermVals = [sub_read_length_sum,subset_cov,arrayfun(@(x) {counts(x),trigscore(x),mapq(x),sh(x)},1:numel(clonotypes),'UniformOutput',false)];
        addmulti(intermKVStore,intermKeys,intermVals);
    end

    function hReducer(intermKey,intermValIter,outKVStore)
        if isequal(intermKey,'read_length')
            read_length = 0;
            while hasnext(intermValIter)
                read_length = read_length + getnext(intermValIter);
            end
            add(outKVStore,intermKey,read_length);
        elseif isequal(intermKey,'subset_cov')
            cov = zeros(numel(vcats)+1,numel(jcats)+1);
            while hasnext(intermValIter)
                cov = cov + getnext(intermValIter);
            end
            cov = array2table(cov);
            cov.Properties.VariableNames = ['NoJGene',jcats];
            cov.Properties.RowNames = ['NoVGene',vcats];
            add(outKVStore,intermKey,cov);
        else
            counts = [];
            trigscore = [];
            mapq = [];
            sh = [];
            while hasnext(intermValIter)
                value = getnext(intermValIter);
                counts = [counts;value{1}];
                trigscore = [trigscore;value{2}];
                mapq = [mapq;value{3}];
                sh = [sh;value{4}];
            end
            add(outKVStore,intermKey,[sum(counts),mean(trigscore),mean(mapq),sum(sh)/sum(counts)]);
        end
    end
    
    function str=trigseqSamDataExtractor(samBioIndexedFile,samBioIndexedFileKey)
        s = struct2table(read(samBioIndexedFile,samBioIndexedFileKey));
        str.ReadLength = length(s.Sequence{1});
        RefCig = unique([s(:,3),s(:,6)]);
        tags = unique(struct2table(s.Tags));
        if ~(isequal(width(RefCig),height(RefCig)) || isequal(height(RefCig),1)) || ~isequal(height(tags),1) || ~gt(width(tags),0)
            error('Checksum error in TRIGSeq sam data extraction; sam file may be corrupted');
        end
        vind = cell2mat(cellfun(@(x) ~isempty(x),strfind(RefCig.ReferenceName,'V'),'UniformOutput',false));
        if sum(vind)>0
            v = RefCig.ReferenceName{vind};
            str.VGene = v(1:end-3);
            if isempty(str.VGene)
                str.VGene = '*';
            end
            j = RefCig.ReferenceName{~vind};
            str.JGene = j(1:end-3);
            if isempty(str.JGene)
                str.JGene = '*';
            end
            str.VCigar = RefCig.CigarString{vind}; 
            str.JCigar = RefCig.CigarString{~vind};
        else
            str.VGene = '*'; str.JGene = '*'; str.VCigar = '*'; str.JCigar = '*';
        end
        if ismember('cl',tags.Properties.VariableNames)
            % clonotype descriptor is present
            str.clonotype = tags.cl{1};
        else
            str.clonotype = '';
        end
        if ismember('cs',tags.Properties.VariableNames)
            % cdr3 sequence start position is present
            str.cdr3start = str2num(tags.cs{1});
        else
            str.cdr3start = 0;
        end
        if ismember('ce',tags.Properties.VariableNames)
            % cdr3 sequence end position is present
            str.cdr3end = str2num(tags.ce{1});
        else
            str.cdr3end = 0;
        end
        if ismember('ts',tags.Properties.VariableNames)
            % TRIGScore is present
            str.trigscore = str2num(tags.ts{1});
        else
            str.trigscore = -1;
        end
        if ismember('mq',tags.Properties.VariableNames)
            % combined MAPQ score is present
            str.mapq = str2num(tags.mq{1});
        else
            str.mapq = -1;
        end
    end
    
    function [tally,leng]=variation_from_reference(INPUT_CIGAR)
        delim = {'M','I','D','N','S','H','P','=','X'};
        vars = {'I','D','X'};
        alvars = {'M','I','D','=','X'};
        cigar = INPUT_CIGAR;
        tally = 0;
        leng = 0;
        while ~isempty(cigar)
            [token,remain]=strtok(cigar,delim);
            if ~isempty(remain)
                let = remain(1);
                if ismember(let,vars)
                    tally = tally + str2num(token);
                end
                if ismember(let,alvars)
                    leng = leng + str2num(token);
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
    end
end