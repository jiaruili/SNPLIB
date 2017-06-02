classdef SNPLIB < handle
    properties (Dependent = true)
        nSNPs
        nSamples
    end
    properties
        nAxes = 1
        nThreads 
    end
    properties (SetAccess = private)
        %% Basic Information
        CHR
        RSID
        POS
        A1
        A2
        IID
        FID
        PID
        MID
        Sex
        GENO
        AF = []
        MAF = []
        CALLRATE = []
        GRM = []
        IBS = []
        HammingMatrix = []
    end
    properties (SetAccess = private)
        loadings_pca = [];
        loadings_spectralibs = [];
    end
    properties (SetAccess = private, Transient = true)
        scores_pca = [];
        scores_spectralibs = [];
        scores_mds = [];
    end
    methods
        function obj = SNPLIB(nThreads)
            if nargin<1
                obj.nThreads = feature('numcores');
                return;
            end
            obj.nThreads = nThreads;
        end
    end
    methods % Gets & Sets
        function out = get.nSNPs(obj)
            out = length(obj.RSID);
        end
        function out = get.nSamples(obj)
            out = length(obj.FID);
        end
        function updateAF(obj)
            obj.AF = calc_af_(obj.GENO,obj.nSamples,obj.nThreads);
        end
        function updateGRM(obj)
            obj.GRM = calc_grm_(obj.GENO, obj.AF, obj.nSamples, obj.nThreads);
        end
        function updateIBS(obj)
            obj.IBS = calc_ibs_(obj.GENO, obj.nSamples, obj.nThreads);
        end
        function updateHammingMatrix(obj)
            obj.HammingMatrix = calc_hamming_(obj.GENO, obj.nSamples, obj.nThreads);
        end
        function updatePcaScores(obj)
            if isempty(obj.GRM)
                updateAF(obj);
                updateGRM(obj);
            end
            [vectors, values] = eig(obj.GRM,'vector');
            [~,ind] = sort(values,'descend');
            obj.scores_pca = vectors(:,ind(1:obj.nAxes));
        end
        function updateSpectralIbsScores(obj)
            if isempty(obj.IBS)
                updateIBS(obj);
            end
            D = sum(obj.IBS);
            D = diag(D.^(-0.5));
            I = D*obj.IBS*D;
            [vectors, values] = eig(I,'vector');
            [~,ind] = sort(values,'descend');
            obj.scores_spectralibs = D*vectors(:,ind(2:obj.nAxes+1));
        end        
        function updateMdsScores(obj)
            if isempty(obj.IBS)
                updateIBS(obj);
            end
            D = (1-obj.IBS).^2;
            n = size(obj.IBS,1);
            B = bsxfun(@plus, bsxfun(@minus, bsxfun(@minus, D, sum(D,1)/n),sum(D,2)/n), sum(D(:))/(n^2))*(-0.5);
            [V,E] = eigs((B+B')./2,obj.nAxes,'LA');
            [e,i] = sort(diag(E)); e = flipud(e); i = flipud(i);
            keep = find(e > max(abs(e)) * eps(class(e))^(3/4));
            obj.scores_mds = bsxfun(@times, V(:,i(keep)), sqrt(e(keep))');
        end
        function updatePcaLoadings(obj)
            if isempty(obj.AF)
                updateAF(obj);
            end
            obj.loadings_pca = calc_pca_loadings_(obj.GENO, obj.AF, obj.nSamples, obj.nAxes, obj.nThreads);
        end
        function updateSpectralIbsLoadings(obj)
            obj.loadings_spectralibs = calc_spectral_ibs_loadings_(obj.GENO, obj.nSamples,obj.nAxes, obj.nThreads);
        end
    end
    methods
        function createSIM(obj,geno,N)
            obj.GENO = geno;
            V = size(geno,2);
            obj.CHR = ones(V,1);
            obj.RSID = (1:V)';
            obj.POS = (1:V)';
            obj.A1 = repmat({'A'},[V,1]);
            obj.A2 = repmat({'B'},[V,1]);
            obj.FID = (1:N)';
            obj.IID = (1:N)';
            obj.PID = zeros(N,1);
            obj.MID = zeros(N,1);
            obj.Sex = ones(N,1);
            obj.Sex(1:floor(N/2)) = 2;
        end
        function importPLINKDATA(obj, bfile)
            % IMPORTPLINKDATA Import Plink binary data
            %   IMPORTPLINKDATA(obj, bfile) import plink binary data
            %   (bfile.bed, bfile.bim, bfile.fam) into obj
            filename = [bfile,'.bim'];
            formatSpec = '%s%s%f%f%s%s';
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, 'EmptyValue' ,NaN, 'ReturnOnError', false);
            fclose(fileID);
            chr = dataArray{1};
            obj.CHR = str2double(chr);
            obj.CHR(strcmp(chr,'X')) = 23;
            obj.CHR(strcmp(chr,'Y')) = 24;
            obj.CHR(strcmp(chr,'XY')) = 25;
            obj.CHR(strcmp(chr,'MT')) = 26;
            obj.RSID = dataArray{2};
            obj.POS = dataArray{4};
            obj.A1 = dataArray{5};
            obj.A2 = dataArray{6};
            filename = [bfile,'.fam'];
            formatSpec = '%s%s%s%s%f%s';
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
            fclose(fileID);
            obj.FID = dataArray{1};
            obj.IID = dataArray{2};
            obj.PID = dataArray{3};
            obj.MID = dataArray{4};
            obj.Sex = dataArray{5};
            fid = fopen([bfile '.bed'],'rb');
            bin = fread(fid,inf,'uint8=>uint8');
            fclose(fid);
            L = ceil(obj.nSamples/4);
            obj.GENO = reshape(bin(4:end),[L obj.nSNPs]);
        end
        function exportPLINKDATA(obj, bfile)
            % EXPORTPLINKDATA Export to Plink binary data
            %   EXPORTPLINKDATA(obj, bfile) export to plink binary data
            %   (bfile.bed, bfile.bim, bfile.fam) from obj
            T = table(obj.CHR,obj.RSID,zeros(obj.nrSNPs,1),obj.POS,obj.A1, obj.A2);
            writetable(T, [bfile,'.txt'], 'WriteVariableNames', false, 'Delimiter', '\t');
            movefile([bfile,'.txt'],[bfile,'.bim']);
            T = table(obj.FID,obj.IID,obj.PID,obj.MID, obj.Sex, zeros(obj.nrSamples,1));
            writetable(T, [bfile,'.txt'], 'WriteVariableNames', false, 'Delimiter', '\t');
            movefile([bfile,'.txt'],[bfile,'.fam']);
            fid = fopen([bfile '.bed'],'wb');
            fwrite(fid,[108, 27,1],'uint8');
            fwrite(fid,obj.GENO,'uint8');
            fclose(fid);
        end
        function keep_snp(obj, ind_s)
            % KEEP_SNP Keep selected SNPs
            %   KEEP_SNP(obj, ind_s)
            obj.GENO = obj.GENO(:,ind_s);
            obj.CHR = obj.CHR(ind_s);
            obj.RSID = obj.RSID(ind_s);
            obj.POS = obj.POS(ind_s);
            obj.A1 = obj.A1(ind_s);
            obj.A2 = obj.A2(ind_s);
        end
    end
    methods
        function out = projectPca(obj, obj_dest)
            if class(obj_dest)~='SNPLIB'
                out = [];
                return;
            end
            if obj_dest.nSNPs ~= obj.nSNPs
                out = [];
                return;
            end
            if isempty(obj.loadings_pca)
                updatePcaLoadings(obj);
            end
            out = project_pca_(obj_dest.GENO,obj.AF,obj.loadings_pca,obj_dest.nSamples,obj.nAxes,obj.nThreads); 
        end
        function out = projectSpectralIBS(obj, obj_dest)
            if class(obj_dest)~='SNPLIB'
                out = [];
                return;
            end
            if obj_dest.nSNPs ~= obj.nSNPs
                out = [];
                return;
            end
            if isempty(obj.loadings_spectralibs)
                updateSpectralIbsLoadings(obj)
            end
            out = project_spectral_ibs_(obj_dest.GENO, obj.GENO,obj.loadings_spectralibs, obj.nSamples, obj_dest.nSamples, obj.nThreads);
        end
    end
end