%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Test case using WCAWE                             %
%                                                                         %
%                             March 2020                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
% Init main program
%--------------------------------------------------------------------------


warning('off', 'MATLAB:nearlySingularMatrix');

%--------------------------------------------------------------------------
% Folders
%--------------------------------------------------------------------------

% Add folders for Mumps, WCAWE, and Mesh functions/files
global Meshfold WCAWEfold Mumps DataMap Derivatives
Meshfold = 'Matrices';
WCAWEfold = 'WCAWE';
Mumps = 'Mumps';
DataMap = 'DataMap';
Derivatives = 'Derivatives';



addpath(genpath(strcat(pwd,'/',Meshfold)));
addpath(genpath(strcat(pwd,'/',WCAWEfold)));
addpath(genpath(strcat(pwd,'/',Mumps)));
addpath(genpath(strcat(pwd,'/',DataMap)));
addpath(genpath(strcat(pwd,'/',Derivatives)));




%--------------------------------------------------------------------------
% Input data for the problem
%--------------------------------------------------------------------------

% Input parameters for Matlab calculation
flag.rerun = 1; % to recalculate FreeFem++ matrices
flag.recalculated = 1; % allow WCAWE and/or FE recalculation
flag.calculateFE = 1;  % calculate FE solution
flag.calculateWCAWE = 0; % calculate WCAWE solution

flag.plotcomparison = 0; % plot comparison between FE and WCAWE
flag.comparisonMULTI = 0;
flag.calculatePrad = 0;

flag.converge = 0;
flag.convert2VTK = 0; % convert SOLFE.mat into a .vkt file
flag.plotMQP = 0;

flag.getmatrices = 1;

if flag.converge || flag.plotMQP || flag.convert2VTK
    flag.rerun = 0;
    flag.recalculated = 0;
    flag.getmatrices = 0;
end



% Input files for mesh and FE matrices
mesh.file = 'Plate';
sizemesh = load('sizemesh.txt');
sizemesh = sizemesh(end);

% define timing struct to access to time calculation of each method                                                    
timing.freefem = 0;
timing.WCAWE = 0;
timing.computeFE = 0;                                                    

% Source
source = [0,0,0];  % coordinates of the source 
S_amp = 0.1;

% Material parameters
param.rho = 1.2;
param.rhoS = 1200;
param.c0 = 340;



% Frequency range
param.fmin = 50; % 300
param.fmax = 50; % 6000
param.f_range = [param.fmin param.fmax];
param.freqincr = 1; % 20
param.freq = param.fmin : param.freqincr : param.fmax; % frequency range
param.nfreq = length(param.freq);

% those frequencies are the frequencies point for Padé expension
param.freqref = [50];%[75 125 175 225 275 325 375 425 475 525 575];
param.nfreqref = length(param.freqref);

% Input data for the loop over expansion orders. Note that for each
% frequency sweep the number of vectors from the WCAWE basis will depend on
% the number of point for Padé expension. For instance, if we have 2 points
% for expansion, and nvecfreq=5 (order of expansion), we will have 15
% vectors in the basis, 5 by intervals.
param.nvecfreqmin = 25;
param.nvecfreqmax = 25;
param.incrvec = 20;
param.vecfreqrange = param.nvecfreqmin : param.incrvec : param.nvecfreqmax;


% path to store results in Matrices/ folder
param.idData = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']'];


% generation of the different folder to store data if they don't already exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genfolders(mesh,param);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%build sub_range with the different freqref
index_freqref = zeros(1,param.nfreqref);
for ii=1:param.nfreqref
    [~,index_freqref_ii] = min(abs(param.freq-param.freqref(ii)));
    index_freqref(ii) = index_freqref_ii;
end

caracteristic_index = [1 index_freqref param.nfreq];
param.sub_range = cell(1,param.nfreqref);
param.n_sub_range = length(param.sub_range);
if param.nfreqref==1
    param.sub_range{1} = param.freq;
else
    for ii=1:param.nfreqref
        if ii==1
            left_edge = 1;
            right_edge = floor(mean([caracteristic_index(ii+1) caracteristic_index(ii+2)]));
            param.sub_range{ii} = param.freq(left_edge) : param.freqincr : param.freq(right_edge);
        elseif ii==param.nfreqref
            left_edge = floor(mean([caracteristic_index(ii) caracteristic_index(ii+1)])) + 1;
            right_edge = param.nfreq;
            param.sub_range{ii} = param.freq(left_edge) : param.freqincr : param.freq(right_edge);
        else
            left_edge = floor(mean([caracteristic_index(ii) caracteristic_index(ii+1)])) + 1;
            right_edge = floor(mean([caracteristic_index(ii+1) caracteristic_index(ii+2)]));
            param.sub_range{ii} = param.freq(left_edge) : param.freqincr : param.freq(right_edge);
        end
    end
end


%--------------------------------------------------------------------------
%matrices calculated with Freefem++
%--------------------------------------------------------------------------
if flag.getmatrices
    matrix_names = ["M.txt","K.txt",...
                    "H.txt","Q.txt",...
                    "C.txt"];

    [FEmatrices,ndof,timing,flag] = get_matrices(timing,flag,mesh,matrix_names,param);
    Nodes = FEmatrices.Nodes;
    RHS = FEmatrices.RHS;
    LHS = FEmatrices.LHS;
    nLHS = length(LHS);
    
    coeff_LHS = {@(f) 1,@(f) -(2*pi*f)^2};
    coeff_RHS = @(f) 1;

    % size of the "reduced" system < 4*ndof
    size_system = size(LHS{1},1);

    % get coordinates of closest node to source_location
    [~,idmin] = min(sqrt((Nodes(:,1)-source(1)).^2+(Nodes(:,2)-...
                       source(2)).^2+(Nodes(:,3)-source(3)).^2));
end
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag.recalculated
    
%--------------------------------------------------------------------------
% Calculate reference finite element frequecy sweep
%--------------------------------------------------------------------------

t_0 = cputime;
if flag.calculateFE == 1
       disp('#######################');
       disp('Recalculate FE solution');
       disp('#######################');
       SOLFE = zeros(size_system,param.nfreq); %size ndof,nfreq
       % Frequency loop calculation
       for ii=1:param.nfreq
           tic;
           disp(['[FE] Frequency : ',num2str(param.freq(ii))]);
           id = initmumps;
           id = zmumps(id);
           id.JOB = 6;
           Aglob = sparse(size(LHS{1},1),size(LHS{1},2));
           for kk = 1:nLHS
              Aglob = Aglob + coeff_LHS{kk}(param.freq(ii))*LHS{kk};
           end %kk
           id.RHS = RHS;
           id = zmumps(id,Aglob);
           resp_P = id.SOL;
           id.JOB = -2; id = zmumps(id);
           SOLFE(:,ii) = resp_P;
           toc;
       end % ii
       timing.computeFE = cputime-t_0;
end

%-----------------------------------------------------------------------
% Calculate WCAWE solution
%-----------------------------------------------------------------------

for nvecfreq=param.vecfreqrange
               
    % Initialize array of RHSderiv, Cell array linear comb coefficients derivative functions
    deriv_deg = [param.nvecfreqmax];

    if exist('Derivatives/derivative_orders.mat','file') ~= 2
        disp('#################################');
        disp('Recalculate all cross derivatives');
        disp('#################################');
        create_cross_derivatives(LHS,coeff_LHS,...
                                 coeff_RHS,deriv_deg,'f');
    end
    load('Derivatives/derivative_orders.mat');
    if ~isempty(find(derivative_orders-deriv_deg<0))
        disp('#################################');
        disp('Recalculate all cross derivatives');
        disp('#################################');
        create_cross_derivatives(LHS,coeff_LHS,...
                                 coeff_RHS,deriv_deg,'f');
    end 

    
    
    coeff_deriv_fun = cell(nLHS,nvecfreq+1);
    RHScoeffderiv = cell(1,nvecfreq+1);
    [coeff_deriv_fun,RHScoeffderiv] = get_coeff_deriv_matrices(...
                            coeff_deriv_fun,RHScoeffderiv,nvecfreq,nLHS);
                        
    coeff_derivgen_fun = @(freq) cellfun(@(cellfunc) cellfunc(freq),coeff_deriv_fun);

    % RHS

    RHSderivmulti = cell(1,param.nfreqref);
    for ii=1:param.nfreqref
        RHSderiv = cell(1,nvecfreq+1);
        for kk=1:nvecfreq+1
            RHSderiv{kk} = RHScoeffderiv{kk}(param.freqref(ii))*RHS;    %derivatives of RHS at freq0
        end
        RHSderivmulti{ii}=RHSderiv;
    end
    
    

    % Fill Cell array of linear comb coefficients derivatives at
    % freqref(ii)
    coeff_deriv_multi = cell(1,param.nfreqref);
    for ii=1:param.nfreqref
        coeff_deriv_multi{ii}=coeff_derivgen_fun(param.freqref(ii));
        % Fix coeff_deriv by replacing all NaN values by 0
        tmpidxnan = find(isnan(coeff_deriv_multi{ii}));
        coeff_deriv_multi{ii}(tmpidxnan) = 0;
        tmpidxinf = find(isinf(coeff_deriv_multi{ii}));
        coeff_deriv_multi{ii}(tmpidxinf) = 0;
    end

%-----------------------------------------------------------------------
%Recalculation of WCAWE basis
%-----------------------------------------------------------------------

   if flag.calculateWCAWE == 1
       %if exist(strcat('Matrices/',mesh.file,'/','[',num2str(param.f_range),']/',...
       %                'SOLWCAWE_',mesh.file,'_nvec',num2str(param.n_sub_range*nvecfreq),...
       %                '_[',num2str(param.freqref),']','.mat'),'file')~=2
       Wtrans = [];
       for ii=1:param.nfreqref
          [Wtranstmp,Ucoeff,timing] = WCAWE_basis(LHS,coeff_deriv_multi{ii},RHSderivmulti{ii},nvecfreq,timing);
           Wtrans = [Wtrans Wtranstmp];
       end
       [uu,vv,ww] = svd(Wtrans,0);
       iiselect = find(diag(vv)>vv(1,1)*1e-15);
       Wtranssvd = uu(:,iiselect);
       nsvd = size(Wtranssvd,2);
       output = sprintf("[SVD:Info] Number of selcted vector %d/%d",nsvd,size(Wtrans,2));
       disp(output);
       [SOLWCAWE] = Solve_WCAWE(LHS,coeff_deriv_fun,RHS,Wtranssvd,param.freq);
       %end
   end
end
end
%--------------------------------------------------------------------------
% Saves
%--------------------------------------------------------------------------
% FE solution
if flag.recalculated
    if flag.calculateFE
        save(['Matrices/',mesh.file,'/',param.idData,'/SOLFE_',mesh.file,'_meshsize_',num2str(sizemesh),'.mat'],'SOLFE');
    end

    if flag.calculateWCAWE
        save(['Matrices/',mesh.file,'/',param.idData,'/SOLWCAWE_',mesh.file,'_nvec',num2str(param.n_sub_range*nvecfreq),'_[',num2str(param.freqref),']_sizemesh_',num2str(sizemesh),'.mat'],'SOLWCAWE');
    end
    
    % save FEmatrices which contains all the data of the simulation for each
    % mesh
    save(['Matrices/',mesh.file,'/',param.idData,'/DATA_sizemesh_',num2str(sizemesh),'.mat'],'FEmatrices','param');
end


%--------------------------------------------------------------------------
% Post processing
%--------------------------------------------------------------------------

if flag.converge
    clear FEmatrices SOLFE SOLWCAWE;
    sizemesh_ARRAY = load('sizemesh.txt');
    fid1 = fopen('converge1.txt','wt');
    if flag.calculateFE == 1
        for ii=1:length(sizemesh_ARRAY)
            sizemeshtmp = sizemesh_ARRAY(ii);    
            SOL = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/SOLFE_',mesh.file,'_meshsize_',num2str(sizemesh),'.mat']));
            SOL = SOL{1};
        end
    elseif flag.calculateWCAWE
        for ii=1:length(sizemesh_ARRAY)
            sizemeshtmp = sizemesh_ARRAY(ii);
            SOL = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/SOLWCAWE_',mesh.file,'_nvec',num2str(param.n_sub_range*param.nvecfreqmin),'_[',num2str(param.freqref),']_sizemesh_',num2str(sizemeshtmp),'.mat']));
            SOL = SOL{1};
        end
    end
    DATA = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/DATA_sizemesh_',num2str(sizemeshtmp),'.mat']));
    FEmatrices = DATA{1};
    param = DATA{2};
    ndof = size(FEmatrices.Nodes,1);
    ndof_acoustic = length(FEmatrices.indexp);
    %ndof_acoustic = length(FEmatrices.acoustic_nodes);
    acoustic_volume = load('acoustic_volume.txt');
    Pressure = zeros(ndof,param.nfreq);
    Pressure(FEmatrices.acoustic_nodes,:) = real(SOL(FEmatrices.indexp,:));
    MQP1 = Pressure(FEmatrices.acoustic_nodes,:)'*FEmatrices.Q*Pressure(FEmatrices.acoustic_nodes,:)/((4e-10)*acoustic_volume);
    MQP1 = 10*log10(MQP1(1:81,1:81));
    fprintf(fid1,strcat(num2str(ndof),'\t',num2str(real(diag(MQP1)')),'\n'));
    fclose(fid1);
end

if flag.plotMQP
    fid = fopen('converge.txt','rt');
    while true
        line = fgets(fid);
        if line == -1
            break;
        else
            line = str2num(strtrim(line));
            plot(param.freq,line(2:end),'DisplayName',['ndof = ' num2str(line(1))]);
            hold on
        end
    end
    
    xlabel('Frequency (Hz)');
    ylabel('Mean quadratic pressure (dB)');
    legend();
    hold off
    
    fclose(fid);
end


%--------------------------------------------------------------------------
% convert
%--------------------------------------------------------------------------
if flag.convert2VTK
    sizemesh_ARRAY = load('sizemesh.txt');
    sizemesh = sizemesh_ARRAY(end);
    if flag.calculateWCAWE
        SOLWCAWE = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/SOLWCAWE_',mesh.file,'_nvec',num2str(param.n_sub_range*param.nvecfreqmin),'_[',num2str(param.freqref),']_sizemesh_',num2str(sizemeshtmp),'.mat']));
        SOLWCAWE = SOLWCAWE{1};
        convertGEO2VTK(FEmatrices,mesh,Nodes,SOLWCAWE,param,1:1:param.nfreq)
    elseif flag.calculateFE
        SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/SOLFE_',mesh.file,'_meshsize_',num2str(sizemesh),'.mat']));
        SOLFE = SOLFE{1};
        convertGEO2VTK(FEmatrices,mesh,Nodes,SOLFE,param,1:1:param.nfreq)
    end
end

if flag.calculatePrad
    clear FEmatrices SOLFE SOLWCAWE;
    sizemeshtmp = load('sizemesh.txt');
    DATA = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/DATA_sizemesh_',num2str(sizemeshtmp),'.mat']));
    FEmatrices = DATA{1};
    param = DATA{2};
    SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/SOLFE_',mesh.file,'_meshsize_',num2str(sizemesh),'.mat']));
    SOLFE = SOLFE{1};
    Un = zeros(size(FEmatrices.Nodes,1),param.nfreq);%Un=U1, indeed, the plate is orthogonal to x component
    Un(FEmatrices.elas_nodes,:) = SOLFE(FEmatrices.indexu1,:);
    Un_red = Un(FEmatrices.coupling_nodes,:);
    Pc = zeros(size(FEmatrices.Nodes,1),param.nfreq);
    Pc(FEmatrices.acoustic_nodes,:) = SOLFE(FEmatrices.indexp,:);
    Pc_red = Pc(FEmatrices.coupling_nodes,:);
    % calculation of the radiated power
    Prad = zeros(param.nfreq,1);
    Pradtmp = FEmatrices.Radiated*Pc_red;
    for ii=1:param.nfreq
        Prad(ii) = 0.5*(1i*param.freq(ii)*Un_red(:,ii))'*Pradtmp(:,ii);
    end
end



















