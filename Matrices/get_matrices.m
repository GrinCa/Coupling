function [FEmatrices,ndof,timing,flag] = get_matrices(timing,flag,mesh,matrix_names,param)

isupdated = 0; % isupdated=1 : recalculate freefrem++ script

for ii=1:length(matrix_names)
    matrix_names{ii} = strcat('Matrices/',mesh.file,'/',matrix_names{ii});
    if exist(matrix_names(ii),"file")
    else
        isupdated=1;
    end
end


%--------------------------------------------------------------------------
%Run FreeFem++ script, IF NEEDED
%--------------------------------------------------------------------------

if (isupdated||flag.rerun) % EDP updated and/or mesh updated
   t_0 = cputime;
   disp('************************');
   disp('*Rerun FreeFem++ script*');
   disp('************************');
   edpcommand = strcat('FreeFem++'," ",mesh.file,'.edp');
   system(edpcommand);
   timing.freefem = cputime-t_0;
   disp('*********************************************************');
   output = sprintf('[Get_matrices:infos] CPUtime for building of matrices %.4f s',timing.freefem);
   disp(output);
   disp('*********************************************************');
end

%--------------------------------------------------------------------------
% Get matrices from the files
%--------------------------------------------------------------------------

listLHS = cell(1,length(matrix_names)); 
size_matrices = zeros(length(matrix_names),3);

% Matrices of the FE problem

for ii=1:length(matrix_names)
    fid = fopen(matrix_names(ii),'rt');
    for jj=1:3
        line = fgets(fid);
        if jj==3
           values = str2num(strtrim(line));
           size_matrices(ii,1) = values(1);
           size_matrices(ii,2) = values(2);
        end
    end
    fclose(fid);
    matrix_data = importdata(matrix_names(ii),' ',3);
    matrix_data = matrix_data.data;
    listLHS{ii} = sparse([matrix_data(:,1);size_matrices(ii,1)]+1,[matrix_data(:,2);size_matrices(ii,1)]+1,[matrix_data(:,3);0]);
end

% RHS
RHSdata = importdata(strcat(mesh.file,'/',"RHS.txt")," ",3);
RHSdata = RHSdata.data;
RHS = diag(sparse(RHSdata(:,1)+1,RHSdata(:,2)+1,RHSdata(:,3)));


% Nodes
Nodes = load(strcat(mesh.file,'/',"Nodes.txt"));
ndof = size(Nodes,1);



%--------------------------------------------------------------------------
% return
%--------------------------------------------------------------------------

FEmatrices.Nodes = Nodes; 

%plot3(FEmatrices.Nodes(find(RHS),1),FEmatrices.Nodes(find(RHS),2),FEmatrices.Nodes(find(RHS),3),'+');

FEmatrices = build_global(FEmatrices,listLHS,RHS,param,mesh.file);



end




function FEmatrices = build_global(FEmatrices,listLHS,RHStmp,param,FILENAME)

ndof = size(FEmatrices.Nodes,1);
%definition of the regions (for tetrahedra) and labels( for triangles) in 3D
acoustic_region = 1;
plate_region = 2;
coupling_label = 4;
ext_pressure = 5;


M = listLHS{1};
K = listLHS{2};
H = listLHS{3};
Q = listLHS{4};
%RBC = listLHS{5};
Ctmp = listLHS{5};


% get the arrays of the nodes belonging to a given region/label
tab_region = get_region([acoustic_region,plate_region],ndof,FILENAME);
FEmatrices.acoustic_nodes = find(tab_region(:,1));
FEmatrices.elas_nodes = find(tab_region(:,2));

labels_cell = get_labels([coupling_label,ext_pressure],...
                         ndof,...
                         FILENAME);
FEmatrices.coupling_nodes = find(labels_cell{1});
FEmatrices.extPressure_nodes = find(labels_cell{2});

tab_elas = [];
for ii=1:length(FEmatrices.elas_nodes)
    for jj=1:3
        tab_elas = [tab_elas 3*(FEmatrices.elas_nodes(ii)-1)+jj];
    end
end

FEmatrices.indexu1 = 1:3:length(tab_elas);
FEmatrices.indexu2 = 2:3:length(tab_elas);
FEmatrices.indexu3 = 3:3:length(tab_elas);
FEmatrices.indexp = length(tab_elas)+1:1:(length(tab_elas)+length(FEmatrices.acoustic_nodes));

M = M(tab_elas,tab_elas);
K = K(tab_elas,tab_elas);
%RBC = RBC(tab_acoustic,tab_acoustic);
H = H(FEmatrices.acoustic_nodes,FEmatrices.acoustic_nodes);% + 1i*RBC;
Q = Q(FEmatrices.acoustic_nodes,FEmatrices.acoustic_nodes);
C = sparse(size(M,1),size(H,2));
C(FEmatrices.indexu1,:) = Ctmp(FEmatrices.elas_nodes,FEmatrices.acoustic_nodes);
%FEmatrices.Radiated = Ctmp(FEmatrices.coupling_nodes,FEmatrices.coupling_nodes);

Kglob = sparse([K -C;sparse(size(C',1),size(C',2)) H]);
Mglob = sparse([M sparse(size(C,1),size(C,2));param.rho*C' Q/param.c0^2]);

FEmatrices.RHS = zeros(size(Mglob,1),1);
FEmatrices.RHS(FEmatrices.indexu1) = RHStmp(FEmatrices.elas_nodes);
FEmatrices.LHS = {Kglob,Mglob};

end


function tab_region = get_region(region_array,ndof,FILENAME)

connectivity_table = load(['Matrices/',FILENAME,'/connectivity_table.txt']);
region_element = load(['Matrices/',FILENAME,'/regions.txt']);

tab_region = zeros(ndof,length(region_array));

for ii=region_array
    id_elements = find(region_element == region_array(ii));
    for jj=1:length(id_elements)
        tab_region(connectivity_table(id_elements(jj),:)+1,ii) = 1;
    end
end

end

function labels_cell = get_labels(label_number,ndof,FILENAME)

tab_labels = load(['Matrices/',FILENAME,'/labels.txt']);
labels_cell = cell(1,length(label_number));

for ii=1:length(label_number)
    jj = find(tab_labels(1,:) == label_number(ii));
    labels_cell{ii} = tab_labels(2:end,jj);
end

end


