function convertGEO2VTK(FEmatrices,mesh,Nodes,SOL,param,index)

%##########################################################################
%Please read the following lines for further informations
%##########################################################################
%this function aims to convert "SOL" array calculated with FE or WCAWE
%method into a .vtk file which Paraview can read. To understand how this
%scripts work, it is important to know what is the structure of .vtk file.
%You may find good info in the "vtk_fil_documentation.pdf" in the
%Documentation folder.
%We firstly need the Nodes file in order to get the coordinates of each
%node. Then we need SOL array (size = (i*ndof,nfreq), i=number of function). 
%It is possible thatwith the time, version compatibilities are no longer 
%working. The main probleme of this type of file(.vtk) is that Paraview 
%won't give you accurate information if it fails to read.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To sum up .vtk file contains(in this particular case):
% -header: gives the version of the vtk file
% -DATA_SET type : in our case type=UNSCTRUCTURED_GRID which enables us to
% choose as we want triangle elements for 2D or tetrahedral elements for 3D
% -POINTS : content of the Nodes file (coordinates)
% -CELLS : which refers to the .msh file. It contains the ID of the nodes
% of each elements
% -CELLS_TYPE : contains the ID of the element to use, according .vtk
% files. For instance triangles ID=5, and tetra ID=10.
% -POINT_DATA : can take several data, from scalar for pressure to tensor
% for stress field. In our case it is SCALAR id.


% FILENAME refers to .msh file and it enables us to create the path to
% store .vtk files
FILENAME = mesh.file;
connectivity_table = load(['Matrices/',FILENAME,'/connectivity_table.txt']);
connectivity_table = connectivity_table(:,[1 2 3 4 5 8 6 7 9 10]);
text_field = get_text_field(Nodes,connectivity_table,FILENAME);
convertVTK3(FEmatrices,text_field,Nodes,SOL,FILENAME,param,index)

end


function convertVTK3(FEmatrices,text_field,Nodes,SOL,FILENAME,param,index)
% -This function create as much .vtk file as we have frequencies.
% -text field contains all the text data that every .vtk file needs,
% and which is the same no matter the frequency considerate. It was
% therefore useful to get it once and for all, print it in each
% file, instead of "recalculating" it every loop iteration...
% -index is given by the user at the call of convertGEO2VTK and contains
% the IDs of the frequencies (regarding freq array) that we want to write 
% in memory.
ndof = size(Nodes,1);
disp('***Converting DATA to 3D VTK file***');

U1 = real(SOL(FEmatrices.indexu1,:));
U2 = real(SOL(FEmatrices.indexu2,:));
U3 = real(SOL(FEmatrices.indexu3,:));
U = cell(1,length(index));

for ii=index
    U{ii} = zeros(ndof,3);
    U{ii}(FEmatrices.elas_nodes,:) = [U1(:,ii) U2(:,ii) U3(:,ii)];
end

Pressure = zeros(ndof,length(index));
Pressure(FEmatrices.acoustic_nodes,:) = real(SOL(FEmatrices.indexp,:));

for ii=index
    file_name = strcat('DataMap/',FILENAME,'/',FILENAME,'_freq_',num2str(param.freq(ii)),'Hz.vtk');
    %if exist(file_name,'file') ~= 2
        fileID = fopen(file_name,'wt');
        fprintf(fileID,text_field);
        text_data = [];
        for jj=1:ndof
            text_data = [text_data [num2str(Pressure(jj,ii)) '\n']];
        end
        %------------------------------------------------------------------
        %if U ~= 0 % if there is no displacement, the convention is U = 0
            text_data = [text_data 'VECTORS DISPLACEMENT float\n'];
            %text_data = [text_data 'LOOKUP_TABLE default\n'];
            for jj=1:ndof
                text_data = [text_data num2str(U{ii}(jj,1)) ' '...
                                       num2str(U{ii}(jj,2)) ' '...
                                       num2str(U{ii}(jj,3)) '\n'];
            end
        %end
        %------------------------------------------------------------------
        fprintf(fileID,text_data);
        fclose(fileID);
    %end
end
end


function text_field = get_text_field(Nodes,connectivity_table,FILENAME)
ndof = size(Nodes,1);
text_field = [];
% -text field contains all the text data that every .vtk file needs,
% and which is the same no matter the frequency considerate. It was
% therefore useful to get it once and for all, print it in each
% file, instead of "recalculating" it every loop iteration...
text_field = [text_field ['# vtk DataFile Version 2.0\n',FILENAME,'\nASCII\n']];
% above is the header refering to the version of vtk. It might have
% changed...
text_field = [text_field 'DATASET UNSTRUCTURED_GRID\n'];       %refer to doc
text_field = [text_field ['POINTS ',num2str(ndof),' float\n']];%refer to doc
% wrinting coordinates of each node
for ii=1:ndof
    text_field = [text_field [num2str(Nodes(ii,1)),' ',...
                              num2str(Nodes(ii,2)),' ',...
                              num2str(Nodes(ii,3)),'\n']];
end

disp('***Initialize conversion 3D***');
text_field = [text_field ['CELLS ',num2str(size(connectivity_table,1)),' ',...
                                   num2str(11*size(connectivity_table,1)),'\n']];
for ii=1:size(connectivity_table,1)
    text_field = [text_field ['10 ',num2str(connectivity_table(ii,:))],'\n'];
end 
text_field = [text_field ['CELL_TYPES ',num2str(size(connectivity_table,1)),'\n']];
for ii=1:size(connectivity_table,1)
    text_field = [text_field '24\n'];
end

text_field = [text_field ['POINT_DATA ' num2str(ndof) '\n']];
text_field = [text_field 'SCALARS PRESSURE float 1\n'];
text_field = [text_field 'LOOKUP_TABLE default\n'];
end


